function result = train_ozone_cnn(options)
% Train a 1-D CNN on paired lidar signals and ozone profiles.

arguments
    options.DataFolder (1,1) string = ""
    options.OutputFolder (1,1) string = ""
    options.AltitudeBounds_m (1,2) double = [1000 2500]
    options.NumBins (1,1) double = 200
    options.FirstBinAltitude_m (1,1) double = 3.75
    options.BinSpacing_m (1,1) double = 3.75
    options.MaxEpochs (1,1) double = 400
    options.MiniBatchSize (1,1) double = 8
    options.InitialLearnRate (1,1) double = 3e-4
    options.ValidationPatience (1,1) double = 80
    options.FeatureSet (1,1) string {mustBeMember(options.FeatureSet,["basic","dial"])} = "dial"
    options.Seed (1,1) double = 7
    options.ShowPlots (1,1) logical = true
end

%% Settings
script_folder = string(fileparts(mfilename("fullpath")));
if options.DataFolder == ""
    options.DataFolder = fullfile(fileparts(script_folder),"profiles");
end
if options.OutputFolder == ""
    options.OutputFolder = fullfile(script_folder,"outputs","correct_altitude");
end

% Stored profiles start at 3.75 m, with 3.75 m between bins.
rng(options.Seed);
altitude_m = linspace(options.AltitudeBounds_m(1),options.AltitudeBounds_m(2),options.NumBins)';
n_features = 2;
if options.FeatureSet == "dial"
    n_features = 6;
end

%% Load profiles
files = dir(fullfile(options.DataFolder,"o3_processed_*_prof*.mat"));
X = zeros(options.NumBins,n_features,0);
Y = zeros(options.NumBins,1,0);
day = strings(0,1);
paths = strings(0,1);
rejectedFile = strings(0,1);
rejectedReason = strings(0,1);

for i = 1:numel(files)
    profile_path = string(fullfile(files(i).folder,files(i).name));
    try
        date_token = regexp(files(i).name,'^o3_processed_(\d{8})_prof\d+\.mat$','tokens','once');
        data = load(profile_path,"mg_on_fr","mg_off_fr","N_o3_fr_mg");
        assert(isvector(data.mg_on_fr) && isvector(data.mg_off_fr) && isvector(data.N_o3_fr_mg),"Each file must contain one profile.");
        z = options.FirstBinAltitude_m + (0:numel(data.mg_on_fr)-1)' * options.BinSpacing_m;
        input = prepare_lidar_input(data.mg_on_fr,data.mg_off_fr,z,altitude_m,options.FeatureSet);
        target = double(data.N_o3_fr_mg(:));
        assert(numel(target) == numel(z),"Ozone and signals have different grids.");

        first = find(z <= altitude_m(1),1,"last");
        last = find(z >= altitude_m(end),1,"first");
        assert(all(isfinite(target(first:last))) && all(target(first:last) > 0),"Ozone labels must be finite and positive in the retrieval interval.");
        target = log(interp1(z(first:last),target(first:last),altitude_m));
        profile_day = string(date_token{1});

        X(:,:,end+1) = input; %#ok<AGROW>
        Y(:,:,end+1) = target; %#ok<AGROW>
        day(end+1,1) = profile_day; %#ok<AGROW>
        paths(end+1,1) = profile_path; %#ok<AGROW>
    catch exception
        rejectedFile(end+1,1) = profile_path; %#ok<AGROW>
        rejectedReason(end+1,1) = string(exception.message); %#ok<AGROW>
    end
end

rejections = table(rejectedFile,rejectedReason);
fprintf("Loaded %d profiles; rejected %d.\n",numel(day),height(rejections));
if ~isempty(rejections)
    disp(rejections);
end

%% Training, validation and test days
days = unique(day);
assert(numel(days) >= 3,"Need at least three valid measurement days.");
days = days(randperm(numel(days)));
n_train_days = min(numel(days)-2,max(1,floor(0.70*numel(days))));
n_val_days = min(numel(days)-n_train_days-1,max(1,floor(0.15*numel(days))));

is_train = ismember(day,days(1:n_train_days));
is_val = ismember(day,days(n_train_days+(1:n_val_days)));
is_test = ~(is_train | is_val);
split = repmat("test",numel(day),1);
split(is_train) = "train";
split(is_val) = "validation";
manifest = table(paths,day,split);
fprintf("Profiles: %d train, %d validation, %d test.\n",sum(is_train),sum(is_val),sum(is_test));

%% Normalize
% Arrays are altitude bins x features x profiles. Learn scaling from training data only.
scaling.inputMean = mean(X(:,:,is_train),[1 3]);
scaling.inputStd = max(std(X(:,:,is_train),0,[1 3]),1e-6);
scaling.targetMean = mean(Y(:,:,is_train),"all");
scaling.targetStd = max(std(Y(:,:,is_train),0,"all"),1e-6);
X_scaled = single((X - scaling.inputMean) ./ scaling.inputStd);
Y_scaled = single((Y - scaling.targetMean) ./ scaling.targetStd);

%% Network
layers = [
    sequenceInputLayer(n_features,Normalization="none")

    convolution1dLayer(5,24,DilationFactor=1,Padding="same")
    reluLayer

    convolution1dLayer(5,24,DilationFactor=2,Padding="same")
    reluLayer

    convolution1dLayer(5,24,DilationFactor=4,Padding="same")
    reluLayer

    convolution1dLayer(5,24,DilationFactor=8,Padding="same")
    reluLayer

    convolution1dLayer(1,1)
];

%% Training options
plots = "none";
if options.ShowPlots
    plots = "training-progress";
end
steps_per_epoch = max(1,floor(sum(is_train)/options.MiniBatchSize));

training = trainingOptions("adam", ...
    InitialLearnRate=options.InitialLearnRate, ...
    MaxEpochs=options.MaxEpochs, ...
    MiniBatchSize=options.MiniBatchSize, ...
    Shuffle="every-epoch", ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropPeriod=60, ...
    LearnRateDropFactor=0.5, ...
    GradientThreshold=1, ...
    ValidationData={X_scaled(:,:,is_val),Y_scaled(:,:,is_val)}, ...
    ValidationFrequency=steps_per_epoch, ...
    ValidationPatience=options.ValidationPatience, ...
    OutputNetwork="best-validation", ...
    InputDataFormats="TCB", ...
    TargetDataFormats="TCB", ...
    ExecutionEnvironment="cpu", ...
    Plots=plots);

%% Train
% MSE is calculated on standardized log ozone density.
[net,history] = trainnet(X_scaled(:,:,is_train),Y_scaled(:,:,is_train),layers,"mse",training);

%% Validation
validation_scaled = minibatchpredict(net,X_scaled(:,:,is_val),InputDataFormats="TCB",OutputDataFormats="TCB");
validationPredicted = exp(double(validation_scaled)*scaling.targetStd + scaling.targetMean);
validationObserved = exp(Y(:,:,is_val));
validationLogMSE = mean((double(validation_scaled)-double(Y_scaled(:,:,is_val))).^2,"all");

predicted_difference = diff(validationPredicted,1,1);
observed_difference = diff(validationObserved,1,1);
validationDifferenceRMSE = sqrt(mean((predicted_difference-observed_difference).^2,"all"));
validationStructureRatio = sqrt(mean(predicted_difference.^2,"all")) / max(sqrt(mean(observed_difference.^2,"all")),eps);
fprintf("Validation log MSE: %.4f; adjacent-bin structure ratio: %.3f\n",validationLogMSE,validationStructureRatio);

%% Test set evaluation
prediction_scaled = minibatchpredict(net,X_scaled(:,:,is_test),InputDataFormats="TCB",OutputDataFormats="TCB");
predicted = exp(double(prediction_scaled)*scaling.targetStd + scaling.targetMean);
observed = exp(Y(:,:,is_test));
assert(all(isfinite(predicted),"all"),"Nonfinite predictions: inspect training.");

% Compare with the mean training ozone profile.
baseline = repmat(mean(exp(Y(:,:,is_train)),3),1,1,sum(is_test));
model = ["CNN"; "Training mean profile"];
rmse = [sqrt(mean((predicted-observed).^2,"all")); sqrt(mean((baseline-observed).^2,"all"))];
mae = [mean(abs(predicted-observed),"all"); mean(abs(baseline-observed),"all")];
bias = [mean(predicted-observed,"all"); mean(baseline-observed,"all")];
metrics = table(model,rmse,mae,bias);
disp("Held-out metrics in ozone number-density units (m^-3):");
disp(metrics);

%% Save results
result.net = net;
result.scaling = scaling;
result.altitude_m = altitude_m;
result.options = options;
result.history = history;
result.metrics = metrics;
result.manifest = manifest;
result.rejections = rejections;
result.observed = observed;
result.predicted = predicted;
result.baseline = baseline;
result.validationObserved = validationObserved;
result.validationPredicted = validationPredicted;
result.validationLogMSE = validationLogMSE;
result.validationDifferenceRMSE = validationDifferenceRMSE;
result.validationStructureRatio = validationStructureRatio;

if ~isfolder(options.OutputFolder)
    mkdir(options.OutputFolder);
end
save(fullfile(options.OutputFolder,"ozone_cnn.mat"),"-struct","result");
writetable(manifest,fullfile(options.OutputFolder,"split_manifest.csv"));
writetable(rejections,fullfile(options.OutputFolder,"rejected_profiles.csv"));
writetable(metrics,fullfile(options.OutputFolder,"test_metrics.csv"));

if options.ShowPlots
    plot_test_profiles(result,OutputFolder=options.OutputFolder);
end
fprintf("Saved model and evaluation to %s\n",options.OutputFolder);
end
