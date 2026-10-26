function result = train_ozone_cnn(options)
%TRAIN_OZONE_CNN A small, supervised ozone profile retrieval, built from scratch.
% Run: result = train_ozone_cnn;
% Quick pipeline check: result = train_ozone_cnn(MaxEpochs=2, ShowPlots=false);
% Requires MATLAB R2024a+ and Deep Learning Toolbox (trainnet).
% No pretrained network or functions from the older training scripts are used.

arguments
    options.DataFolder (1,1) string = ""
    options.OutputFolder (1,1) string = ""
    options.AltitudeBounds_m (1,2) double = [500 2500]
    options.NumBins (1,1) double {mustBeInteger,mustBeGreaterThan(options.NumBins,1)} = 200
    options.FirstBinAltitude_m (1,1) double = 1
    options.BinSpacing_m (1,1) double {mustBePositive} = 1
    options.MaxEpochs (1,1) double {mustBeInteger,mustBePositive} = 400
    options.MiniBatchSize (1,1) double {mustBeInteger,mustBePositive} = 8
    options.InitialLearnRate (1,1) double {mustBePositive} = 3e-4
    options.ValidationPatience (1,1) double {mustBeInteger,mustBePositive} = 80
    options.FeatureSet (1,1) string {mustBeMember(options.FeatureSet,["basic","dial"])} = "dial"
    options.Seed (1,1) double = 7
    options.ShowPlots (1,1) logical = true
end

%% 1. Configuration and data locations
% Resolve paths relative to THIS file, so the working directory does not matter.
here = string(fileparts(mfilename('fullpath')));
if options.DataFolder == ""
    options.DataFolder = fullfile(fileparts(here), 'profiles');
end
if options.OutputFolder == ""
    options.OutputFolder = fullfile(here, 'outputs', 'structured');
end
assert(all(isfinite(options.AltitudeBounds_m)) && ...
    options.AltitudeBounds_m(2) > options.AltitudeBounds_m(1), ...
    'AltitudeBounds_m must contain increasing finite bounds.');
assert(isfinite(options.FirstBinAltitude_m) && isfinite(options.BinSpacing_m), ...
    'The input altitude grid must be finite.');
rng(options.Seed);
altitude_m = linspace(options.AltitudeBounds_m(1), ...
    options.AltitudeBounds_m(2), options.NumBins)';

%% 2. Load paired signals and labels; reject unusable profiles explicitly
% Each file has mg_on_fr, mg_off_fr, N_o3_fr_mg, each a column vector.
% Existing project convention: array index 500 corresponds to 500 metres.
% Change FirstBinAltitude_m/BinSpacing_m if your actual grid differs.
files = dir(fullfile(options.DataFolder, 'o3_processed_*_prof*.mat'));
assert(~isempty(files), 'No profile MAT files found in %s.', options.DataFolder);
nChannels = 2;
if options.FeatureSet == "dial", nChannels = 6; end
X = zeros(options.NumBins, nChannels, 0);
Y = zeros(options.NumBins, 1, 0);
day = strings(0,1);
paths = strings(0,1);
rejectedFile = strings(0,1);
rejectedReason = strings(0,1);
for k = 1:numel(files)
    path = string(fullfile(files(k).folder, files(k).name));
    try
        token = regexp(files(k).name, ...
            '^o3_processed_(\d{8})_prof\d+\.mat$', 'tokens', 'once');
        assert(~isempty(token), 'Filename must identify its measurement day.');
        s = load(path, 'mg_on_fr', 'mg_off_fr', 'N_o3_fr_mg');
        assert(isfield(s,'mg_on_fr') && isfield(s,'mg_off_fr') && ...
            isfield(s,'N_o3_fr_mg'), 'Missing signal or ozone variable.');
        assert(isvector(s.mg_on_fr) && isvector(s.mg_off_fr) && ...
            isvector(s.N_o3_fr_mg), 'Each file must contain one profile.');
        z = options.FirstBinAltitude_m + ...
            (0:numel(s.mg_on_fr)-1)' * options.BinSpacing_m;
        input = prepare_lidar_input(s.mg_on_fr, s.mg_off_fr, z, altitude_m, options.FeatureSet);
        target = double(s.N_o3_fr_mg(:));
        assert(numel(target) == numel(z), 'Ozone and signals have different grids.');
        first = find(z <= altitude_m(1), 1, 'last');
        last = find(z >= altitude_m(end), 1, 'first');
        assert(all(isfinite(target(first:last))) && all(target(first:last) > 0), ...
            'Ozone labels must be finite and positive in the retrieval interval.');
        % Interpolate physical density BEFORE taking its logarithm.
        target = log(interp1(z(first:last), target(first:last), altitude_m));
        X(:,:,end+1) = input; %#ok<AGROW>
        Y(:,:,end+1) = target; %#ok<AGROW>
        day(end+1,1) = string(token{1}); %#ok<AGROW>
        paths(end+1,1) = path; %#ok<AGROW>
    catch exception
        rejectedFile(end+1,1) = path; %#ok<AGROW>
        rejectedReason(end+1,1) = string(exception.message); %#ok<AGROW>
    end
end
rejections = table(rejectedFile, rejectedReason);
fprintf('Loaded %d profiles; rejected %d.\n', numel(day), height(rejections));
if ~isempty(rejections), disp(rejections); end

%% 3. Split by DAY, before learning any scaling parameters
% Profiles from the same day are correlated. Keep each day in only one set.
days = unique(day);
assert(numel(days) >= 3, 'Need at least three valid measurement days.');
days = days(randperm(numel(days)));
nTrainDays = min(numel(days)-2, max(1, floor(0.70*numel(days))));
nValidationDays = min(numel(days)-nTrainDays-1, max(1, floor(0.15*numel(days))));
isTrain = ismember(day, days(1:nTrainDays));
isValidation = ismember(day, days(nTrainDays+(1:nValidationDays)));
isTest = ~(isTrain | isValidation);
split = repmat("test", numel(day), 1);
split(isTrain) = "train";
split(isValidation) = "validation";
manifest = table(paths, day, split);
fprintf('Profiles: %d train, %d validation, %d test.\n', ...
    sum(isTrain), sum(isValidation), sum(isTest));

%% 4. Standardize using TRAINING observations only
% X shape is T x C x B: altitude (time axis), channels, profile batch.
% One mean/std per signal channel, pooled across training bins and profiles.
scaling.inputMean = mean(X(:,:,isTrain), [1 3]);
scaling.inputStd = max(std(X(:,:,isTrain), 0, [1 3]), 1e-6);
scaling.targetMean = mean(Y(:,:,isTrain), 'all');
scaling.targetStd = max(std(Y(:,:,isTrain), 0, 'all'), 1e-6);
Xscaled = single((X - scaling.inputMean) ./ scaling.inputStd);
Yscaled = single((Y - scaling.targetMean) ./ scaling.targetStd);

%% 5. Define a small multi-scale 1-D convolutional network
% A filter slides along altitude, combining neighboring bins and both signals.
% 'same' padding preserves length. No pooling: we need an output at EVERY bin.
% Dilation spaces filter samples farther apart: 1, 2, 4, then 8 bins.
% Four 5-bin convolutions give a 61-bin receptive field (~603 m by default).
% No pooling or extra smoothing removes the output structure.
% The final 1-bin convolution maps 24 learned features to one ozone value.
layers = [
    sequenceInputLayer(nChannels, Normalization='none', Name='signals')
    convolution1dLayer(5, 24, Padding='same', Name='conv1')
    reluLayer(Name='relu1')
    convolution1dLayer(5, 24, DilationFactor=2, Padding='same', Name='conv2')
    reluLayer(Name='relu2')
    convolution1dLayer(5, 24, DilationFactor=4, Padding='same', Name='conv3')
    reluLayer(Name='relu3')
    convolution1dLayer(5, 24, DilationFactor=8, Padding='same', Name='conv4')
    reluLayer(Name='relu4')
    convolution1dLayer(1, 1, Name='ozone')
];

%% 6. Train on standardized log density, select by validation loss
% MSE in log space emphasizes relative errors; it is not physical-density MSE.
plots = 'none';
if options.ShowPlots, plots = 'training-progress'; end
% trainnet discards a final incomplete training batch. Count actual updates.
stepsPerEpoch = max(1, floor(sum(isTrain)/options.MiniBatchSize));
training = trainingOptions('adam', ...
    InitialLearnRate=options.InitialLearnRate, MaxEpochs=options.MaxEpochs, ...
    LearnRateSchedule='piecewise', LearnRateDropPeriod=60, LearnRateDropFactor=0.5, ...
    GradientThreshold=1, L2Regularization=1e-4, ...
    MiniBatchSize=options.MiniBatchSize, Shuffle='every-epoch', ...
    ValidationData={Xscaled(:,:,isValidation), Yscaled(:,:,isValidation)}, ...
    ValidationFrequency=stepsPerEpoch, ValidationPatience=options.ValidationPatience, ...
    OutputNetwork='best-validation', ...
    InputDataFormats='TCB', TargetDataFormats='TCB', ...
    ExecutionEnvironment='cpu', Plots=plots, Verbose=true);
[net, history] = trainnet(Xscaled(:,:,isTrain), Yscaled(:,:,isTrain), ...
    layers, 'mse', training);

%% 7. Evaluate once on held-out days, in the original number-density units
predictionScaled = minibatchpredict(net, Xscaled(:,:,isTest), ...
    InputDataFormats='TCB', OutputDataFormats='TCB');
predicted = exp(double(predictionScaled)*scaling.targetStd + scaling.targetMean);
observed = exp(Y(:,:,isTest));
% Record validation performance/structure as well as aggregate loss, so that
% future architecture choices need not be based on the test profiles.
validationScaled = minibatchpredict(net,Xscaled(:,:,isValidation), ...
    InputDataFormats='TCB',OutputDataFormats='TCB');
validationPredicted = exp(double(validationScaled)*scaling.targetStd+scaling.targetMean);
validationObserved = exp(Y(:,:,isValidation));
validationLogMSE = mean((double(validationScaled)-double(Yscaled(:,:,isValidation))).^2,'all');
% Compare adjacent-bin changes: small values can reveal an over-smooth model.
validationDifferenceRMSE = sqrt(mean((diff(validationPredicted,1,1)- ...
    diff(validationObserved,1,1)).^2,'all'));
validationStructureRatio = sqrt(mean(diff(validationPredicted,1,1).^2,'all')) / ...
    max(sqrt(mean(diff(validationObserved,1,1).^2,'all')),eps);
fprintf('Validation log MSE: %.4f; adjacent-bin structure ratio: %.3f\n', ...
    validationLogMSE,validationStructureRatio);
assert(all(isfinite(predicted), 'all'), 'Nonfinite predictions: inspect training.');
% Baseline: the mean TRAINING ozone profile, repeated for all test profiles.
baseline = repmat(mean(exp(Y(:,:,isTrain)), 3), 1, 1, sum(isTest));
model = ["CNN"; "Training mean profile"];
rmse = [sqrt(mean((predicted-observed).^2,'all')); ...
    sqrt(mean((baseline-observed).^2,'all'))];
mae = [mean(abs(predicted-observed),'all'); mean(abs(baseline-observed),'all')];
bias = [mean(predicted-observed,'all'); mean(baseline-observed,'all')];
metrics = table(model, rmse, mae, bias);
disp('Held-out metrics in ozone number-density units (project labels: m^-3):');
disp(metrics);

%% 8. Save everything needed to reproduce preprocessing and make predictions
if ~isfolder(options.OutputFolder), mkdir(options.OutputFolder); end
result = struct('net',net, 'scaling',scaling, 'altitude_m',altitude_m, ...
    'options',options, 'history',history, 'metrics',metrics, ...
    'manifest',manifest, 'rejections',rejections, ...
    'observed',observed, 'predicted',predicted, 'baseline',baseline, ...
    'validationObserved',validationObserved,'validationPredicted',validationPredicted, ...
    'validationLogMSE',validationLogMSE, ...
    'validationDifferenceRMSE',validationDifferenceRMSE, ...
    'validationStructureRatio',validationStructureRatio);
save(fullfile(options.OutputFolder,'ozone_cnn.mat'), '-struct', 'result');
writetable(manifest, fullfile(options.OutputFolder,'split_manifest.csv'));
writetable(rejections, fullfile(options.OutputFolder,'rejected_profiles.csv'));
writetable(metrics, fullfile(options.OutputFolder,'test_metrics.csv'));
if options.ShowPlots
    plot_test_profiles(result, OutputFolder=options.OutputFolder);
end
fprintf('Saved model and evaluation to %s\n', options.OutputFolder);
end
