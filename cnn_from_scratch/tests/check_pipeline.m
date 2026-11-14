function check_pipeline
% Two-epoch check of training, day separation and saved-model prediction.
script_folder = string(fileparts(fileparts(mfilename("fullpath"))));
addpath(script_folder);
files = ["train_ozone_cnn.m","prepare_lidar_input.m","predict_ozone_cnn.m","plot_test_profiles.m"];
for file = files
    checkcode(fullfile(script_folder,file));
end

result = train_ozone_cnn(MaxEpochs=2,ShowPlots=false,OutputFolder=fullfile(script_folder,"outputs","smoke"));
train_days = result.manifest.day(result.manifest.split == "train");
val_days = result.manifest.day(result.manifest.split == "validation");
test_days = result.manifest.day(result.manifest.split == "test");
assert(isempty(intersect(train_days,val_days)));
assert(isempty(intersect(train_days,test_days)));
assert(isempty(intersect(val_days,test_days)));

data = load(result.manifest.paths(find(result.manifest.split == "test",1)));
z = result.options.FirstBinAltitude_m + (0:numel(data.mg_on_fr)-1)'*result.options.BinSpacing_m;
[ozone,altitude_m] = predict_ozone_cnn(fullfile(result.options.OutputFolder,"ozone_cnn.mat"),data.mg_on_fr,data.mg_off_fr,z);
assert(isequal(altitude_m,result.altitude_m));
assert(all(isfinite(ozone)) && all(ozone > 0));
assert(max(abs(ozone-result.predicted(:,1,1))./max(ozone,1)) < 1e-4,"Saved inference disagrees with evaluation.");

bad_on = data.mg_on_fr;
[~,i] = min(abs(z-mean(result.options.AltitudeBounds_m)));
bad_on(i) = NaN;
rejected = false;
try
    prepare_lidar_input(bad_on,data.mg_off_fr,z,altitude_m);
catch
    rejected = true;
end
assert(rejected,"Invalid signal was not rejected.");
disp("Pipeline check passed.");
end
