function check_pipeline
%CHECK_PIPELINE Two-epoch integration check, not an accuracy experiment.
here = fileparts(mfilename('fullpath'));
for name = ["train_ozone_cnn.m", "prepare_lidar_input.m", "predict_ozone_cnn.m"]
    checkcode(fullfile(here,name));
end
r = train_ozone_cnn(MaxEpochs=2, ShowPlots=false, ...
    OutputFolder=fullfile(here,'outputs','smoke'));
trainDays = r.manifest.day(r.manifest.split == "train");
valDays = r.manifest.day(r.manifest.split == "validation");
testDays = r.manifest.day(r.manifest.split == "test");
assert(isempty(intersect(trainDays,valDays)));
assert(isempty(intersect(trainDays,testDays)));
assert(isempty(intersect(valDays,testDays)));
s = load(r.manifest.paths(find(r.manifest.split == "test",1)));
z = r.options.FirstBinAltitude_m + ...
    (0:numel(s.mg_on_fr)-1)'*r.options.BinSpacing_m;
[n,zOut] = predict_ozone_cnn(fullfile(r.options.OutputFolder,'ozone_cnn.mat'), ...
    s.mg_on_fr,s.mg_off_fr,z);
assert(isequal(zOut,r.altitude_m));
assert(all(isfinite(n)) && all(n > 0));
assert(max(abs(n-r.predicted(:,1,1))./max(n,1)) < 1e-4, ...
    'Single-profile inference disagrees with batched evaluation.');
% Invalid retrieval-interval data must fail instead of being silently repaired.
badOn = s.mg_on_fr;
[~,i] = min(abs(z-mean(r.options.AltitudeBounds_m)));
badOn(i) = NaN;
didReject = false;
try
    prepare_lidar_input(badOn,s.mg_off_fr,z,zOut);
catch
    didReject = true;
end
assert(didReject,'Invalid signal was not rejected.');
disp('Pipeline check passed: training, day separation, saved inference, invalid input.');
end
