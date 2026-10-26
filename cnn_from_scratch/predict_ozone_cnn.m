function [ozone, altitude_m] = predict_ozone_cnn(modelFile, on, off, inputAltitude_m)
%PREDICT_OZONE_CNN Retrieve ONE new profile using saved training preprocessing.
% on/off must have the same units, channel order, and processing as training.
% Always supply the actual input altitude grid in metres. No labels required.
% Example:
%   s = load('new_profile.mat');
%   z = (1:numel(s.mg_on_fr))'; % ONLY if this is the actual grid
%   [N,zout] = predict_ozone_cnn('outputs/ozone_cnn.mat', ...
%       s.mg_on_fr, s.mg_off_fr, z);

model = load(modelFile, 'net', 'scaling', 'altitude_m', 'options');
altitude_m = model.altitude_m;
% Older two-channel models did not save a FeatureSet setting.
featureSet = "basic";
if isfield(model.options,'FeatureSet'), featureSet = model.options.FeatureSet; end
X = prepare_lidar_input(on, off, inputAltitude_m, altitude_m, featureSet);
X = single((X-model.scaling.inputMean) ./ model.scaling.inputStd);
Y = minibatchpredict(model.net, X, ...
    InputDataFormats='TCB', OutputDataFormats='TCB');
ozone = exp(double(Y(:))*model.scaling.targetStd + model.scaling.targetMean);
assert(all(isfinite(ozone)), 'Model returned nonfinite ozone values.');
end
