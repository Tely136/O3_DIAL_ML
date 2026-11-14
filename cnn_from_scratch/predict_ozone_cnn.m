function [ozone,altitude_m] = predict_ozone_cnn(model_file,on,off,input_altitude_m)
% Retrieve one profile using the saved model and training normalization.
model = load(model_file,"net","scaling","altitude_m","options");
altitude_m = model.altitude_m;
feature_set = "basic";
if isfield(model.options,"FeatureSet")
    feature_set = model.options.FeatureSet;
end

X = prepare_lidar_input(on,off,input_altitude_m,altitude_m,feature_set);
X = single((X - model.scaling.inputMean) ./ model.scaling.inputStd);
Y = minibatchpredict(model.net,X,InputDataFormats="TCB",OutputDataFormats="TCB");
ozone = exp(double(Y(:))*model.scaling.targetStd + model.scaling.targetMean);
assert(all(isfinite(ozone)),"Model returned nonfinite ozone values.");
end
