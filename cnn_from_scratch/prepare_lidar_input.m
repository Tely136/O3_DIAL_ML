function X = prepare_lidar_input(on,off,altitude_m,output_altitude_m,feature_set)
% Prepare one profile using the same features for training and prediction.
if nargin < 5
    feature_set = "basic";
end
assert(ismember(string(feature_set),["basic","dial"]),"Unknown feature set.");
on = double(on(:));
off = double(off(:));
altitude_m = double(altitude_m(:));
output_altitude_m = double(output_altitude_m(:));
assert(numel(on) == numel(off) && numel(on) == numel(altitude_m),"Signals and altitude must have the same length.");
assert(all(isfinite(altitude_m)) && all(diff(altitude_m) > 0),"Altitude must be finite and increasing.");
assert(all(isfinite(output_altitude_m)) && all(diff(output_altitude_m) > 0),"Output altitude must be finite and increasing.");
assert(output_altitude_m(1) >= altitude_m(1) && output_altitude_m(end) <= altitude_m(end),"Requested altitude range lies outside this measurement.");

% Include the bins immediately outside the requested interval for interpolation.
first = find(altitude_m <= output_altitude_m(1),1,"last");
last = find(altitude_m >= output_altitude_m(end),1,"first");
z = altitude_m(first:last);
signals = [on(first:last),off(first:last)];
assert(all(isfinite(signals),"all") && all(signals > 0,"all"),"Signals must be finite and positive throughout the retrieval interval.");
X = log(interp1(z,signals,output_altitude_m,"linear"));

if string(feature_set) == "dial"
    log_ratio = X(:,2) - X(:,1);
    slope_fine = gradient(movmean(log_ratio,5),output_altitude_m);
    slope_broad = gradient(movmean(log_ratio,15),output_altitude_m);
    X = [X,log_ratio,slope_fine,slope_broad,output_altitude_m/1000];
end
end
