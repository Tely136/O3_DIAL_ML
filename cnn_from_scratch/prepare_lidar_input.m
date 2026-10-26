function X = prepare_lidar_input(on, off, altitude_m, outputAltitude_m, featureSet)
%PREPARE_LIDAR_INPUT Shared preprocessing for training and new measurements.
% Inputs are vectors for ONE profile, on the same altitude grid.
% Output: [altitude bins x channels]. Legacy default: log(on), log(off).
% 'dial' adds log(off/on), two smoothed slopes, and altitude in kilometres.
% No per-profile scaling: it could remove information useful for retrieval.

if nargin < 5, featureSet = "basic"; end
assert(ismember(string(featureSet),["basic","dial"]),'Unknown feature set.');
on = double(on(:));
off = double(off(:));
altitude_m = double(altitude_m(:));
outputAltitude_m = double(outputAltitude_m(:));
assert(numel(on) == numel(off) && numel(on) == numel(altitude_m), ...
    'Signals and altitude must have the same length.');
assert(all(isfinite(altitude_m)) && all(diff(altitude_m) > 0), ...
    'Altitude must be finite and strictly increasing.');
assert(all(isfinite(outputAltitude_m)) && all(diff(outputAltitude_m) > 0), ...
    'Output altitude must be finite and strictly increasing.');
assert(outputAltitude_m(1) >= altitude_m(1) && ...
    outputAltitude_m(end) <= altitude_m(end), ...
    'Requested altitude range lies outside this measurement.');

% Include the interpolation neighbors just outside the requested interval.
first = find(altitude_m <= outputAltitude_m(1), 1, 'last');
last = find(altitude_m >= outputAltitude_m(end), 1, 'first');
z = altitude_m(first:last);
signals = [on(first:last), off(first:last)];
assert(all(isfinite(signals), 'all') && all(signals > 0, 'all'), ...
    'Signals must be finite and positive throughout the retrieval interval.');
X = log(interp1(z, signals, outputAltitude_m, 'linear'));
if string(featureSet) == "dial"
    % Subtract logs rather than dividing tiny signals. All features use only
    % measured signals/altitude, never ozone labels. Two slope scales expose
    % local differential structure without asking the CNN to discover it
    % by subtracting two much larger, correlated signal channels.
    logRatio = X(:,2)-X(:,1);
    slopeFine = gradient(movmean(logRatio,5),outputAltitude_m);
    slopeBroad = gradient(movmean(logRatio,15),outputAltitude_m);
    X = [X, logRatio, slopeFine, slopeBroad, outputAltitude_m/1000];
end
end
