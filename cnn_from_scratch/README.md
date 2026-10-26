# Learnable CNN ozone retrieval

This is an independent MATLAB starting point: two lidar signal profiles go in,
and one ozone number-density profile comes out. The CNN starts with random
weights. It learns to reproduce `N_o3_fr_mg`, your existing retrieval product.
Agreement with those labels does **not** establish agreement with independent
atmospheric truth.

## Run it

Requires MATLAB R2024a or newer and Deep Learning Toolbox. No Python, external
project paths, datastore setup, or pretrained model is required.

Open this folder in MATLAB, then run:

```matlab
% First, check that the entire pipeline runs (not an accuracy experiment).
check = train_ozone_cnn(MaxEpochs=2, ShowPlots=false, ...
    OutputFolder=fullfile(pwd,'outputs','smoke'));

% Alternatively, also verify day separation, saved inference, and rejection
% of an invalid input (this performs its own two-epoch run).
check_pipeline;

% Then train with validation-based early stopping (up to 400 epochs).
result = train_ozone_cnn;
```

The default input folder is `../profiles`, resolved relative to the script.
The files should be named `o3_processed_YYYYMMDD_profN.mat` and each contain
one vector for each of `mg_on_fr`, `mg_off_fr`, and `N_o3_fr_mg`.
The input units and processing must be consistent across all files.

**Verify the altitude convention before interpreting results.** The default
matches `ccny_train.m`: index 1 means 1 metre, with 1 metre between indices.
The training grid spans 500–2500 m, resampled to 200 bins. For another grid:

```matlab
result = train_ozone_cnn(FirstBinAltitude_m=0, BinSpacing_m=3.75, ...
    AltitudeBounds_m=[500 2500], NumBins=200);
```

Those settings describe the **stored processed vectors**, which need not have
the same spacing as the raw instrument bins. Changing the settings cannot
correct an incorrectly identified physical coordinate.

## Read the code in this order

1. `train_ozone_cnn.m`: numbered sections show the complete experiment.
2. `prepare_lidar_input.m`: shared signal checks, interpolation, and logarithms.
3. `predict_ozone_cnn.m`: how to apply the saved experiment to a new measurement.

The main arrays have shape **altitude bins × channels × profiles** (`TCB` in
MATLAB). By default there are six input channels and one target channel. The sequence
dimension here is altitude, not chronological time.

The six features are log(on), log(off), log(off/on), slopes of that log ratio
after 5-bin and 15-bin moving averages, and altitude in kilometres. The ratio
and slopes make small differential features easier to access than learning
them solely from two correlated signal channels. At the default grid, those
smoothing windows cover roughly 50 and 151 metres. These are input features,
not a calibrated physical DIAL inversion. They contain no ozone labels.
`FeatureSet="basic"` uses only the two log signals for an input-feature ablation.
All channels are standardized using the training set. Prediction automatically
uses the feature set stored with the model; older models remain supported.

Signals must be finite and positive within the interpolation interval. We
reject invalid profiles instead of silently filling gaps or clipping them.
We interpolate the signals and density first, then take natural logarithms.
Logs compress their dynamic range. Training-set means and standard deviations
then put the inputs and targets on manageable numerical scales. Validation,
test, and new measurements reuse those exact training statistics.

The four convolution layers learn patterns across nearby heights and both
wavelength channels. ReLU introduces nonlinearity. The last convolution
produces one value per height; no pooling removes altitude bins. With the
default grid the convolutional receptive field is 61 bins, spanning about 603 m
(the smoothed input features add further context). Dilation factors 1, 2, 4,
and 8 increase context without downsampling. This is a
context window, not a claim of physical retrieval resolution. Zero padding
means predictions near the interval boundaries have less measured context.

The loss is mean squared error of standardized **log density**. Undoing scaling
and exponentiating gives positive density in the original label units (the
project identifies these as m^-3). This objective does not directly minimize
physical-density RMSE; exponentiation also does not guarantee an unbiased
estimate of the conditional arithmetic mean.

## How to judge it

Whole measurement days are assigned to training (~70%), validation (~15%), and
test (remaining days), with at least one day in each set. Random profile splits
can give overly optimistic scores when adjacent profiles look alike. The day
split is repeatable with `Seed`, but a small number of days still means large
uncertainty. This split measures performance on held-out days, not necessarily
on a future season or a different instrument.

Validation loss selects the saved network and stops training after 80 checks
without improvement, with one check per epoch. The initial learning rate is
3e-4, halved every 60 epochs; batches contain 8 profiles. Gradient clipping
limits large updates. These changes allow more time for learning, but longer
patience alone cannot fix missing input information or overfitting.
The test set is used afterward. Reported RMSE, MAE, and
signed bias pool all test bins equally in physical density units. Compare the
CNN with the saved **training mean profile** baseline: a useful learned
retrieval should improve on simply predicting a typical profile. A two-epoch
check is not expected to do so.

`outputs/structured/` contains the revised model, normalization, settings, training history,
held-out references/predictions, split manifest, rejection report, metrics,
and plots of all test profiles when plots are enabled. Reusing an output folder replaces
its previous run, so choose a different `OutputFolder` for each experiment.
The original `outputs/ozone_cnn.mat` remains available for comparison.

The result also saves validation reference/prediction arrays, standardized
log MSE, adjacent-bin difference RMSE, and a structure ratio (RMS predicted
adjacent-bin changes divided by RMS reference changes). A ratio near zero
signals flat predictions; a ratio near one alone does not establish accuracy:
noise or misplaced layers can also increase it. Use both validation accuracy
and structure diagnostics. The loss remains log-density MSE; no gradient loss
forces the network to copy potentially noisy label structure.

You can adjust the exposed training settings without editing the function:

```matlab
result = train_ozone_cnn(MaxEpochs=600, ValidationPatience=100, ...
    InitialLearnRate=1e-4, MiniBatchSize=8, ...
    OutputFolder='outputs/experiment_02');
```

On the existing 80 profiles with seed 7 and the same day split, the revision
reduced validation standardized log MSE from 0.581 to 0.355 and validation
density RMSE from 1.71e17 to 1.21e17 m^-3. Validation structure ratio increased
from 0.280 to 0.854. The best revised checkpoint was around epoch 21; training
stopped at epoch 101 after later validation scores did not improve. Longer
training alone was not the solution. After this change was fixed, test RMSE
was 1.32e17 versus the original 2.05e17 m^-3 on the same 22 profiles. These are
one-split results, not evidence of robustness across seasons or instruments.
`outputs/structured/validation_comparison.csv` records the validation comparison.

The existing `profiles` folder was already quality filtered upstream. These
scripts check numerical validity, but cannot recreate unavailable merge-fit
flags or independently verify cloud screening. Results apply to that selected
population. Use the same preprocessing and quality screening for deployment.
Interpolation to 200 bins does not add new physical information.

Change one setting at a time, using validation data to guide decisions. After
repeated tuning, reserve fresh days for a final assessment; repeatedly choosing
models from test scores makes that test set part of model development.

## Compare multiple test profiles without retraining

From this folder, display all test profiles in pages of six panels:

```matlab
plot_test_profiles;  % Prefers outputs/structured/ozone_cnn.mat; falls back to original

% Or use a result already in your workspace:
plot_test_profiles(result);

% Choose particular test profiles and save plots plus per-profile metrics:
metrics = plot_test_profiles('outputs/structured/ozone_cnn.mat', ...
    ProfileIndices=[1 3 5 8], OutputFolder='outputs/comparison');
```

Every panel shows reference and predicted density against altitude, identifies
the source profile, and reports its RMSE. All panels share axis limits.
New training runs also save and show the training-mean baseline. Older saved
models still work, showing reference and CNN only. The returned table includes
RMSE, MAE, and bias for every test profile; profile indices follow the test rows
in `split_manifest.csv`. With `ShowPlots=true`, training automatically displays
and saves every page. Set `ProfilesPerFigure` to change the page size.

## Retrieve a new profile

```matlab
s = load('path_to_new_profile.mat');
z = (1:numel(s.mg_on_fr))'; % Replace with the actual stored altitude grid!
[N, zOut] = predict_ozone_cnn('outputs/structured/ozone_cnn.mat', ...
    s.mg_on_fr, s.mg_off_fr, z);
plot(N, zOut);
xlabel('Ozone number density (m^{-3})'); ylabel('Altitude (m)');
```

No ozone label is supplied at prediction time. Keep wavelength ordering,
signal units, averaging, background subtraction, and channel merging consistent
with the training data. This first model has no explicit DIAL physics loss,
uncertainty estimate, or out-of-distribution detector.
