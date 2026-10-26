function [profileMetrics, figures] = plot_test_profiles(resultOrFile, options)
%PLOT_TEST_PROFILES Compare reference and CNN retrievals on held-out profiles.
% No retraining: use either a result struct or a saved ozone_cnn.mat file.
%   plot_test_profiles(result)
%   plot_test_profiles('outputs/ozone_cnn.mat')
%   plot_test_profiles(result, ProfileIndices=[1 4 7], OutputFolder='plots')
% Test indices refer to the rows of manifest with split == "test".

arguments
    resultOrFile = ""
    options.ProfileIndices (1,:) double {mustBeInteger,mustBePositive} = []
    options.ProfilesPerFigure (1,1) double {mustBeInteger,mustBePositive} = 6
    options.OutputFolder (1,1) string = ""
    options.Visible (1,1) logical = true
end

if ~isstruct(resultOrFile) && string(resultOrFile) == ""
    outputRoot = fullfile(fileparts(mfilename('fullpath')),'outputs');
    resultOrFile = fullfile(outputRoot,'structured','ozone_cnn.mat');
    if ~isfile(resultOrFile), resultOrFile = fullfile(outputRoot,'ozone_cnn.mat'); end
end
if isstruct(resultOrFile)
    r = resultOrFile;
else
    % Only load the evaluation arrays; the trained network is not needed.
    r = load(resultOrFile, 'observed','predicted','altitude_m','manifest');
    savedVariables = whos('-file',resultOrFile);
    if any(strcmp({savedVariables.name},'baseline'))
        savedBaseline = load(resultOrFile,'baseline');
        r.baseline = savedBaseline.baseline;
    end
end
assert(isequal(size(r.observed),size(r.predicted)), ...
    'Reference and prediction arrays must have the same shape.');
nTest = size(r.observed,3);
assert(nTest > 0, 'No held-out profiles available.');
testRows = r.manifest(r.manifest.split == "test",:);
assert(height(testRows) == nTest, 'Manifest does not match the test arrays.');
indices = options.ProfileIndices;
if isempty(indices), indices = 1:nTest; end
assert(all(indices <= nTest), 'ProfileIndices exceeds the number of test profiles.');

% One metric row per test profile, including the source file for traceability.
error = r.predicted-r.observed;
testIndex = (1:nTest)';
rmse = reshape(sqrt(mean(error.^2,1)),[],1);
mae = reshape(mean(abs(error),1),[],1);
bias = reshape(mean(error,1),[],1);
profileMetrics = table(testIndex,testRows.paths,testRows.day,rmse,mae,bias, ...
    VariableNames={'TestIndex','SourceFile','Day','RMSE','MAE','Bias'});
disp(profileMetrics(indices,:));

% Use the same axes across panels so the visual comparison is fair.
values = [r.observed(:,:,indices),r.predicted(:,:,indices)];
hasBaseline = isfield(r,'baseline');
if hasBaseline
    values = [values,r.baseline(:,:,indices)];
end
xLimits = [min(values,[],'all'),max(values,[],'all')];
margin = max(0.05*diff(xLimits),eps(max(abs(xLimits))));
xLimits = xLimits+[-margin margin];
visibility = 'off';
if options.Visible, visibility = 'on'; end
if options.OutputFolder ~= "" && ~isfolder(options.OutputFolder)
    mkdir(options.OutputFolder);
end
nPages = ceil(numel(indices)/options.ProfilesPerFigure);
figures = gobjects(nPages,1);
for page = 1:nPages
    selected = indices((page-1)*options.ProfilesPerFigure+1: ...
        min(page*options.ProfilesPerFigure,numel(indices)));
    figures(page) = figure('Name',sprintf('Test profiles: page %d of %d',page,nPages), ...
        'Visible',visibility,'Color','w','Position',[100 100 1200 750]);
    nColumns = min(3,numel(selected));
    layout = tiledlayout(figures(page),ceil(numel(selected)/nColumns),nColumns, ...
        'TileSpacing','compact','Padding','compact');
    for i = selected
        ax = nexttile(layout);
        plot(ax,r.observed(:,1,i),r.altitude_m,'k-', ...
            r.predicted(:,1,i),r.altitude_m,'r--','LineWidth',1.2);
        labels = {'Reference retrieval','CNN'};
        if hasBaseline
            hold(ax,'on');
            plot(ax,r.baseline(:,1,i),r.altitude_m,'b:','LineWidth',1);
            labels{end+1} = 'Training mean';
        end
        xlim(ax,xLimits); ylim(ax,[min(r.altitude_m),max(r.altitude_m)]);
        grid(ax,'on');
        set(ax,'Color','w','XColor','k','YColor','k');
        [~,source] = fileparts(testRows.paths(i));
        title(ax,sprintf('Test %d | %s\nRMSE = %.3g',i,source,rmse(i)), ...
            'Interpreter','none','Color','k');
        xlabel(ax,'Ozone number density (m^{-3})','Color','k');
        ylabel(ax,'Altitude (m)','Color','k');
        if i == selected(1)
            legend(ax,labels,'Location','best','Color','w','TextColor','k');
        end
    end
    title(layout,sprintf('Held-out ozone profiles | page %d of %d',page,nPages),'Color','k');
    if options.OutputFolder ~= ""
        exportgraphics(figures(page),fullfile(options.OutputFolder, ...
            sprintf('test_profiles_%02d.png',page)),'Resolution',150,'BackgroundColor','white');
    end
end
if options.OutputFolder ~= ""
    writetable(profileMetrics,fullfile(options.OutputFolder,'test_profile_metrics.csv'));
end
end
