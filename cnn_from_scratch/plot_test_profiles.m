function [profile_metrics,figures] = plot_test_profiles(result_or_file,options)
% Plot held-out profiles from a result struct or a saved model file.
arguments
    result_or_file = ""
    options.ProfileIndices (1,:) double = []
    options.ProfilesPerFigure (1,1) double = 6
    options.OutputFolder (1,1) string = ""
    options.Visible (1,1) logical = true
end

%% Load results
if ~isstruct(result_or_file) && string(result_or_file) == ""
    output_folder = fullfile(fileparts(mfilename("fullpath")),"outputs");
    candidates = [fullfile(output_folder,"correct_altitude","ozone_cnn.mat"),fullfile(output_folder,"structured","ozone_cnn.mat"),fullfile(output_folder,"ozone_cnn.mat")];
    result_or_file = candidates(find(isfile(candidates),1));
end
if isstruct(result_or_file)
    result = result_or_file;
else
    result = load(result_or_file,"observed","predicted","altitude_m","manifest");
    saved_variables = whos("-file",result_or_file);
    if any(strcmp({saved_variables.name},"baseline"))
        saved = load(result_or_file,"baseline");
        result.baseline = saved.baseline;
    end
end

test_rows = result.manifest(result.manifest.split == "test",:);
n_test = size(result.observed,3);
indices = options.ProfileIndices;
if isempty(indices)
    indices = 1:n_test;
end

%% Metrics
residual = result.predicted - result.observed;
test_index = (1:n_test)';
rmse = reshape(sqrt(mean(residual.^2,1)),[],1);
mae = reshape(mean(abs(residual),1),[],1);
bias = reshape(mean(residual,1),[],1);
profile_metrics = table(test_index,test_rows.paths,test_rows.day,rmse,mae,bias,VariableNames=["TestIndex","SourceFile","Day","RMSE","MAE","Bias"]);
disp(profile_metrics(indices,:));

%% Shared axes
has_baseline = isfield(result,"baseline");
labels = ["Reference retrieval","CNN"];
values = [result.observed(:,:,indices),result.predicted(:,:,indices)];
if has_baseline
    labels = ["Reference retrieval","CNN","Training mean"];
    values = [values,result.baseline(:,:,indices)];
end
x_limits = [min(values,[],"all"),max(values,[],"all")];
margin = max(0.05*diff(x_limits),eps(max(abs(x_limits))));
x_limits = x_limits + [-margin margin];
visibility = "off";
if options.Visible
    visibility = "on";
end
if options.OutputFolder ~= "" && ~isfolder(options.OutputFolder)
    mkdir(options.OutputFolder);
end

%% Profile plots
n_pages = ceil(numel(indices)/options.ProfilesPerFigure);
figures = gobjects(n_pages,1);
for page = 1:n_pages
    first = (page-1)*options.ProfilesPerFigure + 1;
    last = min(page*options.ProfilesPerFigure,numel(indices));
    selected = indices(first:last);
    figures(page) = figure(Name=sprintf("Test profiles: page %d of %d",page,n_pages),Visible=visibility,Color="w",Position=[100 100 1200 750]);
    n_columns = min(3,numel(selected));
    layout = tiledlayout(figures(page),ceil(numel(selected)/n_columns),n_columns,TileSpacing="compact",Padding="compact");

    for i = selected
        ax = nexttile(layout);
        plot(ax,result.observed(:,1,i),result.altitude_m,"k-",result.predicted(:,1,i),result.altitude_m,"r--",LineWidth=1.2);
        if has_baseline
            hold(ax,"on");
            plot(ax,result.baseline(:,1,i),result.altitude_m,"b:",LineWidth=1);
        end
        xlim(ax,x_limits);
        ylim(ax,[min(result.altitude_m),max(result.altitude_m)]);
        grid(ax,"on");
        set(ax,Color="w",XColor="k",YColor="k");
        [~,source] = fileparts(test_rows.paths(i));
        title(ax,sprintf("Test %d | %s\nRMSE = %.3g",i,source,rmse(i)),Interpreter="none",Color="k");
        xlabel(ax,"Ozone number density (m^{-3})",Color="k");
        ylabel(ax,"Altitude (m)",Color="k");
        if i == selected(1)
            legend(ax,labels,Location="best",Color="w",TextColor="k");
        end
    end

    title(layout,sprintf("Held-out ozone profiles | page %d of %d",page,n_pages),Color="k");
    if options.OutputFolder ~= ""
        exportgraphics(figures(page),fullfile(options.OutputFolder,sprintf("test_profiles_%02d.png",page)),Resolution=150,BackgroundColor="white");
    end
end
if options.OutputFolder ~= ""
    writetable(profile_metrics,fullfile(options.OutputFolder,"test_profile_metrics.csv"));
end
end
