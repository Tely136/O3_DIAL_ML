%[text] reads GLASS file table produced by "produce\_glass\_database"
%[text] matches ozone profiles in GLASS files to lidar signals in signal\_folder
%[text] produces signal\_database stored in specified location which contains the training data for emulating GLASS
clc; clearvars;
addpath(genpath("C:\Users\tely1\OneDrive\Matlab Onedrive\O3_Lidar\proc_functions")); %[control:filebrowser:8db7]{"position":[17,82]}
addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\functions")); %[control:filebrowser:11da]{"position":[17,84]}
addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\functions")); %[control:filebrowser:5527]{"position":[17,79]}

save_name = "GLASS_training_data.mat"; %[control:editfield:1f32]{"position":[13,38]}
save_path = "C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\dir"; %[control:filebrowser:58c5]{"position":[13,69]}
%%
signal_folder = dir("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\signal_data"); %[control:filebrowser:4d22]{"position":[21,85]}
glass_folder = dir("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\GLASS"); %[control:filebrowser:155f]{"position":[20,78]}

signal_folder = signal_folder(~ismember({signal_folder.name}, {'.','..'}));
%%
signal_database = table;

n_sig_f = length(signal_folder);
for i = 1:n_sig_f
    temp_folder = fullfile(signal_folder(i).folder, signal_folder(i).name);

    temp_files = dir(temp_folder);
    temp_files = temp_files(not([temp_files.isdir]) & not(endsWith({temp_files.name},{'.dat','.zip'})));
    n_good = size(temp_files,1);

    names = string({temp_files.name});
    times = parseLicelFilename(names);

    paths = fullfile(string({temp_files.folder}), string({temp_files.name}));

    signal_database = [signal_database; table(times', paths')];
end

signal_database.Properties.VariableNames = ["Datetime", "Path"];
%%
load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\dir\glass_profiles_cldrem.mat"); %[control:filebrowser:13b6]{"position":[6,88]}

n_ds = size(glass_dataset,1);
total_n_prof = 0;
for i = 1:n_ds
    total_n_prof = total_n_prof + length(cell2mat(glass_dataset.id_good(i)));
end 

alt = (1e3:30:5e3)';

signals = NaN(length(alt),1,2,total_n_prof);
profiles = NaN(length(alt),total_n_prof);

f = waitbar(0,"Processing");

count = 1;
for i = 1:n_ds
    temp_ds = glass_dataset(i,:);
    n_prof = length(cell2mat(temp_ds.id_good));

    temp_glass_data = read_glass(temp_ds.file);
    glass_alt = temp_glass_data.alt;
    glass_N = temp_glass_data.N';
    glass_N(glass_N<0) = NaN;

    id_good = cell2mat(temp_ds.id_good);

    for j = 1:n_prof
        temp_idx = id_good(j);
        profiles(:,count) = interp1(glass_alt, glass_N(:,temp_idx), alt);

        temp_start = temp_glass_data.dt_start(temp_idx);
        temp_stop = temp_glass_data.dt_stop(temp_idx);

        % Filter the signal database based on the current dataset time range
        timeMask = (signal_database.Datetime >= temp_start) & (signal_database.Datetime <= temp_stop);
        matched_signals = signal_database(timeMask, :);

        if isempty(matched_signals)
            continue
        end

        [scaled_data,~,~] = read_lidar_datafiles(matched_signals.Path,"list","");
        on = mean(scaled_data(:,:,1,2), 2);
        off = mean(scaled_data(:,:,2,2), 2);

        alt_sig = (1:length(on)) .* 3.75;

        on_interp = interp1(alt_sig, on, alt);
        off_interp = interp1(alt_sig, off, alt);

        signals(:,:,1,count) = on_interp';
        signals(:,:,2,count) = off_interp';

        count = count+1;
    end

    waitbar(i./n_ds,f);
end
close(f)
%%
f = any(isnan(profiles),1);
profiles(:,f) = [];
signals(:,:,:,f) = [];

save(fullfile(save_path,save_name), "alt","signals","profiles");
%%
function dts = parseLicelFilename(filename)
    dts = NaT(1,length(filename));

    for i = 1:length(filename)

        % Remove extension if present
        [~, name, j] = fileparts(filename(i));
        
        % Extract portion after "ol"
        % s = char(erase(name,{'ol','te','2d','al','he','xy','xx','oz'}));
        s = char(name);
        s = s(3:end);

        j = char(erase(j,'.'));
        
        year  = 2000 + str2double(s(1:2));
        month = str2double(s(3));
        day   = str2double(s(4:5));
        hour  = str2double(s(6:7));
        minute = str2double(j(1:2));
        
        dts(i) = datetime(year,month,day,hour,minute,0);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[control:filebrowser:8db7]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:11da]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:5527]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:editfield:1f32]
%   data: {"defaultValue":"\"\"","label":"save_name","run":"Nothing","valueType":"String"}
%---
%[control:filebrowser:58c5]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Nothing"}
%---
%[control:filebrowser:4d22]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:155f]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:13b6]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"Nothing"}
%---
