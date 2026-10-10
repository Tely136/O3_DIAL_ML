%[text] Script for testing reading licel signal files, reading simulated lidar signals, and comparing them
clc; clearvars;
%%
addpath(genpath("C:\Users\tely1\OneDrive\Matlab Onedrive\O3_Lidar\proc_functions")); %[control:filebrowser:4f47]{"position":[17,82]}

folder = "C:\Users\tely1\OneDrive\Desktop\EEJ9903\ai\signal_data\groundbased_lidar.o3_ccny001_raw_new.york.ny_20260608t000000z_20260608t235959z_001"; %[control:filebrowser:1683]{"position":[10,149]}

[scaled_data,times,configs] = read_lidar_datafiles(folder,"bin","ol");
%%
sim_prof = load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\lidar_model\sim_prof.mat"); %[control:filebrowser:536b]{"position":[17,83]}
sim_alt = sim_prof.alt_ext./1000; %[output:3a4b1a28]

N_sim = [sim_prof.N_287(:,1) sim_prof.N_299(:,1)];
N_sim(sim_alt < .75,:) = NaN;

figure;
semilogy(sim_alt, N_sim(:,1),'DisplayName','287')
hold on
plot(sim_alt, N_sim(:,2),'DisplayName','299')
hold off
legend
grid on
%%
dt = 1 / (40e6);

alt = (1:size(scaled_data,1)) .* 3.75 ./1000;

data = squeeze(scaled_data(:,1,[1,2],2));
data = data./(1-data/270);
b = mean(data(5000:7999,:),1);
r = data - b;
N = r .* dt;

data_an = squeeze(scaled_data(:,1,[1,2],1));
b_an = mean(data_an(5000:7999,:),1);

v = data_an - b_an;

nan_id = v<0 | r<0;
v(nan_id) = NaN;
r(nan_id) = NaN;

[r_glued_287, ~, ~, ~] = glue_an_pc(v(:,1),r(:,1),3,60); 
[r_glued_299, ~, ~, ~] = glue_an_pc(v(:,2),r(:,2),3,60); 

r = r .* 1e6;
r_glued = [r_glued_287 r_glued_299] .* 1e6;

N_glued = r_glued .* dt;
N_glued_interp = [interp1(alt', N_glued(:,1), sim_alt') interp1(alt', N_glued(:,2), sim_alt')];
%%
figure;
semilogy(alt,r(:,1),'DisplayName','287')
hold on
plot(alt,r(:,2),'DisplayName','299')
legend
grid on
title("Photon counting")
figure;
semilogy(alt,v(:,1),'DisplayName','287')
hold on
plot(alt,v(:,2),'DisplayName','299')
legend
grid on
title("Analog")
figure;
semilogy(alt,r_glued(:,1),'DisplayName','287')
hold on
plot(alt,r_glued(:,2),'DisplayName','299')
legend
grid on
title("Glued Photon Counting")
figure;
semilogy(alt,N_glued(:,1),'DisplayName','287')
hold on
plot(alt,N_glued(:,2),'DisplayName','299')
legend
grid on
xlim([0 3])
title("Measured Photon Counts per Range Bin - Glued analog and photon counting")
figure;
semilogy(sim_alt,N_glued_interp(:,1),'DisplayName','287')
hold on
plot(sim_alt,N_glued_interp(:,2),'DisplayName','299')
legend
grid on
xlim([0 3])
title("Measured Photon Counts per Range Bin - Glued analog and photon counting - Interpolated")
%%
factor = N_sim ./ N_glued_interp;

figure;
plot(sim_alt,factor(:,1),'DisplayName','287')
hold on
plot(sim_alt, factor(:,2),'DisplayName','299')
hold off
legend
grid on
ids = sim_alt>=2 & sim_alt<=4;

f_287 = mean(factor(ids,1));
f_299 = mean(factor(ids,2));

figure;
semilogy(alt,N_glued(:,1),'DisplayName','287 measured')
hold on
plot(sim_alt, N_sim(:,1),'DisplayName','287 simulated')
plot(alt,N_glued(:,2),'DisplayName','299 measured')
plot(sim_alt, N_sim(:,2),'DisplayName','299 simulated')

legend
grid on
xlim([0 3])
title("Measured Photon Counts per Range Bin - Glued analog and photon counting")

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[control:filebrowser:4f47]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:1683]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:536b]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[output:3a4b1a28]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized field name \"alt_ext\"."}}
%---
