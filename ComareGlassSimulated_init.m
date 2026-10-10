%[text] # Load Experiment Data
%[text] This initialization function template loads the training and validation data for the experiment as 4-D arrays. The training and validation data sets each contain 5000 images of handwritten digits and their digit labels. Each image shows a digit from 0 to 9, rotated by a certain angle. The regression values correspond to the angles of rotation of the digits. In the setup function, use `params.InitializationFunctionOutput` to access the initialization function output.
%[text] ## Input
%[text] - This function does not accept any input arguments. \
%[text] ## Output
%[text] - `output` $-$ Experiment data, returned as a structure. The fields  `XTrain` and `TTrain` contains the training data, and the fields  `XValidation` and `TValidation` contains the validation data. \
%[text] ## Initialization Function
%[text] Initialize the required data for the experiment.
function output = ComareGlassSimulated_init
    
    addpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML"); %[control:filebrowser:8145]{"position":[13,65]}
    % addpath("C:\Users\tely1\OneDrive\Matlab Onedrive\O3_Lidar\proc_functions"); %[control:filebrowser:98f4]{"position":[15,80]}
    addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\functions")); %[control:filebrowser:1eae]{"position":[21,88]}
    
    % [output.signals_train, output.signals_sim_train, output.profiles_train, output.signals_validation, output.signals_sim_validation,output.profiles_validation,output.signals_test,output.signals_sim_test,output.profiles_test,output.net] = init_sim;
    load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\dir\GLASS_training_data.mat"); %[control:filebrowser:6877]{"position":[10,90]}
    network = load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\dir\network.mat"); %[control:filebrowser:879b]{"position":[20,88]}
%%
    % Ozone Absorption Cross Section
    BDM = load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\cross_sections\ozone\bdm\bdm_287_299.mat"); %[control:filebrowser:26ef]{"position":[16,114]}
    
    o3_287 = BDM.BDM_c287(1) .* 1e-4; % use first temp dependent cross section for now
    o3_299 = BDM.BDM_c299(1) .* 1e-4;
    
    % do3 = BDM.BDM_d_287_299(1) .* 1e-4;
    
    % Alt grid for calculating simulated profile
    alt_inst = 100;
    alt_sim = (alt_inst:1:alt(end))';
    
    % Air Number Density
    Nair = airNumberDensity(alt_sim);
    
    % Prepare parameters for signal calculation
    E0_287 = 10.6e-3;
    E0_299 = 6.3e-3;
    
    eta_287 = 0.005; % Efficiency factor
    eta_299 = 0.009; % 
    
    etaQ = 0.7;
    
    t_r = .2;   % Radius of the receiver in meters
    o_r = 0.1;  % Radius of the obstruction
    A = pi.*(t_r^2 - o_r^2);
    
    f = 40e6;
    tau = 1/f; % Laser pulse duration
    
    load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\O_test.mat"); %[control:filebrowser:6514]{"position":[10,78]}
    
    O_interp = interp1(O_grid,O,alt_sim,'linear',1);
    profiles_interp = interp1(alt,profiles,alt_sim,'linear',0);
    
    % Calculate the signal using the defined parameters
    N_287 = signal(alt_sim,287.2e-9,E0_287,A,tau,eta_287,etaQ,Nair,o3_287,profiles_interp,O_interp) .*1e-6./tau; 
    N_299 = signal(alt_sim,299.1e-9,E0_299,A,tau,eta_299,etaQ,Nair,o3_299,profiles_interp,O_interp) .*1e-6./tau; 
    
    N_287 = interp1(alt_sim,N_287,alt);
    N_299 = interp1(alt_sim,N_299,alt);
    
    signals_sim = NaN(length(alt), 1, 2, size(N_287, 2));
    
    for i = 1:size(N_287,2)
        signals_sim(:,1,1,i) = N_287(:,i);
        signals_sim(:,1,2,i) = N_299(:,i);
    end
%%
    n_prof = size(profiles,2);
    
    [idx_train,idx_validation, idx_test] = trainingPartitions(n_prof,[.85, .1, .05]);
    
    output.signals_train = signals(:,:,:,idx_train);
    output.signals_sim_train = signals_sim(:,:,:,idx_train);
    
    output.profiles_train = profiles(:,idx_train);
    
    output.signals_validation = signals(:,:,:,idx_validation);
    output.signals_sim_validation = signals_sim(:,:,:,idx_validation);
    
    output.profiles_validation = profiles(:,idx_validation)';
    
    output.signals_test = signals(:,:,:,idx_test);
    output.signals_sim_test = signals_sim(:,:,:,idx_test);
    
    output.profiles_test = profiles(:,idx_test)';
    
    output.net = network.net_1;
end


%[appendix]{"version":"1.0"}
%---
%[control:filebrowser:8145]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:98f4]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:1eae]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:6877]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"AllSections"}
%---
%[control:filebrowser:879b]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"AllSections"}
%---
%[control:filebrowser:26ef]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:6514]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"Section"}
%---
