%[text] Loads training and neural network data
%[text] Ouputs lidar signals, ozone profiles, and network
function [signals_train, signals_sim_train, profiles_train, signals_validation, signals_sim_validation, profiles_validation, signals_test, signals_sim_test, profiles_test, net] = init_sim
    addpath(genpath("..\signal_analysis"));

    load("dir\GLASS_training_data.mat");
    network = load("dir\network.mat");
%%
    % Ozone Absorption Cross Section
    BDM = load("..\signal_analysis\cross_sections\ozone\bdm\bdm_287_299.mat");
    
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

    % Dead-time constant
    td = 270e-5;
    
    eta_287 = 0.005; % Efficiency factor
    eta_299 = 0.009; % 
    
    etaQ = 0.7;
    
    t_r = .2;   % Radius of the receiver in meters
    o_r = 0.1;  % Radius of the obstruction
    A = pi.*(t_r^2 - o_r^2);
    
    f = 40e6;
    tau = 1/f; % Laser pulse duration

    overlap_data = load("..\signal_analysis\overlap_data.mat");

    O_interp = interp1(overlap_data.overlap_params_save.alt,overlap_data.O,alt_sim,'linear',1);
    profiles_interp = interp1(alt,profiles,alt_sim,'linear',0);
    
    % Calculate the signal using the defined parameters
    N_287 = signal(alt_sim,287.2e-9,E0_287,A,tau,eta_287,etaQ,Nair,o3_287,profiles_interp,O_interp) .*1e-6./tau; 
    N_299 = signal(alt_sim,299.1e-9,E0_299,A,tau,eta_299,etaQ,Nair,o3_299,profiles_interp,O_interp) .*1e-6./tau; 

    N_287 = add_dt(N_287,td);
    N_299 = add_dt(N_299,td);

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
    
    signals_train = signals(:,:,:,idx_train);
    signals_sim_train = signals_sim(:,:,:,idx_train);
    
    profiles_train = profiles(:,idx_train);
    
    signals_validation = signals(:,:,:,idx_validation);
    signals_sim_validation = signals_sim(:,:,:,idx_validation);
    
    profiles_validation = profiles(:,idx_validation)';
    
    signals_test = signals(:,:,:,idx_test);
    signals_sim_test = signals_sim(:,:,:,idx_test);
    
    profiles_test = profiles(:,idx_test)';

    net = network.net_1;
end

function Nmeas = add_dt(Nreal,dt)
    Nmeas = Nreal ./ (1 + Nreal * dt);
end

%[appendix]{"version":"1.0"}
%---
