%[text] Loads training and neural network data
%[text] Ouputs lidar signals, ozone profiles, and network
function [signals_train, profiles_train, signals_validation, profiles_validation, signals_test, profiles_test, net] = init_sim
    addpath(genpath("..\signal_analysis"));

    network = load("dir\network.mat");
%%
    n_prof = size(profiles,2);
    
    [idx_train, idx_validation, idx_test] = trainingPartitions(n_prof,[.85, .1, .05]);
    
    signals_train = signals(:,:,:,idx_train);
    signals_test = signals(:,:,:,idx_test);
    signals_validation = signals(:,:,:,idx_validation);

    profiles_train = profiles(:,idx_train);
    profiles_test = profiles(:,idx_test)';
    profiles_validation = profiles(:,idx_validation)';

    net = network.net_1;
end

%[appendix]{"version":"1.0"}
%---
