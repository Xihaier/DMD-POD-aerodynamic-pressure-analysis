%%

clear;
close all;
clc;

addpath('./data');
addpath('./utils');
addpath('./viz');

%-------------------step 1: load data-------------------
load Windward.mat;

% define parameters
num_snap = 1030;
num_d = 30;

for k1 = 1:numel(num_snap)
    for k2 = 1:numel(num_d)
        
        num.snapshots = num_snap(k1);
        num.truncate = 125;
        num.delay = num_d(k2);
        
        dt = 1/1000;
        
        % assemble a Hankel matrix
        [DMD_data] = Hankel_matrix(data, num);
        
        %-------------------step 2: DMD-------------------
        X1 = DMD_data(:, 1:(end-1));
        X2 = DMD_data(:, 2:end);
        
        % perform decomposition
        [DMD_infor] = DMD_method(X1, X2, num, dt);
        
        %-------------------step 3: results-------------------
        results_summary(DMD_infor, X1, dt);
        
    end
end


