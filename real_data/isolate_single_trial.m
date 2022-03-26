clear all;
close all;
clc;

% Load data
data_filename = 'ForceData_Optim_Format.mat';
load(data_filename);

% Savename
savename = 'ForceData_Optim_Format_Single_Trajectory_No%03d.mat';

% Overwrite previous
overwrite = true;
overwrite_num = 1;

% Pick one trajectory and do IO on it
leg = 1;
speed = 2;
trial = 1;

Q = squeeze(Q_list(leg, trial, speed, :));
phi = squeeze(phi_list(leg, trial, speed, :));
A = squeeze(A_list(leg, trial, speed, :));
b = squeeze(b_list(leg, trial, speed, :));
Aeq = squeeze(Aeq_list(leg, trial, speed, :));
beq = squeeze(beq_list(leg, trial, speed, :));
y = squeeze(y_list(leg, trial, speed, :));

varnames = {'Q', 'phi', 'A', 'b', 'Aeq', 'beq', 'y'};

% If no overwriting save to the next possible number
if ~overwrite
    % Initialize save number
    ii = 1;

    % While there are other saves
    while exist(sprintf(savename, ii), 'file')
        ii = ii + 1;
    end
    
    % Save the variables into the first available savenumber
    save(sprintf(savename, ii), varnames{:});
    
% If overwrite
else
    
    % Save the variables into the overwritable number
    save(sprintf(savename, overwrite_num), varnames{:});
end