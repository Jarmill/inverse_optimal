clear all;
close all;
clc;

% Which file to load
leg = 1;
speed = 8;
trial = 5;
filename = sprintf('force_data_l%d_s%d_t%02d.mat', leg, speed-4, trial);

if ~exist(filename, 'file')
    error('Inexistant file');
end

% Load
load(filename);

% Determine the number of variables
n = size(raw_data.f, 1);
% Determine the number of samples
N = size(raw_data.f, 3);
% Determine the number of cost functions
nf = size(inv_opt.Q, 3);


% The raw data
raw_data;

% The inverse otpimization structure
inv_opt
