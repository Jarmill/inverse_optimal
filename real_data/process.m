clear all;
close all;
clc;

% Load data
num_trial = 1;
data_filename = sprintf('ForceData_Optim_Format_Single_Trajectory_No%03d.mat', num_trial);
load(data_filename);

% Plot forces
f = cell2mat(y.');
Time = linspace(0, 100, size(f, 2));
hfig = figure;
[h_ax, h] = plot_vector_quantities(Time, f, []);
