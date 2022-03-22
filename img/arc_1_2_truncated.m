%EXAMPLE_PLOTTING_LOADED_DATA shows how to use GENERATE_RANDOM_Q, 
%GENERATE_RANDOM_X, PLOT_MATHCAL_G, PLOT_X_STAR, PLOT_LEVEL_SETS and 
%PLOT_ARC, on loaded data.
clear all;
close all;
clc;

% Add plotting functions to path
addpath('../plotting');

% Load the data
% filename = 'opt_3_2d.mat';
filename = 'opt_5_2d.mat';
load(filename);
% Q = Q(1:2);
% x_star = x_star(1:2);

Q = {[0.8, 0.5; 0.5, 0.8], [1.25, 0.2; 0.2, 0.7]};
x_star = {[2; 0.5], [-3; -2]};

% Q = Q(4:5);
% x_star = x_star(4:5);

% Plot
fig = figure(10);
clf
fig.Position(1:2) = fig.Position(1:2) / 2;
fig.Position(3:4) = [720, 540];
set(0,'defaultfigurecolor',[1,1,1])
hold all

% tiledlayout(1,1)
% ax = nexttile;
colormap('copper')


% Sample points along the arc
% randarc = [1 nf-1];\
% n_alpha = 100;

N_sample = 400;
N_view = 4;

alpha_sample = linspace(0, 1, N_sample);
alpha_hess_view = linspace(0, 1, N_view);

[Q_sample, x_sample] = get_arc(Q, x_star, alpha_sample);
[Q_view, x_view] = get_arc(Q, x_star, alpha_hess_view);

x_view_cell = mat2cell(x_view, 2, ones(1, N_view));
% Plot level sets

contour_level = 0.5;

h_lev = plot_level_sets(Q_view, x_view_cell, 1000, contour_level, copper(N_view));
for ii = 1 : length(h_lev)
%     h_lev(ii).HandleVisibility = 'Off';
set(h_lev(ii), 'HandleVisibility', 'off', 'LineStyle', '--', 'LineWidth', 2)
end


%Plot the Arc
% plot(x_sample)
Arc = colormapline(x_sample(1, :), x_sample(2, :), [], copper(128));
set(Arc,'LineWidth',3,'HandleVisibility','off')
% [h_arc, h_hess] = plot_arc(Q, x_star, linspace(0, 1, n_alpha), 1000, 1, copper(n_alpha));

%Plot the center points of the sampled level sets
scatter(x_view(1, :), x_view(2, :), 400, copper(N_view), 'filled')


%% set up the options

xlabel('$x_1$', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
ylabel('$x_2$', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
title(sprintf('Set of global optima $\\mathcal{G}$'), 'Interpreter', 'LaTeX', 'FontSize', 20);
% legend('Interpreter', 'LaTeX', 'Location', 'eastoutside');
cbar = colorbar('east', 'AxisLocation','out');
cbar.Label.String = '\alpha';
cbar.Label.FontSize= 12;
cbar.Ticks = linspace(0, 1, N_view);
cbar.TickLabels = {'0', '1/3', '2/3', '1'};
% grid;
% axis square
pbaspect([diff(xlim), diff(ylim), 1])