%EXAMPLE shows how to use GENERATE_RANDOM_Q, GENERATE_RANDOM_X,
%PLOT_MATHCAL_G, PLOT_X_STAR, PLOT_LEVEL_SETS and PLOT_ARC.
clear all;
close all;
clc;

% Dimension
n = 3;
% Number of random matrices to generate
nf = 5;

% Get random Q's and x's
Q = generate_random_Q(n, nf, [0.5 3]);
x_star = generate_random_x(n, nf, [2 5]);
    
% Plot
fig = figure;
fig.Position(1:2) = fig.Position(1:2) / 2;
fig.Position(3:4) = [720, 540];
hold all

% Plot mesh
n_mesh = 750;
h_combined_sol = plot_mathcal_G(Q, x_star, n_mesh, 2);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G} =  \\{x | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
h_combined_sol.FaceAlpha = 0.5;

% Plot solutions
h_sol = plot_x_star(x_star, 150, winter(nf));
h_sol.DisplayName = sprintf('$x^*_{1:%d}$', length(x_star));

% Plot level sets
h_lev = plot_level_sets(Q, x_star, 1000, 1, winter(nf));
for ii = 1 : length(h_lev)
    h_lev(ii).HandleVisibility = 'Off';
end

% % Plot arc
% randarc = [1 nf];
% alpha_arc = linspace(0, 1, 1000);
% alpha_hess = linspace(0, 1, 5);
% [h_arc, h_hess, x_arc, H_arc, x_hess, H_hess] = plot_arc(Q( randarc ), x_star( randarc ), alpha_arc, alpha_hess, 1000, 1);
% for ii = 1 : length(h_hess)
%     h_hess(ii).HandleVisibility = 'Off';
% end
% Esthetics
if n == 3
    view(3)
elseif n == 2
    view(2)
end
xlabel('$X_1$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
ylabel('$X_2$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
zlabel('$X_3$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
title(sprintf('Set of global optima $\\mathcal{G}$'), 'Interpreter', 'LaTeX');
legend('Interpreter', 'LaTeX', 'Location', 'eastoutside');
grid;
axis square