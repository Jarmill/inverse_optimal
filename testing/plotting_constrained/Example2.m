clear all;
close all;
clc;

% Define dimension
n = 2;
% Define number of cost functions
nf = 3;

% Seed the random processes
rng(351, 'twister');

% Design Q's
V1 = eye(n);
D1 = diag([3 0.5]);
V2 = eye(n);
D2 = diag([0.5 3]);
V3 = eye(n);
D3 = diag([1 1]);

Q{1} = V1 * D1 * V1.';
Q{2} = V2 * D2 * V2.';
Q{3} = V3 * D3 * V3.';

% Design phi's
phi{1} = zeros(2, 1);
phi{2} = zeros(2, 1);
phi{3} = zeros(2, 1);

% Calculate x_stars
x_star{1} = zeros(2, 1);
x_star{2} = zeros(2, 1);
x_star{3} = zeros(2, 1);

% Define lower bounds
% xl = zeros(2, 1);
xl = zeros(2, 1);

% Define upper bounds
xu = 3*ones(2, 1);

% % Define x1 + x2 <= 1 & -x1 + x2 <= 1 & x1 - x2 <= 1 & -x1 - x2 <= 1
% C = ones(4, 2);
% C([2 4], 1) = -1;
% C([3 4], 2) = -1;
% d = ones(4, 1);

% Define no inequalities
C = [];
d = [];

% Define 2*x1 + x2 = 4;
A = [2, 1];
b = 4;


% % No equalities
% A = [];
% b = [];

% Plot
fig = figure;
fig.Position(1:2) = fig.Position(1:2) / 2;
fig.Position(3:4) = [720, 540];
hold all

% Plot unconstrained minima
h_uncons = plot_x_star(x_star, 150, winter(length(x_star)));
h_uncons.DisplayName = sprintf('$x^*_{1:%d} \\in \\mathcal{R}^2$', nf);
% Plot inconstrained global optima set
h_global_uncons = plot_mathcal_G(Q, x_star, 750);
h_global_uncons.DisplayName = sprintf('$\\mathcal{G} = \\{ x \\in \\mathcal{R}^2 | \\alpha_{%d : %d} \\}$', 1, nf);
h_global_uncons.LineStyle = 'None';

% Plot feasible set
h_feas = plot_feasible_set_constrained(A, b, C, d, xl, xu);
h_feas.DisplayName = sprintf('$X \\quad (x \\in X)$');

% Plot the set of all optima
n_mesh = 20;
m_grid = 1500;
h_combined_sol = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, m_grid, n_mesh, false);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G}_{X} =  \\{x \\in X | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
% h_combined_sol.FaceAlpha = 1;


% Plot minima
h_star = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, 150, copper(nf));
h_star(round(length(h_star)/2)).DisplayName = sprintf('$x^*_{1:%d} \\in X$', nf);

% Plot level sets
lev = 0.1;
h_level_sets = plot_level_sets(Q, x_star, 1000, lev, winter(length(x_star)));
for ii = 1 : length(h_level_sets)
    h_level_sets(ii).HandleVisibility = 'Off';
end
h_level_sets(round(length(h_level_sets) / 2)).HandleVisibility = 'On';
h_level_sets(round(length(h_level_sets) / 2)).DisplayName = sprintf('$f_i = %.2f$', lev);

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

% Axis range
% off = .25;
% xlim([xl(1)-off, xu(1)+off])
% ylim([xl(2)-off, xu(2)+off])
axis square