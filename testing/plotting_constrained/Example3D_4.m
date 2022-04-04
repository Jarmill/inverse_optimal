clear all;
close all;
clc;

% Define dimension
n = 3;
% Define number of cost functions
nf = 7;

% Seed the random processes
rng(312, 'twister');

% Design Q's
V1 = eye(n);
D1 = diag([3, 0.5, 1]);
V2 = eye(n);
D2 = diag([3, 1, 0.5]);
V3 = eye(n);
D3 = diag([1, 3, 0.5]);
V4 = eye(n);
D4 = diag([1, 0.5, 3]);
V5 = eye(n);
D5 = diag([0.5, 3, 1]);
V6 = eye(n);
D6 = diag([0.5, 1, 3]);
V7 = eye(n);
D7 = diag([1, 1, 1]);

Q{1} = V1 * D1 * V1.';
Q{2} = V2 * D2 * V2.';
Q{3} = V3 * D3 * V3.';
Q{4} = V4 * D4 * V4.';
Q{5} = V5 * D5 * V5.';
Q{6} = V6 * D6 * V6.';
Q{7} = V7 * D7 * V7.';

% Design phi's
phi{1} = zeros(3, 1);
phi{2} = zeros(3, 1);
phi{3} = zeros(3, 1);
phi{4} = zeros(3, 1);
phi{5} = zeros(3, 1);
phi{6} = zeros(3, 1);
phi{7} = zeros(3, 1);

% Calculate x_stars
x_star{1} = zeros(3, 1);
x_star{2} = zeros(3, 1);
x_star{3} = zeros(3, 1);
x_star{4} = zeros(3, 1);
x_star{5} = zeros(3, 1);
x_star{6} = zeros(3, 1);
x_star{7} = zeros(3, 1);


% Define lower bounds
% xl = zeros(3, 1);
xl = zeros(3, 1);

% Define upper bounds
xu = 6*ones(3, 1);

% No inequalities
C = [];
d = [];

% Define 3*x1 + 2*x2 + 1/2*x3 == 5
A = [3 2 0.5];
b = 5;

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
h_uncons.DisplayName = sprintf('$x^*_{1:%d} \\in \\mathcal{R}^3$', nf);
% % Plot unconstrained global optima set
% h_global_uncons = plot_mathcal_G(Q, x_star, 750);
% h_global_uncons.DisplayName = sprintf('$\\mathcal{G} = \\{ x \\in \\mathcal{R}^3 | \\alpha_{%d : %d} \\}$', 1, nf);
% h_global_uncons.LineStyle = 'None';
% h_global_uncons.FaceAlpha = 0.5;

% Plot feasible set
h_feas = plot_feasible_set_constrained(A, b, C, d, xl, xu);
h_feas.DisplayName = sprintf('$X \\quad (x \\in X)$');


% % Plot the set of all optima
n_mesh = 10;
m_grid = 5000;
h_combined_sol = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, m_grid, n_mesh, true, inf);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G} =  \\{x \\in X | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
% h_combined_sol.FaceAlpha = 1;

% Plot minima
h_star = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, 150, winter(nf));
h_star(round(length(x_star)/2)).DisplayName = sprintf('$x^*_{1:%d} \\in X$', nf);

% Plot level sets
lev = 0.1;
h_level_sets = plot_level_sets(Q, x_star, 1000, lev, winter(length(x_star)));
for ii = 1 : length(h_level_sets)
    h_level_sets(ii).HandleVisibility = 'Off';
    h_level_sets(ii).FaceAlpha = 0.1;
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
% off = .05;
% xlim([xl(1)-off, xu(1)+off])
% ylim([xl(2)-off, xu(2)+off])
axis square