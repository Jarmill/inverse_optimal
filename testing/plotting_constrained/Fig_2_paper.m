clear all;
close all;
clc;

% Define dimension
n = 2;
% Define number of cost functions
nf = 3;

% Generate random cost functions
rng(351, 'twister');
[Q, phi, x_star] = generate_random_Q_and_phi(n, nf, [0.5, 1.5], [1, 1.1], [1 2], 0);

% Define lower bounds
% xl = zeros(2, 1);
xl = -ones(2, 1);

% Define upper bounds
xu = ones(2, 1);

% % Define x1 + x2 <= 1 & -x1 + x2 <= 0
% C = ones(2, 2);
% C(2, 1) = -1;
% d = [1; 0];

% Define x1 + x2 <= 1 & -x1 + x2 <= 1 & x1 - x2 <= 1 & -x1 - x2 <= 1
C = ones(4, 2);
C([2 4], 1) = -1;
C([3 4], 2) = -1;
d = ones(4, 1);

% No equalities
A = [];
b = [];
% % Define x1 = 0.4
% A = [1, 0];
% b = 0.4;

% Plot
fig = figure;
fig.Position(1:2) = zeros(1, 2);
fig.Position(3:4) = [780, 780];
hold all

% Plot inconstrained global optima set
m_grid = 5000;
h_global_uncons = plot_mathcal_G(Q, x_star, m_grid);
h_global_uncons.DisplayName = sprintf('$\\mathcal{G} \\; = \\{ x \\in \\mathcal{R}^2 | \\alpha_{%d : %d} \\}$', 1, nf);
h_global_uncons.LineStyle = 'None';

% Plot black arcs between minima
h_arcs_between = plot_mathcal_G_all_arcs(Q, x_star, 5000);

% Plot unconstrained minima
coluncons = copper(3);
coluncons = repmat(coluncons(2, :), length(x_star), 1);
h_uncons = plot_x_star(x_star, 150, coluncons);
h_uncons.DisplayName = sprintf('$x^f_{1:%d} \\in \\mathcal{R}^2$', nf);

% Plot feasible set
h_feas = plot_feasible_set_constrained(A, b, C, d, xl, xu);
h_feas.DisplayName = sprintf('$X \\quad (x \\in X)$');

% Plot the set of all optima
n_mesh = 20;
m_grid = 50000;
h_combined_sol = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, m_grid, n_mesh, false);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G}^c =  \\{x \\in X | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
% h_combined_sol.FaceAlpha = 1;
% h_combined_sol.FaceColor = zeros(3, 1);


% Plot minima
colcons = copper(nf);
colcons = repmat(colcons(3, :), nf, 1);
h_star = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, 150, colcons);
h_star(round(length(h_star)/2)).DisplayName = sprintf('$x^{f, c}_{1:%d} \\in X$', nf);

% Plot level sets
x_star_cons = mat2cell([ [h_star.XData]; [h_star.YData] ], 2, ones(1, length(h_star)));
lev = evaluate_quadratics(Q, x_star, x_star_cons);  % Make the ellipses tangential
n_ell = 1000;
collevsets = winter(nf);
collevsets = repmat(collevsets(1, :), nf, 1);
h_lev = plot_level_sets(Q, x_star, n_ell, lev, collevsets);
for ii = 1 : length(h_lev)
    h_lev(ii).HandleVisibility = 'Off';
end
h_lev(round(length(h_lev)/3)).HandleVisibility = 'On';
h_lev(round(length(h_lev)/3)).DisplayName = sprintf('$\\{ x_{1:%d} | f_{1:%d}(x) = %s \\}$', nf, nf, sprintf(['[ ', repmat('%.2f, ', 1, nf-1), '%.2f]'], lev) );

% Esthetics
if n == 3
    view(3)
elseif n == 2
    view(2)
end
xlabel('$X_1$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
ylabel('$X_2$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
zlabel('$X_3$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
title(sprintf('Sets of unconstrained and constrained global optima $\\mathcal{G}$ and $\\mathcal{G}^c$'), 'Interpreter', 'LaTeX', 'FontSize', 15);
h_leg = legend('Interpreter', 'LaTeX', 'Location', 'southoutside', 'orientation', 'horizontal', 'numcolumns', 3, 'FontSize', 13);
% Permute 4th and 5th legend element
% grid;

% Axis range
% off = .25;
xlim([-1.31, 1.27])
ylim([-1.15 1.54])
axis square

%%
h_rec = rectangle('Position',[-0.95 0.1 0.95 0.85], 'EdgeColor', [1 0 0], 'HandleVisibility', 'Off');

%%
exportgraphics(gca, '../../img/cons_and_uncons_2d.png', 'ContentType', 'Image', 'Resolution', 400);

%%
xlim([h_rec.Position(1), sum(h_rec.Position([1 3]))]);
ylim([h_rec.Position(2), sum(h_rec.Position([2 4]))]);
legend('toggle')
%%
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'xlabel', [])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ylabel', [])
set(gca, 'title', [])
