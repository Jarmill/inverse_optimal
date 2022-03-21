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

% Get random y
Nrand = 100;
y = generate_random_x(n, Nrand, [6 7]);
y = cell2mat(y);

% Prealocate solutions of IOC

% Using QCQP formulation
x_qcqp = zeros(size(y));
alpha_qcqp = zeros(nf, size(y, 2));
d_qcqp = zeros(1, size(y, 2));
% Using Keshavaraz Relaxation (LSQ formulation)
x_lsq = zeros(size(y));
alpha_lsq = zeros(nf, size(y, 2));
d_lsq = zeros(1, size(y, 2));
% Using Brute Force formulation
x_bf = zeros(size(y));
alpha_bf = zeros(nf, size(y, 2));
d_bf = zeros(1, size(y, 2));
% Combining global Brute Force with local QCQP
x_bfqcqp = zeros(size(y));
alpha_bfqcqp = zeros(nf, size(y, 2));
d_bfqcqp = zeros(1, size(y, 2));
% Combining global Brute Force with local LSQ
x_bflsq = zeros(size(y));
alpha_bflsq = zeros(nf, size(y, 2));
d_bflsq = zeros(1, size(y, 2));

% Do IOC for all points
for ii = 1 : size(y, 2)
    
    % Generate random initial point
    [x0, alpha0] = generate_random_direct_solution(Q, x_star);

    % Perform IOC
    [x_qcqp(:, ii), alpha_qcqp(:, ii), d_qcqp(:, ii)] = QCQP_IO(Q, x_star, y(:, ii), [x0;alpha0]);
    [x_lsq(:, ii), alpha_lsq(:, ii), d_lsq(:, ii)] = LSQ_IO(Q, x_star, y(:, ii), alpha0);
    [x_bf(:, ii), alpha_bf(:, ii), d_bf(:, ii)] = BF_IO(Q, x_star, y(:, ii), 1500);    
    [x_bfqcqp(:, ii), alpha_bfqcqp(:, ii), d_bfqcqp(:, ii)] = QCQP_IO(Q, x_star, y(:, ii), [x_bf(:, ii); alpha_bf(:, ii)]);
    [x_bflsq(:, ii), alpha_bflsq(:, ii), d_bflsq(:, ii)] = LSQ_IO(Q, x_star, y(:, ii), alpha_bf(:, ii));   
end
    
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

% Plot the query point y and the solutions of inverse optimization
cmap = copper(Nrand);
sz = round(150 / ( Nrand + 1)) + 10;
if n == 2
h_query = scatter(y(1, :), y(2, :), sz, cmap, 'x', 'DisplayName', 'Query', 'LineWidth', 2);
h_qcqp = scatter(x_qcqp(1, :), x_qcqp(2, :), sz, cmap, 'o', 'DisplayName', 'QCQP', 'LineWidth', 2);
h_lsq = scatter(x_lsq(1, :), x_lsq(2, :), sz, cmap, '^', 'DisplayName', 'LSQ', 'LineWidth', 2);
h_bf = scatter(x_bf(1, :), x_bf(2, :), sz, cmap, 's', 'DisplayName', 'BF', 'LineWidth', 2);
h_bfqcqp = scatter(x_bfqcqp(1, :), x_bfqcqp(2, :), sz, cmap, 'p', 'DisplayName', 'BF+QCQP', 'LineWidth', 2);
h_bflsq = scatter(x_bflsq(1, :), x_bflsq(2, :), sz, cmap, 'h', 'DisplayName', 'BF+LSQ', 'LineWidth', 2);
end
if n == 3
h_query = scatter3(y(1, :), y(2, :), y(3, :), sz, cmap, 'x', 'DisplayName', 'Query', 'LineWidth', 2);
h_qcqp = scatter3(x_qcqp(1, :), x_qcqp(2, :), x_qcqp(3, :), sz, cmap, 'o', 'DisplayName', 'QCQP', 'LineWidth', 2);
h_lsq = scatter3(x_lsq(1, :), x_lsq(2, :), x_lsq(3, :), sz, cmap, '^', 'DisplayName', 'LSQ', 'LineWidth', 2);
h_bf = scatter3(x_bf(1, :), x_bf(2, :), x_bf(3, :), sz, cmap, 's', 'DisplayName', 'BF', 'LineWidth', 2);
h_bfqcqp = scatter3(x_bfqcqp(1, :), x_bfqcqp(2, :), x_bfqcqp(3, :), sz, cmap, 'p', 'DisplayName', 'BF+QCQP', 'LineWidth', 2);
h_bflsq = scatter3(x_bflsq(1, :), x_bflsq(2, :), x_bflsq(3, :), sz, cmap, 'h', 'DisplayName', 'BF+LSQ', 'LineWidth', 2);
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


%% Plot distances
dist_bars = {'QCQP', 'LSQ', 'BF', 'BF+QCQP', 'BF+LSQ'};
dist_arr = [mean(d_qcqp), mean(d_lsq), mean(d_bf), mean(d_bfqcqp), mean(d_bflsq)];
xtick_arr = 1 : length(dist_arr);

% Sort in ascending manner
[dist_arr, sortind] = sort(dist_arr, 'ascend');
dist_bars = dist_bars(sortind);

figure;
h_bar = bar(xtick_arr, dist_arr, 'FaceColor', 'Flat', 'CData', copper(length(dist_bars)));
xticks(xtick_arr);
xticklabels(dist_bars);
set(gca, 'TickLabelInterpreter', 'latex');
title({sprintf('Mean distance of randomly generated points \n to set of global optima $\\mathcal{G}$ estimated by different standard methods')}, 'interpreter', 'latex');