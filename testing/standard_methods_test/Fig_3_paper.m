%EXAMPLE shows how to use GENERATE_RANDOM_Q, GENERATE_RANDOM_X,
%PLOT_MATHCAL_G, PLOT_X_STAR, PLOT_LEVEL_SETS and PLOT_ARC.
clear all;
close all;
clc;

% The MOSEK solver directory
mosek_dir = input_your_mosek_directory;
% The yalmip directory
yalmip_dir = input_your_yalmip_directory;

% Dimension
n = 2;
% Number of random matrices to generateclose a
nf = 5;

% Set random seed
rng(16, 'twister');

% Get random Q's and x's
Q = generate_random_Q(n, nf, [0.5 3]);
x_star = generate_random_x(n, nf, [2 5]);

% Get random y
Nrand = 100;
% y = generate_random_x(n, Nrand, [6 7]);
% y = cell2mat(y);
y0 = [6; 2;];
y = [y0, y0 + randn(2, Nrand-1)];

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
% Using SDP relaxation
x_sdp = zeros(size(y));
alpha_sdp = zeros(nf, size(y, 2));
d_sdp = zeros(1, size(y, 2));
% Using Manopt formulation
x_mano = zeros(size(y));
alpha_mano = zeros(nf, size(y, 2));
d_mano = zeros(1, size(y, 2));

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

% Do IO for all points using SDP and Manopt
% Add mosek and yalmip
addpath(mosek_dir);
addpath(genpath(yalmip_dir));
% For each point
for ii = 1 : size(y, 2)
    
    % Perform sdp inverse optimization
    [d_sdp(:, ii), info_sdp] = sdp_uncons_solve(y(:, ii), Q, x_star);
    % Extract info
    x_sdp(:, ii) = info_sdp.x_rec;
    alpha_sdp(:, ii) = info_sdp.alpha_rec;
    
    % Perform manopt inverse optimization
    [mano_out, info_mano] = manopt_search_uncons(Q, x_star, y(:, ii));
    % Extract info
    d_mano(:, ii) = mano_out.dist;
    x_mano(:, ii) = mano_out.x;
    alpha_mano(:, ii) = mano_out.alpha;
end
% Remove mosek and yalmip
rmpath(mosek_dir);
rmpath(genpath(yalmip_dir));

    
%% Plot
fig = figure;
fig.Position(1:2) = zeros(1, 2);
fig.Position(3:4) = [780, 780];
hold all

% Plot mesh
n_mesh = 25000;
h_combined_sol = plot_mathcal_G(Q, x_star, n_mesh, 2);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G} =  \\{x | \\forall \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
h_combined_sol.FaceAlpha = 0.5;

% Plot black arcs between minima
h_arcs_between = plot_mathcal_G_all_arcs(Q, x_star, 5000);

% Plot solutions
% h_sol = plot_x_star(x_star, 150, winter(nf));
h_sol = plot_x_star(x_star, 150, zeros(nf, 3));
h_sol.DisplayName = sprintf('$x^f_{1:%d}$', length(x_star));
% h_sol.MarkerFaceColor = 'flat';
% h_sol.MarkerFaceAlpha = 1;

% 
% % Plot the query point y and the solutions of inverse optimization
% % cmap = copper(5+3);
% % cmap = cmap(3+1:end, :);
% cmap = jet(5);
% sz = round(150 / ( Nrand + 1)) + 10;
% % sz = 50;
% if n == 2
% h_query = scatter(y(1, :), y(2, :), sz, cmap(1, :), 'o', 'DisplayName', 'Query', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_qcqp = scatter(x_qcqp(1, :), x_qcqp(2, :), sz, cmap(2, :), 'o', 'DisplayName', 'QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_lsq = scatter(x_lsq(1, :), x_lsq(2, :), sz, cmap(3, :), 'o', 'DisplayName', 'LSQ', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_bfqcqp = scatter(x_bfqcqp(1, :), x_bfqcqp(2, :), sz, cmap(4, :), 'o', 'DisplayName', 'BF+QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_sdp = scatter(x_sdp(1, :), x_sdp(2, :), sz, cmap(5, :), 'o', 'DisplayName', 'SDP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% end
% if n == 3
% h_query = scatter3(y(1, :), y(2, :), y(3, :), sz, cmap(1, :), 'o', 'DisplayName', 'Query', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_qcqp = scatter3(x_qcqp(1, :), x_qcqp(2, :), x_qcqp(3, :), sz, cmap(2, :), 'o', 'DisplayName', 'QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_lsq = scatter3(x_lsq(1, :), x_lsq(2, :), x_lsq(3, :), sz, cmap(3, :), 'o', 'DisplayName', 'LSQ', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_bfqcqp = scatter3(x_bfqcqp(1, :), x_bfqcqp(2, :), x_bfqcqp(3, :), sz, cmap(4, :), 'o', 'DisplayName', 'BF+QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_sdp = scatter3(x_sdp(1, :), x_sdp(2, :), x_sdp(3, :), sz, cmap(5, :), 'o', 'DisplayName', 'SDP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% end


% Plot the query point y and the solutions of inverse optimization
% cmap = copper(5+3);
% cmap = cmap(3+1:end, :);
% cmap = flipud(jet(7));
cmap = linspecer(9);
cmap = cmap(5:end, :);
% sz = round(150 / ( Nrand + 1)) + 10;
sz = 250;
% Choose a random point out of the Nrand ones
% iir = randi([1, Nrand], 1, 1);
iir = 1;

% Plot random index point and the projection lines
if n == 2

% % Random index point
% h_query = scatter(y(1, iir), y(2, iir), sz, cmap(1, :), 'o', 'DisplayName', 'Query', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_bfqcqp = scatter(x_bfqcqp(1, iir), x_bfqcqp(2, iir), sz, cmap(4, :), 'o', 'DisplayName', 'BF+QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_qcqp = scatter(x_qcqp(1, iir), x_qcqp(2, iir), sz, cmap(2, :), 'o', 'DisplayName', 'QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_lsq = scatter(x_lsq(1, iir), x_lsq(2, iir), sz, cmap(3, :), 'o', 'DisplayName', 'LSQ', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_sdp = scatter(x_sdp(1, iir), x_sdp(2, iir), sz, cmap(5, :), 'o', 'DisplayName', 'SDP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% 
% % Plot connecting line
% h_line_bfqcqp = plot([y(1, iir), x_bfqcqp(1, iir)], [y(2, iir), x_bfqcqp(2, iir)], '-', 'Color', cmap(4, :),  'HandleVisibility', 'Off', 'LineWidth', 2);
% h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(2, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
% h_line_lsq = plot([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], ':', 'Color', cmap(3, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
% h_line_sdp = plot([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], '-.', 'Color', cmap(5, :), 'HandleVisibility', 'Off', 'LineWidth', 2);

% % Random index point
% h_query = scatter(y(1, iir), y(2, iir), sz, cmap(1, :), 'o', 'DisplayName', 'y', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_bfqcqp = scatter(x_bfqcqp(1, iir), x_bfqcqp(2, iir), sz, cmap(4, :), 'o', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_qcqp = scatter(x_qcqp(1, iir), x_qcqp(2, iir), sz, cmap(2, :), '+', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_lsq = scatter(x_lsq(1, iir), x_lsq(2, iir), sz, cmap(3, :), 's', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_sdp = scatter(x_sdp(1, iir), x_sdp(2, iir), sz, cmap(5, :), 'x', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% 
% % Plot connecting line
% h_line_bfqcqp = plot([y(1, iir), x_bfqcqp(1, iir)], [y(2, iir), x_bfqcqp(2, iir)], '-', 'Color', cmap(4, :),  'DisplayName', '$\textrm{Proj}_{\mathcal{G}}(y)$', 'LineWidth', 2);
% % h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(2, :), 'DisplayName', '$\alpha^* = \textrm{argmin}_{\alpha} ||y - x^*||_2^2; \; x^* = \textrm{argmin}_{x} f_{\alpha^*}(x)$', 'LineWidth', 2);
% % h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(2, :), 'DisplayName', '$x^* = \pi^x \textrm{argmin}_{(x, \alpha) \in \hat{\mathcal{G}}} ||y-x||_2^2$', 'LineWidth', 2);
% h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(2, :), 'DisplayName', 'ProjGM', 'LineWidth', 2);
% h_line_lsq = plot([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], ':', 'Color', cmap(3, :), 'DisplayName', 'Keshavaraz', 'LineWidth', 2);
% h_line_sdp = plot([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], '-.', 'Color', cmap(5, :), 'DisplayName', 'LMI', 'LineWidth', 2);


% Random index point
h_query = scatter(y(1, iir), y(2, iir), sz, cmap(1, :), 'o', 'DisplayName', 'y', 'LineWidth', 4);
h_bfqcqp = scatter(x_bfqcqp(1, iir), x_bfqcqp(2, iir), sz, cmap(2, :), 'o', 'HandleVisibility', 'Off', 'LineWidth', 4);
h_qcqp = scatter(x_qcqp(1, iir), x_qcqp(2, iir), sz, cmap(3, :), '+', 'HandleVisibility', 'Off', 'LineWidth', 4);
h_lsq = scatter(x_lsq(1, iir), x_lsq(2, iir), sz, cmap(4, :), 's', 'HandleVisibility', 'Off', 'LineWidth', 4);
h_sdp = scatter(x_sdp(1, iir), x_sdp(2, iir), sz, cmap(5, :), 'x', 'HandleVisibility', 'Off', 'LineWidth', 4);

% Plot connecting line
h_line_bfqcqp = plot([y(1, iir), x_bfqcqp(1, iir)], [y(2, iir), x_bfqcqp(2, iir)], '-', 'Color', cmap(2, :),  'DisplayName', '$\textrm{Proj}_{\mathcal{G}}(y)$', 'LineWidth', 2);
% h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(2, :), 'DisplayName', '$\alpha^* = \textrm{argmin}_{\alpha} ||y - x^*||_2^2; \; x^* = \textrm{argmin}_{x} f_{\alpha^*}(x)$', 'LineWidth', 2);
% h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(2, :), 'DisplayName', '$x^* = \pi^x \textrm{argmin}_{(x, \alpha) \in \hat{\mathcal{G}}} ||y-x||_2^2$', 'LineWidth', 2);
h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(3, :), 'DisplayName', 'ProjGM', 'LineWidth', 2);
h_line_lsq = plot([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], ':', 'Color', cmap(4, :), 'DisplayName', 'Keshavaraz', 'LineWidth', 2);
h_line_sdp = plot([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], '-.', 'Color', cmap(5, :), 'DisplayName', 'LMI', 'LineWidth', 2);


end
if n == 3
        
% % Random index point
% h_query = scatter3(y(1, iir), y(2, iir), y(3, iir), sz, cmap(1, :), 'o', 'DisplayName', 'Query', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_qcqp = scatter3(x_qcqp(1, iir), x_qcqp(2, iir), x_qcqp(3, iir), sz, cmap(2, :), 'o', 'DisplayName', 'QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_lsq = scatter3(x_lsq(1, iir), x_lsq(2, iir), x_lsq(3, iir), sz, cmap(3, :), 'o', 'DisplayName', 'LSQ', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_bfqcqp = scatter3(x_bfqcqp(1, iir), x_bfqcqp(2, iir), x_bfqcqp(3, iir), sz, cmap(4, :), 'o', 'DisplayName', 'BF+QCQP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% h_sdp = scatter3(x_sdp(1, iir), x_sdp(2, iir), x_sdp(3, iir), sz, cmap(5, :), 'o', 'DisplayName', 'SDP', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
% 
% % Plot connecting line
% h_line_bfqcqp = plot3([y(1, iir), x_bfqcqp(1, iir)], [y(2, iir), x_bfqcqp(2, iir)], [y(3, iir), x_bfqcqp(3, iir)], '-', 'Color', cmap(4, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
% h_line_qcqp = plot3([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], [y(3, iir), x_qcqp(3, iir)], '--', 'Color', cmap(2, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
% h_line_lsq = plot3([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], [y(3, iir), x_lsq(3, iir)], ':', 'Color', cmap(3, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
% h_line_sdp = plot3([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], [y(3, iir), x_sdp(3, iir)], '-.', 'Color', cmap(5, :), 'HandleVisibility', 'Off', 'LineWidth', 2);

% Random index point
h_query = scatter3(y(1, iir), y(2, iir), y(3, iir), sz, cmap(1, :), 'o', 'DisplayName', 'Query', 'LineWidth', 3, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_bfqcqp = scatter3(x_bfqcqp(1, iir), x_bfqcqp(2, iir), x_bfqcqp(3, iir), sz, cmap(4, :), 'o', 'HandleVisibility', 'Off', 'LineWidth', 3, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_qcqp = scatter3(x_qcqp(1, iir), x_qcqp(2, iir), x_qcqp(3, iir), sz, cmap(2, :), '+', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_lsq = scatter3(x_lsq(1, iir), x_lsq(2, iir), x_lsq(3, iir), sz, cmap(3, :), 's', 'HandleVisibility', 'Off', 'LineWidth', 3, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_sdp = scatter3(x_sdp(1, iir), x_sdp(2, iir), x_sdp(3, iir), sz, cmap(5, :), 'x', 'HandleVisibility', 'Off', 'LineWidth', 3, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);

% Plot connecting line
h_line_bfqcqp = plot3([y(1, iir), x_bfqcqp(1, iir)], [y(2, iir), x_bfqcqp(2, iir)], [y(3, iir), x_bfqcqp(3, iir)], '-', 'Color', cmap(4, :), 'DisplayName', '$\textrm{Proj}_{\mathcal{G}}(y)$', 'LineWidth', 2);
h_line_qcqp = plot3([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], [y(3, iir), x_qcqp(3, iir)], '--', 'Color', cmap(2, :), 'DisplayName', 'ProjGM', 'LineWidth', 2);
h_line_lsq = plot3([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], [y(3, iir), x_lsq(3, iir)], ':', 'Color', cmap(3, :), 'DisplayName', 'Keshavaraz', 'LineWidth', 2);
h_line_sdp = plot3([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], [y(3, iir), x_sdp(3, iir)], '-.', 'Color', cmap(5, :), 'DisplayName', 'LMI', 'LineWidth', 2);

end


% Plot level sets
lev = 1;
lev_col = winter(nf);
h_lev = plot_level_sets(Q, x_star, 25000, lev, lev_col);
for ii = 1 : length(h_lev)
    h_lev(ii).HandleVisibility = 'Off';
end
% h_lev(round(length(h_lev)/3)).HandleVisibility = 'On';
% h_lev(round(length(h_lev)/3)).DisplayName = sprintf('$\\{ x | f(x) = %.2f\\}$', lev);

% Esthetics
if n == 3
    view(3)
    % For rng seed 16
    a = rng;
    if a.Seed == 16
        view(-150, 20);
    end
elseif n == 2
    view(2)
end
xlabel('$X_1$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
ylabel('$X_2$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
zlabel('$X_3$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
title(sprintf('Projection of a random point onto unconstrained set of global optima $\\mathcal{G}$'), 'Interpreter', 'LaTeX', 'FontSize', 15);
legend('Interpreter', 'LaTeX', 'Location', 'southoutside', 'orientation', 'vertical', 'FontSize', 15, 'numcolumns', 4);
grid;
xlim([-5, 6.5]);
ylim([-5, 6]);
% zlim([-6.5, 6.5]);
axis square

% exportgraphics(gca, '../../img/projection_on_uncons.pdf', 'ContentType', 'Vector');
exportgraphics(gca, '../../img/projection_on_uncons.png', 'ContentType', 'Image', 'Resolution', 400);

%% Plot distances
% Remove test point from colormap
cmap = cmap(2:end, :);
% Sort different methods
dist_bars = {'$\textrm{Proj}_{\mathcal{G}}(y)$', 'ProjGM', 'Keshavaraz', 'LMI'};
dist_arr = [mean(d_bfqcqp), mean(d_qcqp), mean(d_lsq), mean(d_sdp)];
dist_mat = [d_sdp; d_qcqp; d_lsq; d_bfqcqp];

xtick_arr = 1 : length(dist_arr);

% Sort in ascending manner
[dist_arr, sortind] = sort(dist_arr, 'ascend');
dist_bars = dist_bars(sortind);
dist_mat = dist_mat(sortind, :);
cmap = cmap(sortind, :);

fig = figure;
fig.Position(1:2) = zeros(1, 2);
fig.Position(3:4) = [700, 540];
hold all
h_bar = bar(xtick_arr, dist_arr, 'FaceColor', 'Flat', 'CData', cmap);
xticks(xtick_arr);
xticklabels(dist_bars);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
title({sprintf('Mean distance of randomly generated points \n to set of global optima $\\mathcal{G}$ estimated by different methods')}, 'interpreter', 'latex', 'FontSize', 15);
ylabel('Mean Distance', 'Interpreter', 'latex');

exportgraphics(gca, '../../img/projection_on_uncons_mean_distance.png', 'ContentType', 'Image', 'Resolution', 400);

% figure;
% hold all;
% for ii = 1 : length(dist_arr)
%     plot(log10(dist_mat(ii, :)), 'DisplayName', dist_bars{ii});
% end
% xlabel('Example No.', 'Interpreter', 'Latex', 'FontSize', 15);
% ylabel('Distance', 'Interpreter', 'Latex', 'FontSize', 15);
% title(sprintf('Projection distance of randomly generated points to set of global optima $\\mathcal{G}$ obtained by different methods.'));