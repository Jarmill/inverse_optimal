clear all;
close all;
clc;

% Define dimension
n = 3;
% Define number of cost functions
nf = 6;

% Seed the random processes
rng(312, 'twister');

% The MOSEK solver directory
mosek_dir = input_your_mosek_directory;
% The yalmip directory
yalmip_dir = input_your_yalmip_directory;


% Design Q's
V1 = eye(n);
D1 = diag([2, 0.5, 1]);
V2 = eye(n);
D2 = diag([2, 1, 0.5]);
V3 = eye(n);
D3 = diag([1, 2, 0.5]);
V4 = eye(n);
D4 = diag([1, 0.5, 2]);
V5 = eye(n);
D5 = diag([0.5, 2, 1]);
V6 = eye(n);
D6 = diag([0.5, 1, 2]);
% V7 = eye(n);
% D7 = diag([1, 1, 1]);

Q{1} = V1 * D1 * V1.';
Q{2} = V2 * D2 * V2.';
Q{3} = V3 * D3 * V3.';
Q{4} = V4 * D4 * V4.';
Q{5} = V5 * D5 * V5.';
Q{6} = V6 * D6 * V6.';
% Q{7} = V7 * D7 * V7.';

% Design phi's
phi{1} = zeros(3, 1);
phi{2} = zeros(3, 1);
phi{3} = zeros(3, 1);
phi{4} = zeros(3, 1);
phi{5} = zeros(3, 1);
phi{6} = zeros(3, 1);
% phi{7} = zeros(3, 1);

% Calculate x_stars
x_star{1} = zeros(3, 1);
x_star{2} = zeros(3, 1);
x_star{3} = zeros(3, 1);
x_star{4} = zeros(3, 1);
x_star{5} = zeros(3, 1);
x_star{6} = zeros(3, 1);
% x_star{7} = zeros(3, 1);


% Define lower bounds
% xl = zeros(3, 1);
xl = zeros(3, 1);

% Define upper bounds
xu = 3*ones(3, 1);

% No inequalities
C = [];
d = [];

% Define x1 + x2 + x3 == 1
% A = ones(1, 3);
A = [.5 .5 .5];
% A = [1 0 0];
b = 2;

    %% Generate random points y from the shape defined by constraints
% Get vertices of feasible polytope
V = get_vertices_of_polytope(A, b, C, d, xl, xu);
% Get number of vertices
nv = size(V, 2);

% Generate random points in the feasible polytope
% Number of points to generate
Nrand = 200;
% Parametrization of polytope vertices
beta = rand(nv, Nrand); % between 0 and 1
beta = beta ./ sum(beta); % all normalized to sum to 1
% Generate random points from the polytope
y = V * beta;

% Prealocate inverse optimization for each one of them
x_bl = zeros(n, Nrand);
alpha_bl = zeros(nf, Nrand);
d_bl = zeros(1, Nrand);

x_bf = zeros(n, Nrand);
alpha_bf = zeros(nf, Nrand);
d_bf = zeros(1, Nrand);

x_qcqp = zeros(n, Nrand);
alpha_qcqp = zeros(nf, Nrand);
d_qcqp = zeros(1, Nrand);

x_lsq = zeros(n, Nrand);
alpha_lsq= zeros(nf, Nrand);
d_lsq = zeros(1, Nrand);

x_sdp = zeros(n, Nrand);
alpha_sdp = zeros(nf, Nrand);
d_sdp = zeros(1, Nrand);

% Perform IO on each sample
for ii = 1 : Nrand
    
    % Initialize at middle of simplex
    Qmix = 0;
    phimix = 0;
    for cf = 1 : nf
        Qmix = Qmix + 1/nf * Q{cf};
        phimix = phimix + 1/nf * phi{cf};
    end
    
    % Prepare the initialization
    [x0, ~, ~, ~, mult0] = quadprog(Qmix, phimix, C, d, A, b, xl, xu);
    alpha0 = 1/nf * ones(nf, 1);
    lambda0 = mult0.eqlin;
    mu0 = mult0.ineqlin;
    nu0_minus = mult0.lower;
    nu0_plus = mult0.upper;
    
    % Initialize
    z0_lsq = [alpha0; lambda0; mu0; nu0_minus; nu0_plus;];
    
    % Perform IO
    [x_bf(:, ii), alpha_bf(:, ii), d_bf(:, ii)] = BF_IO_constrained(Q, phi, A, b, C, d, xl, xu, y(:, ii), 1000, 20);
    [x_lsq(:, ii), alpha_lsq(:, ii), d_lsq(:, ii)] = LSQ_IO_constrained(Q, phi, A, b, C, d, xl, xu, y(:, ii), z0_lsq);

   
    

end

% Add mosek and yalmip
addpath(mosek_dir);
addpath(genpath(yalmip_dir));
    
% For each sample perform specifically the SDP
for ii = 1 : Nrand
    % Get the matrices in corresponding form
    Aeqsdp = A;
    beqsdp = b;
    Asdp = [C; -eye(n); eye(n);];
    bsdp = [d; -xl; xu;];
    % Solve
    [d_sdp(:, ii), info_sdp] = sdp_cons_lin_solve(y(:, ii), Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);
    x_sdp(:, ii) = info_sdp.x_rec;
    alpha_sdp(:, ii) = info_sdp.alpha_rec;
    
    % Solve QCQP using fmincon
    [x_fmincon, dist_fmincon, alpha_con, exitflag, output] = fmincon_cons_search(y(:, ii), Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);

    alpha_qcqp(:, ii) = alpha_con;
    d_qcqp(:, ii) = dist_fmincon;
    
%     % Redo DOC to get x
%     Qmix = 0;
%     phimix = 0;
%     for cf = 1 : nf
%         Qmix = Qmix + alpha_con(cf) * Q{cf};
%         phimix = phimix + alpha_con(cf) * phi{cf};
%     end    
%     
%     [x_qcqp(:, ii), ~, ~, ~, mult0] = quadprog(Qmix, phimix, C, d, A, b, xl, xu);
    x_qcqp(:, ii) = x_fmincon;
end

% Remove mosek and yalmip
rmpath(mosek_dir);
rmpath(yalmip_dir);

%%

% % No equalities
% A = [];
% b = [];

%% Plot
fig = figure;
fig.Position(1:2) = zeros(1, 2);
fig.Position(3:4) = [780, 780];
hold all

% Plot unconstrained minima
x_star_uncons_sz = 150;
% x_star_uncons_col = winter(length(x_star));
x_star_uncons_col = zeros(length(x_star), 3);
h_uncons = plot_x_star(x_star, x_star_uncons_sz, x_star_uncons_col);
h_uncons.DisplayName = sprintf('$x^f_{1:%d} \\in \\mathcal{R}^3$', nf);
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
m_grid = 1000;
h_combined_sol = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, m_grid, n_mesh, true, inf);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G}^c =  \\{x \\in X | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
h_combined_sol.FaceAlpha = 0.7;
h_combined_sol.FaceColor = zeros(3, 1);


% Plot minima
x_star_cons_sz = 250;
% x_star_cons_col = copper(nf);
x_star_cons_col = copper(3);
x_star_cons_col = repmat(x_star_cons_col(3, :), nf, 1);
h_star = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, x_star_cons_sz, x_star_cons_col);
h_star(round(length(x_star)/2)).DisplayName = sprintf('$x^{f, c}_{1:%d} \\in X$', nf);

% Plot level sets
lev = 0.05;
h_level_sets = plot_level_sets(Q, x_star, 1000, lev, winter(length(x_star)));
for ii = 1 : length(h_level_sets)
    h_level_sets(ii).HandleVisibility = 'Off';
    h_level_sets(ii).FaceAlpha = 0.1;
end
h_level_sets(round(length(h_level_sets) / 2)).HandleVisibility = 'On';
h_level_sets(round(length(h_level_sets) / 2)).DisplayName = sprintf('$f_i = %.2f$', lev);


% % Plot optimal points
% cmap_sol = summer(3);
% sz = 100;
% h_sol = scatter3(y(1, :), y(2, :), y(3, :), sz, [1, 0, 0], 'filled', 'DisplayName', '$y$');
% h_bf = scatter3(x_bf(1, :), x_bf(2, :), x_bf(3, :), sz, cmap_sol(2, :), 'filled', 'DisplayName', 'Proj$_{ \\mathcal{G}^c} y$');
% h_lsq = scatter3(x_lsq(1, :), x_lsq(2, :), x_lsq(3, :), sz, cmap_sol(1, :), 'filled', 'DisplayName', '$x_{AKKT}$');
% % h_lsq = scatter3(x_lsq(1, :), x_lsq(2, :), x_lsq(3, :), sz, cmap_sol(3, :));

% Plot the query point y and the solutions of inverse optimization
cmap = linspecer(9);
cmap = cmap(5:end, :);
% sz = round(150 / ( Nrand + 1)) + 10;
sz = 250;
% Choose a random point out of the Nrand ones
iir = randi([1, Nrand], 1, 1);
[~, iir] = max(d_sdp);

% Plot random index point and the projection lines
if n == 2

% Random index point
h_query = scatter(y(1, iir), y(2, iir), sz, cmap(1, :), 'o', 'DisplayName', 'y', 'LineWidth', 2);
h_bf = scatter(x_bf(1, iir), x_bf(2, iir), sz, cmap(2, :), 'o', 'HandleVisibility', 'Off', 'LineWidth', 2);
h_qcqp = scatter(x_qcqp(1, iir), x_qcqp(2, iir), sz, cmap(3, :), '+', 'HandleVisibility', 'Off', 'LineWidth', 2);
h_lsq = scatter(x_lsq(1, iir), x_lsq(2, iir), sz, cmap(4, :), 's', 'HandleVisibility', 'Off', 'LineWidth', 2);
h_sdp = scatter(x_sdp(1, iir), x_sdp(2, iir), sz, cmap(5, :), 'x', 'HandleVisibility', 'Off', 'LineWidth', 2);

% Plot connecting line
h_line_bf = plot([y(1, iir), x_bf(1, iir)], [y(2, iir), x_bf(2, iir)], '-', 'Color', cmap(2, :),  'DisplayName', '$\textrm{Proj}_{\mathcal{G}}(y)$', 'LineWidth', 2);
h_line_qcqp = plot([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], '--', 'Color', cmap(3, :), 'DisplayName', 'ProjGM', 'LineWidth', 2);
h_line_lsq = plot([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], ':', 'Color', cmap(4, :), 'DisplayName', 'Keshavaraz', 'LineWidth', 2);
h_line_sdp = plot([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], '-.', 'Color', cmap(5, :), 'DisplayName', 'LMI', 'LineWidth', 2);

end
if n == 3
    
% Random index point
h_query = scatter3(y(1, iir), y(2, iir), y(3, iir), sz, cmap(1, :), 'o', 'DisplayName', 'y', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_bf = scatter3(x_bf(1, iir), x_bf(2, iir), x_bf(3, iir), sz, cmap(2, :), 'o', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_qcqp = scatter3(x_qcqp(1, iir), x_qcqp(2, iir), x_qcqp(3, iir), sz, cmap(3, :), '+', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_lsq = scatter3(x_lsq(1, iir), x_lsq(2, iir), x_lsq(3, iir), sz, cmap(4, :), 's', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);
h_sdp = scatter3(x_sdp(1, iir), x_sdp(2, iir), x_sdp(3, iir), sz, cmap(5, :), 'x', 'HandleVisibility', 'Off', 'LineWidth', 2, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 1);

% Plot connecting line
h_line_bf = plot3([y(1, iir), x_bf(1, iir)], [y(2, iir), x_bf(2, iir)], [y(3, iir), x_bf(3, iir)], '-', 'Color', cmap(2, :), 'DisplayName', '$\textrm{Proj}_{\mathcal{G}^c}(y)$', 'LineWidth', 2);
h_line_qcqp = plot3([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], [y(3, iir), x_qcqp(3, iir)], '--', 'Color', cmap(3, :), 'DisplayName', 'ProjGM', 'LineWidth', 2);
h_line_lsq = plot3([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], [y(3, iir), x_lsq(3, iir)], ':', 'Color', cmap(4, :), 'DisplayName', 'Keshavaraz', 'LineWidth', 2);
h_line_sdp = plot3([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], [y(3, iir), x_sdp(3, iir)], '-.', 'Color', cmap(5, :), 'DisplayName', 'LMI', 'LineWidth', 2);

end


% Esthetics
if n == 3
    view(3)
elseif n == 2
    view(2)
end
xlabel('$X_1$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
ylabel('$X_2$-axis', 'FontSize', 15, 'FontWeight', 'bol', 'Interpreter', 'LaTeX');
zlabel('$X_3$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
title(sprintf('Projection of a random point onto constrained set of global optima $\\mathcal{G}^c$'), 'Interpreter', 'LaTeX', 'FontSize', 15);
legend('Interpreter', 'LaTeX', 'Location', 'southoutside', 'orientation', 'vertical', 'FontSize', 15, 'numcolumns', 4);
grid;

% Axis range
% off = .05;
% xlim([xl(1)-off, xu(1)+off])
% ylim([xl(2)-off, xu(2)+off])
axis square

 
%% Plot distances
% Remove test point from colormap
cmap = cmap(2:end, :);
% Sort different methods
dist_bars = {'$\textrm{Proj}_{\mathcal{G}^c}(y)$', 'ProjGM', 'Keshavaraz', 'LMI'};
dist_arr = [mean(d_bf), mean(d_qcqp), mean(d_lsq), mean(d_sdp)];
xtick_arr = 1 : length(dist_arr);

% Sort in ascending manner
[dist_arr, sortind] = sort(dist_arr, 'ascend');
dist_bars = dist_bars(sortind);
cmap = cmap(sortind, :);

figure;
h_bar = bar(xtick_arr, dist_arr, 'FaceColor', 'Flat', 'CData', cmap);
xticks(xtick_arr);
xticklabels(dist_bars);
set(gca, 'TickLabelInterpreter', 'latex');
title({sprintf('Mean distance of randomly generated points \n to set of global optima $\\mathcal{G}$ estimated by different standard methods')}, 'interpreter', 'latex');