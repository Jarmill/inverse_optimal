clear all;
close all;
clc;

% Define dimension
n = 3;
% Define number of cost functions
nf = 3;

% Seed the random processes
rng(312, 'twister');

% The MOSEK solver directory
mosek_dir = input_your_mosek_directory;
% The yalmip directory
yalmip_dir = input_your_yalmip_directory;


% Define center of shapes
xc = [3;3;3];
% xc = zeros(2, 1);
r = 2; % "radius" of shapes

% Design Q's
V1 = eye(n);
D1 = diag([2, 0.8, 1.5]);
D1 = eye(n);
V2 = eye(n);
D2 = diag([1, 3, 1.3]);
V3 = eye(n);
D3 = diag([0.5, 1.6, 2]);

Q{1} = V1 * D1 * V1.';
Q{2} = V2 * D2 * V2.';
Q{3} = V3 * D3 * V3.';

% Design true Q's
% Q_true{1} = 

% Design x_stars
% x_star{1} = 2*[xc(1); 0];
x_star{1} = [0; 0; 0];
x_star{2} = [r; 0; 0];
x_star{3} = 2*[0; r; 0];

% Calculate phis
phi{1} = -Q{1} * x_star{1};
phi{2} = -Q{2} * x_star{2};
phi{3} = -Q{3} * x_star{3};

% Define lower bounds
xl = xc - r;

% Define upper bounds
xu = xc + r;
xu(3) = 3;

% 4 inequalities: 3*x1 + 3*x2 <= 2*r && 3*x1 - 3*x2 <= 2*r && 
%                 -3*x1 + 3*x2 <= 2*r && -3*x1 - 3*x2 <= 2*r
% Make them centered at xc: Cx <= d  with variable transf x = x-xc
%   Cx <= d + C xc
C = ones(8, 3);
C([5 6 7 8], 1) = -1;
C([3 4 7 8], 2) = -1;
C([2 4 6 8], 3) = -1;
C = C([[2, 4, 6, 8]], :);
% d = r*ones(8, 1);
d = r*ones(4, 1);
d = d + C*xc;

% No equalities
A = [];
b = [];

%% Generate random points y from the shape defined by constraints
% Number of points to generate
Nrand = 100;
% % Parametrization of polytope vertices
% beta = rand(nv, Nrand); % between 0 and 1
% beta = beta ./ sum(beta); % all normalized to sum to 1
% % Generate random points from the polytope
% y = V * beta;

% Choose point
% y0 = [6; 6;];
% y0 = [4; 4;];
y0 = [3.5; 3; 2];
y = [y0, randn(n, Nrand-1)];

% Prealocate inverse optimization for each one of them
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
    x_qcqp(:, ii) = x_fmincon;
    
%     % Redo DOC to get x
%     Qmix = 0;
%     phimix = 0;
%     for cf = 1 : nf
%         Qmix = Qmix + alpha_con(cf) * Q{cf};
%         phimix = phimix + alpha_con(cf) * phi{cf};
%     end
%     
%     % 
%     [x_qcqp(:, ii), ~, ~, ~, mult0] = quadprog(Qmix, phimix, C, d, A, b, xl, xu);
end

% Remove mosek and yalmip
rmpath(mosek_dir);
rmpath(yalmip_dir);

%% Plot
fig = figure;
fig.Position(1:2) = zeros(1, 2);
fig.Position(3:4) = [780, 780];
hold all

% % Plot unconstrained minima
% x_star_uncons_sz = 200;
% % x_star_uncons_col = winter(length(x_star));
% x_star_uncons_col = zeros(length(x_star), 3);
% h_uncons = plot_x_star(x_star, x_star_uncons_sz, x_star_uncons_col);
% h_uncons.DisplayName = sprintf('$x^*_{1:%d} \\in \\mathcal{R}^3$', nf);
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
m_grid = 10000;
h_combined_sol = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, m_grid, n_mesh, false, inf);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G}^c =  \\{x \\in X | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
% h_combined_sol.FaceAlpha = 0.7;
% h_combined_sol.FaceColor = zeros(3, 1);


% Plot minima
x_star_cons_sz = 350;
% x_star_cons_col = copper(nf);
x_star_cons_col = copper(3);
x_star_cons_col = repmat(x_star_cons_col(3, :), nf, 1);
h_star = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, x_star_cons_sz, x_star_cons_col);
h_star(round(length(x_star)/2)).DisplayName = sprintf('$x^{f,c}_{1:%d} \\in X$', nf);

% % Plot level sets
% x_star_cons = mat2cell([ [h_star.XData]; [h_star.YData]; [h_star.ZData] ], 3, ones(1, length(h_star)));
% lev = evaluate_quadratics(Q, x_star, x_star_cons);  % Make the ellipses tangential
% n_ell = 1000;
% collevsets = winter(nf);
% collevsets = repmat(collevsets(1, :), nf, 1);
% h_level_sets = plot_level_sets(Q, x_star, n_ell, lev, collevsets);
% for ii = 1 : length(h_level_sets)
%     h_level_sets(ii).HandleVisibility = 'Off';
% %     h_level_sets(ii).FaceAlpha = 0.1;
% end
% % h_lev(round(length(h_lev)/3)).HandleVisibility = 'On';
% % h_lev(round(length(h_lev)/3)).DisplayName = sprintf('$\\{ x_{1:%d} | f_{1:%d}(x) = %s \\}$', nf, nf, sprintf(['[ ', repmat('%.2f, ', 1, nf-1), '%.2f]'], lev) );


% Plot optimal points
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
sz = 250;% Choose a random point out of the Nrand ones
iir = 1;


% Plot random index point and the projection lines
if n == 2

% Random index point
h_query = scatter(y(1, iir), y(2, iir), sz, cmap(1, :), 'o', 'DisplayName', 'y', 'LineWidth', 2);
h_bf = scatter(x_bf(1, iir), x_bf(2, iir), sz, cmap(2, :), 'o', 'HandleVisibility', 'Off', 'LineWidth', 2);
h_qcqp = scatter(x_qcqp(1, iir), x_qcqp(2, iir), sz, cmap(3, :), '+', 'HandleVisibility', 'Off', 'LineWidth', 2);
h_lsq = scatter(x_lsq(1, iir), x_lsq(2, iir), sz, cmap(4, :), 's', 'HandleVisibility', 'Off', 'LineWidth', 2);
h_sdp = scatter(x_sdp(1, iir), x_sdp(2, iir), sz, cmap(5, :), 'x', 'HandleVisibility', 'Off', 'LineWidth', 2);

% Plot connecting line
h_line_bf = plot([y(1, iir), x_bf(1, iir)], [y(2, iir), x_bf(2, iir)], '-', 'Color', cmap(2, :),  'DisplayName', '$\textrm{Proj}_{\mathcal{G}^c}(y)$', 'LineWidth', 2);
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
h_line_bf = plot3([y(1, iir), x_bf(1, iir)], [y(2, iir), x_bf(2, iir)], [y(3, iir), x_bf(3, iir)], '-', 'Color', cmap(2, :), 'DisplayName', '$\textrm{Proj}_{\mathcal{G}}(y)$', 'LineWidth', 2);
h_line_qcqp = plot3([y(1, iir), x_qcqp(1, iir)], [y(2, iir), x_qcqp(2, iir)], [y(3, iir), x_qcqp(3, iir)], '--', 'Color', cmap(3, :), 'DisplayName', 'ProjGM', 'LineWidth', 2);
h_line_lsq = plot3([y(1, iir), x_lsq(1, iir)], [y(2, iir), x_lsq(2, iir)], [y(3, iir), x_lsq(3, iir)], ':', 'Color', cmap(4, :), 'DisplayName', 'Keshavaraz', 'LineWidth', 2);
h_line_sdp = plot3([y(1, iir), x_sdp(1, iir)], [y(2, iir), x_sdp(2, iir)], [y(3, iir), x_sdp(3, iir)], '-.', 'Color', cmap(5, :), 'DisplayName', 'LMI', 'LineWidth', 2);

end


% Esthetics
if n == 3
    view(3)
    view(76, 4)
elseif n == 2
    view(2)
end
xlabel('$X_1$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
ylabel('$X_2$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
zlabel('$X_3$-axis', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'LaTeX');
title(sprintf('Projection of a random point onto constrained set of global optima $\\mathcal{G}^c$'), 'Interpreter', 'LaTeX', 'FontSize', 15);
legend('Interpreter', 'LaTeX', 'Location', 'southoutside', 'orientation', 'vertical', 'FontSize', 15, 'numcolumns', 5);
grid;

% Axis range
% off = .05;
% xlim([xl(1)-off, xu(1)+off])
% ylim([xl(2)-off, xu(2)+off])
% xlim([-0.25, 12.9])
% ylim([-0.25, 12.2])
axis square
exportgraphics(gca, '../../img/projection_on_cons.png', 'ContentType', 'Image', 'Resolution', 400);


 
%% Plot distances
% Remove test point from colormap
cmap = cmap(2:end, :);
% Sort different methods
dist_bars = {'$\textrm{Proj}_{\mathcal{G}}(y)$', 'ProjGM', 'Keshavaraz', 'LMI'};
dist_arr = [mean(d_bf), mean(d_qcqp), mean(d_lsq), mean(d_sdp)];
dist_mat = [d_sdp; d_qcqp; d_lsq; d_bf];

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
