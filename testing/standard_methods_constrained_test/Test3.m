clear all;
close all;
clc;

% Define dimension
n = 3;
% Define number of cost functions
nf = 4;

% The MOSEK solver directory
mosek_dir = input_your_mosek_directory;
% The yalmip directory
yalmip_dir = input_your_yalmip_directory;

% Generate random cost functions
rng(312, 'twister');
[Q, phi, x_star] = generate_random_Q_and_phi(n, nf, [0.5, 1.5], [1, 1.1], [1 2], 0);


% Define lower bounds
% xl = zeros(2, 1);
xl = -1*ones(3, 1);

% Define upper bounds
xu = 1*ones(3, 1);

% % Define x1 + x2 <= 1 & -x1 + x2 <= 0
% C = ones(2, 2);
% C(2, 1) = -1;
% d = [1; 0];

% Define x1 + x2 <= 1 & -x1 + x2 <= 1 & x1 - x2 <= 1 & -x1 - x2 <= 1
C = ones(8, 3);
C([5 6 7 8], 1) = -1;
C([3 4 7 8], 2) = -1;
C([2 4 6 8], 3) = -1;
d = ones(8, 1);

% No equalities
A = [];
b = [];
% % Define x1 = 0.4
% A = [1, 0];
% b = 0.4;


% Generate random test point and do IO
% y = generate_random_x(n, 3, [0, sqrt(2)]);
% y = y{3};
% y = [-1.5; -0.49];
y = [2; -2; 2];
% Guess Z0 by solving random direct problem
alpha0 = rand(nf, 1);
alpha0 = alpha0 ./ sum(alpha0);
Qmix = 0;
phimix = 0;
for ii = 1 : nf
    Qmix = Qmix + alpha0(ii) * Q{ii};
    phimix = phimix + alpha0(ii) * phi{ii};
end
[x0, ~, ~, ~, multipliers] = quadprog(Qmix, phimix, C, d, A, b, xl, xu);

lambda0 = multipliers.eqlin;
mu0 = multipliers.ineqlin;
nu_minus0 = multipliers.lower;
nu_plus0 = multipliers.upper;
z0_lsq = [alpha0; lambda0; mu0; nu_minus0; nu_plus0];
z0_qcqp = [x0; alpha0; lambda0; mu0; nu_minus0; nu_plus0];
% [x_retrieval, alpha_retrieval, d_retrieval] = QCQP_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, z0);
[x_bf, alpha_bf, d_bf] = BF_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, 10000, 20);
[x_lsq, alpha_lsq, d_lsq] = LSQ_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, z0_lsq);
[x_qcqp, alpha_qcqp, d_qcqp, lambda_qcqp, mu_qcqp, nu_minus_qcqp, nu_plus_qcqp] = QCQP_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, z0_qcqp);

% Add mosek and yalmip
addpath(mosek_dir);
addpath(genpath(yalmip_dir));

% Use the SDP relaxation
Aeqsdp = A;
beqsdp = b;
Asdp = [C; -eye(n); eye(n);];
bsdp = [d; -xl; xu;];
[d_sdp, info_sdp] = sdp_cons_lin_solve(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);
x_sdp = info_sdp.x_rec;
alpha_sdp = info_sdp.alpha_rec;
% d_sdp = info_sdp.dist_rec;


% Use manopt
[out_rec, info] = manopt_search_cons(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);
x_mano = out_rec.x;
d_mano = out_rec.dist;
alpha_mano = out_rec.alpha;

% Remove mosek and yalmip
rmpath(mosek_dir);
rmpath(genpath(yalmip_dir));

%% Plot
fig = figure;
fig.Position(1:2) = fig.Position(1:2) / 2;
fig.Position(3:4) = [720, 540];
hold all

% Plot feasible set
h_feas = plot_feasible_set_constrained(A, b, C, d, xl, xu);
h_feas.DisplayName = sprintf('$X \\quad (x \\in X)$');

% % Plot unconstrained minima
% h_star_uncons = plot_x_star(x_star, 150, copper(nf));
% h_star_uncons.DisplayName = sprintf('$x^*_{1:%d \\; \\rm Unconstrained}$', nf);

% Plot the set of all optima
n_mesh = 10;
m_grid = 750;
h_combined_sol = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, m_grid, n_mesh, true, 1);
h_combined_sol.DisplayName = sprintf('$\\mathcal{G} =  \\{x | \\alpha_{%d : %d} \\} $', 1, nf);
h_combined_sol.LineStyle = 'none';
h_combined_sol.FaceAlpha = 0.3;


% Plot minima
h_star = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, 150, copper(nf));

% Plot test point
h_testpoint = scatter3(y(1), y(2), y(3), 150, 'r', 'filled');
h_testpoint.DisplayName = '$y$';
h_bfsol = scatter3(x_bf(1), x_bf(2), x_bf(3), 150, 'g', 'filled');
h_bfsol.DisplayName = '$x_{\rm BF}$';
h_lsqsol = scatter3(x_lsq(1), x_lsq(2), x_lsq(3), 150, 'b', 'filled');
h_lsqsol.DisplayName = '$x_{\rm LSQ}$';
h_lsqsol = scatter3(x_qcqp(1), x_qcqp(2), x_qcqp(3), 150, 'm', 'filled');
h_lsqsol.DisplayName = '$x_{\rm QCQP}$';
h_sdpsol = scatter3(x_sdp(1), x_sdp(2), x_sdp(3), 150, 'c', 'filled');
h_sdpsol.DisplayName = '$x_{\rm SDP}$';
h_manosol = scatter3(x_mano(1), x_mano(2), x_mano(3), 150, 'y', 'filled');
h_manosol.DisplayName = '$x_{\rm MANO}$';

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
% xlim([x_bf(1), y(1)])
% ylim([x_bf(2), y(2)])
axis square