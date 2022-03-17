load("sdp_opt_5.mat")
load('opt_5_2d.mat')
n = length(x_star{1});
m = length(x_star);

% opt_round = struct('n', n, 'm', m, 'x_opt', @(alpha) x_opt_uncons(alpha, Q, x_star));

y = [2; 1];

% [alpha_rec, x_rec] = round_alpha_matrix(Z, opt_round);

rrPar = struct;
rrPar.Q = Q;
rrPar.x_star = x_star;
rrPar.y = y;
% rrPar.m = m;


% try rounding alone
[Xmat_round, pobjround_round, info_round] = local_search_uncons(Z, [], rrPar, [], true);


% % try calling manopt
[Xmat, pobjround, info] = local_search_uncons(Z, [], rrPar, [], false);