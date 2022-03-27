%run the unconstrained lower bound code

load('opt_5_2d.mat')

% SOLVE = 1;

y_out = [2; 1];

 [dist_rec, dist_info] = sdp_uncons_solve(y_out, Q, x_star);