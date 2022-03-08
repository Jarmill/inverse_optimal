% load('opt_3_2d.mat')
load('opt_5_2d.mat')

%% test if a point is inside the feasible set
% this is a linear programming problem

y_in  = [-1.5; 1];
y_out = [2; 2];
% y_test = [1.3086; 1.0162];
% y_test = [1.3; 1.0];
% y_test = [-1.2765; 0.2654];
% y_test = [1.76223692633813;0.679691280061691];

alpha_in  = feas_lp_uncons(y_in,  Q, x_star);
alpha_out = feas_lp_uncons(y_out, Q, x_star);