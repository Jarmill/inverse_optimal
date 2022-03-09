load('opt_3_2d.mat')
% load('opt_5_2d.mat')

m = length(x_star);

%% test if a point is inside the feasible set
% this is a linear programming problem

y_in  = [-1.5; 1];
y_out = [2; 2];
% y_test = [1.3086; 1.0162];
% y_test = [1.3; 1.0];
% y_test = [-1.2765; 0.2654];
% y_test = [1.76223692633813;0.679691280061691];

alpha_in_1  = feas_lp_uncons(y_in,  Q, x_star);
alpha_in_2 = feas_lp_uncons(y_in, Q, x_star, [zeros(m-2); 2;0]);

x_opt_1 = x_opt_uncons(alpha_in_1, Q, x_star);
x_opt_2 = x_opt_uncons(alpha_in_2, Q, x_star);

alpha_out = feas_lp_uncons(y_out, Q, x_star);