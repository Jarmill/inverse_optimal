rng(33, 'twister')
load('opt_5_2d.mat')


alpha = [0.2; 0.1; 0.3; 0.15; 0.25];
z = sqrt(alpha);
y_out = [2; 1];


x_opt = x_opt_uncons(alpha, Q, x_star);

grad = hproj_egrad(alpha, Q, x_star, y_out);
gz = hproj_egrad_hadamard(z, Q, x_star, y_out);
