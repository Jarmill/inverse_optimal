load('opt_5_2d.mat');

y = [2; 1];
[out_rec, info] = manopt_search_uncons(Q, x_star, y);