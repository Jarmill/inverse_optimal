load('test3_data.mat');

%% test out sdp 

[d_sdp, info_sdp] = sdp_cons_lin_solve(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);
x_sdp = info_sdp.x_rec;
alpha_sdp = info_sdp.alpha_rec; 

%% Use manopt
[out_rec, info] = manopt_search_cons(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);
x_mano = out_rec.x;
d_mano = out_rec.dist;
alpha_mano = out_rec.alpha;

%% test out feasible point solver
[alpha_feas, lambda_feas, mu_feas, exitflag] = ...
    feas_lp_cons(x_mano, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);

%% define optimizer objects and debug
[P_func,P_selector] = cons_local_optimizers(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);

[dist_rec, func_min, x_rec] = cost_cons_alpha(alpha_mano, P_func, P_selector, y)