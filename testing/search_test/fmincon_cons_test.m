load('test3_data.mat');

n = length(y);
m = length(phi);



[P_func,P_selector] = cons_local_optimizers(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);


alpha0 = ones(m, 1)/m;
f0 = f_in(alpha0, P_func, y);

%% fmincon command
Aeq= ones(1, m);
beq = [1];

A = -speye(m);
b = sparse(m, 1);

options = optimset;

[alpha_con, fval, exitflag, output, lambda, grad, hessian]...
    = fmincon(@(alpha) f_in(alpha, P_func, y), alpha0, A, b, Aeq, beq, [], [], [], options);


dist_0 = sqrt(f0);
dist_fmincon = sqrt(fval);

%% test out sdp 

[dist_sdp, info_sdp] = sdp_cons_lin_solve(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);
x_sdp = info_sdp.x_rec;
alpha_sdp = info_sdp.alpha_rec; 


function f_out = f_in(alpha, P_func, y)
%cost function to minimize the mixed quadratic
funcdata = P_func(alpha);
% f_out = xdata(1);
f_out = norm(y - funcdata(2:end))^2;

end