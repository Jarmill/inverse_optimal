load('opt_5_2d.mat')

y = [2; 2];

%define constraints
n = length(x_star{1});
m = length(x_star);

Aeq= ones(1, m);
beq = [1];

A = -speye(m);
b = sparse(m, 1);

%initial guess
alpha0 = ones(m, 1)/m;

cost_fun = @(alpha) cost_uncons_alpha(alpha, Q, x_star, y);
grad_fun = @(alpha) hproj_egrad(alpha, Q, x_star, y, 0);
grad_proj_fun = @(alpha) hproj_egrad(alpha, Q, x_star, y, 0);
grad_fun_s = @(alpha, s) hproj_egrad(alpha, Q, x_star, y, s);
%set options
options = optimset;

CUSTOM_GRAD = false;

if CUSTOM_GRAD
    %custom gradient debugging from hproj_egrad
    
end

%run
[alpha_con, fval, exitflag, output, lambda, grad, hessian]...
    = fmincon(cost_fun, alpha0, A, b, Aeq, beq, [], [], [], options);


x_opt = x_opt_uncons(alpha_con, Q, x_star);

dist0 = sqrt(2*cost_fun(alpha0));
dist = sqrt(2*fval);

grad_con = grad_fun(alpha_con);
norm(grad_con - grad)