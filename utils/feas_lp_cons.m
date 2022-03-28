function [alpha, lambda, mu, exitflag] = feas_lp_cons(y, Q, phi, A, b, Aeq, beq)
%FEAS_UNCONS: is the point y a global optimum of (Q, x_star)?
% does there exist simplex-weights alpha such that:
% y = x_opt_uncons(alpha, Q, x_star)?
%
%INPUTS (may change around)
%   y:      test point outside the GO set
%   Q:      Dictionary of hessians
%   x_star: Corresponding optimal points to quadratic functions
%   f:      (optional) cost function of LP when alpha is nonunique
%Output:
n = size(Q{1}, 1);
m = length(Q);
n_ineq = length(b);
n_eq = length(beq);

alpha = sdpvar(m, 1);
lambda = sdpvar(n_eq, 1);
mu = sdpvar(n_ineq, 1);


% stationarity = Q*y + phi
stat = zeros(n, 1, 'like', sdpvar);
for j = 1:m
    stat = stat + alpha(j)*Q{j}*y + alpha(j)*phi{j};
end

if n_eq
    stat = stat + Aeq'*lambda;
end

if n_ineq
    stat = stat + A'*mu;
end

cons = [(sum(alpha)==1):'simplex sum'; (alpha >= 0):'simplex nonneg'; (mu >= 0):'dual nonneg'; (sum(mu'*(A*y - b)) == 0):'comp slack'; (stat==0):'stationarity'];

opts = sdpsettings('solver', 'linprog');
sol = optimize(cons, [], opts);

exitflag = sol.problem;
if sol.problem == 0
    alpha = value(alpha);
    lambda = value(lambda);
    mu = value(mu);
else
    alpha = [];
    lambda = [];
    mu = [];
end

% if nargin < 4
%     f = sparse(m, 1);
% end

% Aeq = zeros(n+1, m);
% %find orthogonal projector to sum_j alpha_j Q_j(xstar-x0_j) = 0
% for j = 1:m
%     Aeq(1:(end-1), j) = Q{j}*(y - x_star{j});
% end
% Aeq(end, :) = 1;
% 
% beq = [zeros(n, 1); 1];
% 
% A = -speye(m);
% b = sparse(m, 1);
% 
% % options = optimset('linprog');
% % options.Display = 'on';
% 
% options = optimoptions('linprog');
% % options.linprog.algorithm = 'interior-point';
% [alpha, fval, exitflag, output] = linprog(f, A, b, Aeq, beq, [], [],options);

end

