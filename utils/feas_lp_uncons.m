function [alpha, exitflag] = feas_lp_uncons(y, Q, x_star, f)
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
n = length(x_star{1});
m = length(x_star);

if nargin < 4
    f = sparse(m, 1);
end

Aeq = zeros(n+1, m);
%find orthogonal projector to sum_j alpha_j Q_j(xstar-x0_j) = 0
for j = 1:m
    Aeq(1:(end-1), j) = Q{j}*(y - x_star{j});
end
Aeq(end, :) = 1;

beq = [zeros(n, 1); 1];

A = -speye(m);
b = sparse(m, 1);

% options = optimset('linprog');
% options.Display = 'on';

options = optimoptions('linprog');
% options.linprog.algorithm = 'interior-point';
[alpha, fval, exitflag, output] = linprog(f, A, b, Aeq, beq, [], [],options);

end

