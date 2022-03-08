function [alpha] = feas_lp_uncons(y, Q, x_star)
%FEAS_UNCONS: is the point y a global optimum of (Q, x_star)?
% does there exist simplex-weights alpha such that:
% y = x_opt_uncons(alpha, Q, x_star)?
n = length(x_star{1});
m = length(x_star);

Aeq = zeros(n+1, m);
%find orthogonal projector to sum_j alpha_j Q_j(xstar-x0_j) = 0
for j = 1:m
    Aeq(1:(end-1), j) = Q{j}*(y - x_star{j});
end
Aeq(end, :) = 1;

beq = [zeros(n, 1); 1];

A = -speye(m);
b = sparse(m, 1);
f = sparse(m, 1);
options = optimset('linprog');
options.Display = 'off';

[alpha, fval, exitflag, output] = linprog(f, A, b, Aeq, beq, [], [],options);

end

