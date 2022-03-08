function grad_perp = hproj_egrad(alpha, Q, x_star, y, perp_proj)

if nargin == 4
    perp_proj = true;
end

%derived from http://www.matrixcalculus.org/
%norm2(inv(a1*Q1+a2*Q2+a3*Q3)*(a1*x1+a2*x2+a3*x3)-y)^2
% alpha = z.*z;

%to determine the hessian, use matrixcalculus.org on:
%2*(inv(a1*Q1+a2*Q2+a3*Q3)*(a1*x1+a2*x2+a3*x3)-y)'*inv(a1*Q1+a2*Q2+a3*Q3)*x1-2*(inv(a1*Q1+a2*Q2+a3*Q3)*(a1*x1+a2*x2+a3*x3)-y)'*inv(a1*Q1+a2*Q2+a3*Q3)*Q1*inv(a1*Q1+a2*Q2+a3*Q3)*(a1*x1+a2*x2+a3*x3) 

%% compute the gradient

% x_star_mat = cell2mat(x_star);
n = length(x_star{1});

m = length(alpha);
Q_accum = 0;
x_accum = 0;
for j = 1:m
    Q_accum = Q_accum +  alpha(j) * Q{j};
    x_accum = x_accum + alpha(j) * x_star{j};
end

t2 = x_accum - y;

grad = zeros(m, 1);

for j = 1:m
    x_inv = Q_accum \ x_star{j};
    Q_inv_term = Q_accum \ (Q{j}*x_inv);
    
    grad(j) = 2*t2'*x_inv - 2*t2'*Q_inv_term;
end

%% perform the orthogonal projection away from linear constraints
if perp_proj
%TODO: validate this, fix in case m <= n+1 (?)
Qx_star = zeros(n+1, m);
x_opt = x_opt_uncons(alpha, Q, x_star);
%find orthogonal projector to sum_j alpha_j Q_j(xstar-x0_j) = 0
for j = 1:m
    Qx_star(1:(end-1), j) = Q{j}*(x_opt - x_star{j});
end
Qx_star(end, :) = 1;

% proj_matrix = Qx_star'*((Qx_star*Qx_star') \ Qx_star);
grad_proj = Qx_star'*((Qx_star*Qx_star') \ (Qx_star * grad));
grad_perp = grad - grad_proj;
else
    grad_perp = grad;
end
