function grad_perp = hproj_egrad(alpha, Q, x_star, y, perp_proj)
%evaluate the gradient of the unconstrained GO distance function
if nargin == 4
    perp_proj = false;
end

%derived from http://www.matrixcalculus.org/
%norm2(inv(a1*Q1+a2*Q2+a3*Q3)*(a1*Q1*x1+a2*Q2*x2+a3*Q3*x3)-y)^2
% alpha = z.*z;

%to determine the hessian, use matrixcalculus.org on:
%2⋅t⊤3⋅T0⋅t1−2⋅t⊤3⋅T0⋅Q1⋅t2where
%T0=inv(a1⋅Q1+a2⋅Q2+a3⋅Q3) [inv(Q_accum)]
%t1=Q1⋅x1
%t2=T0⋅(a1⋅t1+a2⋅Q2⋅x2+a3⋅Q3⋅x3)   [x_opt]
%t3=t2−y 


%% compute the gradient

% x_star_mat = cell2mat(x_star);
n = length(x_star{1});

m = length(alpha);
Q_accum = 0;
x_accum = 0;
for j = 1:m
    Q_accum = Q_accum +  alpha(j) * Q{j};
    x_accum = x_accum + alpha(j) * Q{j}* x_star{j};
end

x_opt = Q_accum \ x_accum;
t3 = x_opt - y;

grad = zeros(m, 1);

for j = 1:m
    x_inv_j = Q_accum \ (Q{j}*x_star{j});
    Q_inv_term = Q_accum \ (Q{j}*x_opt);
    
    grad(j) = 2*t3'*(x_inv_j -Q_inv_term);
end

%% perform the orthogonal projection away from linear constraints

%this section may be bugged and inaccurate.
%should probably drop.
%the set of alphas describing a valid point is a polytope

if perp_proj == 3
    grad_perp = grad - sum(grad)/m;
elseif perp_proj == 2
    %TODO: validate this, fix in case m <= n+1 (?)
    Qx_star = zeros(n+1, m);
    % x_opt = x_opt_uncons(alpha, Q, x_star);
    %find orthogonal projector to sum_j alpha_j Q_j(xstar-x0_j) = 0
    for j = 1:m
        Qx_star(1:(end-1), j) = Q{j}*(x_opt - x_star{j});
    end
    Qx_star(end, :) = 1;

    % proj_matrix = Qx_star'*((Qx_star*Qx_star') \ Qx_star);
    grad_proj = Qx_star'*((Qx_star*Qx_star') \ (Qx_star * grad));
    grad_perp = grad - grad_proj;
elseif perp_proj == 1
    %do not perform projection against  simplex constraint
    %TODO: validate this, fix in case m <= n+1 (?)
    Qx_star = zeros(n, m);
    % x_opt = x_opt_uncons(alpha, Q, x_star);
    %find orthogonal projector to sum_j alpha_j Q_j(xstar-x0_j) = 0
    for j = 1:m
        Qx_star(1:(end), j) = Q{j}*(x_opt - x_star{j});
    end

    % proj_matrix = Qx_star'*((Qx_star*Qx_star') \ Qx_star);
    grad_proj = Qx_star'*((Qx_star*Qx_star') \ (Qx_star * grad));
    grad_perp = grad - grad_proj;
else
    grad_perp = grad;
end
