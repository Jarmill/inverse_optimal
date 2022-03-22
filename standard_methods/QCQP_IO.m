function [x, alpha, d] = QCQP_IO(Q, x_star, y, z0)
%QCQP_IO performs the QCQP formulation of the inverse optimization problem
%with quadratic functions.
%   [x, alpha, d] = QCQP_IO(Q, x_star, y, z0)
%   
%   Inputs:
%   Q ~ cell array of symmetric positive definite matrices of the
%   individual quadratics
%   x_star ~ cell array of vectors representing the optima of the
%   individual quadratics
%   y ~ input to the inverse optimization
%   z0 ~ initial guess for solution (preferably feasible) (z0 =
%   [x0;alpha0])
%   
%   Outputs:
%   x ~ point from the global optima set that is closest to y
%   alpha ~ corresponding parametrization of the global optima closest to y
%   d ~ distance from y to the point x

% Optimization over x and alpha concatenated: z = (x, alpha)
n = numel(y);   % Dimensionality of variable
m = numel(Q);   % Dimensionality of parametrization

% Create the equality and inequality matrices
Aeq = [zeros(1, n), ones(1, m)];
beq = 1;
A = [zeros(m, n), -eye(m)];
b = zeros(m, 1);

% Create the options
qcqp_op = optimoptions(@fmincon,...
                       'Display', 'Iter-Detailed',...
                       'Algorithm', 'sqp',...
                       'CheckGradients', true);

% Call the optimization process
[z_opt, d2] = fmincon(@(z)QCQP_COST(z,y,n,m), z0, A, b, Aeq, beq, [], [], @(z)QCQP_CONSTRAINT(z,Q,x_star,n,m), qcqp_op);

% Extract result
x = z_opt(1:n);
alpha = z_opt(n+1 : n+m);
d = sqrt(d2);
end

function [f, df] = QCQP_COST(z,y,n,m)
%QCQP_COST implements the cost function of the QCQP inverse optimization
%formulation.

% Separate into x and alpha
x = z(1:n);
alpha = z(n+1 : n+m);

% Cost function = || y-x ||^2
f = (y - x).' * (y - x);

% Gradient = 2 * (x - y)
df = zeros(1, n+m);
df(1:n) = 2 * (x-y).';
end

function [c, ceq, dc, dceq] = QCQP_CONSTRAINT(z,Q,x_star,n,m)
%QCQP_CONSTRAINT implements the constraint function of the QCQP inverse optimization
%formulation.

% Separate into x and alpha
x = z(1:n);
alpha = z(n+1 : n+m);

% Nonlinear inequalities none
c = [];
dc = [];

% Nonlinear equality constraint (sum alphai Qi)*x - (sum alphai Qi xi)
sumQ = 0;
sumQx = 0;
Qdx = zeros(n, m);
for ii = 1 : m
    sumQ = sumQ + alpha(ii) * Q{ii};
    sumQx = sumQx + alpha(ii) * Q{ii} * x_star{ii};
    Qdx(:, ii) = Q{ii} * (x - x_star{ii});
end
ceq = sumQ * x - sumQx;
dceq = [sumQ, Qdx];
end