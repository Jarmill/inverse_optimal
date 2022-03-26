function [x, alpha, d] = LSQ_IO(Q, x_star, y, alpha0)
%LSQ_IO performs the LSQ formulation of the inverse optimization problem
%with quadratic functions.
%   [x, alpha, d] = LSQ_IO(Q, x_star, y, alpha0)
%   
%   Inputs:
%   Q ~ cell array of symmetric positive definite matrices of the
%   individual quadratics
%   x_star ~ cell array of vectors representing the optima of the
%   individual quadratics
%   y ~ input to the inverse optimization
%   alpha0 ~ initial guess for solution (preferably feasible)
%   
%   Outputs:
%   x ~ point from the global optima set that is closest to y according to
%   the "minimum residual" principle
%   alpha ~ corresponding parametrization of the global optima closest to y
%   according to the "minimum residual" principle
%   d ~ distance from y to the point x

% Optimization alpha
n = numel(y);   % Dimensionality of variable
m = numel(Q);   % Dimensionality of parametrization

% Create the regressor
Qdy = zeros(n, m);
for ii = 1 : m
    Qdy(:, ii) = Q{ii} * (y - x_star{ii});
end

% Create the vector to which we regress
d = zeros(n, 1);

% Create the equality and inequality matrices
Aeq = ones(1, m);
beq = 1;
A = -eye(m);
b = zeros(m, 1);

% Create the options
lsq_op = optimoptions(@lsqlin,...
                       'Display', 'Iter-Detailed',...
                       'Algorithm', 'active-set');

% Call the optimization process
[alpha_opt, ~] = lsqlin(Qdy, d, A, b, Aeq, beq, [], [], alpha0, lsq_op);

% Calculate additional quantities for calculating global minima
sumQ = 0;
sumQx = 0;
for ii = 1 : m
    sumQ = sumQ + alpha_opt(ii) * Q{ii};
    sumQx = sumQx + alpha_opt(ii) * Q{ii} * x_star{ii};
end

% Extract result
x = inv(sumQ) * sumQx;
alpha = alpha_opt;
d = norm(y-x);
end