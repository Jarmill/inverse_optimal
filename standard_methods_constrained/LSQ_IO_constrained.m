function [x, alpha, d] = LSQ_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, z0)
%LSQ_IO_CONSTRAINED performs the LSQ formulation of the inverse 
%optimization on a problem with quadratic cost functions and linear
%constraints.
%
%   min_x    1/2*x'*Q*x + phi'*x
%   s.t.     xl <= x <= xu
%            A*x == b
%            C*x <= d
%
%   [x, alpha, d] = LSQ_IO(Q, phi, A, b, C, d, xl, xu, y, z0)
%   
%   Inputs:
%   Q ~ cell array of symmetric positive definite matrices of the
%   individual quadratics
%   phi ~ cell array of vectors representing the linear part of the
%   quadratic cost functions
%   A ~ matrix of linear equality constraints of the quadratic
%   programming direct problem
%   b ~ vector of linear equality constraints of the quadratic
%   programming direct problem
%   C ~ matrix of linear inequality constraints of the quadratic
%   programming direct problem
%   d ~ vector of linear inequality constraints of the quadratic
%   programming direct problem
%   xl ~ vector of lower bounds on the variables in the direct problem
%   xu ~ vector of upper bounds on the variables in the direct problem
%   y ~ input to the inverse optimization
%   z0 ~ initial guess for solution (preferably feasible); A concatenation
%   [alpha0; lambda0; mu0];
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
p = size(A, 1);     % Dimensionality of equality dual variables
q = size(C, 1) + length(xl) + length(xu);   % Dimensionality of inequality dual variables


% Create the regressor
% Derivatives of cost functions
dF = zeros(n, m);
for ii = 1 : m
    dF(:, ii) = Q{ii} * y + phi{ii};
end
% Derivatives of equality constraint functions
dH = A.';
% Derivatives of inequality constraint functions
dG = C.';
% Derivatives of bound inequality constraint functions
dGlb = -eye(n);
dGub = eye(n);
% Create the derivative-of-lagrangian part of regressor
dL = [dF, dH, dG, dGlb, dGub];
% Calculate the inequality constraints at input point
if isempty(C)
    G = [];
else
    G = C * y - d;
end
Glb = -y + xl;
Gub = y - xu;
% Create the complementarity-slackness part of regressor
Cs = [zeros(q, m+p), diag([G;Glb;Gub])];

% Whole regressor
R_lsq = [dL; Cs];

% Create the vector to which we regress
d_lsq = zeros(n+q, 1);

% Create the equality matrices (sum of cost function parametrization is 1)
Aeq_lsq = zeros(1, size(R_lsq, 2));
Aeq_lsq(1:m) = 1;
beq_lsq = 1;
% Create inequality matrices (Cost function parametrization and inequality
% lagrangian multipliers are bigger than 0)
A_lsq = zeros(m+q, m+p+q);
A_lsq(1:m, 1:m) = -eye(m);
A_lsq(end-q+1:end, end-q+1:end) = -eye(q);
b_lsq = zeros(m+q, 1);

% % Create the options
% lsq_op = optimoptions(@lsqlin,...
%                        'Display', 'Iter-Detailed',...
%                        'Algorithm', 'active-set');
% Create the options
lsq_op = optimoptions(@lsqlin,...
                       'Display', 'None',...
                       'Algorithm', 'active-set');

% Call the optimization process
[z_opt, ~] = lsqlin(R_lsq, d_lsq, A_lsq, b_lsq, Aeq_lsq, beq_lsq, [], [], z0, lsq_op);

% Extract multipliers
alpha = z_opt(1:m);         % Cost function parametrization
lambda = z_opt(m+1:m+p);    % Equality constraint multipliers
mu = z_opt(m+p+1:m+p+q);    % Inequality constraint multipliers

% Divide inequality multipliers 
nu_minus = mu(size(C, 1)+1 : size(C, 1)+length(xl));                        % Lower bound constraint multipliers  
nu_plus= mu(size(C, 1)+length(xl)+1 : size(C, 1)+length(xl)+length(xu));    % Upper bound constraint multipliers
mu = mu(1:size(C,1));                                                       % General inequality constraint multipliers


% Calculate additional quantities for calculating global minima
Qmix = 0;
phimix = 0;
for ii = 1 : m
    Qmix = Qmix + alpha(ii) * Q{ii};
    phimix = phimix + alpha(ii) * phi{ii};
end

% Redo an lsqlin to find the 
% % Extract result
x = quadprog(Qmix, phimix, C, d, A, b, xl, xu, []);
% alpha = z_opt;
d = norm(y-x);
end