function [x, alpha, d, varargout] = QCQP_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, z0)
%QCQP_IO_CONSTRAINED performs the QCQP formulation of the inverse 
%optimization on a problem with quadratic cost functions and linear
%constraints.
%
%   min_x    1/2*x'*Q*x + phi'*x
%   s.t.     xl <= x <= xu
%            A*x == b
%            C*x <= d
%
%   [x, alpha, d] = QCQP_IO_CONSTRAINED(Q, phi, A, b, C, d, xl, xu, y, z0)
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
%   [x0; alpha0; lambda0; mu0];
%   
%   Outputs:
%   x ~ point from the global optima set that is closest to y
%   alpha ~ corresponding parametrization of the global optima closest to y
%   d ~ distance from y to the point x

% Optimization over x and alpha concatenated: z = (x, alpha)
n = numel(y);   % Dimensionality of variable
m = numel(Q);   % Dimensionality of parametrization
p = size(A, 1);     % Dimensionality of equality dual variables
q = size(C, 1) + length(xl) + length(xu);   % Dimensionality of inequality dual variables

% Create the equality matrices (Direct optimization constraints,
% sum of cost function parametrization is 1)
Aeq_do = zeros(p, n+m+p+q);
Aeq_do(1:p, 1:n) = A;
beq_do = b;
Aeq_io = [zeros(1, n), ones(1, m), zeros(1, p), zeros(1, q)];
beq_io = 1;
Aeq_qcqp = [Aeq_do; Aeq_io];
beq_qcqp = [beq_do; beq_io];
% Create the inequality matrices (Direct optimization constraints, cost 
% functions parametrization and inequality lagrangian multipliers are 
% bigger than 0)
A_do = zeros(q, n+m+p+q);
A_do(1:size(C,1), 1:n) = C;
A_do(size(C,1)+1 : size(C,1)+length(xl), 1:n) = -eye(n);
A_do(size(C,1)+length(xl)+1 : size(C,1)+length(xl)+length(xu), 1:n) = eye(n);
b_do = [d; -xl; xu];
A_io = zeros(m+q, n+m+p+q);
A_io(n+1:n+m, n+1:n+m) = -eye(m);          % Cost function parametrization
A_io(end-q+1:end, end-q+1:end) = -eye(q);  % Inequality Lagrangian Multipliers
b_io = zeros(m+q, 1);
A_qcqp = [A_do; A_io];
b_qcqp = [b_do; b_io];

% Before hessian
% % Create the options
% qcqp_op = optimoptions(@fmincon,...
%                        'Display', 'Iter-Detailed',...
%                        'Algorithm', 'active-set',...
%                        'StepTolerance', 1e-15,...
%                        'OptimalityTolerance', 1e-10, ...
%                        'PlotFcn', {'optimplotx', 'optimplotfval', 'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotstepsize', 'optimplotfirstorderopt'}, ...
%                        'SpecifyObjectiveGradient', true,...
%                        'SpecifyConstraintGradient', true,...
%                        'FiniteDifferenceStepSize', 1e-5,...
%                        'CheckGradients', false);
% 
% % Call the optimization process
% [z_opt, d2] = fmincon(@(z)QCQP_COST(z,y,n,m,p,q), z0, A_qcqp, b_qcqp, Aeq_qcqp, beq_qcqp, [], [], @(z)QCQP_CONSTRAINT(z,Q,phi,A,b,C,d,xl,xu,n,m,p,q), qcqp_op);

% % WITH HESSIAN
% % Create the options
% qcqp_op = optimoptions(@fmincon,...
%                        'Display', 'Iter-Detailed',...
%                        'Algorithm', 'interior-point',...
%                        'StepTolerance', 1e-15,...
%                        'PlotFcn', {'optimplotx', 'optimplotfval', 'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotstepsize', 'optimplotfirstorderopt'}, ...
%                        'OptimalityTolerance', 1e-10, ...
%                        'SpecifyObjectiveGradient', true,...
%                        'SpecifyConstraintGradient', true,...
%                        'HessianFcn', @(~, lagr_mult) QCQP_HESSIAN(Q,C,n,m,p,q,lagr_mult), ...
%                        'CheckGradients', false);
% Create the options
qcqp_op = optimoptions(@fmincon,...
                       'Display', 'None',...
                       'Algorithm', 'interior-point',...
                       'StepTolerance', 1e-15,...
                       'OptimalityTolerance', 1e-10, ...
                       'SpecifyObjectiveGradient', true,...
                       'SpecifyConstraintGradient', true,...
                       'HessianFcn', @(~, lagr_mult) QCQP_HESSIAN(Q,C,n,m,p,q,lagr_mult), ...
                       'CheckGradients', false);

% Call the optimization process
[z_opt, d2, ef, ~] = fmincon(@(z)QCQP_COST(z,y,n,m,p,q), z0, A_qcqp, b_qcqp, Aeq_qcqp, beq_qcqp, [], [], @(z)QCQP_CONSTRAINT(z,Q,phi,A,b,C,d,xl,xu,n,m,p,q), qcqp_op);

% Extract result
x = z_opt(1:n);
alpha = z_opt(n+1 : n+m);
d = sqrt(d2);

% Out
if nargout > 3
    lambda = z_opt(n+m+1:n+m+p);
    varargout{1} = lambda;
    if nargout > 4
        mu = z_opt(n+m+p+1:n+m+p+size(C, 1));
        varargout{2} = mu;
        if nargout > 5
            nu_minus = z_opt(n+m+p+size(C, 1)+1 : n+m+p+size(C, 1)+length(xl));
            varargout{3} = nu_minus;
            if nargout > 6
                nu_plus = z_opt(n+m+p+size(C, 1)+length(xl)+1 : n+m+p+q);
                varargout{4} = nu_plus;
            end
        end
    end
end

end

function [f, df] = QCQP_COST(z,y,n,m,p,q)
%QCQP_COST implements the cost function of the QCQP inverse optimization
%formulation.

% Separate into x, alpha, lambda, mu, 
x = z(1:n);
% alpha = z(n+1 : n+m);
% lambda = z(n+m+1:n+m+p);
% mu = z(n+m+p+1:n+m+p+q);

% Cost function = || y-x ||^2
f = (y - x).' * (y - x);

% Gradient = 2 * (x - y)
df = zeros(1, n+m+p+q);
df(1:n) = 2 * (x-y).';
end

function [c, ceq, dc, dceq] = QCQP_CONSTRAINT(z,Q,phi,A,b,C,d,xl,xu,n,m,p,q)
%QCQP_CONSTRAINT implements the constraint function of the QCQP inverse optimization
%formulation.

% Separate into x, alpha, lambda, mu
x = z(1:n);
alpha = z(n+1 : n+m);
lambda = z(n+m+1:n+m+p);
mu = z(n+m+p+1:n+m+p+q);

% Separate mu into nu_minus, nu+plus and mu
nu_minus = mu(size(C,1)+1 : size(C,1)+length(xl));
nu_plus = mu(size(C,1)+length(xl)+1 : size(C,1)+length(xl)+length(xu));
mu = mu(1 : size(C,1));

% Nonlinear inequalities none
c = [];
dc = [];

% Nonlinear equality constraint: Stationarity Condition
%(sum alphai Qi)*x + (sum alphai phii) + A'*lambda + C'*mu - nu_minus + nu_plus == 0
Qmix = 0;
phimix = 0;
dAlpha = zeros(n, m);
for ii = 1 : m
    Qmix = Qmix + alpha(ii) * Q{ii};
    phimix = phimix + alpha(ii) * phi{ii};
    dAlpha(:, ii) =  Q{ii} * x + phi{ii};
end

% If neither A nor C are empty
if ~isempty(A) && ~isempty(C)
    cstatcond = Qmix*x + phimix + A.' * lambda + C.' * mu - nu_minus + nu_plus;
    dcstatcond = [Qmix, dAlpha, A.', C.', -eye(length(xl)), eye(length(xu))];    
% If both are empty
elseif isempty(A) && isempty(C)
    cstatcond = Qmix*x + phimix - nu_minus + nu_plus;
    dcstatcond = [Qmix, dAlpha, -eye(length(xl)), eye(length(xu))];
% If only A is empty
elseif isempty(A)
    cstatcond = Qmix*x + phimix + C.' * mu - nu_minus + nu_plus;
    dcstatcond = [Qmix, dAlpha, C.', -eye(length(xl)), eye(length(xu))];
% If only C is empty
elseif isempty(C)
    cstatcond = Qmix*x + phimix + A.' * lambda - nu_minus + nu_plus;
    dcstatcond = [Qmix, dAlpha, A.', -eye(length(xl)), eye(length(xu))];
end
    
% Nonlinear equality constraint: Complementarity Slackness for General
% Inequalities
% mu .* (C*x - d) == 0

% If C is nonempty
if ~isempty(C)
    ccompslackgen = mu .* (C*x - d);
    dccompslackgen = [diag(mu)*C, zeros(size(C,1), m+p), diag(C*x-d), zeros(size(C,1), length(xl)+length(xu))];
else
    ccompslackgen = [];
    dccompslackgen = [];
end

% Nonlinear equality constraint: Complementarity Slackness for Lower Bounds
% nu_minus.* (-x + xl) == 0
ccompslacklb = nu_minus .* (-x+xl);
dccompslacklb = [-diag(nu_minus), zeros(length(xl), m+p+size(C, 1)), diag(-x+xl), zeros(length(xl), length(xu))];

% Nonlinear equality constraint: Complementarity Slackness for Upper Bounds
% nu_plus.* (x + xu) == 0
ccompslackub = nu_plus .* (x-xu);
dccompslackub = [diag(nu_plus), zeros(length(xu), m+p+size(C, 1)), zeros(length(xu), length(xl)), diag(x-xu)];

ceq = [cstatcond; ccompslackgen; ccompslacklb; ccompslackub];
dceq = [dcstatcond; dccompslackgen; dccompslacklb; dccompslackub].';
% ceq = [cstatcond];
% dceq = [dcstatcond].';
end

function H = QCQP_HESSIAN(Q,C,n,m,p,q,lagr_mult)
%QCQP_HESSIAN calculates the QCQP problem hessian.


% To avoid storing too much individual hessians, just add them to a common
% sum as soon as you generate them

% Initialize with the cost function hessian, as it's 2*eye(n) for the DO variable
% part.
H = zeros(n+m+p+q, n+m+p+q);
H(1:n, 1:n) = 2*eye(n);

% Stationarity equality constraint ( there are n of them )
%(sum alphai Qi)*x + (sum alphai phii) + A'*lambda + C'*mu - nu_minus + nu_plus == 0
% Calculate the hessian of each of them and add them to the sum
% Extract lagrangian multipliers
lagr_mult_stat = lagr_mult.eqnonlin(1:n);
for ii = 1 : n
    % Stationarity hessian is only w.r.t. x and alpha
    Hx = zeros(n+m+p+q, n);
    Halpha = zeros(n+m+p+q, m);
    for jj = 1 : m
        Hx(n+jj, :) = Q{jj}(ii, :);
        Halpha(1:n, jj) = Q{jj}(:, ii);
    end
    % Do not create full matrix but add only the part that's calculated,
    % multiplied by the lagrangian multipliers
    H(:, 1:n) = H(:, 1:n) + lagr_mult_stat(ii) * Hx;
    H(:, n+1:n+m) = H(:, n+1:n+m) + lagr_mult_stat(ii) * Halpha;
end

% General complementarity slackness equality constraint ( there are size(C, 1) of them )
% mu .* (C*x - d) == 0
% Calculate the hessian of each of them and add them to the sum
% Extract lagrangian multipliers
lagr_mult_compgen = lagr_mult.eqnonlin(n+1:n+size(C,1));
for ii = 1 : size(C, 1)
    % General complementarity hessian is only w.r.t. x and mu
    Hx = zeros(n+m+p+q, n);
    Hmu = zeros(n+m+p+q, size(C, 1));
    % Only one row and one column
    Hx(n+m+p+ii, :) = C(ii, :);
    Hmu(1:n, ii) = C(ii, :).';
    % Do not create full matrix but add only the part that's calculated,
    % multiplied by the lagrangian multipliers
    H(:, 1:n) = H(:, 1:n) + lagr_mult_compgen(ii) * Hx;
    H(:, n+m+p+1:n+m+p+size(C, 1)) = H(:, n+m+p+1:n+m+p+size(C, 1)) + lagr_mult_compgen(ii) * Hmu;
end

% Complementarity Slackness for Lower Bounds
% nu_minus.* (-x + xl) == 0
% Calculate the hessian of each of them and add them to the sum
% Extract lagrangian multipliers
lagr_mult_complb = lagr_mult.eqnonlin(n+size(C,1)+1:n+size(C,1)+n);
for ii = 1 : n
    % General complementarity hessian is only w.r.t. x and nu_minus
    Hx = zeros(n+m+p+q, n);
    Hnu_minus = zeros(n+m+p+q, n);
    % Only one row and one column
    Hx(n+m+p+size(C, 1)+ii, ii) = -1;
    Hnu_minus(ii, ii) = -1;
    % Do not create full matrix but add only the part that's calculated,
    % multiplied by the lagrangian multipliers
    H(:, 1:n) = H(:, 1:n) + lagr_mult_complb(ii) * Hx;
    H(:, n+m+p+size(C, 1)+1:n+m+p+size(C, 1)+n) = H(:, n+m+p+size(C, 1)+1:n+m+p+size(C, 1)+n) + lagr_mult_complb(ii) * Hnu_minus;
end

% Complementarity Slackness for Upper Bounds
% nu_plus.* (x - xu) == 0
% Calculate the hessian of each of them and add them to the sum
% Extract lagrangian multipliers
lagr_mult_compub = lagr_mult.eqnonlin(n+size(C,1)+n+1:n+size(C,1)+2*n);
for ii = 1 : n
    % General complementarity hessian is only w.r.t. x and nu_minus
    Hx = zeros(n+m+p+q, n);
    Hnu_plus = zeros(n+m+p+q, n);
    % Only one row and one column
    Hx(n+m+p+size(C, 1)+n+ii, ii) = 1;
    Hnu_plus(ii, ii) = 1;
    % Do not create full matrix but add only the part that's calculated,
    % multiplied by the lagrangian multipliers
    H(:, 1:n) = H(:, 1:n) + lagr_mult_compub(ii) * Hx;
    H(:, n+m+p+size(C, 1)+1:n+m+p+size(C, 1)+n) = H(:, n+m+p+size(C, 1)+1:n+m+p+size(C, 1)+n) + lagr_mult_compub(ii) * Hnu_plus;
end
end