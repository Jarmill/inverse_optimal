function [x, alpha, d] = BF_IO_constrained(Q, phi, A, b, C, d, xl, xu, y, Ngrid, Nmesh)
%BF_IO_CONSTRAINED performs the brute force formulation of the inverse 
%optimization on a problem with quadratic cost functions and linear
%constraints.
%
%   min_x    1/2*x'*Q*x + phi'*x
%   s.t.     xl <= x <= xu
%            A*x == b
%            C*x <= d
%
%   [x, alpha, d] = BF_IO_CONSTRAINED(Q, phi, A, b, C, d, xl, xu, y, Ngrid, Nmesh)
%   
%   Inputs:
%   Q ~ cell array of symmetric positive semidefinite matrices of the
%   individual quadratic cost functions
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
%   Ngrid ~ total number of regularly spaced grid points on the simplex of 
%   cost function parametrizations
%   Nmesh ~ total number of regularly spaced points in the polyhedron of
%   solutions when we have a cost function that has multiple solutions
%   
%   Outputs:
%   x ~ point from the global optima set that is closest to y according to
%   the brute force computation
%   alpha ~ corresponding grid point on the simplex which produces the 
%   global optima closest to y
%   d ~ distance from y to the point x

% Optimization alpha
n = numel(y);   % Dimensionality of variable
m = numel(Q);   % Dimensionality of parametrization

% Determine the number of regular gridpoints to use
ngrid = determine_simplex_grid_partition(m-1, Ngrid);

% Get alphas from the simplex
alpha_list = prob_simplex_ndim(m-1, ngrid);

% Initialize minimum distance
dmin = inf;
alphamin = [];
xmin = [];

% Perform a search over the alphas
for ii = 1 : size(alpha_list, 1)
    
    % Get current alpha
    alpha = alpha_list(ii, :);
    
    % Calculate the mixture of Q's and phi's
    Qmix = 0;
    phimix = 0;
    for cf = 1 : m
        Qmix = Qmix + alpha(cf) * Q{cf};
        phimix = phimix + alpha(cf) * phi{cf};
    end
    
    % Calculate the new set of global optima
    xcurr = global_optima_set_qp(Qmix, phimix, A, b, C, d, xl, xu, Nmesh);
    
    % If a single optima
    if size(xcurr, 2) == 1    
        % Calculate distance to query point
        dcurr = norm(xcurr - y);

        % Check if distance is optimal and if so update
        if dcurr < dmin
            dmin = dcurr;
            alphamin = alpha;
            xmin = xcurr;
        end
    % If mutliple optima
    else
        % Calculate the distance to y for each xcurr
        for jj = 1 : size(xcurr, 2)
            % Calculate distance to query point
            dcurr = norm(xcurr(:, jj) - y);
            
            % Check if distance is optimal and if so, update
            if dcurr < dmin
                dmin = dcurr;
                alphamin = alpha;
                xmin = xcurr;
            end
        end
    end
end

% Return global minima
x = xmin;
alpha = alphamin;
d = dmin;
end