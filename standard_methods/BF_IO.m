function [x, alpha, d] = BF_IO(Q, x_star, y, Ngrid)
%BF_IO performs the brute force formulation of the inverse optimization 
%problem with quadratic functions.
%   [x, alpha, d] = BF_IO(Q, x_star, y, Ngrid)
%   
%   Inputs:
%   Q ~ cell array of symmetric positive definite matrices of the
%   individual quadratics
%   x_star ~ cell array of vectors representing the optima of the
%   individual quadratics
%   y ~ input to the inverse optimization
%   Ngrid ~ total number of regularly spaced grid points on the simplex of 
%   cost function parametrizations
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
    
    % Calculate the mixture of Q's and x's
    Qmix = 0;
    xmix = 0;
    for cf = 1 : m
        Qmix = Qmix + alpha(cf) * Q{cf};
        xmix = xmix + alpha(cf) * Q{cf} * x_star{cf};
    end
    
    % Calculate the new global optima
    xcurr = Qmix \ xmix;
    
    % Calculate distance to query point
    dcurr = norm(xcurr - y);
    
    % Check if distance is optimal and if so update
    if dcurr < dmin
        dmin = dcurr;
        alphamin = alpha;
        xmin = xcurr;
    end
end

% Return global minima
x = xmin;
alpha = alphamin;
d = dmin;
end