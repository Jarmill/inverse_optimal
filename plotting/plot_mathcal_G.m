function h  = plot_mathcal_G(Q, x_star, n_grid, varargin)
%PLOT_MATHCAL_G plots the set of all global optima of the convex combination
%of a library of quadratic cost functions, given that library.
%   
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_MATHCAL_G(Q, x_star, npts)
%   h = PLOT_MATHCAL_G(Q, x_star, npts, dots_vs_shape)
%   h = PLOT_MATHCAL_G(Q, x_star, npts, dots_vs_shape, smoothness)
%
%   Inputs:
%   Q ~ cell array of positive semidefinite matrices defining the quadratic
%   term
%   x_star ~ cell array of vectors corresponding to global minima of each
%   quadratic
%   n_grid ~ upper bound on the number of points that should be calculated 
%   to generate the mesh
%   dots_vs_shape ~ flag indicating whether we represent the set using a
%   cloud of dots or an alphashape (false for dots, true for alphashape,
%   default is alphashape).
%   smoothness ~ the scalar alpha parameter of the alphaShape representing 
%   the sets mesh
%   
%   Outputs:
%   h ~ handle of the plotted object

% Get the number of cost functions in the library
nf = length(Q);
% Get the size of the optimization variables
n = length(x_star{1});

% Check the size of input
if n ~= 2 && n ~= 3
    error("Dimensionality of the problem is not 2D nor 3D. Please supply 2D or 3D x_stars.");
end

% Treat inputs:
% First argument is passed
if length(varargin) >= 1
    % Check validity
    if ~islogical(varargin{1}) && ~isscalar(varargin{1})
        error('dots_vs_shape parameter must be of a numeric scalar or logical type.');
    end
    % If it is numeric or logical assign its value
    shape = varargin{1};
% First argument is not passed
else
    % Default value is true (corresponds to alphashape representation)
    shape = true;
end

% Second argument is passed
if length(varargin) >= 2
    % Check validity
    if ~isscalar(varargin{2})
        error('smoothness parameter must be of a numeric scalar type.');
    end
    % If it is logical assign its value
    smoothness = varargin{2};
% First argument is not passed
else
    % Default value is false (doesn't set smoothness)
    smoothness = false;
end

% Initialize the affine dimension of the simplex of parametrizations
n_simplex = nf - 1;
% Define the upper bound on the number of gridpoints that you want
m_tot_grid = n_grid;
% Approximate the number of grid points along each dimension of the simplex
m_grid = floor(nthroot(factorial(n_simplex) * m_tot_grid, n_simplex));

% Get the parametrization
alpha = prob_simplex_ndim(n_simplex, m_grid);

% Compute set of minimums of convex combinations
% Prealocate
xcomb = zeros(n, size(alpha, 1));

% Calculate optimum for different combinations
for cc = 1 : size(alpha, 1)
    
    % Find linear combinations of Q matrices and F matrics
    Qcomb = 0;
    x_star_comb = 0;
    for ii = 1 : nf
        Qcomb = Qcomb + alpha(cc, ii) * Q{ii};
        x_star_comb = x_star_comb + alpha(cc, ii) * Q{ii} * x_star{ii};
    end
    
    % Find optimal solution
    xcomb(:, cc) = Qcomb \ x_star_comb;    
end

% Get the set shape differently for 2D and 3D
if n == 2
    % Differentiate between black dot and shape representation
    % If shape
    if shape    
        % Get the 2D shape differently if smoothness parameter is given
        if smoothness
            set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).', smoothness);
        else
            set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).');
        end
        % Get the shape plot
        h = plot(set_shape);
    % If dots
    else
        % Get the dots
        h = plot(xcomb(1, :), xcomb(2, :), 'k.', 'MarkerSize', 20, 'LineStyle', 'None', 'LineWidth', 2);
    end
    
elseif n == 3
    % If shape
    if shape    
        % Get the 3D shape differently if smoothness parameter is given
        if smoothness
            set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).', xcomb(3, :).', smoothness);
        else
            set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).', xcomb(3, :).');
        end
        % Get the shape plot
        h = plot(set_shape);
    % If dots
    else
        % Get the dots
        h = plot3(xcomb(1, :), xcomb(2, :), xcomb(3, :), 'k.', 'MarkerSize', 20, 'LineStyle', 'None', 'LineWidth', 2);
    end
end


end