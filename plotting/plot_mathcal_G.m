function h  = plot_mathcal_G(Q, x_star, n_grid, varargin)
%PLOT_MATHCAL_G plots the set of all global optima of the convex combination
%of a library of quadratic cost functions, given that library.
%   
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_MATHCAL_G(Q, x_star, npts)
%   h = PLOT_MATHCAL_G(Q, x_star, npts, smoothness)
%
%   Inputs:
%   Q ~ cell array of positive semidefinite matrices defining the quadratic
%   term
%   x_star ~ cell array of vectors corresponding to global minima of each
%   quadratic
%   n_grid ~ upper bound on the number of points that should be calculated 
%   to generate the mesh
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
    % Get the 2D shape differently if smoothness parameter is given
    if ~isempty(varargin)
        set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).', varargin{1});
    else
        set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).');
    end
    % Get the plot
    h = plot(set_shape);
elseif n == 3
    % Get the 3D shape differently if smoothness parameter is given
    if ~isempty(varargin)
        set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).', xcomb(3, :).', varargin{1});
    else
        set_shape = alphaShape(xcomb(1, :).', xcomb(2, :).', xcomb(3, :).');
    end
    % Get the plot
    h = plot(set_shape);
end


end