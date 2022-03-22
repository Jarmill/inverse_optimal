function h  = plot_mathcal_G_constrained(Q, phi, A, b, C, d, xl, xu, n_grid, n_mesh, varargin)
%PLOT_MATHCAL_G_CONSTRAINED plots the set of all global optima of the 
%convex combination of a library of quadratic cost functions, given that 
%library, given the bounds on the decision variables and the constraints.
%   
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_MATHCAL_G(Q, phi, A, b, C, d, xl, xu, n_grid, n_mesh)
%   h = PLOT_MATHCAL_G(Q, phi, A, b, C, d, xl, xu, n_grid, n_mesh, smoothness)
%
%   Inputs:
%   Q ~ cell array of positive semidefinite matrices defining the quadratic
%   term
%   x_star ~ cell array of vectors corresponding to global minima of each
%   quadratic
%   n_grid ~ upper bound on the number of grid-points on the simplex that 
%   should be calculated to generate the mesh
%   n_mesh ~ number of points that should be, at most, used to represent
%   solutions of quadratics that have mutliple solutions
%   smoothness ~ the scalar alpha parameter of the alphaShape representing 
%   the sets mesh
%   
%   Outputs:
%   h ~ handle of the plotted object

% Get the number of cost functions in the library
nf = length(Q);
% Get the size of the optimization variables
n = size(Q{1}, 1);

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
xcomb = cell(1, size(alpha, 1));

% Calculate optimum for different combinations
for cc = 1 : size(alpha, 1)
    
    % Find linear combinations of Q matrices and F matrics
    Qmix = 0;
    phimix = 0;
    for ii = 1 : nf
        Qmix = Qmix + alpha(cc, ii) * Q{ii};
        phimix = phimix + alpha(cc, ii) * phi{ii};
    end
    
    % Find optimal solution(s)
    x_opt = global_optima_set_qp(Qmix, phimix, A, b, C, d, xl, xu, n_mesh);
    
    % Add the optimal solution(s) to the whole set
    xcomb{cc} = x_opt;
end

% Transform the whole set to a matrix
xcomb = cell2mat(xcomb);
size(xcomb)

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