function h  = plot_level_sets(Q, x_star, n_ell, varargin)
%PLOT_LEVEL_SETS plots one level set of all individual functions within a 
%library of quadratic cost functions, given that library.
%   
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_LEVEL_SETS(Q, x_star, n_ell)
%   h = PLOT_LEVEL_SETS(Q, x_star, n_ell, levels)
%   h = PLOT_LEVEL_SETS(Q, x_star, n_ell, levels, cmap)
%
%   Inputs:
%   Q ~ cell array of positive semidefinite matrices defining the quadratic
%   term
%   x_star ~ cell array of vectors corresponding to global minima of each
%   quadratic
%   n_ell ~ maximum number of mesh points to generate the level sets
%   levels ~ array of scalars corresponding to the level we should plot for
%   each quadratic function (can be a single scalar)
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

% Check for inputs: levels
if length(varargin) >= 1
    if isscalar(varargin{1})
        levels = varargin{1} * ones(nf, 1);
    elseif ~isequal(length(varargin{1}), nf)
        error('levels must be either a scalar or an array with the same number of elements as cell array Q.');
    else
        levels = varargin{1};
    end
% Default if no inputs for levels
else
    levels = ones(nf, 1);
end

% Check for inputs: cmap
if length(varargin) >= 2
    if isequal(size(varargin{2}), [1, 3])
        cmap = repmat(varargin{2}, nf, 1);
    elseif ~isequal(size(varargin{2}), [nf, 3])
        error('cmap must be either an RGB-triplet or a matrix of RGB-triplets with the same number of rows as cell array Q has elements.');
    else
        cmap = varargin{2};
    end
% Default if no inputs for cmap
else
    cmap = cool(nf);
end


% For two dimensions of input
if n == 2
    % Number of angle discretizations
    nphi = n_ell;
    % Discretization of angle
    phi = linspace(0, 2*pi, nphi);
    % Transform angles
    Cphi = cos(phi);
    Sphi = sin(phi);
    
    % Initialize plot objects for each ellipsoid
    h = gobjects(nf, 1);
    
    % For each function
    for ii = 1 : nf
        
        % Decompose hessian
        [V, D] = eig(Q{ii});
        
        % Extract level
        l = levels(ii);
        
        % Get the mesh of the level set
        x = V * sqrt(l * D^(-1)) * [Cphi; Sphi] + x_star{ii};
        
        % Get the plot of the level set
        h(ii) = plot(x(1, :), x(2, :), 'Color', cmap(ii, :));
    end
    
% For three dimensions of input
elseif n == 3
    % Number of angle discretizations
    nphi = round(n_ell^(2/3));
    ntheta = round(n_ell / nphi);
    
    % Angle discretizations
    phi = linspace(0, 2*pi, nphi);
    theta = linspace(-pi/2, pi/2, ntheta);
    
    % Transform angles
    CphiCtheta = reshape(cos(phi).' * cos(theta), 1, []);
    SphiCtheta = reshape(sin(phi).' * cos(theta), 1, []);
    Stheta = reshape(ones(nphi, 1) * sin(theta), 1, []);
    
    % Initialize plot objects for each ellipsoid
    h = gobjects(nf, 1);
    
    % For each function
    for ii = 1 : nf
        
        % Decompose hessian
        [V, D] = eig(Q{ii});
        
        % Extract level
        l = levels(ii);
        
        % Get the mesh points of the level set
        x = V * sqrt(l * D^(-1)) * [CphiCtheta; SphiCtheta; Stheta] + x_star{ii};
        
        % Get the mesh of the convex hull
        CH = convhull(x(1, :).', x(2, :).', x(3, :).');
        
        % Get the plot of the level set
        h(ii) = trisurf(CH, x(1, :).', x(2, :).', x(3, :).', 'LineStyle', 'none', 'FaceAlpha', 0.5, 'FaceVertexCData', random_colormap(size(CH, 1), cmap(ii, :)));
    end
end

end