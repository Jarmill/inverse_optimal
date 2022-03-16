function h  = plot_x_star(x_star, varargin)
%PLOT_X_STAR plots the optima of a library of strictly convex quadratic 
%cost functions.
%   
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_X_STAR(x_star)
%   h = PLOT_X_STAR(x_star, sz)
%   h = PLOT_X_STAR(x_star, sz, cmap)
%
%   Inputs:
%   x_star ~ cell array of vectors corresponding to global minima of each
%   quadratic
%   sz ~ marker size on the plot
%   cmap ~ colormap to use for the plotting (optional, default: cool)
%   
%   Outputs:
%   h ~ handle of the plotted object

% Get the number of cost functions in the library
nf = length(x_star);
% Get the size of the optimization variables
n = length(x_star{1});

% Check the size of input
if n ~= 2 && n ~= 3
    error("Dimensionality of the problem is not 2D nor 3D. Please supply 2D or 3D x_stars.");
end

% If cmap argument is passed
if length(varargin) >= 2    
    cmap = varargin{2};
    % Sizecheck
    if ~isequal(size(cmap, 2), 3)
        error('The colormap must be a matrix with 3 columns.');
    elseif size(cmap, 1) < nf
        error('The colormap must be a matrix with at least the same amount of rows as there are "x_star"s.');
    end
% Use default
else
    cmap = cool(nf);
end

% If size argument is passed
if length(varargin) >= 1
    sz = varargin{1};
% Use default
else
    sz = 75;
end

% Transform points into matrix
x_star = cell2mat(x_star);

% Reshape if needed
if size(x_star, 1) ~= n
    reshape(x_star.', n, []);
end

% Get the set shape differently for 2D and 3D
if n == 2
    h = scatter(x_star(1, :), x_star(2, :), sz, cmap, 'filled');
elseif n == 3
    h = scatter3(x_star(1, :), x_star(2, :), x_star(3, :), sz, cmap, 'filled');
end


end