function [h_arc, h_hess, varargout]  = plot_arc(Q, x_star, alpha, n_ell, levels, varargin)
%PLOT_ARC plots the arc of global optima which are a mixture of two
%quadratic cost functions, and plots the level sets of the functions along
%the way.
%   
%   Notes:
%   - Call figure before this function
%
%   [h_arc, h_hess] = PLOT_ARC(Q, x_star, alpha, n_ell, levels)
%   [h_arc, h_hess] = PLOT_ARC(Q, x_star, alpha, n_ell, levels, cmap)
%   [h_arc, h_hess, x, H] = PLOT_ARC(...)
%
%   Inputs:
%   Q ~ length-2 cell array of positive semidefinite matrices defining the \
%   quadratic terms
%   x_star ~ length-2 cell array of vectors corresponding to global minima 
%   of each quadratic
%   alpha ~ array of scalars between 0 and 1 which will be used to combine
%   the two given quadratic functions
%   n_ell ~ scalar indicating number of points used to generate the
%   levels ~ scalar or array of the same length as alpha
%   
%   Outputs:
%   h_arc ~ handles of the plotted piece-wise arc
%   x ~ cell array of 

% Get the number of cost functions in the library
nf = length(Q);

% Check the number of cost functions
if nf ~= 2
    error("There may be only 2 cost functions in the library (length(Q) = 2; length(x_star) = 2).");
end

% Get the size of the optimization variables
n = length(x_star{1});

% Check the size of input
if n ~= 2 && n ~= 3
    error("Dimensionality of the problem is not 2D nor 3D. Please supply 2D or 3D x_stars.");
end

% Check alpha_list 
if ~all(0 <= alpha & alpha <= 1)
    error('alpha_list must have elements in the range between 0 and 1.');
end

% Check for inputs: cmap
if length(varargin) >= 1
    if isequal(size(varargin{1}), [1, 3])
        cmap = repmat(varargin{1}, length(alpha), 1);
    elseif ~isequal(size(varargin{1}), [length(alpha), 3])
        error('cmap must be either an RGB-triplet or a matrix of RGB-triplets with an element less than alpha_list.');
    else
        cmap = varargin{1};
    end
% Default if no inputs for cmap
else
    cmap = jet(length(alpha));
end

% Prepare matrix of optimal locations
x = zeros(n, length(alpha));

% For each alpha along the list calculate the optimum location and hessian
H = cell(1, length(alpha) - 1);
for ii = 1 : length(alpha)
    
    % Combine hessians
    Qmix = alpha(ii) * Q{1} + (1-alpha(ii)) * Q{2};
    
    % Store combined hessian
    H{ii} = Qmix;
    
    % Combine minima
    xmix = alpha(ii) * Q{1} * x_star{1} + (1-alpha(ii)) * Q{2} * x_star{2};
    
    % Find compound global minima
    x(:, ii) = Qmix \ xmix;
end

% Get an array of graphical objects
h_arc = gobjects(length(alpha) - 1, 1);


% Plot each part of the line in different colors
for ii = 1 : length(alpha) - 1
    
    % For 2D vectors
    if n == 2
        h_arc(ii) = plot(x(1, ii:ii+1), x(2, ii:ii+1), 'Color', cmap(ii, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
    end
    
    % For 3D vectors
    if n == 3
        h_arc(ii) = plot3(x(1, ii:ii+1), x(2, ii:ii+1), x(3, ii:ii+1), 'Color', cmap(ii, :), 'HandleVisibility', 'Off', 'LineWidth', 2);
    end
end

% Get the x's into cell array format
x = mat2cell(x, n, ones(1, length(alpha)));

% Plot the level sets of the intermediate points in the list
h_hess = plot_level_sets(H, x, n_ell, levels, cmap);

% If additionnal arguments are requested
if nargout > 2
    % Give minima along the way as outputs in a cell array
    varargout{1} = x;
    % Give hessians along the way as outputs
    varargout{2} = H;    
end

end