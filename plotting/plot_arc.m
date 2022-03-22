function [h_arc, h_hess, varargout]  = plot_arc(Q, x_star, alpha_arc, alpha_hess, n_ell, levels, varargin)
%PLOT_ARC takes in two arrays of weights. It plots the arc of global optima
%which are a mixture of two quadratic cost functions given by one array of 
%weights, and plots the level sets of the mixtures given by another set of 
%weights.
%   
%   Notes:
%   - Call figure before this function
%
%   [h_arc, h_hess] = PLOT_ARC(Q, x_star, alpha_arc, alpha_hes, n_ell, levels)
%   [h_arc, h_hess] = PLOT_ARC(Q, x_star, alpha_arc, alpha_hes, n_ell, levels, cmap_arc, cmap_hess)
%   [h_arc, h_hess, x_arc, H_arc] = PLOT_ARC(...)
%   [h_arc, h_hess, x_arc, H_arc, x_hess, H_hess] = PLOT_ARC(...)
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
%   cmap ~ colormap (length(alpha)x3)
%   
%   Outputs:
%   h_arc ~ handles of the plotted arc
%   x ~ cell array of optimal points of the mixture
%   H ~ cell array of hessians of the mixture

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

% Check alpha_arc
if ~all(0 <= alpha_arc & alpha_arc <= 1)
    error('alpha_arc must have elements in the range between 0 and 1.');
end

% Check alpha_hess 
if ~all(0 <= alpha_hess & alpha_hess <= 1)
    error('alpha_hess must have elements in the range between 0 and 1.');
end

% Sort alphas in descending order
alpha_arc = sort(alpha_arc, 'descend');
alpha_hess = sort(alpha_hess, 'descend');

% Check for inputs: cmap
if length(varargin) >= 1
    
    % Check cmap_arc
    if isequal(size(varargin{1}), [1, 3])
        cmap_arc = repmat(varargin{1}, length(alpha_arc), 1);
    elseif ~isequal(size(varargin{1}), [length(alpha_arc), 3])
        error('cmap_arc must be either an RGB-triplet or a matrix of RGB-triplets with an element less than alpha_list.');
    else
        cmap_arc = varargin{1};
    end
    
    % Check cmap_hess
    if isequal(size(varargin{2}), [1, 3])
        cmap_hess = repmat(varargin{2}, length(alpha_hess), 1);
    elseif ~isequal(size(varargin{2}), [length(alpha_hess), 3])
        error('cmap_arc must be either an RGB-triplet or a matrix of RGB-triplets with an element less than alpha_list.');
    else
        cmap_hess = varargin{2};
    end
% Default if no inputs for either cmap
else
    cmap_arc = copper(length(alpha_arc));
    cmap_hess = copper(length(alpha_hess));
end

% Prepare matrix of optimal locations for the arc
x_arc = zeros(n, length(alpha_arc));
% Prepare matrix of optimal locations for the hessians along the way
x_hess = zeros(n, length(alpha_hess));

% For each alpha along the arc calculate the hessian
H_arc = cell(1, length(alpha_arc) - 1);
% For each alpha along the way calculate hessian
H_hess = cell(1, length(alpha_hess) - 1);

% Compute the arc related 
for ii = 1 : length(alpha_arc)
    
    % Combine hessians
    Qmix = alpha_arc(ii) * Q{1} + (1-alpha_arc(ii)) * Q{2};
    
    % Store combined hessian
    H_arc{ii} = Qmix;
    
    % Combine minima
    xmix = alpha_arc(ii) * Q{1} * x_star{1} + (1-alpha_arc(ii)) * Q{2} * x_star{2};
    
    % Find compound global minima
    x_arc(:, ii) = Qmix \ xmix;
end

% Compute the hessian-level-set related info
for ii = 1 : length(alpha_hess)
    
    % Combine hessians
    Qmix = alpha_hess(ii) * Q{1} + (1-alpha_hess(ii)) * Q{2};
    
    % Store combined hessian
    H_hess{ii} = Qmix;
    
    % Combine minima
    xmix = alpha_hess(ii) * Q{1} * x_star{1} + (1-alpha_hess(ii)) * Q{2} * x_star{2};
    
    % Find compound global minima
    x_hess(:, ii) = Qmix \ xmix;
end

% Get x_hess in cell array format
x_hess = mat2cell(x_hess, n, ones(1, length(alpha_hess)));
% Plot the level sets of the intermediate points in the list
h_hess = plot_level_sets(H_hess, x_hess, n_ell, levels, cmap_hess);

% Get an array of graphical objects for the arc
h_arc = gobjects(length(alpha_arc) - 1, 1);
% Plot each part of the line in different colors
for ii = 1 : length(alpha_arc) - 1
    
    % For 2D vectors
    if n == 2
        h_arc(ii) = plot(x_arc(1, ii:ii+1), x_arc(2, ii:ii+1), 'Color', cmap_arc(ii, :), 'HandleVisibility', 'Off', 'LineWidth', 3);
    end
    
    % For 3D vectors
    if n == 3
        h_arc(ii) = plot3(x_arc(1, ii:ii+1), x_arc(2, ii:ii+1), x_arc(3, ii:ii+1), 'Color', cmap_arc(ii, :), 'HandleVisibility', 'Off', 'LineWidth', 3);
    end
end
% Get the x's into cell array format
x_arc = mat2cell(x_arc, n, ones(1, length(alpha_arc)));

% If additionnal arguments are requested
if nargout > 2
    % Give minima along the arc as outputs in a cell array
    varargout{1} = x_arc;
    % Give hessians along the arc as outputs
    varargout{2} = H_arc;
end
% If additionnal arguments are requested
if nargout > 4
    % Give minima of the functions with plotted level sets
    varargout{3} = x_hess;
    % Give hessians of the functions with plotted level sets
    varargout{4} = H_hess;
end


end