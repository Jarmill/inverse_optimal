function h  = plot_x_star_constrained(Q, phi, A, b, C, d, xl, xu, varargin)
%PLOT_X_STAR_CONSTRAINED plots the optima of a library of convex quadratic 
%programs, differing only by their cost functions.
%
%   min_x    1/2*x'*Q*x + phi'*x
%   s.t.     xl <= x <= xu
%            A*x == b
%            C*x <= d
%
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_X_STAR_CONSTRAINED(Q, phi, A, b, C, d, xl, xu)
%   h = PLOT_X_STAR_CONSTRAINED(Q, phi, A, b, C, d, xl, xu, sz)
%   h = PLOT_X_STAR_CONSTRAINED(Q, phi, A, b, C, d, xl, xu, sz, cmap)
%
%   Inputs:
%   Q ~ cell array of quadratic-term matrices for the library of cost
%   functions
%   phi ~ cell array of linear-term vectors for the library of cost
%   functions
%   A,b,C,d,xu,xl ~ constraint parameters (Linear equality matrices, linear
%   inequality matrices, bound constraints vectors)
%   cmap ~ colormap to use for the plotting (optional, default: cool)
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

% Number of mesh points for multiple minima
Nmesh = 20;

% Prealocate minima cell array (since the output isn't of uniform size)
x_star = cell(1, nf);
multiplicity_flags = zeros(1, nf, 'logical');

% Calculate the minima
for ii = 1 : nf
    ii
    % Get the matrix of minima
    x_opt = global_optima_set_qp(Q{ii}, phi{ii}, A, b, C, d, xl, xu, Nmesh);
    
    % Set the mutliplicity flag to true when there are multiple minima
    if size(x_opt, 2) > 1
        multiplicity_flags(ii) = true;
    end
    
    % Store the matrix of minima
    x_star{ii} = x_opt;
end


% Prepare to plot the minima separately
h = gobjects(nf, 1);

% Get the set shape differently for 2D and 3D
if n == 2
    
    % Plot each minima separately
    for ii = 1 : nf
        % Show legend for only one case
        if ii == round(nf/2)
            h(ii) = scatter(x_star{ii}(1, :), x_star{ii}(2, :), sz, cmap(ii, :), 'filled', 'DisplayName', sprintf('$x^*_{1:%d}$', nf));
        else
            h(ii) = scatter(x_star{ii}(1, :), x_star{ii}(2, :), sz, cmap(ii, :), 'filled', 'HandleVisibility', 'Off');
        end
    
    end
elseif n == 3
    
    % Plot each minima separately
    for ii = 1 : nf
        
        % Show legend for only one case
        if ii == round(nf/2)
            h(ii) = scatter3(x_star{ii}(1, :), x_star{ii}(2, :), x_star{ii}(3, :), sz, cmap(ii, :), 'filled', 'DisplayName', sprintf('$x^*_{1:%d}$', nf));
        else
            h(ii) = scatter3(x_star{ii}(1, :), x_star{ii}(2, :), x_star{ii}(3, :), sz, cmap(ii, :), 'filled', 'HandleVisibility', 'Off');
        end
    end
end


end