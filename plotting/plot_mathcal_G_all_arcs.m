function h  = plot_mathcal_G_all_arcs(Q, x_star, n_grid)
%PLOT_MATHCAL_G_ALL_ARCS plots the arcs between pairs of all global optima 
%of the convex combination of a library of quadratic cost functions, given 
%that library.
%   
%   Notes:
%   - Call figure before this function
%
%   h = PLOT_MATHCAL_G_ALL_ARCS(Q, x_star, npts)
%
%   Inputs:
%   Q ~ cell array of positive semidefinite matrices defining the quadratic
%   term
%   x_star ~ cell array of vectors corresponding to global minima of each
%   quadratic
%   n_grid ~ number of points that should be calculated To generate the mesh
%   
%   Outputs:
%   h ~ handleS of the plotted objectS

% Get the number of cost functions in the library
nf = length(Q);
% Get the size of the optimization variables
n = length(x_star{1});

% Check the size of input
if n ~= 2 && n ~= 3
    error("Dimensionality of the problem is not 2D nor 3D. Please supply 2D or 3D x_stars.");
end

% Initialize a parametrization
alpha = linspace(0, 1, n_grid).';

% Compute set of minimums of convex combinations
% Prealocate
xcomb = zeros(n, size(alpha, 1), nf*(nf-1)/2);

% Initialize current pair number cuz i'm too lazy/stupid to deduce the formula
% from ii and jj
curr_pair = 1;
% For each pair of cost functions
for ii = 1 : nf - 1
    for jj = ii + 1 : nf
    
    % Calculate optimum for different combinations of the current pair
    for cc = 1 : size(alpha, 1)

        % Find linear combinations of Q matrices and F matrics
        Qcomb = alpha(cc) * Q{ii} + (1 - alpha(cc)) * Q{jj};
        x_star_comb = alpha(cc) * Q{ii} * x_star{ii} + (1 - alpha(cc)) * Q{jj} * x_star{jj};

        % Find optimal solution
        xcomb(:, cc, curr_pair) = Qcomb \ x_star_comb;    
    end
    
    % Update current pair
    curr_pair = curr_pair + 1;
    
    end
end

% Prealocate handle array
h = gobjects(nf * (nf - 1) / 2, 1);

% Get the set shape differently for 2D and 3D
if n == 2
    
    % For each pair
    for ii = 1 : nf*(nf-1)/2
        h(ii) = plot(squeeze(xcomb(1, :, ii)), squeeze(xcomb(2, :, ii)), 'Color', [0, 0, 0], 'LineWidth', 1, 'HandleVisibility', 'Off');
    end
    
elseif n == 3

    % For each pair
    for ii = 1 : nf*(nf-1)/2
        h(ii) = plot3(squeeze(xcomb(1, :, ii)), squeeze(xcomb(2, :, ii)), squeeze(xcomb(3, :, ii)), 'Color', [0, 0, 0], 'LineWidth', 1, 'HandleVisibility', 'Off');
    end
    
end


end