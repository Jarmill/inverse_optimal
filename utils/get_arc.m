function [Q_mix, x_mix_star] = get_arc(Q, x_star, alpha)
%GET_ARC get points along the arc of global minima from (Q{1}, x_star{1})
%  to (Q{2}, x_star{2}). Only two cost functions are allowed
%
%   Inputs:
%   Q ~ length-2 cell array of positive semidefinite matrices defining the \
%   quadratic terms
%   x_star ~ length-2 cell array of vectors corresponding to global minima 
%   of each quadratic
%   alpha ~ array of scalars between 0 and 1 which will be used to combine
%   the two given quadratic functions
%
%   Q_mix ~ cell array of mixed quadratic costs (hessians)
%   x_mix ~ cell array of mixed quadratic optima


nf = length(Q);

% Check the number of cost functions
if nf ~= 2
    error("There may be only 2 cost functions in the library (length(Q) = 2; length(x_star) = 2).");
end

% Get the size of the optimization variables
n = length(x_star{1});


% Prepare matrix of optimal locations
x_mix_star = zeros(n, length(alpha));

% For each alpha along the list calculate the optimum location and hessian
Q_mix = cell(1, length(alpha) - 1);
for ii = 1 : length(alpha)
    
    % Combine hessians
    Qmix = alpha(ii) * Q{1} + (1-alpha(ii)) * Q{2};
    
    % Store combined hessian
    Q_mix{ii} = Qmix;
    
    % Combine minima
    xmix = alpha(ii) * Q{1} * x_star{1} + (1-alpha(ii)) * Q{2} * x_star{2};
    
    % Find compound global minima
    x_mix_star(:, ii) = Qmix \ xmix;
end

% x_mix_star = x_mix_star, n, ones(1, length(alpha)));

end

