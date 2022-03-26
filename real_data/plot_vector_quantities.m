function [h_ax, h] = plot_vector_quantities(t, f, h_ax, varargin)
%PLOT_VECTOR_QUANTITIES automatically determines the size of the subgrid in
%which to plot the quantites, and performs the plotting.
%
%   
%   Takes in plot options which are common to all plots.
%   e.g.
%   [h_ax, h] = PLOT_VECTOR_QUANTITIES(t, f, h_ax, 'LineWidth', 2);

% Size of quantities to plot 
% n - number of quantities
% N - number of samples per quantity
[n, N] = size(f);

% Round
figrows = ceil(sqrt(n));
figcols = ceil(sqrt(n));

% Try to reduce number of columns
if (figcols-1) * figrows >= n
    figcols = figcols - 1;
end

% Create an array of graphical objects
h = gobjects(figrows, figcols);
if isempty(h_ax)
    h_ax = gobjects(figrows, figcols);
end

% For each row and column
for ii = 1 : figrows
    for jj = 1 : figcols
        
        % Determine the current order
        curr = (ii-1)*figcols + jj;
        % If current order surpasses the number of elements to plot
        if curr > n
            break
        end
        
        % Create current subplot
        h_ax(ii, jj) = subplot(figrows, figcols, curr);
        % Hold
        hold on;
        
        % If arguments are passed plot with them
        if ~isempty(varargin)
            
            h(ii, jj) = plot(t, f(curr, :), varargin{:});
            
        % If not dont
        else
           
            h(ii, jj) = plot(t, f(curr, :));
            
        end
    end
end


end