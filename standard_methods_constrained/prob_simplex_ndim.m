function arr = prob_simplex_ndim(n_simplex, m_grid)
% PROB_SIMPLEX_NDIM is an iteratively implemented function for generation of
% a grid of points on a probability simplex embedded in n_simplex+1 dimensions
% (of affine dimension n_simplex).
%   arr = PROB_SIMPLEX_NDIM(n_simplex,m_grid)
%
% Inputs:
%   n_simplex : is the affine dimension of the probability simplex (thus the 
%   dimension in which it is embedded minus 1, or in other words: the number
%   of its vertices minus 1).
%   m_grid : is the number of grid points along each dimension of the simplex.
%
% Outputs:
%   arr: a matrix which will contain nchoosek(m_grid + n_simplex - 1, n_simplex) rows
%        which represent the (n_simplex+1)-dimensional gridpoints.

% Initialize according to prob_simplex_ndim_iter
d = 1;
index_arr = zeros(1, n_simplex+1);
sum = 0;
out.array = zeros(nchoosek(n_simplex+m_grid-1, n_simplex), n_simplex+1);
out.ind = 0;

% Call the iterative implementation
out = prob_simplex_ndim_iter(d, n_simplex, m_grid, index_arr, sum, out);

% The return value is the array within out
arr = out.array;

end