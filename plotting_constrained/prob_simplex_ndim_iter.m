function out = prob_simplex_ndim_iter(d, n, m, index_arr, sum, out)
%PROB_SIMPLEX_NDIM_ITER is an iteratively implemented function for generation of
%a grid of points on a probability simplex embedded in n+1 dimensions
%(of affine dimension n).
% out = PROB_SIMPLEX_NDIM_ITER(d,n,m,index_arr,sum,out)
% Inputs:
%
%   d : is the currently visited dimension.
%   n : is the affine dimension of the probability simplex (thus the 
%   dimension in which it is embedded minus 1, or in other words: the number
%   of its vertices minus 1).
%   m : is the number of grid points along each dimension of the simplex, which
%   can be indexed from 0 to m-1.
%   index_arr : is an array with (n+1) elements whose elements all correspond to
%   a single dimension, and contain the indices (ordinal numbers) of the current
%   grid point along the corresponding dimension.
%   sum : is the accumulated sum of indices that have been attributed until
%   now.
%   out : is the structure in which the output should be stored.
%   
% Outputs:
%   
%   out : a structure containing two fields
%           - array : a matrix which will contain the (n+1)-dimensional
%           gridpoints as rows.
%           - ind : a scalar which should contain the index within the
%           array of the last assigned gridpoint.
%
% Should be initialized as follows
% d = 1
% n_simplex = # of vertices - 1
% m_grid = # of evenly space points along each dimension
% index_arr = zeros(1, n+1);
% sum = 0;
% out.array = zeros(nchoosek(n+m-1, n), n+1);
% out.ind = 0;

% If current dimension is the last, it means all indices have been assigned
% and one should call the function
if d > n
    % Add a last element to it
    index_arr(d) = (m-1)-sum;
    
    % Add to output
    out.array(out.ind + 1, :) = index_arr ./ (m-1);
    out.ind = out.ind + 1;
    return
end

% Otherwise, visit each element of the current dimension
for ii = 0 : ((m-1) - sum)
    
    % Set the current element index
    index_arr(d) = ii;
    
    % Add the current index to the sum
    sumii = sum + ii;
        
    % Call the for loop for the next dimension
    out = prob_simplex_ndim_iter(d+1, n, m, index_arr, sumii, out);
end

end