function d = affine_dimension(P, tol)
%AFFINE_DIMENSION find the affine dimension of a set of points.
%   d = AFFINE_DIMENSION(P, tol)
%
%   Inputs: 
%   P ~ matrix of input points, points are column vectors
%   tol ~ tolerance of the rank calculation for the affine dimension
%
%   Outputs:
%   d ~ affine dimension

% Take a reference point and subtract it away from the others
Pminusref = P(:, 1:end-1) - P(:, end);

% Determine the dimension as the rank of the matrix
d = rank(Pminusref, tol);
end