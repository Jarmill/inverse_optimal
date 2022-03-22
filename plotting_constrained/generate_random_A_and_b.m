function [A, b] = generate_random_A_and_b(n, na, range_b)
%GENERATE_RANDOM_A_AND_B generates two arrays, one of which is a matrix
%containing in its rows unit vectors from the n-sphere, the other one
%being a uni-dimensional array whose elements are uniformly sampled from a
%given range.
%
%The outputs can be interpreted as random linear inequality or equality 
%constraints.
%   
%   [A, b] = GENERATE_RANDOM_A_AND_B(n, na, range_b)
%   Inputs:
%   n ~ dimensionality of rows of A
%   na ~ number of random rows of A, and elements of b to generate.
%   range_b ~ the range from which elements of b should be uniformly
%   sampled.

% Define vector normalization
rows_normalize = @(M) M ./ sqrt(sum(M.^2, 2));

% Prealocate output A as a rotationally invariant random vector
A = randn(na, n);

% Normalize columns of A so they are on the unit sphere
A = rows_normalize(A);

% Generate b from a random range
b = diff(range_b) * rand(na, 1) + range_b(1);

end