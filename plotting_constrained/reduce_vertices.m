function V2keep = reduce_vertices(V, tol)
%REDUCE_VERTICES reduces the number of vertices of a polyhedron by
%eliminating vertices which are a distance tol from each other.
%   
%   Notes: Not very efficient :)
%
%   V2keep = REDUCE_VERTICES(V, tol)
%
%   Inputs:
%   V ~ matrix of vertices, each column representing a point in
%   n-dimensional space
%   tol ~ minimum distance allowed between two vertices
%   
%   Outputs:
%   V2keep ~ matrix of kept vertices, with vertices being less than tol 
%   away from the kept vertices being removed.

% Number of vertices
[n, nv] = size(V);

% Distance matrix
% D = zeros(nv, nv);
D = diag((1+tol) * ones(1, nv)); % Initialize so diagonal elements are greater than tolerance

% Look through pairs of vertices and calculate distances
for ii = 1 : nv - 1
    for jj = ii + 1 : nv
        D(ii, jj) = norm(V(:, ii) - V(:, jj));
        D(jj, ii) = D(ii, jj);
    end
end

% Discarding vertices: Keep the ones which have the most vertices close to
% them, and discard the ones that are closest to the ones that have the
% most
V2keep = V;

% While there are vertices close enough for discarding
while any(D <= tol, 'all')
    
    % For each vertex count the number of vertices closer than tol
    NumCV = sum(D <= tol, 1);
    
    % Sort
    [~, indsort] = sort(NumCV, 'descend');
    
    % Pick the vertex with the largest amount of closest vertices closest
    % to it
    ind2keep = indsort(1);
    
    % Find all vertices at a distance less than tol from said vertex and
    % throw them away
    ind2discard = find(D(:, ind2keep) <= tol);
    
    % Discard those vertices
    D(ind2discard, :) = [];
    D(:, ind2discard) = [];
    V2keep(:, ind2discard) = [];
end


end