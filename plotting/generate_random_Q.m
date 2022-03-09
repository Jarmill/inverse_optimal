function Q = generate_random_Q(n, nq, range_ev)
%GENERATE_RANDOM_Q generates a cell array of random symmetric positive
%definite matrices of size n.
%The algorithm works by iteratively choosing eigenvectors from the
%intersection of the unit sphere with the orthogonal complement of the
%subspace spanned by currently chosen eigenvectors.
%Each eigenvector is associated with an eigenvalue which is uniformly
%sampled from a provided range.
%   
%   Q = GENERATE_RANDOM_Q(n, nq, range_ev)
%   Inputs:
%   n ~ size of square matrix Q.
%   nq ~ number of random Q's to generate.
%   range_ev ~ two-element vector describing the range from which random
%   eigenvalues can be sampled uniformly.

% Define vector normalization
vector_normalize = @(V) V ./ norm(V);

% Prealocate output cell array
Q = cell(1, nq);

% For each requested Q
for ii = 1 : nq
    
    % Prealocate matrix of eigenvectors
    V = zeros(n, n);
    
    % Sample first vector
    V(:, 1) = vector_normalize(randn(n, 1));
    
    % For each subsequent vector, sample it from the intersection of the
    % unit sphere and the orthogonal complement of the subspace of the
    % already chosen vectors
    for jj = 2 : n
        
        % Basis of the subspace spanned by currently chosen vectors
        B = V(:, 1:jj-1);
        
        % Find an (orthonormal) basis of the orthogonal complement of 
        % subspace spanned by currently chosen vectors        
        N = null(B.');
        
        % Generate a random vector from the unit sphere of dimension
        % n-jj+1
        V_prime = vector_normalize(randn(n-jj+1, 1));
        
        % Use previous unit vector to generate a linear combination of the
        % orthonormal basis spanning the orthogonal complement
        %
        % Equivalent to generating a random vector from the intersection of
        % the orthogonal complement and the unit spher in n dimensions
        V(:, jj) = N * V_prime;
    end
    
    % Generate uniformly distributed eigenvalues from the range provided
    s = diff(range_ev) * rand(1, n) + range_ev(1);
    
    % Diagonalize
    S = diag(s);
    
    % Create matrix Q from spectral decomposition
    Q{ii} = V * S * V.';
end

end