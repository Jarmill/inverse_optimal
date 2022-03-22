function [Q, phi, varargout] = generate_random_Q_and_phi(n, nq,range_ev, range_ring, range_ring_null, rankdef_prob)
%GENERATE_RANDOM_Q_AND_PHI generates two cell arrays, one of which contains
%random symmetric positive definite matrices of size n, and another which
%contains random vectors from a ring about the origin.
%
%The algorithm works by iteratively choosing eigenvectors from the
%intersection of the unit sphere with the orthogonal complement of the
%subspace spanned by currently chosen eigenvectors.
%Each eigenvector is associated with an eigenvalue which is uniformly
%sampled from a provided range.
%
%   
%   [Q, phi] = GENERATE_RANDOM_Q_AND_PHI(n, nq, range_ev, range_ring, rankdef_prob)
%   [Q, phi, x] = GENERATE_RANDOM_Q_AND_PHI(n, nq, range_ev, range_ring, rankdef_prob)
%   
%   Definition:
%   Quadratic invariant line : Affine space where the quadratic-term is
%   constant. Also the affine space for which that constant is the lowest
%   out of all the affine spaces that satisfy that constancy.
%
%   Inputs:
%   n ~ size of square matrix Q and of vector phi.
%   nq ~ number of random Qs and phis to generate.
%   range_ev ~ two-element vector describing the range from which random
%   eigenvalues can be sampled uniformly.
%   range_ring ~ two-element vector describing the radii of the spheres
%   which delimit the ring which must intersect the quadratic-invariant line
%   range_ring_null ~ two-element vector describing the radii of the spheres
%   which delimit the ring from which we choose the part of phi in the
%   nullspace of Q (in the direction of the quadratic-invariant-line).
%   rankdef_prob ~ probability of a column of the matrix Q being rank
%   deficient
%
%   Outputs:
%   Q ~ cell array of positive semidefinite matrices representing the
%   quadratic terms in our library of functions
%   phi ~ cell array of vectors representing the linear terms in our
%   library of cost functions
%   x ~ cell array of vectors which are contained in both the ring and in
%   the quadratic invariant line

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
    
    
    % Generate a vector of binary variables which are 0 with probability
    % rankdef_prob
    rankdef_col = rand(1, n);
    rankdef_col = rankdef_col >= rankdef_prob;
    
    % Generate uniformly distributed eigenvalues from the range provided
    s = diff(range_ev) * rand(1, n) + range_ev(1);
        
    % Multiply the random eigenvalue vector with the rankdef_prob
    s = s .* rankdef_col; 
    
    % Diagonalize
    S = diag(s);    
    
    % Create matrix Q from spectral decomposition
    Q{ii} = V * S * V.';
    % Numerically insure symmetricity
    Q{ii} = (Q{ii} + Q{ii}.')/2;
end

% Define vector normalization
vector_normalize = @(V) V ./ norm(V);

% Prealocate output cell array
x = cell(1, nq);

% Generate nx random n-vectors
V = randn(n, nq);

% Generate random nx scaling factors from range given by input range
s = diff(range_ring) * rand(1, nq) + range_ring(1); 

% For each requested x
for ii = 1 : nq
    
    % Set output vectors
    x{ii} = s(ii) * vector_normalize(V(:, ii));
    
end

% Prealocate output cell array
phi = cell(1, nq);

% For each cost function
for ii = 1 : nq
    
    % The part of phi that is in the range of the matrix Q should be
    % consistent with the quadratic-invariant line condition
    phiR = - Q{ii} * x{ii};
        
    % The part of phi that is in the nullspace of the matrix Q should be
    % chosen at random from the nullspace
    NQ = null(Q{ii});
    
    % Choose a random vector from a ring
    s = diff(range_ring_null) * rand(size(NQ, 2), 1) + range_ring_null(1);
    
    % Choose the part of the vector phi that is in the nullspace of the
    % matrix Q from a random ring lying in the nullspace of Q
    phiN = NQ * s;
    
    % The final vector phi is the sum of its components in the nullspace
    % and in the range of matrix Q
    % Except if nullspace vector is empty
    if ~isempty(phiN)
        phi{ii} = phiR + phiN;
    else
        phi{ii} = phiR;
    end
end

% If multiple outputs requested
if nargout > 2
    varargout{1} = x;
end

end