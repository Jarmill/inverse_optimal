function x = generate_random_x(n, nx, range_ring)
%GENERATE_RANDOM_X generates a cell array of random vectors from a
%n-dimensional sphere ring.
%Generates an n-dimensional point on an n-sphere and then scales it by a
%number sampled uniformly from a supplied interval (radii of the rings).
%   
%   X = GENERATE_RANDOM_X(n, nqx, range_ring)
%   Inputs:
%   n ~ size of points x.
%   nx ~ number of random x's to generate.
%   range_ring ~ two-element vector describing the range from which random
%   ring radii can be sampled uniformly.

% Define vector normalization
vector_normalize = @(V) V ./ norm(V);

% Prealocate output cell array
x = cell(1, nx);

% Generate nx random n-vectors
V = randn(n, nx);

% Generate random nx scaling factors from range given by input range
s = diff(range_ring) * rand(1, nx) + range_ring(1); 

% For each requested x
for ii = 1 : nx
    
    % Set output vectors
    x{ii} = s(ii) * vector_normalize(V(:, ii));
end

end