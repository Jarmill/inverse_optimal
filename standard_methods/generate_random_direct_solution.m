function [x0, alpha0] = generate_random_direct_solution(Q, x_star)
%GENERATE_RANDOM_INITIAL_SOLUTION generates a random initial solution for
%of the direct problem.
%   
%   [x0, alpha0] = generate_random_direct_solution(Q, x_star)


% Extract numbers
n = numel(x_star{1});   % Dimensionality of variable
nf = numel(Q);   % Dimensionality of parametrization

% Generate random alpha
alpha0 = rand(nf, 1);

% Generate corresponding solution
Qmix = 0;
xmix = 0;
for ii = 1 : nf
    Qmix = Qmix + alpha0(ii) * Q{ii};
    xmix = xmix + alpha0(ii) * Q{ii} * x_star{ii};
end

x0 = Qmix \ xmix;
end