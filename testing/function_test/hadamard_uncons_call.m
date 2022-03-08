rng(33, 'twister')
load('opt_5_2d.mat')

y = [2; 1];

% [SDP_out,recoverdata,diagnostic,interfacedata] = sdp_full_uncons(y, Q, x_star);

%% testing for the manopt

n = length(x_star{1});
m = length(x_star);

% x = euclideanfactory(n, 1);
%hadamard map f(alpha) = f(z .* z)
z = spherefactory(m, 1);

elements = struct('A', z);
% manifold = productmani
% manifold    = productmanifold(elements);

problem.M   = z;
problem.cost = @(z) hproj_cost(z, Q, x_star, y);
problem.egrad = @(z) hproj_egrad_hadamard(z, Q, x_star, y, 0);

checkgradient(problem);
 
[x, xcost, info, options] = trustregions(problem);

figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

% problem.cost = @(xin) norm(xin.A-y);

% 
% z = sphere_sample(1, 3);
% fz = hproj_cost(z, Q, x_star, y);



function f = hproj_cost(z, Q, x_star, y)
    alpha = z.*z;
    f = cost_uncons_alpha(alpha, Q, x_star, y);
end

