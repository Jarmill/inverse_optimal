function [out_rec, info] = manopt_search_uncons(Q, x_star, y, alpha0)
%MANOPT_SEARCH_UNCONS perform a local search through manopt for the
%optimizing alpha
%   Detailed explanation goes here
n = length(x_star{1});
m = length(x_star);

% x = euclideanfactory(n, 1);
%hadamard map f(alpha) = f(z .* z)
z = spherefactory(m, 1);

elements = struct('A', z);
% manifold = productmani
% manifold    = productmanifold(elements);

problem.M   = z;
problem.cost = @(z) cost_uncons_alpha(z.^2, Q, x_star, y);
problem.egrad = @(z) hproj_egrad_hadamard(z, Q, x_star, y, 0);

if nargin == 5
    problem.X0 = sqrt(alpha0);
end

% checkgradient(problem);
 
[z_rec, xcost, info, options] = trustregions(problem);




out_rec = struct;
out_rec.z = z_rec;
out_rec.alpha = z_rec.^2;
out_rec.x = x_opt_uncons(out_rec.alpha, Q, x_star);
out_rec.cost = xcost;
out_rec.dist = sqrt(2*xcost);

end

