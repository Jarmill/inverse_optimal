function [out_rec, info] = manopt_search_cons(y, Q, f, A, b, Aeq, beq, alpha0)
%MANOPT_SEARCH_CONS perform a local search through manopt for the
%optimizing alpha
%   Detailed explanation goes here
n = length(y);
m = length(f);

% x = euclideanfactory(n, 1);
%hadamard map f(alpha) = f(z .* z)
z = spherefactory(m, 1);

elements = struct('A', z);
% manifold = productmani
% manifold    = productmanifold(elements);

%get the yalmip optimizers for quadratic programming
[P_func,P_selector] = cons_local_optimizers(y, Q, f, A, b, Aeq, beq);

problem.M   = z;

%TODO: fix minimum distance selector code
problem.cost = @(z) cost_cons_alpha(z.^2, P_func, P_selector, y);
% problem.egrad = @(z) hproj_egrad_hadamard(z, Q, x_star, y, 0);

options = struct('tolgradnorm', 1e-4);

if nargin == 8
%     problem.X0 = sqrt(alpha0);
    z0 = sqrt(alpha0);
    [z_rec, xcost, info, options] = trustregions(problem, z0, options);
else
    [z_rec, xcost, info, options] = trustregions(problem, [], options);

end

% checkgradient(problem);



out_rec = struct;
out_rec.z = z_rec;
out_rec.alpha = z_rec.^2;

%slightly ineffecient to do one more QP call
[out_rec.norm2, out_rec.obj, out_rec.x] = cost_cons_alpha(out_rec.alpha, P_func, P_selector);
% out_rec.x = x_opt_uncons(out_rec.alpha, Q, x_star);
% out_rec.cost = xcost;
out_rec.dist = sqrt(out_rec.norm2);

end

