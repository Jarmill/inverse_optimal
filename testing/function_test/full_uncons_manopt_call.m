load('opt_3_2d.mat')

y_out = [2; 1];

% [SDP_out,recoverdata,diagnostic,interfacedata] = sdp_full_uncons(y_out, Q, x_star);

%% testing for the manopt

n = length(x_star{1});
m = length(x_star);

x = euclideanfactory(n, 1);
alpha = euclideanfactory(m, 1);

elements = struct('A', x, 'B', alpha);
% manifold = productmani
manifold    = productmanifold(elements);

problem.M   = manifold;
problem.cost = @(xin) norm(xin.A-y);