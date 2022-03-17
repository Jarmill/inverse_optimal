function [alpha_rec, x_rec] = round_alpha_matrix(Z, opt)
%ROUND_ALPHA_MATRIX round the moment matrix (top eigenvector) in order to
%acquire a point (alpha_rec, x_rec) that is (hopefully) a global minimizer.
%
%extract the alpha from the eigenvector, and then find x_star given that
%alpha.
%
%   Detailed explanation goes here

%yalmip may add other variables to the SDP
%the semidefinite block is last (blk.s)
M = Z{end};

% x_index = 1 + (1:n);
alpha_index = (opt.n+1) + (1:(opt.m));
aug_index = [1, alpha_index];
M_alpha = M(aug_index, aug_index);

%extract the top eigenvector from the moment matrix of alpha
%then project on to the simplex
[v_rec, lam_rec] = eigs(M_alpha, 1);

alpha_rec_infeas = v_rec/v_rec(1);
alpha_rec_feas = simplex_project(alpha_rec_infeas(2:end)')';

%an extremely rough rounding method based on the pseudomoments
%TODO: replace this with top eigenvalue
alpha_rec_mom = M(alpha_index, 1);

cost_feas = opt.cost(alpha_rec_feas);
cost_mom = opt.cost(alpha_rec_mom);
alpha_rec_mom(alpha_rec_mom <= 0) = 0;

%take the rounded alpha with the smaller cost.
if cost_feas <= cost_mom
    alpha_rec = alpha_rec_feas;
else
    alpha_rec = alpha_rec_mom;
end

if nargout == 2
%now get the x point associated with the alpha (may not be necessary)
x_rec = opt.x_opt(alpha_rec);
end


end

