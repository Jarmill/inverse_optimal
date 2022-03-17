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
% M_alpha = M(aug_index, aug_index);

%extract the top eigenvector from the moment matrix of alpha
% [alpha_rec, lam_rec] = eigs(M_alpha, 1);

%an extremely rough rounding method based on the pseudomoments
%TODO: replace this with top eigenvalue
alpha_rec = M(alpha_index, 1);

if nargout == 2
%now get the x point associated with the alpha (may not be necessary)
x_rec = opt.x_opt(alpha_rec);
end


end

