function [dist_rec, info] = sdp_cons_lin_solve(y, Q, f, A, b, Aeq, beq, DUAL)
%SDP_CONS_LIN_SOLVE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 8
    DUAL = 0;
end
[Fd, objd, indexer] = sdp_full_cons_lin_yalmip(y, Q, f, A, b, Aeq, beq, DUAL);

opts = sdpsettings('solver', 'mosek');

% opts = sd\psettings('solver', 'mosek', 'savesolveroutput', true, 'savesolverinput', true);
sol = optimize(Fd, objd, opts);

info = struct;
info.status = sol.problem;

if sol.problem == 0
    M_rec = value(indexer.M);
    M_rec_clean = M_rec;
    M_rec_clean(abs(M_rec_clean) <= 1e-6) = 0;
    e_rec = sort(eig(M_rec), 'descend');
    dist_rec = sqrt(value(indexer.dist));
    
    info.M_rec = M_rec_clean;
    info.e_rec = e_rec;
    info.alpha_rec = M_rec(indexer.a, 1);
    info.x_rec = M_rec(indexer.a, 1);
    info.mu_rec = M_rec(indexer.mu, 1);
    info.lambda = value(indexer.lambda);
    
    
    %round M to get the approximate alpha, x
%     rrPar = struct;
%     rrPar.Q = Q;
%     rrPar.x_star = x_star;
%     rrPar.y = y;
%     [info.alpha_rec, info.x_rec] = local_search_uncons({info.M_rec}, [], rrPar, [], true);
else
    dist_rec = NaN;
end

    info.dist_rec = dist_rec;



end

