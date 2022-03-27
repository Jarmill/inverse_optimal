function [dist_rec, info] = sdp_uncons_solve(y, Q, x_star)
%UNTITLED Summary of this function goes here
%SDP_UNCONS_SOLVE get alower bound for the distance to  the
%unconstrained global minimizers.
[Fd, objd, indexer] = sdp_full_uncons_yalmip(y, Q, x_star);

opts = sdpsettings('solver', 'mosek');
% opts = sdpsettings('solver', 'mosek', 'savesolveroutput', true, 'savesolverinput', true);
sol = optimize(Fd, objd, opts);

%export the result of the sdp
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
    info.x_rec = M_rec(indexer.x, 1);
    
    
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


% info.status = sol.problem;

end

