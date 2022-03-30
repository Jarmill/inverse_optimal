function [dist_rec, info] = reweight_single(input, bound, rank_opts)
%REWEIGHT_SINGLE Reweighted trace heuristic based rank minimization of the
%projgm problem. 
%
% bisection on the distance 'bound', with objective bound^2.
%
%reweight_single_uncons 

if input.cons == 1
    [Fd, objd, indexer] = sdp_full_cons_lin_yalmip(input.y, input.Q, input.f, ...
        input.A, input.b, input.Aeq, input.beq, 0);
else
    [Fd, objd, indexer] = sdp_full_uncons_yalmip(input.y, input.Q, input.x_star);
end

Fd_bound = [Fd; objd<=bound^2];

if nargin == 2
    rank_opts = struct('delta', 1e-3, 'rank1tol', 1e-2, 'N_iter', 10);
end

W = eye(indexer.nm);

opts = sdpsettings('solver', 'mosek');


%% start the for loop
% rank_opts.N_iter = 1;
for i = 1:rank_opts.N_iter

    %solve the problem

% opts = sdpsettings('solver', 'mosek', 'savesolveroutput', true, 'savesolverinput', true);
cost_bound = trace(W*indexer.M);
sol = optimize(Fd_bound, cost_bound, opts);

%export the result of the sdp
info = struct;
info.status = sol.problem;
info.rank1ness = [];
% 
if sol.problem == 0
    M_rec = value(indexer.M);
    M_rec_clean = M_rec;
    M_rec_clean(abs(M_rec_clean) <= 1e-6) = 0;
    e_rec = sort(eig(M_rec), 'descend');
    rank1_rec = (e_rec(1)/e_rec(2));
    
%     M_x = M_rec(
    
    dist_rec = sqrt(value(indexer.dist));
%     
    info.M_rec = M_rec_clean;
    info.e_rec = e_rec;
    info.alpha_rec = M_rec(indexer.a, 1);
    info.x_rec = M_rec(indexer.x, 1);
    info.rank1ness = [info.rank1ness; rank1_rec];
    
    if input.cons
        info.lambda_rec = value(indexer.lambda);
        info.mu_rec = value(M_rec(indexer.mu, 1));
    end   
    
    if (rank1_rec >= 1/rank_opts.rank1tol) && (e_rec(2) <= rank_opts.rank1tol)
        break
    end
else
    dist_rec = NaN;
    break;
    
end

% update the reweighting
W = inv(M_rec + rank_opts.delta*eye(size(M_rec, 1)));


% 
%     info.dist_rec = dist_rec;
% 

% info.status = sol.problem;

end

end

