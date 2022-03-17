function [Xmat,pobjround,info] = local_search_uncons(Zmat,C,rrPar,opt,roundonly)
% Round and refine for Global Optimum Projection
% (structure copied from local_search_quasar in STRIDE example)
% Input:
% Zmat: SDP primal iterate
% C: cost matrix of the SDP
% rrPar: a structure where you can pass any data necessary for local search
% opt: indices of eigenvectors for rounding
% roundonly: if true, then only round a solution without NLP
% Output:
% Xmat: a rank-one SDP iterate
% pobjround: primal SDP cost attained by Xmat

%
x_star = rrPar.x_star;
Q = rrPar.Q;
y = rrPar.y;

n = length(x_star{1});
m = length(x_star);

if nargin < 5
    roundonly = false;
end
if nargin < 4
    % default round the first two eigenvectors
    opt = [1,2]; 
end

%extract the variate information
opt_round = struct('n', n, 'm', m, 'x_opt', @(alpha) x_opt_uncons(alpha, Q, x_star));

[alpha_rec, x_rec] = round_alpha_matrix(Zmat, opt_round);


% if roundonly
    %perform only rounding
    
%     v_rec = [1; alpha_rec; x_rec];
%     Xmat = v_rec * v_rec';
info = struct;
if ~roundonly
    
% else
    %perform local search through manopt
    [out_rec, info_manopt] = manopt_search_uncons(Q, x_star, y, alpha_rec);

    pobjround = info_manopt(end).cost;
    if pobjround == inf
        nlpsuccess   = false;
    else
        nlpsuccess = true;
        alpha_rec = out_rec.alpha;
        x_rec = out_rec.x;
    end

    info.manopt = info_manopt;
end

v_rec = [1; x_rec; alpha_rec];
Xmat = v_rec * v_rec';

% Xmat= 0;
pobjround = norm(x_rec - y, 2)^2;
% info = 0;

info.pobjround = pobjround;

% if nargout > 2
%     info = 
% end

end