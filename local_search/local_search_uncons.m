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

Xmat= 0;
pobjround = 0;
info = 0;
end