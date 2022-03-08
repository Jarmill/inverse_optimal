function X =  simplex_project(Y)
%pass in Y as a row vector (matrix)
%project Y onto the simplex sum(X)=1, X>=0

%https://arxiv.org/pdf/1309.1541.pdf
%use Algorithm 1 from the above paper

%The following vectorized Matlab code implements algorithm 1. It projects each row vector in the N × D
%matrix Y onto the probability simplex in D dimensions.


[N,D] = size(Y);
X = sort(Y,2,'descend');
Xtmp = (cumsum(X,2)-1)*diag(sparse(1./(1:D)));
X = max(bsxfun(@minus,Y,Xtmp(sub2ind([N,D],(1:N)',sum(X>Xtmp,2)))),0);

end