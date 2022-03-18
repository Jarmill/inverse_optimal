%% figure out the function sequence of calling the mexsvec function, which
%  is currently crashing



n = 4;
% A = magic(n);
% X = A'*A;
X = ones(n);
blk = {'s', n};
blk2 = {'s', n; 's', n};

% vX = mexsvec(blk, X)
vX2= svec(blk2, {X; X})