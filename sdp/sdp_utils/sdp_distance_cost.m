function [C] = sdp_distance_cost(y_in, n, m)
%SDP_COST generate the distance cost for the Uncons SDP
%   distance from in moment matrix from x to y.
%1/2 norm(y-x, 2)^2 = (sum(x.^2) + sum(y.^2) - sum(2 x.*y))

M_dim = 1+n+m;

x_index = 1 + (1:n);
% alpha_index = (n+1) + (1:(m));


%constant corner entry (sum(y.^2))
C_i = 1;
C_j = 1;
C_v = sum(y.^2);

%affine - sum(2 x.*y)
C_i = [C_i, ones(1, n), x_index];
C_j = [C_j, x_index, ones(1, n)];
C_v = [C_v, -y', -y'];
% C_v = [C_v, -2*y'];

%squared term sum(x.^2)
C_i = [C_i, x_index];
C_j = [C_j, x_index];
C_v = [C_v, ones(1, n)];
% C_i = [C_i; ]
C = sparse(C_i, C_j, C_v, M_dim, M_dim);
end

