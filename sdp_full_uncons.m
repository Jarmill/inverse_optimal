function [model,recoverdata,diagnostic,interfacedata] = sdp_full_uncons(y, Q, x_star)
%SDP_FULL_UNCONS: create an SDP in SDPT3 format for the distance to inverse
%optimal control problem
%
% f_i(x) = (x-x_star{i})'*Q{i}*(x-x_star{i})
% f_alpha = sum_i f_i(x)
%
%INPUTS (may change around)
%   y:      test point outside the GO set
%   Q:      Dictionary of hessians
%   x_star: Corresponding optimal points to quadratic functions

%% create indexers and define variables
%form and index out the moment matrix. remember that alpha_n is eliminated
n = length(x_star{1});
m = length(x_star);

x_index = 1 + (1:n);
alpha_index = (n+1) + (1:(m));

M = sdpvar(1+n+m);

%variables
x = M(x_index, 1);
alpha = M(alpha_index, 1);
x2 = diag(M(x_index, x_index));

xa = M(x_index, alpha_index);

%% start developing constraints
%start the constraints
cons = [(M>=0):'PSD moments', (M(1,1) == 1):'1=1'; ...
    (alpha >= 0):'alpha >= 0', (sum(alpha)==1):'sum(alpha)==1'];

%now do the gradient constraints
% grad_term = Q{m}*(x - x_star{m});
grad_term = zeros(n, 1, 'like', sdpvar);
for j = 1:m
    curr_term = Q{j}*xa(:, j) - Q{j}*x_star{j}*alpha(j);
%     curr_term = (Q{j} - Q{m})*xa(:, j) - (Q{j}*x_star{j} - Q{m}*x_star{m})*alpha(j);
    grad_term = grad_term + curr_term;
end

cons = [cons; (grad_term==0):'KKT stationarity'];

%% redundant constraints (valid inequalities)
con_alpha_diag = (sum(M(alpha_index, alpha_index), 2) == alpha);

cons = [cons; con_alpha_diag:'alpha(i) sum(alpha) = alpha(i)'];
con_alpha2 = [];
for i = 1:(m)
    curr_i = alpha_index(i);

    %off diagonal
    for j = 1:(i-1)        
        curr_j = alpha_index(j);
        con_alpha2 = [con_alpha2; M(curr_i, curr_j) <= 0.25; ...
            M(curr_i, curr_j) >= 0; M(curr_i, curr_j) <= [alpha(i); alpha(j)]];
    end
    
    
%     con_alpha2 = [con_alpha2; alpha(i) 
end
cons = [cons; con_alpha2:'off-diagonal constraints'];

dist = sum(y.^2 - 2*(x.*y) + x2);
cons = [cons; (dist>=0):'distance is nonnegative'];

opts = sdpsettings('solver', 'sdpt3');

% SDP
[model,recoverdata,diagnostic,interfacedata] = export(cons, dist, opts);
end