function [Fd, objd, indexer] = sdp_full_cons_lin_yalmip(y, Q, f, A, b, Aeq, beq, DUAL)
%SDP_FULL_UNCONS: create an SDP in YALMIP format for the distance to inverse
%optimal control problem with constraints
%
% f_i(x) = 1/2 x' Qi x - fi' x
% f_alpha = sum_i f_i(x)
%
%Aeq x <= beq
%A x = b
%
%
%This problem very specifically has only equality constraints
%so there is no multiplication between lambda and x anywhere.
%
%INPUTS (may change around)
%   y:      test point outside the GO set
%   Q:      Dictionary of hessians
%   f:      vectors in cost functions    
%   A:      Inequality constraint matrix (format copied from linprog)
%   b:      Inequality constraint vector
%   Aeq     Equality constraint matrix (format copied from linprog)
%   beq     Equality constraint vector
%
%Outputs 
%   Fd:     YALMIP constraints
%   objd:   YALMIP objective for distance (primalized)
%   indexer:index and output structure

%% create indexers and define variables
%form and index out the moment matrix. remember that alpha_n is eliminated

if nargin < 8
    %dualize the SDP (put it into primal form. Fewer variables needed to
    %represent the problem)
    DUAL = 1;
end

n = length(f{1});
m = length(f);

%linear and affine constraints
n_ineq = length(b);
n_eq = length(beq);


%[1, x, alpha, lambda, mu]
%[1, x, simplex, eq, ineq]
x_index = 1 + (1:n);
alpha_index = (n+1) + (1:m);
mu_index = (1+n+m) + (1:n_ineq);

%technically there is correlative sparsity:
%(1,x,alpha), (1,x,lambda), (1,x,mu).
%alpha and mu are both nonnegative, so the cross-terms mu*alpha may be
%imposed to be nonnegative as well.

M = sdpvar(1+n+m+n_ineq);

if n_eq
    lambda= sdpvar(n_eq, 1);
else
    lambda = [];
end

%variables
x = M(x_index, 1);
alpha = M(alpha_index, 1);
mu = M(mu_index, 1);

x2 = diag(M(x_index, x_index));

xa = M(x_index, alpha_index);

%% start developing constraints
%start the constraints
cons = [(M>=0):'PSD moments', (M(1,1) == 1):'1=1'; ...
    (alpha >= 0):'alpha >= 0', (sum(alpha)==1):'sum(alpha)==1'];

%now do the gradient constraints
grad_term = zeros(n, 1, 'like', sdpvar);
for j = 1:m
    curr_term = Q{j}*xa(:, j) + f{j}*alpha(j);
%     curr_term = (Q{j} - Q{m})*xa(:, j) - (Q{j}*x_star{j} - Q{m}*x_star{m})*alpha(j);
    grad_term = grad_term + curr_term;
end

if n_eq
    grad_term = grad_term + Aeq'*lambda;
end

if n_ineq
    grad_term = grad_term + A'*mu;
end

cons = [cons; (grad_term==0):'KKT stationarity'];

%% Primal and dual feasibility

% con_pd
if n_eq
    cons = [cons; (Aeq*x == beq):'Aeq x = beq'];
end

if n_ineq
    cons = [cons; (A*x <= b):'A x <= b'; (mu>=0):'mu>=0'];
end


%% Complementary slackness

%TODO
%mu.*(Ax - b) = 0
%(Ax mu(j) - b mu(j)) = 0

if n_ineq
    bmu = b .* mu;

    comp_slack = -bmu;
    for j = 1:n_ineq
        A_curr = A(j, :);
        xmu_curr = M(x_index, mu_index(j));

        comp_slack(j) = comp_slack(j) + A_curr*xmu_curr;
    end

    cons = [cons; (comp_slack==0):'comp slack mu(Ax-b)=0'];
%     xmu = M(mu_index, x_index);

%     Axmu = A*xmu;
    
end




%% Valid Inequalities (alpha)
%these redundant constraints arise from the preordering on inequality and
%equality constraints, and other properties of symmetric polynomials over
%the probability simplex

con_alpha_diag_sum = (sum(M(alpha_index, alpha_index), 2) == alpha);

con_alpha_diag_cmp = diag(M(alpha_index, alpha_index)) <= alpha;
cons = [cons; con_alpha_diag_cmp:'alpha(i)^2 <= alpha(i)'; ...
    con_alpha_diag_sum:'alpha(i) sum(alpha) = alpha(i)'];
% con_alpha2 = [];


%mixed term alpha(i)alpha(j) <= [alpha(i); alpha(j); 0.25],
%alpha(i)alpha(j) >= 0
con_alpha_mix_cmp_lhs = [];
con_alpha_mix_cmp_rhs = [];

% con_alpha_nonneg_off;
for i = 1:(m)
    curr_i = alpha_index(i);
    %off diagonal
    for j = 1:(i-1)        
        curr_j = alpha_index(j);
        con_alpha_mix_cmp_lhs = [con_alpha_mix_cmp_lhs; [1;1;1]*M(curr_i, curr_j)];
        con_alpha_mix_cmp_rhs = [con_alpha_mix_cmp_rhs;  [alpha(i); alpha(j); 0.25]];
    end       
end
cons = [cons; (con_alpha_mix_cmp_lhs <= con_alpha_mix_cmp_rhs):'mixed alpha cmp'];

%% Valid Inequalities (nonnegativity, alpha and mu only)

%each alpha and mu is nonnegative, so their mixed products are also
%nonnegative

M_amu = M([alpha_index, mu_index], [alpha_index, mu_index]);

con_nonneg_mix = zeros((m+n_ineq)*(m+n_ineq+1)/2, 1, 'like', sdpvar);
count_mix = 1;
for i = 1:(m + n_ineq)
    for j = 1:i
        con_nonneg_mix(count_mix) = M_amu(i, j);
        count_mix = count_mix+1;
    end
end

cons = [cons; (con_nonneg_mix>=0):'alpha & mu nonneg'];

%% Valid inequalities (mu alpha smaller than mu)

M_amu2 = M(alpha_index, mu_index);
con_mua_lt_mu = zeros(m*n_ineq, 1, 'like', sdpvar);
count_mix = 1;
for i = 1:m
    for j = 1:n_ineq
        con_mua_lt_mu(count_mix) = mu(j) - M_amu2(i, j);
        count_mix = count_mix+1;
    end
end

cons = [cons; (con_mua_lt_mu>=0): 'alpha * mu less than mu'];

%% Norm of mu ^2
maxval = 2;

cons_mu_lt_thresh = zeros(n_ineq, 1, 'like', sdpvar);
count_mix = 1;
for ii = 1 : n_ineq
    cons_mu_lt_thresh(count_mix) = maxval - mu(ii);
    count_mix = count_mix+1;
end

cons = [cons; (cons_mu_lt_thresh>=0): 'mu bounded'];
%% Norm of mu ^2
maxval = maxval^2;
M_mu2 = M(mu_index, mu_index);

cons_mu2_lt_thresh = zeros(n_ineq * (n_ineq+1)/2, 1, 'like', sdpvar);
count_mix = 1;
for ii = 1 : n_ineq
    for jj = ii : n_ineq
        cons_mu2_lt_thresh(count_mix) = maxval - M_mu2(ii, jj);
        count_mix = count_mix+1;
    end
end

cons = [cons; (cons_mu2_lt_thresh>=0): 'mu squared bounded'];

%% get the objectives
dist = sum(M(1,1)*y.^2 - 2*(x.*y) + x2);
cons = [cons; (dist>=0):'distance is nonnegative'];

objective = dist;

%dualize the yalmip model to get the problem into primal form (constraints
%in moment matrix and scalar inequalities)


if DUAL
    [Fd, objd] = dualize(cons, objective);
    objd = -objd;
else
    Fd = cons;
    objd = objective;
end
indexer = struct('x', x_index, 'a', alpha_index, 'mu', mu_index, 'const', 1, 'M', M, 'lambda', lambda, 'dist', dist, 'nm', size(M, 1));

% [model,recoverdata,diagnostic,interfacedata] = export(Fd, -objd, opts);

% SDP
% [model,recoverdata,diagnostic,interfacedata] = export(cons, dist, opts);
end