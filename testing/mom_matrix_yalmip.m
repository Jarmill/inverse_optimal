%% testing script to ensure that the moment matrix indexing functions correctly

load('opt_3_2d.mat')

y_in  = [-1.5; 1];
y_out = [1; 0];
y_out = y_in;
% y_out = [-1; 0];

n = length(x_star{1});
m = length(x_star);


%form and index out the moment matrix. remember that alpha_n is eliminated
x_index = 1 + (1:n);
alpha_index = (n+1) + (1:(m-1));

M = sdpvar(n+m);

%variables
x = M(x_index, 1);
alpha = [M(alpha_index, 1); 1 - sum(M(1, alpha_index))];
x2 = diag(M(x_index, x_index));

xa = M(x_index, alpha_index);

%start the constraints
cons = [(M>=0):'PSD moments', (M(1,1) == 1):'1=1'; (alpha >= 0):'alpha >= 0'];

%now do the gradient constraints
grad_term = Q{m}*(x - x_star{m});
for j = 1:m-1
    curr_term = (Q{j} - Q{m})*xa(:, j) - (Q{j}*x_star{j} - Q{m}*x_star{m})*alpha(j);
    grad_term = grad_term + curr_term;
end

cons = [cons; (grad_term==0):'KKT stationarity'];

%add redundant constraints to the entries of alpha^2
%alphai^2 <= 1
%alphai*alphaj <= 0.25 if i~=j
con_alpha2 = [];
con_alpha2_m = [];
% for i = 1:(m-1)
%     sum(
for i = 1:(m-1)
    curr_i = alpha_index(i);
    %on diagonal
%     con_alpha2 = [con_alpha2; M(curr_i, curr_i) <= 1];
    con_alpha2 = [con_alpha2; M(curr_i, curr_i) <= alpha(i)];
    
    m_term_off =alpha(i) - M(curr_i, alpha_index);
    con_alpha2_m = [con_alpha2_m; m_term_off <= 0.25; m_term_off >= 0];
    %off diagonal
    for j = 1:(i-1)        
        curr_j = alpha_index(j);
        con_alpha2 = [con_alpha2; M(curr_i, curr_j) <= 0.25; ...
            M(curr_i, curr_j) >= 0; M(curr_i, curr_j) <= alpha(i)];
    end
end
%alpham^2 <= 1
all_m_term = 1 - sum(alpha) + sum(sum(M(alpha_index, alpha_index)));
con_alpha2_m = [con_alpha2_m; all_m_term <= 1; all_m_term >= 0];
%now deal with the last alpha(m)

dist = sum(y_out.^2 - 2*(x.*y_out) + x2);

cons = [cons; con_alpha2; con_alpha2_m; dist>=0];




opts = sdpsettings('solver', 'mosek');

% objective = -alpha(2);

opt = optimize(cons, dist, opts)

M_rec = value(M)
x_rec = value(x);
alpha_rec = [value(alpha); 1-sum(value(alpha))];
Mx = value(M(1:(n+1), 1:(n+1)));
dist_rec = sqrt(value(dist))