%% testing script to ensure that the moment matrix indexing functions correctly

load('opt_3_2d.mat')

y_in  = [-1.5; 1];
y_out = [2; 2];

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


dist = sum(y_out - 2*(x.*y_out) + x2);

opts = sdpsettings('solver', 'mosek');

% objective = -alpha(2);

opt = optimize(cons, dist, opts)

M_rec = value(M);
Mx = value(M(x_index, x_index))