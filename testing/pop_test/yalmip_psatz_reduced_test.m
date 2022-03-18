load('opt_5_2d.mat');

% y = [2; 1];
y = [2; 3];

n = length(x_star{1});
m = length(x_star);


%% define variables
x = sdpvar(n, 1);
alpha = sdpvar(m, 1);

vars = [x; alpha];

gamma = sdpvar(1);

Q_accum = zeros(n, 1, 'like', sdpvar);
for j = 1:m
    Q_accum = Q_accum + alpha(j)*Q{j}*(x - x_star{j});
end


%reduced form
alpha_red = alpha(1:end-1);
vars_red = [x; alpha_red];
Q_accum_red = replace(Q_accum, alpha(end), 1 - sum(alpha_red));

% X = struct('ineq', [alpha], 'eq', [Q_accum; 1-sum(alpha)]);


%add in valid inequalities (preordering)


X_red = struct('ineq', [alpha_red; 1-sum(alpha_red)], 'eq', Q_accum_red);

dist2 = norm(y-x, 2)^2;
p = dist2 - gamma;



%% perform psatz
order = 2;
d = 2*order;


% [p_out, cons, Gram] = psatz(dist2-gamma, X, order, vars, n+m);

% [p_out, cons, coeff_list] = constraint_psatz(p, X, vars, d);

[p_out, cons, coeff_list] = constraint_psatz(p, X_red, vars_red, d);

objective = -gamma;

opts = sdpsettings('solver', 'mosek');

%% solve the program
% sol = optimize(cons, objective, opts);

[sol,u,Q] = solvesos(cons, objective, opts, [gamma; coeff_list]);

dist_rec = sqrt(value(gamma))

