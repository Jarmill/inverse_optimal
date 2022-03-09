load('opt_5_2d.mat')


%% test if a point is inside the feasible set
% this is a linear programming problem

y_in  = [-1.5; 1];
% y_corner = [3.96153846153846;-2.03846153846154];
y_corner = x_star{1};
% y_corner = [3.1;-2.0384];
y_out = [2; 2];
% y_test = [1.3086; 1.0162];
% y_test = [1.3; 1.0];
% y_test = [-1.2765; 0.2654];
% y_test = [1.76223692633813;0.679691280061691];

[alpha_rec, opt] = GO_contain_full(y_in, Q, x_star);
% [alpha_rec, opt] = GO_contain_full(y_corner, Q, x_star);
% [alpha_rec, opt] = GO_contain(y_out, Q, x_star);
% [alpha_rec, opt] = GO_contain(y_test, Q, x_star);

x_rec =  x_opt_uncons(alpha_rec, Q, x_star);

disp(alpha_rec)

function [alpha_rec, opt] = GO_contain_full(y,  Q, x_star)
%GO_contain: does the set of global optima associated with the quadratic
%functions (Q, x_star) contasin the point y? (convex combinations of
%objectives)
%
%Solve this problem through linear programming.

n = length(x_star{1});
m = length(x_star);


alpha = sdpvar(m, 1);

grad_term = zeros(n, 1, 'like', sdpvar);

Q_accum = zeros(n, n,'like', sdpvar);
x_accum = zeros(n, 1,'like', sdpvar);
obj_term = zeros(1,1,'like', sdpvar);
% m_grad_term = Q{m}*(y - x_star{m});
% grad_term = m_grad_term;
for j = 1:m
    curr_grad_term = alpha(j) * Q{j}*(y - x_star{j});
    
    %for debugging purposes (why does x_rec ~= y when y is feasible?)
    obj_term = obj_term + (y - x_star{j})'*curr_grad_term/2;
    
    Q_accum = Q_accum + alpha(j) * Q{j};
    x_accum = x_accum + alpha(j) * Q{j} * x_star{j};
    
    grad_term = grad_term + curr_grad_term;
end
% grad_term = grad_term + m_grad_term;

cons = [grad_term == 0; alpha >= 0; sum(alpha)==1];

opts = sdpsettings('verbose', false);

objective = 0;

opt = optimize(cons, objective, opts);

if opt.problem == 0
    alpha_rec = value(alpha);
    x_rec = value(Q_accum) \ value(x_accum);
else
    alpha_rec = nan(m, 1);
end


end

function [alpha_rec, opt] = GO_contain(y,  Q, x_star)
%GO_contain: does the set of global optima associated with the quadratic
%functions (Q, x_star) contasin the point y? (convex combinations of
%objectives)
%
%Solve this problem through linear programming.

n = length(x_star{1});
m = length(x_star);


alpha = sdpvar(m-1, 1);

grad_term = zeros(n, 1, 'like', sdpvar);
m_grad_term = Q{m}*(y - x_star{m});
% grad_term = m_grad_term;
for j = 1:m-1
    curr_grad_term = alpha(j) * (Q{j}*(y - x_star{j}) - m_grad_term);
    grad_term = grad_term + curr_grad_term;    
end
grad_term = grad_term + m_grad_term;

cons = [grad_term == 0; alpha >= 0];

opts = sdpsettings('verbose', false);

objective = -alpha(2);

opt = optimize(cons, objective, opts);

if opt.problem == 0
    alpha_rec = [value(alpha); 1-sum(value(alpha))];
else
    alpha_rec = nan(m, 1);
end


end