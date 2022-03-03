load('opt_3_2d.mat')


%% test if a point is inside the feasible set
% this is a linear programming problem

y_in  = [-1.5; 1];
y_out = [2; 2];
% y_test = [1.3086; 1.0162];
y_test = [1.3; 1.0];

% [alpha_rec, opt] = GO_contain(y_in, Q, x_star);
% [alpha_rec, opt] = GO_contain(y_out, Q, x_star);
[alpha_rec, opt] = GO_contain(y_test, Q, x_star);

disp(alpha_rec)

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