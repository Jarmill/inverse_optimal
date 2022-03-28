function [func_bound, selector, vars] = cons_local_yalmip_model(y, Q, f, A, b, Aeq, beq)
%CONS_YALMIP_MODEL generate yalmip models for the constrained projection
%problem
%   Outputs:
%       func_bound: minimize the quadratic sum_j alphaj (x'Qj x/2  + fj'x)
%       selector:   minimize norm(x-y)^2 when x has a cost from func_bound
%       vars:       store sdpvars for x and alpha

n = length(y);
m = length(f);

n_ineq = length(b);
n_eq = length(beq);

x = sdpvar(n, 1);
alpha = sdpvar(m, 1);
func_min = sdpvar(1, 1);

vars = struct('x', x, 'alpha', alpha, 'func_min', func_min);

% X = [A*x <= b; Aeq*x == beq];
X = [];
if n_ineq
    X = [X; A*x <= b];
end
if n_eq
    X = [X; Aeq*x == beq];
end

obj_func = zeros(1,1,'like', sdpvar);
for j = 1:m
    obj_func = obj_func + alpha(j)*(0.5*x'*Q{j}*x + f{j}'*x);
end

obj_selector = norm(x-y)^2;


X_selector = [X; obj_func <= func_min];

func_bound = struct('cons', X, 'obj', obj_func);
selector = struct('cons', X_selector, 'obj', obj_selector);

end

