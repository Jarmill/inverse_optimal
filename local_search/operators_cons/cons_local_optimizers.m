function [P_func,P_selector] = cons_local_optimizers(y, Q, f, A, b, Aeq, beq)
%CONS_LOCAL_OPTIMIZERS Create yalmip optimizer objects for the constrained
%projection problem
%   Detailed explanation goes here



[func_bound, selector, vars] = cons_local_yalmip_model(y, Q, f, A, b, Aeq, beq);

opts = sdpsettings('solver', 'mosek');


P_func = optimizer(func_bound.cons, func_bound.obj, opts, vars.alpha, [func_bound.obj; vars.x]);

P_selector = optimizer(selector.cons, selector.obj, opts, [vars.alpha; vars.func_min], [selector.obj; vars.x]);

% P = optimizer(Con,Obj,Options,Parameters,WantedVariables)

end

