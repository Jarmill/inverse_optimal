%% formulate the constrained optimum sdp
Q = {[1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1]};
f = {[0;0], [0;1], [1; 1], [0.5; 0]};

% y = [2; 0.5];
y = [1; -0.5];


EQ = false;
INEQ = true;

if EQ
%     Aeq=[1 0];
%     beq=[0.5];

    Aeq = [1 1];
    beq = [1];

else
    Aeq = [];
    beq = [];
end


if INEQ
    A=[1 0; 1 1];
    b=[0.6; 1.5];
else
    A=[];
    b=[];
end

%% generate the yalmip model

[out_rec, info] = manopt_search_cons(y, Q, f, A, b, Aeq, beq)

% [func_bound, selector, vars] = cons_local_yalmip_model(y, Q, f, A, b, Aeq, beq);

% [P_func,P_selector] = cons_local_optimizers(y, Q, f, A, b, Aeq, beq);
% 
% 
% alpha_test = [0.2; 0.3; 0; 0.5];
%  [dist_rec, obj_rec, x_rec] = cost_cons(alpha_test, P_func, P_selector);