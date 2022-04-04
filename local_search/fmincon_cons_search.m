function [x_fmincon,dist_fmincon,alpha_con, exitflag, output] = fmincon_cons_search(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp)
%FMINCON_CONS_SEARCH Perform a local search of the constrained ProjGM
%function using fmincon on alpha
%   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

%define the operators
    n = length(y);
    m = length(phi);
    [P_func,P_selector] = cons_local_optimizers(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);

    %% fmincon command
    Aeq= ones(1, m);
    beq = [1];

    A = -speye(m);
    b = sparse(m, 1);

    options = optimset;
    
    alpha0 = ones(m, 1)/m; %initial point, likely not the greatest choice

    [alpha_con, fval, exitflag, output, lambda, grad, hessian]...
        = fmincon(@(alpha) f_in(alpha, P_func, y), alpha0, A, b, Aeq, beq, [], [], [], options);

    dist_fmincon = sqrt(fval);
    % Out data
    x_fmincon = P_func(alpha_con);
    x_fmincon = x_fmincon(2:end);
    
    function f_out = f_in(alpha, P_func, y)
        %cost function to minimize the mixed quadratic
        funcdata = P_func(alpha);
        % f_out = xdata(1);
        f_out = norm(y - funcdata(2:end))^2;

    end

end

