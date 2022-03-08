function [x_star] = x_opt(Q, x_star, alpha)
%X_OPT Get the global optimum solution of convex quadratics
m = length(Q);
Qall = 0;
hall = 0;
for i = 1:m
    Qall = Qall + alpha(i)*Q{i};
    hall = hall + alpha(i)*Q{i}*x_star{i};
end

x_star = Qall \ hall;
end

