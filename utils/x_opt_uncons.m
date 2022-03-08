function [x_opt] = x_opt_uncons(alpha, Q, x_star)
%X_OPT_UNCONS Summary of this function goes here
%   Detailed explanation goes here
m = length(alpha);
Q_accum = 0;
x_accum = 0;
for j = 1:m
    Q_accum = Q_accum +  alpha(j) * Q{j};
    x_accum = x_accum + alpha(j) * x_star{j};
end

x_opt = Q_accum \ x_accum;

end

