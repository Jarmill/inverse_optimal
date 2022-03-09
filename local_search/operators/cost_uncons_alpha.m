function [f] = cost_uncons_alpha(alpha, Q, x_star, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x_opt = x_opt_uncons(alpha, Q, x_star);
 
f = 0.5*norm(x_opt - y)^2;


end

