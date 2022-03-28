function [dist_rec, func_min, x_rec] = cost_cons_alpha(alpha, P_func, P_selector)
%COST_CONS Cost of the constrained selection problem
%   Detailed explanation goes here


%get the minimum value of the quadratic parameterized by alpha
out_func = P_func(alpha);

func_min = out_func(1);

%find the x point with a quadratic value of func_min that is as close to y
%as possible
out_selector = P_selector([alpha; func_min]);

dist_rec = out_selector(1);
x_rec = out_selector(2:end);

end

