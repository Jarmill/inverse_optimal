function v  = evaluate_quadratics(Q, x_star, x)

v = zeros(1, length(x));

for ii = 1 : length(x)
    v(ii) = 1/2 * (x{ii} - x_star{ii}).' * Q{ii} * (x{ii} - x_star{ii});
end

end