grad_mano = zeros(size(x_sdp));

for jj = 1 : size(x_sdp, 2)
    
   % Total gradient vector 
    total_grad = zeros(size(x_sdp, 1), 1);
    
    % Add each cost function
    for cf = 1 : nf
        total_grad = total_grad + alpha_mano(cf,jj) * Q{cf} * x_mano(:, jj) - alpha_mano(cf, jj) * Q{cf} * x_star{cf};
    end
    
    % Add to check sdp
    grad_mano(:, jj) = total_grad;
end

% Check that gradients are zero
all(grad_mano <= 1e-6, 'all')