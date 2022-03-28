grad_sdp = zeros(size(x_sdp));

for jj = 1 : size(x_sdp, 2)
    
   % Total gradient vector 
    total_grad = zeros(size(x_sdp, 1), 1);
    
    % Add each cost function
    for cf = 1 : nf
        total_grad = total_grad + alpha_sdp(cf,jj) * Q{cf} * x_sdp(:, jj) - alpha_sdp(cf, jj) * Q{cf} * x_star{cf};
    end
    
    % Add to check sdp
    grad_sdp(:, jj) = total_grad;
end

% Check
all(grad_sdp <= 1e-6, 'all')


% Check that it is a lower bound
all(d_sdp <= d_bf, 'all')
all(d_sdp <= d_lsq, 'all')
all(d_sdp <= d_qcqp, 'all')
all(d_sdp <= d_bflsq, 'all')
all(d_sdp <= d_bfqcqp, 'all')
all(d_sdp <= d_mano, 'all')