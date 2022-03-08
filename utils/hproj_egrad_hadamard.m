function grad_perp = hproj_egrad_hadamard(z, Q, x_star, y, perp_proj)
%HPROJ_EGRAD_HADAMARD Summary of this function goes here
%   Detailed explanation goes here
if nargin == 4
    perp_proj = true;
end
alpha = z.*z;
grad_perp_alpha = hproj_egrad(alpha, Q, x_star, y, perp_proj);

grad_z = 2*z.*grad_perp_alpha;

grad_perp = grad_z - (z'*grad_z)*z;

end
