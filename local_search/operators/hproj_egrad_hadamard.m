function grad_perp = hproj_egrad_hadamard(z, Q, x_star, y, perp_proj)
%HPROJ_EGRAD_HADAMARD gradient of the GO function in the hadamard
%parameterization of the simplex (z.*z = alpha with z in a sphere)

if nargin == 4
    perp_proj = false;
end
alpha = z.*z;
grad_perp_alpha = hproj_egrad(alpha, Q, x_star, y, perp_proj);

grad_z = 2*z.*grad_perp_alpha;

grad_perp = grad_z;

%should manopt already be performing this projection onto the sphere
%tangent space? That should likely be occuring in the SphereFactory routine
% grad_perp = grad_z - (z'*grad_z)*z;

end

