function hess_z = hproj_ehess_hadamard(z, Q, x_star, y)
%HPROJ_HESS_HADAMARD hessian of the GO function in the hadamard
%parameterization of the simplex (z.*z = alpha with z in a sphere)

alpha = z.*z;
grad = hproj_grad(alpha, Q, x_star);
hess = hproj_hess(alpha, Q, x_star, y);

%https://arxiv.org/pdf/2112.05273.pdf

hess_z = 2*diag(grad) + 4*diag(z)*hess*diag(z);

%let spherefactory.ehess2rhess() take care of the rest of the conversion
%routines to get the Riemannian hessian

% hess_perp = hess;
end

