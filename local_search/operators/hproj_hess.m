function hess = hproj_hess(alpha, Q, x_star, y)
%HPROJ_HESS Hessian of the uncons GO distance function
%   may be too complicated and expensive to implement

%use matrixcalculus.org to implement this
%this is the grad expression for a1:
%grad = 2*(inv(a1*Q1+a2*Q2+a3*Q3)*(a1*Q1*x1+a2*Q2*x2+a3*Q3*x3)-y)'*inv(a1*Q1+a2*Q2+a3*Q3)*Q1*x1-2*(inv(a1*Q1+a2*Q2+a3*Q3)*(a1*Q1*x1+a2*Q2*x2+a3*Q3*x3)-y)'*inv(a1*Q1+a2*Q2+a3*Q3)*Q1*inv(a1*Q1+a2*Q2+a3*Q3)*(a1*Q1*x1+a2*Q2*x2+a3*Q3*x3) 
%
%dgrad/da1 : main diagonal
%dgrad/da2 : off diagonal
%dgrad/da3 : elements not involved in the hessian diferentiation
%(equivalent to dgrad/da2)

hess = 0;
end

