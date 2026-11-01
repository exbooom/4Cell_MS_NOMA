function [f,gradient,hessian] = objective(z,problem)
%OBJECTIVE Scaled sum(c_i/p_i), with exact first and second derivatives.
c = problem.scaledC;
f = sum(c./z);
if nargout>1, gradient = -c./z.^2; end
if nargout>2, hessian = diag(2*c./z.^3); end
end
