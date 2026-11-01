function [value,gradient,z,hessian] = dualValue(lambda,problem,rows,fixed)
%DUALVALUE Negative Lagrange dual on the explicit positive power box.
% Complete coupling: price=A'*lambda, including all incoming constraints.
if nargin<3, rows=(1:size(problem.scaledA,1))'; end
if nargin<4, fixed=zeros(size(problem.scaledB)); end
allLambda=fixed;
allLambda(rows)=lambda;
assert(all(allLambda>=0) && all(isfinite(allLambda)));
price=problem.scaledA'*allLambda;
c=problem.scaledC;
lo=problem.scaledLb;
hi=problem.scaledUb;
z=hi;
positive=price>0 & c>0;
z(positive)=sqrt(c(positive)./price(positive));
z=max(lo,min(hi,z));
% For c=0 the scalar objective is linear; choose a box endpoint exactly.
z(c==0 & price>=0)=lo(c==0 & price>=0);
value=problem.scaledB'*allLambda-sum(c./z+price.*z);
A=problem.scaledA(rows,:);
gradient=problem.scaledB(rows)-A*z;
if nargout>3
    curvature=zeros(size(z));
    free=c>0 & z>lo & z<hi;
    curvature(free)=z(free).^3./(2*c(free));
    hessian=(A.*curvature')*A';
    hessian=(hessian+hessian')/2;
end
end
