function out = checkFeasibility(problem,opts)
%CHECKFEASIBILITY Analytic rejection followed by a maximum-positive-margin LP.
out = struct('feasible',false,'status','not_checked','powerW',[], ...
    'scaledPoint',[],'positiveMargin',NaN,'exitflag',NaN,'message','');
if any(problem.bqos<0)
    out.status = 'noise_alone_violates_qos';
    out.message = 'At least one communication link has a negative interference allowance.';
    return;
end
if problem.q>1
    out.status = 'all_pair_hearability_impossible';
    out.message = 'q>1 makes mutually positive all-pair received powers impossible.';
    return;
end
n = numel(problem.c);
Alp = [problem.scaledA,zeros(size(problem.scaledA,1),1);-eye(n),ones(n,1)];
blp = [problem.scaledB;zeros(n,1)];
options = optimoptions('linprog','Display','none');
[x,~,flag,output] = linprog([zeros(n,1);-1],Alp,blp,[],[], ...
    [problem.scaledLb;0],[problem.scaledUb;max(problem.scaledUb)],options);
out.exitflag = flag;
out.message = output.message;
if flag<=0 || isempty(x)
    out.status = 'linear_constraints_infeasible';
    return;
end
z = x(1:n);
residual = max([0;problem.scaledA*z-problem.scaledB;problem.scaledLb-z]);
out.positiveMargin = x(end);
if residual>opts.feasibilityTolerance
    out.status = 'linear_solver_residual_failed';
    return;
end
% Homogeneous hearability permits scaling toward zero; bounds still respected.
interior = max(problem.scaledLb,z/2);
if max(problem.scaledA*interior-problem.scaledB)<=opts.feasibilityTolerance
    z = interior;
end
out.feasible = true;
out.status = 'feasible_positive_domain_with_numeric_lower_bound';
out.scaledPoint = z;
out.powerW = reshape(z*problem.powerScaleW,problem.K,problem.M);
end
