function result = estimateRanges(bsPosM,rangesM,valid,opts)
%ESTIMATERANGES Known-clock 2D range least squares, without a truth input.
% Three noncollinear anchors are required for an unambiguous first fix.
validateattributes(bsPosM,{'numeric'},{'2d','ncols',2,'real','finite'});
rangesM=rangesM(:); valid=logical(valid(:));
assert(numel(rangesM)==size(bsPosM,1) && numel(valid)==numel(rangesM));
ids=find(valid & isfinite(rangesM) & rangesM>=0);
result=struct('status','insufficient_ranges','valid',false,'positionM',[NaN,NaN], ...
    'usedBs',ids,'iterations',0,'residualM',nan(numel(ids),1), ...
    'residualRmsM',NaN,'conditionNumber',Inf,'source','estimated_ranges_no_truth');
if numel(ids)<3, return; end
b=bsPosM(ids,:); r=rangesM(ids);
A=2*(b(2:end,:)-b(1,:));
if rank(A,opts.rankTolerance*max(1,norm(A,2)))<2
    result.status='anchor_rank_deficient'; return;
end
rhs=r(1)^2-r(2:end).^2+sum(b(2:end,:).^2,2)-sum(b(1,:).^2);
x=(A\rhs)'; damping=1e-6; converged=false;
for iteration=1:opts.maxIterations
    offset=x-b; d=sqrt(sum(offset.^2,2));
    if any(d<1e-10)
        x=x+[1e-6,-1e-6]; offset=x-b; d=sqrt(sum(offset.^2,2));
    end
    J=offset./d; residual=d-r; gradient=J'*residual;
    if norm(gradient,inf)<=opts.gradientToleranceM
        converged=true; break;
    end
    accepted=false;
    for trial=1:20
        step=-(J'*J+damping*eye(2))\gradient;
        next=x+step'; nextResidual=sqrt(sum((next-b).^2,2))-r;
        if sum(nextResidual.^2)<sum(residual.^2)
            accepted=true; x=next; damping=max(1e-12,damping/3); break;
        end
        damping=damping*10;
    end
    if ~accepted, break; end
    if norm(step)<=opts.positionToleranceM
        converged=true; break;
    end
end
offset=x-b; d=sqrt(sum(offset.^2,2));
if any(d==0), result.status='zero_anchor_distance'; return; end
J=offset./d; singular=svd(J,0);
result.iterations=iteration; result.positionM=x;
result.residualM=d-r; result.residualRmsM=sqrt(mean(result.residualM.^2));
result.conditionNumber=singular(1)/singular(end);
if singular(end)<=opts.rankTolerance*singular(1)
    result.status='solution_rank_deficient';
elseif ~converged
    result.status='iteration_or_line_search_failure';
elseif result.residualRmsM>opts.maximumResidualRmsM
    result.status='range_residual_exceeded';
elseif any(x'<opts.domainM(:,1) | x'>opts.domainM(:,2))
    result.status='outside_known_domain';
else
    result.status='converged'; result.valid=true;
end
end
