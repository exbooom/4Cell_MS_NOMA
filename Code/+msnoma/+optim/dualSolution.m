function solution = dualSolution(problem,opts,state,algorithm)
%DUALSOLUTION Build and independently certify a dual-recovered power vector.
solution.algorithm=algorithm;
solution.powerW=reshape(state.z*problem.powerScaleW,problem.K,problem.M);
solution.objectiveM2=sum(problem.c./solution.powerW(:));
solution.exitflag=state.exitflag;
solution.solverOutput=state.output;
solution.history=state.history;
solution.historyColumns=state.historyColumns;
solution.dual=state;
solution.dualBoundM2=-state.negativeDual*problem.objectiveScaleM2;
solution.dualityGapM2=solution.objectiveM2-solution.dualBoundM2;
solution.relativeDualityGap=solution.dualityGapM2/max(solution.objectiveM2,realmin);
lambda=state.lambda;
[~,gradient]=msnoma.optim.objective(state.z,problem);
stationarity=gradient+problem.scaledA'*lambda;
lower=zeros(size(state.z)); upper=lower;
atLower=state.z<=problem.scaledLb*(1+1e-12);
atUpper=state.z>=problem.scaledUb*(1-1e-12);
lower(atLower)=max(0,stationarity(atLower));
upper(atUpper)=max(0,-stationarity(atUpper));
solution.multipliers=struct('ineqlin',lambda,'lower',lower,'upper',upper);
solution.diagnostics=msnoma.optim.diagnostics(problem,solution);
d=solution.diagnostics;
hearTolerance=opts.feasibilityTolerance;
if isfield(problem,'hearabilityComparisonTolerance')
    hearTolerance=problem.hearabilityComparisonTolerance;
end
if d.maxNormalizedViolation<=opts.feasibilityTolerance ...
        && d.maxRelativeHearabilityViolation<=hearTolerance ...
        && d.kktStationarity<=opts.kktTolerance ...
        && d.kktComplementarity<=opts.kktTolerance ...
        && abs(solution.relativeDualityGap)<=opts.dual.relativeGapTolerance ...
        && d.minimumMultiplier>=0
    solution.status='converged_verified';
elseif d.maxNormalizedViolation<=opts.feasibilityTolerance
    solution.status='feasible_not_verified_optimal';
else
    solution.status='dual_recovery_infeasible';
end
end
