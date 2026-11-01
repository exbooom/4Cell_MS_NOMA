function d = diagnostics(problem,solution)
%DIAGNOSTICS Recompute constraints in SI and scaled units independently.
p = solution.powerW(:);
z = p/problem.powerScaleW;
d.qosResidualWPerHz = problem.Aqos*p-problem.bqos;
d.budgetResidualW = problem.Bbudget*p-problem.budgetW;
d.hearabilityResidualW = problem.Chear*p;
d.scaledResidual = problem.scaledA*z-problem.scaledB;
d.maxNormalizedViolation = max([0;d.scaledResidual;problem.scaledLb-z;z-problem.scaledUb]);
d.maxQosViolationWPerHz = max([0;d.qosResidualWPerHz]);
d.maxBudgetViolationW = max([0;d.budgetResidualW]);
d.maxHearabilityViolationW = max([0;d.hearabilityResidualW]);
hearDenominator=max(problem.Chear,0)*p;
d.maxRelativeHearabilityViolation=max([0;d.hearabilityResidualW./max(hearDenominator,realmin)]);
d.objectiveM2 = sum(problem.c./p);
d.kktStationarity = NaN;
d.kktComplementarity = NaN;
d.minimumMultiplier = NaN;
if isfield(solution,'multipliers') && ~isempty(solution.multipliers)
    l = solution.multipliers;
    [~,gradient] = msnoma.optim.objective(z,problem);
    stationarity = gradient+problem.scaledA'*l.ineqlin-l.lower+l.upper;
    d.kktStationarity = norm(stationarity,inf)/max(1,norm(gradient,inf));
    complementarity = [l.ineqlin.*d.scaledResidual; ...
        l.lower.*(problem.scaledLb-z);l.upper.*(z-problem.scaledUb)];
    d.kktComplementarity = norm(complementarity,inf);
    d.minimumMultiplier = min([l.ineqlin;l.lower;l.upper]);
end
end
