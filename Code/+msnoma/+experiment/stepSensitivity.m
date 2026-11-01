function study = stepSensitivity(problem,cfg,feasibility)
%STEPSENSITIVITY Engineering initial step trials; full certificates decide.
o=cfg.experiment.stepSensitivity;
study=struct('status','disabled','authorStepsKnown',false,'options',o,'trials',[]);
if ~o.enabled, return; end
if ~feasibility.feasible, study.status=feasibility.status; return; end
trials=struct('multiplier',{},'innerInitialStep',{},'outerInitialStep',{},'solution',{});
for i=1:numel(o.multipliers)
    factor=o.multipliers(i); numerical=cfg.numerics;
    numerical.subgradient.initialStep=cfg.numerics.subgradient.initialStep*factor;
    numerical.subgradient.outerInitialStep=cfg.numerics.subgradient.outerInitialStep*factor;
    numerical.subgradient.diminishingInitialStep=cfg.numerics.subgradient.diminishingInitialStep*factor;
    numerical.subgradient.maximumStep=max(numerical.subgradient.maximumStep, ...
        numerical.subgradient.initialStep);
    solution=msnoma.optim.solveSubgradient(problem,numerical,feasibility);
    trials(end+1)=struct('multiplier',factor, ...
        'innerInitialStep',numerical.subgradient.initialStep, ...
        'outerInitialStep',numerical.subgradient.outerInitialStep,'solution',solution); %#ok<AGROW>
end
study.status='completed'; study.trials=trials;
study.interpretation='initial_steps_in_scaled_dual_coordinates_backtracking_changes_accepted_steps';
values=arrayfun(@(x) x.solution.objectiveM2,trials);
study.verifiedCount=nnz(arrayfun(@(x) strcmp(x.solution.status,'converged_verified'),trials));
study.objectiveRelativeSpread=NaN;
if all(isfinite(values)) && min(values)>0
    study.objectiveRelativeSpread=max(values)/min(values)-1;
end
end
