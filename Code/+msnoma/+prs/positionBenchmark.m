function out = positionBenchmark(scene,ch,num,averagePowerW,noisePsdWPerHz,spec,geometrySpec)
%POSITIONBENCHMARK Average-power-matched broadcast PRS, using shared geometry.
allocation=msnoma.prs.allocation(spec,num.BHz);
budget=msnoma.prs.powerBudget(averagePowerW,spec.measurementPeriodS, ...
    spec.peakLimitW,allocation,spec);
out.allocation=allocation;
out.budget=budget;
out.referencePowerPolicy=spec.referencePowerPolicy;
out.label='ASSUMED average-positioning-power-matched broadcast PRS CRB';
if ~all(budget.feasible)
    out.status='observation_energy_or_peak_limit_infeasible';
    out.position=[];
    return;
end
cn0=ch.gP.*budget.instantaneousPowerW/noisePsdWPerHz;
crb=msnoma.prs.delayCrb(allocation,cn0,budget.symbolCount,num.c0Mps);
geom=msnoma.geometry.weights(scene,true(num.K,num.M),geometrySpec);
out.crb=crb;
out.geometry=geom;
out.position=msnoma.geometry.positionAccuracy(crb.rangeVarianceM2,geom);
out.status='computed_assumed_AWGN_broadcast_bound';
out.detectionCoverageEvaluated=false;
out.observationCaveat='finite_PRS_measurement_vs_steady_state_MS_DLL_not_identical_estimators';
end
