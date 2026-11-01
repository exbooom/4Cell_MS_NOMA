function out = comparison(scene,ch,num,benchmark,reference,cfg,berModel)
%COMPARISON Audit partial matches; never claim a fully fair receiver comparison.
out.status='prs_observation_schedule_infeasible';
out.referenceAlgorithm=reference.id;
out.fullyMatched=false; out.observationMatched=false;
out.detectionCoverageEvaluated=false;
out.positioningBroadcastEnergyMatched=false;
out.geometryConventionMatched=false;
out.policies=struct('id',{},'crb',{},'tdoa',{},'service',{});
out.observationCaveat='finite_PRS_known_signal_CRB_vs_MS_steady_state_DLL_variance_no_time_conversion';
out.msCovarianceCaveat='independent_MS_TOA_error_surrogate_not_a_finite_observation_CRB';
if ~isfield(benchmark,'position') || isempty(benchmark.position), return; end
spec=cfg.models.prs.comparison;
mask=true(num.K,num.M);
out.commonBsMask=mask;
out.msTdoaCommon=msnoma.geometry.tdoaAccuracy(scene,reference.evaluation.rangeVarianceM2,spec.tdoa,mask);
out.msTdoaHearable=msnoma.geometry.tdoaAccuracy(scene,reference.evaluation.rangeVarianceM2, ...
    spec.tdoa,reference.evaluation.hearability.mask);
out.msMinimumBsForTdoa=3;
out.kernel=msnoma.prs.interferenceKernels(ch,num,benchmark.allocation,spec);
ids={'global_muting','other_bs_communication'};
for i=1:numel(ids)
    if i==1
        crb=benchmark.crb;
    else
        crb=msnoma.prs.delayCrbInterference(benchmark.allocation,ch.gP,benchmark.budget, ...
            cfg.assumed.noisePsdWPerHz,out.kernel,num.c0Mps);
    end
    tdoa=msnoma.geometry.tdoaAccuracy(scene,crb.rangeVarianceM2,spec.tdoa,mask);
    service=msnoma.prs.serviceAudit(ch,num,benchmark.budget,out.kernel, ...
        cfg.assumed.noisePsdWPerHz,berModel,cfg.paper.berThreshold,spec,ids{i},reference.evaluation.ber);
    out.policies(i)=struct('id',ids{i},'crb',crb,'tdoa',tdoa,'service',service);
end
out.positioningBroadcastEnergyMatched=true;
out.geometryConventionMatched=true;
out.status='computed_partial_match_service_and_observation_audit';
out.maskCaveat='all_BS_oracle_comparison_mask_not_PRS_detection_coverage';
out.clockCaveat='synchronized_BS_common_receiver_clock_eliminated_local_frozen_covariance';
out.interferenceCaveat='finite_OFDM_PSD_diagonal_Gaussian_noise_approximation_no_C_to_C';
out.receiverAudit=struct('status','disabled');
if spec.receiverAuditEnabled
    out.receiverAudit=msnoma.prs.receiverAudit(ch,num,benchmark.allocation,benchmark.budget, ...
        cfg.assumed.noisePsdWPerHz,berModel,cfg.paper.berThreshold,spec,out.kernel,reference.evaluation.ber);
end
end
