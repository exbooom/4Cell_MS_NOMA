function validation = validate(cfg)
%VALIDATE Fail before expensive work when a requested model is unsupported.
validateattributes(cfg.paper.K,{'numeric'},{'scalar','integer','>=',2});
validateattributes(cfg.paper.M,{'numeric'},{'scalar','integer','positive'});
assert(isequal(size(cfg.paper.bsPosM),[cfg.paper.K,2]),'BS positions must be K-by-2.');
validateattributes(cfg.assumed.noisePsdWPerHz,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.assumed.minimumDistanceM,{'numeric'},{'scalar','positive'});
validateattributes(cfg.paper.positioningBudgetW,{'numeric'},{'scalar','positive'});
validateattributes(cfg.models.hearability.ratio,{'numeric'},{'scalar','nonnegative','finite'});
assert(cfg.models.geometry.dimension==2 && strcmp(cfg.models.geometry.source,'oracle'), ...
    'Only oracle 2D geometry is implemented in this stage.');
assert(strcmp(cfg.models.geometry.clockModel,'paper_no_clock_state'), ...
    'The optimization geometry must retain the no-clock paper convention; TDoA is a separate post-evaluation.');
assert(strcmp(cfg.assumed.communicationDistribution,'uniform_nearest_bs_region') ...
    && strcmp(cfg.assumed.positioningDistribution,'uniform_square'), ...
    'Requested user distribution is not implemented.');
assert(cfg.assumed.minimumDistanceM==cfg.assumed.propagation.minimumDistanceM, ...
    'Sampling and propagation minimum distances must agree.');
assert(strcmp(cfg.models.bandwidth,'integer_nearest'),'Unsupported bandwidth mode.');
assert(cfg.models.hearability.minimumBs>=2 && cfg.models.hearability.minimumBs<=cfg.paper.K);
assert(all(cfg.assumed.domainM(:,2)>cfg.assumed.domainM(:,1)));
assert(license('test','Optimization_Toolbox') && ~isempty(which('fmincon')) ...
    && ~isempty(which('linprog')),'Optimization Toolbox is required.');
assert(cfg.numerics.minimumPowerW>0 && ...
    cfg.numerics.minimumPowerW*cfg.paper.M<cfg.paper.positioningBudgetW);
if cfg.models.rangingIntegral.enabled || cfg.models.dllBias.enabled
    q=cfg.models.rangingIntegral;
    assert(any(strcmp(q.frontendMode,{'global_frontend','local_symmetric'})));
    assert(q.chunkSize>=1 && q.chunkSize==fix(q.chunkSize));
    assert(all(q.userIds>=1 & q.userIds<=cfg.paper.M & q.userIds==fix(q.userIds)));
end
if cfg.models.dllBias.enabled
    ids=cfg.models.dllBias.userIds;
    assert(all(ids>=1 & ids<=cfg.paper.M & ids==fix(ids)) && numel(ids)==numel(unique(ids)));
end
if cfg.models.sampledReceiver.enabled || cfg.experiment.positionFeedback.enabled ...
        || cfg.experiment.interferenceAblation.enabled
    r=cfg.models.sampledReceiver;
    validateattributes(r.oversampling,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(r.codeChips,{'numeric'},{'scalar','integer','>=',31});
    validateattributes(r.maximumRangeM,{'numeric'},{'scalar','positive','finite'});
    validateattributes(r.noiseMultiplier,{'numeric'},{'scalar','nonnegative','finite'});
    validateattributes(r.falseAlarmProbability,{'numeric'},{'scalar','>',0,'<',1});
    validateattributes(r.delayToleranceSamples,{'numeric'},{'scalar','positive','<',0.1});
    validateattributes(r.coarseOversampling,{'numeric'},{'scalar','integer','>=',1,'<=',8});
    for seed=[r.sequenceSeed,r.noiseSeed]
        validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<',2^32-1});
    end
    assert(r.carrierHz==cfg.paper.carrierHz);
    assert(any(strcmp(r.whitening,{'none','known_average_C_psd','known_average_CP_psd'})));
    assert(r.whiteningFloorFraction>0 && r.whiteningFloorFraction<1);
    assert(any(strcmp(r.codeMode,{'author_L1CP','seeded_independent_bpsk'})));
    assert(all(r.userIds>=1 & r.userIds<=cfg.paper.M & r.userIds==fix(r.userIds)) ...
        && numel(unique(r.userIds))==numel(r.userIds));
    assert(r.oversampling>=cfg.paper.frontendMultiplier);
    assert(isequal(r.position.domainM,cfg.assumed.domainM));
    validateattributes(r.position.maxIterations,{'numeric'},{'scalar','positive','integer'});
    validateattributes(r.position.maximumResidualRmsM,{'numeric'},{'scalar','positive'});
    a=r.interferenceCancellation;
    validateattributes(a.maxIterations,{'numeric'},{'scalar','positive','integer'});
    validateattributes(a.maxBacktracks,{'numeric'},{'scalar','positive','integer'});
    validateattributes(a.maximumStepSamples,{'numeric'},{'scalar','positive','finite'});
    validateattributes(a.delayToleranceSamples,{'numeric'},{'scalar','positive','finite'});
    validateattributes(a.rankTolerance,{'numeric'},{'scalar','positive','<',1});
    validateattributes(a.maximumDelayStdSamples,{'numeric'},{'scalar','positive','finite'});
    validateattributes(a.maximumRelativeGainStd,{'numeric'},{'scalar','positive','finite'});
    validateattributes(a.maximumRelativeResidualEnergy,{'numeric'},{'scalar','>',0,'<=',1,'finite'});
end
if cfg.experiment.positionFeedback.enabled
    f=cfg.experiment.positionFeedback;
    assert(any(strcmp(f.solverId,{'central','dual','cp4a_subgradient_candidate'})));
    validateattributes(f.maxFixes,{'numeric'},{'scalar','integer','>=',2});
    assert(f.minimumFixes>=2 && f.minimumFixes<=f.maxFixes);
    assert(f.positionToleranceM>0 && f.relativePowerTolerance>0);
end
if cfg.experiment.stepSensitivity.enabled
    validateattributes(cfg.experiment.stepSensitivity.multipliers,{'numeric'}, ...
        {'vector','positive','finite','nonempty'});
end
if cfg.experiment.interferenceAblation.enabled
    ids=cfg.experiment.interferenceAblation.userIds;
    validateattributes(ids,{'numeric'},{'vector','positive','integer','nonempty','<=',cfg.paper.M});
    assert(numel(unique(ids))==numel(ids));
end
if cfg.models.prs.enabled
    p=cfg.models.prs;
    assert(strcmp(p.referencePowerPolicy,'match_reference_BS_total_positioning_power'));
    assert(any(strcmp(p.mode,{'generalized_ofdm','lte_single_symbol_reference'})));
    assert(p.combSpacing>=2 && p.combSpacing==fix(p.combSpacing));
    assert(p.combOffset>=0 && p.combOffset<p.combSpacing && p.combOffset==fix(p.combOffset));
    assert(p.measurementPeriodS>0 && p.peakLimitW>0);
    if p.comparison.enabled
        s=p.comparison;
        assert(strcmp(cfg.models.ber.id,'qpsk_bit_energy'),'PRS service audit currently supports explicit QPSK.');
        validateattributes(s.chunkSize,{'numeric'},{'scalar','positive','integer'});
        validateattributes(s.minimumAttemptedRateRatio,{'numeric'},{'scalar','positive','<=',1});
        validateattributes(s.communicationPeakLimitW,{'numeric'},{'scalar','positive','finite'});
        assert(any(strcmp(s.tdoa.referencePolicy,{'minimum_variance','fixed'})));
        validateattributes(s.tdoa.rankTolerance,{'numeric'},{'scalar','positive','<',1});
        if strcmp(s.tdoa.referencePolicy,'fixed')
            validateattributes(s.tdoa.referenceBs,{'numeric'},{'scalar','integer','>=',1,'<=',cfg.paper.K});
        end
    end
end
validation = struct('canRun',true,'message','Explicitly assumed research profile validated.');
if any(ismember(cfg.experiment.algorithmIds,{'dual','cp4a_consistent_candidate','cp4a_subgradient_candidate'}))
    d=cfg.numerics.dual;
    assert(d.maxPasses>=1 && d.rowsPerPass>=1 && d.maxIterations>=1);
    assert(d.constraintTolerance>0 && d.optimalityTolerance>0 && d.relativeGapTolerance>0);
end
if any(strcmp(cfg.experiment.algorithmIds,'cp4a_consistent_candidate'))
    assert(strcmp(cfg.numerics.cp4a.mode,'full_coupling_hierarchical_candidate'));
end
if any(strcmp(cfg.experiment.algorithmIds,'cp4a_subgradient_candidate'))
    s=cfg.numerics.subgradient;
    assert(any(strcmp(s.stepRule,{'backtracking','diminishing'})));
    assert(s.initialStep>0 && s.outerInitialStep>0 && s.maximumStep>=s.initialStep);
    assert(s.diminishingExponent>0.5 && s.diminishingExponent<=1);
    assert(s.innerFeasibilityTolerance>0 && s.innerRelativeGapTolerance>0);
    assert(s.maxInnerIterations>=0 && s.maxOuterIterations>=0);
end
end
