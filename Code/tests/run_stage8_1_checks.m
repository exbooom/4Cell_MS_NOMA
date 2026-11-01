function verification = run_stage8_1_checks()
%RUN_STAGE8_1_CHECKS Legacy negative control and residual-gated paired mains.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
reports={run_tests(),run_precision_tests(),run_dual_tests(),run_subgradient_tests(), ...
    run_prs_comparison_tests(),run_bias_tests(),run_sampled_tests(), ...
    run_cancellation_tests(),run_robust_cancellation_tests()};
files=dir(fullfile(root,'**','*.m')); sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep]), continue; end
    diagnostics=checkcode(path,'-id');
    if ~isempty(diagnostics), disp(path); disp(diagnostics); end
    assert(isempty(diagnostics),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,6); metrics=cell(1,6); acquisitions=cell(1,6); feedback=cell(1,6);
for b=1:6
    c=msnoma.config.base; c.io.showFigures=false;
    c.models.sampledReceiver.enabled=true; c.experiment.positionFeedback.enabled=true;
    c.models.sampledReceiver.interferenceCancellation.enabled=b~=1;
    c.experiment.name='stage8_1_residual_dominance_policy';
    if b<=3
        c.paper.requestedBandwidthHz=20e6; c.paper.positioningBudgetW=1e-3;
        c.numerics.powerScaleW=c.paper.positioningBudgetW/c.paper.M;
    elseif b==5
        c.paper.requestedBandwidthHz=20e6;
    elseif b==6
        c.models.sampledReceiver.includeCommunication=false;
        c.models.sampledReceiver.noiseMultiplier=0;
        c.models.sampledReceiver.codeMode='seeded_independent_bpsk';
    end
    if b==1, c.experiment.name='stage8_1_low_budget_original_receiver'; end
    if b==2
        c.models.sampledReceiver.interferenceCancellation.maximumRelativeResidualEnergy=1;
        c.experiment.name='stage8_1_low_budget_legacy_negative_control';
    end
    folder=main(c); loaded=load(fullfile(folder,'run.mat'),'run'); r=loaded.run;
    s=r.sampledReceiver; f=r.positionFeedback;
    assert(strcmp(r.schemaVersion,'8.1'));
    central=r.algorithms(strcmp({r.algorithms.id},'central'));
    for name={'central','dual','cp4a_consistent_candidate','cp4a_subgradient_candidate'}
        a=r.algorithms(strcmp({r.algorithms.id},name{1})); sol=a.solution;
        assert(strcmp(sol.status,'converged_verified') && all(a.evaluation.valid));
        assert(a.evaluation.communicationQosSatisfied);
        assert(abs(sol.objectiveM2/central.solution.objectiveM2-1)<1e-6);
        assert(sol.diagnostics.maxNormalizedViolation<=c.numerics.feasibilityTolerance);
        assert(sol.diagnostics.maxRelativeHearabilityViolation<=c.models.hearability.comparisonTolerance);
        assert(sol.diagnostics.kktStationarity<=c.numerics.kktTolerance);
        assert(sol.diagnostics.kktComplementarity<=c.numerics.kktTolerance);
        if ~strcmp(name{1},'central')
            assert(~sol.usedCentralFallback && abs(sol.relativeDualityGap)<=c.numerics.dual.relativeGapTolerance);
        end
    end
    assert(~s.truthUsedByEstimator && ~s.communicationDataUsedByEstimator);
    assert(~s.estimatedCancellationUncertaintyPropagated && ~s.steadyStateDllImplemented);
    assert(~s.fullyMatched && ~r.prsComparison.fullyMatched);
    assert(all(isnan(s.positionErrorM(~s.validPosition))));
    eligible=0; residualRejected=0;
    if b>1
        t=readtable(fullfile(folder,'sampled_estimated_P_channels.csv'));
        assert(height(t)==80);
        for u=1:20
            fit=s.estimatedPChannels{u};
            assert(~fit.truthUsed && ~fit.communicationDataUsed);
            expected=fit.valid && fit.precisionGatePassed && fit.residualGatePassed;
            assert(expected==fit.cancellationEligible);
            assert(all(cellfun(@(x) x.pCancellationApplied,s.receivers(:,u)))==expected);
        end
        eligible=nnz(cellfun(@(x) x.cancellationEligible,s.estimatedPChannels));
        residualRejected=nnz(cellfun(@(x) x.valid && x.precisionGatePassed ...
            && ~x.residualGatePassed,s.estimatedPChannels));
    end
    if b==2 || b==3
        old=acquisitions{1}; assert(isequal(s.powerW,old.powerW));
        assert(isequal(s.packet.chips,old.packet.chips));
        assert(isequal(s.packet.communicationSymbols,old.packet.communicationSymbols));
        assert(s.options.noiseSeed==old.options.noiseSeed);
    end
    if b==3
        assert(residualRejected>0 && eligible==0);
        assert(isequaln(s.rangeM,acquisitions{1}.rangeM));
        assert(isequaln(s.estimatedPositionM,acquisitions{1}.estimatedPositionM));
        assert(isequal(s.validRanges,acquisitions{1}.validRanges));
        assert(isequaln(f.history(1).acquisition.rangeM,feedback{1}.history(1).acquisition.rangeM));
    end
    if b>=4, assert(all(s.validRanges,'all') && all(s.validPosition)); end
    if b==4 || b==6, assert(f.converged); end
    assert(f.firstFixUnitLambda && ~f.usedOracleFallback && ~f.truthGeometryUsed);
    assert(all(f.history(1).lambda==1,'all'));
    for j=2:f.iterations
        previous=f.history(j-1).acquisition;
        scene=struct('bsPosM',r.scene.bsPosM,'pUserPosM',previous.estimatedPositionM);
        spec=c.models.geometry; spec.source='estimated_sampled_ranges';
        g=msnoma.geometry.weights(scene,true(4,20),spec);
        assert(isequal(g.lambda,f.history(j).lambda));
    end
    figures=13+double(b>1); assert(numel(dir(fullfile(folder,'figure_*.png')))==figures);
    last=f.history(end).acquisition;
    metrics{b}=struct('name',c.experiment.name,'bandwidthHz',r.num.BHz, ...
        'budgetW',c.paper.positioningBudgetW,'residualLimit', ...
        c.models.sampledReceiver.interferenceCancellation.maximumRelativeResidualEnergy, ...
        'validRanges',nnz(s.validRanges),'totalRanges',80, ...
        'validPositions',nnz(s.validPosition),'totalPositions',20, ...
        'conditionalMeanActualErrorM',mean(s.positionErrorM,'omitnan'), ...
        'eligibleUsers',eligible,'residualRejectedAfterPrecisionPassed',residualRejected, ...
        'feedbackStatus',f.status,'feedbackIterations',f.iterations, ...
        'feedbackFinalValidPositions',nnz(last.validPosition), ...
        'feedbackFinalConditionalMeanErrorM',mean(last.positionErrorM,'omitnan'),'figures',figures);
    runs{b}=folder; acquisitions{b}=s; feedback{b}=f;
    if b==3, sensitivity=checkSensitivity(r); end
    fprintf('Verified stage8.1 main %d: %s\n',b,folder);
end
% Keep and check the negative control rather than silently replacing it.
assert(nnz(acquisitions{2}.validRanges)<nnz(acquisitions{1}.validRanges));
hashes=msnoma.report.sourceHashes(root);
for b=1:6
    manifest=jsondecode(fileread(fullfile(runs{b},'manifest.json')));
    assert(isequal({hashes.file},{manifest.sourceHashes.file}));
    assert(isequal({hashes.sha256},{manifest.sourceHashes.sha256}));
end
counts=cellfun(@height,reports);
verification=struct('date','2026-11-01','schemaVersion','8.1','testCounts',counts, ...
    'totalTests',sum(counts),'sourceFiles',sourceCount,'codeAnalyzerDiagnostics',0, ...
    'sourceHashesMatchedAllRuns',true,'runs',{runs},'scenarios',{metrics}, ...
    'lowBudgetThresholdSensitivity',{sensitivity},'accuracyGuarantee',false, ...
    'cancellationUncertaintyPropagated',false);
assert(isscalar(verification) && numel(verification.lowBudgetThresholdSensitivity)==3);
save(fullfile(root,'docs','STAGE8_1_VERIFICATION.mat'),'verification','reports','hashes');
fid=fopen(fullfile(root,'docs','STAGE8_1_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true));
disp(verification);
end

function rows = checkSensitivity(r)
% Reuse the same samples/powers, not true errors, to decide cancellation.
s=r.sampledReceiver; n=r.num; limits=[.25,.5,.75]; rows=cell(1,numel(limits));
o=s.options; packet=msnoma.waveform.makePacket(n,o);
for i=1:numel(limits)
    o.interferenceCancellation.maximumRelativeResidualEnergy=limits(i);
    trial=msnoma.experiment.sampledReceiver(r.scene,r.ch,n,s.powerW,o, ...
        r.cfg.assumed.noisePsdWPerHz,packet);
    assert(isequaln(trial.rangeM,s.rangeM));
    rows{i}=struct('limit',limits(i),'eligibleUsers', ...
        nnz(cellfun(@(x) x.cancellationEligible,trial.estimatedPChannels)), ...
        'rangesExactlyMatchDefault',true,'validRanges',nnz(trial.validRanges), ...
        'validPositions',nnz(trial.validPosition));
end
end
