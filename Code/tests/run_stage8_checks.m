function verification = run_stage8_checks()
%RUN_STAGE8_CHECKS Paired original/estimated-P cancellation and full feedback.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
reports={run_tests(),run_precision_tests(),run_dual_tests(),run_subgradient_tests(), ...
    run_prs_comparison_tests(),run_bias_tests(),run_sampled_tests(),run_cancellation_tests()};
files=dir(fullfile(root,'**','*.m')); sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep]), continue; end
    diagnostics=checkcode(path,'-id');
    if ~isempty(diagnostics), disp(path); disp(diagnostics); end
    assert(isempty(diagnostics),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,7); metrics=cell(1,7); baseline=cell(1,3);
for b=1:7
    c=msnoma.config.base; c.io.showFigures=false;
    c.models.sampledReceiver.enabled=true;
    c.experiment.positionFeedback.enabled=true;
    scenario=mod(b-1,3)+1;
    if scenario>=2 && b<7, c.paper.requestedBandwidthHz=20e6; end
    if scenario==3 && b<7
        c.paper.positioningBudgetW=1e-3;
        c.numerics.powerScaleW=c.paper.positioningBudgetW/c.paper.M;
    end
    cancellation=b>=4;
    c.models.sampledReceiver.interferenceCancellation.enabled=cancellation;
    % Explicit legacy policy, so historical stage8 scenarios remain replayable.
    c.models.sampledReceiver.interferenceCancellation.maximumRelativeResidualEnergy=1;
    c.experiment.name='stage8_paired_baseline_CP_PSD';
    if cancellation, c.experiment.name='stage8_estimated_P_cancellation_Gaussian_precision_gate'; end
    if b==1, c.experiment.interferenceAblation.enabled=true; end
    if b==7
        c.models.sampledReceiver.includeCommunication=false;
        c.models.sampledReceiver.noiseMultiplier=0;
        c.models.sampledReceiver.codeMode='seeded_independent_bpsk';
        c.experiment.name='stage8_original_domain_no_C_no_AWGN_all_P_retained';
    end
    folder=main(c); loaded=load(fullfile(folder,'run.mat'),'run'); r=loaded.run;
    assert(strcmp(r.schemaVersion,c.schemaVersion));
    central=r.algorithms(strcmp({r.algorithms.id},'central'));
    for name={'central','dual','cp4a_consistent_candidate','cp4a_subgradient_candidate'}
        a=r.algorithms(strcmp({r.algorithms.id},name{1})); s=a.solution;
        assert(strcmp(s.status,'converged_verified') && all(a.evaluation.valid));
        assert(a.evaluation.communicationQosSatisfied);
        assert(abs(s.objectiveM2/central.solution.objectiveM2-1)<1e-6);
        assert(s.diagnostics.maxNormalizedViolation<=c.numerics.feasibilityTolerance);
        assert(s.diagnostics.maxRelativeHearabilityViolation<=c.models.hearability.comparisonTolerance);
        assert(s.diagnostics.kktStationarity<=c.numerics.kktTolerance);
        assert(s.diagnostics.kktComplementarity<=c.numerics.kktTolerance);
        if ~strcmp(name{1},'central')
            assert(~s.usedCentralFallback && abs(s.relativeDualityGap)<=c.numerics.dual.relativeGapTolerance);
        end
    end
    s=r.sampledReceiver; f=r.positionFeedback;
    assert(~s.truthUsedByEstimator && ~s.communicationDataUsedByEstimator);
    assert(~s.steadyStateDllImplemented && ~s.fullyMatched && ~r.prsComparison.fullyMatched);
    assert(all(isnan(s.positionErrorM(~s.validPosition))));
    assert(height(readtable(fullfile(folder,'sampled_ranges.csv')))==80);
    assert(height(readtable(fullfile(folder,'sampled_positions.csv')))==20);
    if b<=3
        baseline{scenario}=s;
    elseif b<=6
        old=baseline{scenario};
        assert(isequal(old.packet.chips,s.packet.chips));
        assert(isequal(old.packet.communicationSymbols,s.packet.communicationSymbols));
        assert(isequal(old.options.noiseSeed,s.options.noiseSeed) && isequal(old.powerW,s.powerW));
    end
    fitsValid=0; eligible=0;
    if cancellation
        fitsValid=nnz(cellfun(@(x) x.valid,s.estimatedPChannels));
        eligible=nnz(cellfun(@(x) x.cancellationEligible,s.estimatedPChannels));
        assert(~s.estimatedCancellationUncertaintyPropagated);
        assert(height(readtable(fullfile(folder,'sampled_estimated_P_channels.csv')))==80);
        for u=1:20
            fit=s.estimatedPChannels{u};
            assert(~fit.truthUsed && ~fit.communicationDataUsed);
            assert(all(diff(fit.history(:,2))<=0));
            applied=cellfun(@(x) x.pCancellationApplied,s.receivers(:,u));
            assert(all(applied)==(fit.valid && fit.cancellationEligible));
            if b<=6 && ~fit.cancellationEligible
                assert(isequaln(s.rangeM(:,u),baseline{scenario}.rangeM(:,u)));
            end
        end
    end
    assert(f.firstFixUnitLambda && ~f.usedOracleFallback && ~f.truthGeometryUsed);
    assert(all(f.history(1).lambda==1,'all'));
    for j=2:f.iterations
        previous=f.history(j-1).acquisition;
        scene=struct('bsPosM',r.scene.bsPosM,'pUserPosM',previous.estimatedPositionM);
        spec=c.models.geometry; spec.source='estimated_sampled_ranges';
        expected=msnoma.geometry.weights(scene,true(4,20),spec);
        assert(isequal(expected.lambda,f.history(j).lambda));
    end
    if strcmp(f.status,'acquisition_or_position_failure')
        last=f.history(end).acquisition;
        assert(~all(last.validPosition) || ~all(last.validRanges,'all'));
    end
    if b==7
        assert(all(s.validRanges,'all') && all(s.validPosition));
        assert(f.iterations>=2 && f.converged);
    end
    figures=13+double(cancellation)+double(b==1);
    assert(numel(dir(fullfile(folder,'figure_*.png')))==figures);
    if b==1
        assert(height(readtable(fullfile(folder,'interference_ablation.csv')))==24);
        audit=r.interferenceAblation.rows;
        only=audit(audit.transmitterMask=="desired_only",:);
        assert(all(only.valid) && max(abs(only.rangeErrorForEvaluationM))<.01);
    end
    last=f.history(end).acquisition;
    metrics{b}=struct('name',c.experiment.name,'bandwidthHz',r.num.BHz, ...
        'budgetW',c.paper.positioningBudgetW,'cancellationEnabled',cancellation, ...
        'communicationPresent',c.models.sampledReceiver.includeCommunication, ...
        'noiseMultiplier',c.models.sampledReceiver.noiseMultiplier, ...
        'codeMode',s.packet.codeMode,'validRanges',nnz(s.validRanges), ...
        'validPositions',nnz(s.validPosition),'totalPositions',20, ...
        'conditionalMeanActualErrorM',mean(s.positionErrorM,'omitnan'), ...
        'conditionalMaxActualErrorM',max(s.positionErrorM,[],'omitnan'), ...
        'channelFitsConverged',fitsValid,'cancellationEligibleUsers',eligible, ...
        'feedbackStatus',f.status,'feedbackIterations',f.iterations, ...
        'feedbackConverged',f.converged,'feedbackFinalValidPositions',nnz(last.validPosition), ...
        'feedbackFinalMeanErrorM',mean(last.positionErrorM,'omitnan'),'figures',figures);
    runs{b}=folder; fprintf('Verified stage8 main %d: %s\n',b,folder);
end
hashes=msnoma.report.sourceHashes(root);
for b=1:7
    manifest=jsondecode(fileread(fullfile(runs{b},'manifest.json')));
    assert(numel(manifest.sourceHashes)==numel(hashes));
    for i=1:numel(hashes)
        assert(strcmp(hashes(i).file,manifest.sourceHashes(i).file));
        assert(strcmp(hashes(i).sha256,manifest.sourceHashes(i).sha256));
    end
end
counts=cellfun(@height,reports);
verification=struct('date','2026-11-01','timezone','Asia/Shanghai', ...
    'testCounts',counts,'totalTests',sum(counts),'cancellationTests',counts(end), ...
    'sourceFiles',sourceCount,'codeAnalyzerDiagnostics',0,'sourceHashesMatchedAllRuns',true, ...
    'runs',{runs},'scenarios',{metrics});
save(fullfile(root,'docs','STAGE8_VERIFICATION.mat'),'verification','reports','hashes');
fid=fopen(fullfile(root,'docs','STAGE8_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true)); disp(verification);
end
