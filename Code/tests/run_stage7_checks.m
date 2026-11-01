function verification = run_stage7_checks()
%RUN_STAGE7_CHECKS Existing invariants, sampled paths and full-size first fixes.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
basic=run_tests(); precision=run_precision_tests(); dual=run_dual_tests();
projected=run_subgradient_tests(); comparison=run_prs_comparison_tests();
bias=run_bias_tests(); sampled=run_sampled_tests();
files=dir(fullfile(root,'**','*.m')); sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep]), continue; end
    diagnostics=checkcode(path,'-id');
    if ~isempty(diagnostics), disp(path); disp(diagnostics); end
    assert(isempty(diagnostics),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,5); metrics=cell(1,5);
for b=1:5
    c=msnoma.config.base; c.io.showFigures=false;
    c.models.sampledReceiver.enabled=true;
    c.experiment.positionFeedback.enabled=true;
    if b==2 || b==3, c.paper.requestedBandwidthHz=20e6; end
    if b==3
        c.paper.positioningBudgetW=1e-3;
        c.numerics.powerScaleW=c.paper.positioningBudgetW/c.paper.M;
        c.experiment.stepSensitivity.enabled=true;
    end
    if b>=4
        c.models.sampledReceiver.includeCommunication=false;
        c.models.sampledReceiver.noiseMultiplier=0;
        c.models.sampledReceiver.codeMode='seeded_independent_bpsk';
        c.experiment.name='stage7_no_C_no_AWGN_original_domain_P_interference_retained';
        if b==5
            % Positive full-size chain control, without the near-BS users.
            % The original-domain no-C/no-AWGN case remains a separate run.
            c.assumed.domainM=[50,150;50,150];
            c.models.sampledReceiver.position.domainM=c.assumed.domainM;
            c.experiment.name='stage7_interior_domain_no_C_no_AWGN_full_feedback';
        end
    else
        c.experiment.name='stage7_finite_sampled_C_P_AWGN_acquisition_feedback';
    end
    folder=main(c); loaded=load(fullfile(folder,'run.mat'),'run'); r=loaded.run;
    assert(strcmp(r.schemaVersion,c.schemaVersion));
    central=r.algorithms(strcmp({r.algorithms.id},'central'));
    for name={'central','dual','cp4a_consistent_candidate','cp4a_subgradient_candidate'}
        a=r.algorithms(strcmp({r.algorithms.id},name{1})); s=a.solution;
        assert(strcmp(s.status,'converged_verified') && all(a.evaluation.valid));
        assert(a.evaluation.communicationQosSatisfied);
        assert(abs(s.objectiveM2/central.solution.objectiveM2-1)<1e-6);
        assert(abs(a.evaluation.objectiveAllBsM2/s.objectiveM2-1)<1e-12);
        assert(s.diagnostics.maxNormalizedViolation<=c.numerics.feasibilityTolerance);
        assert(s.diagnostics.maxRelativeHearabilityViolation<=c.models.hearability.comparisonTolerance);
        assert(s.diagnostics.kktStationarity<=c.numerics.kktTolerance);
        assert(s.diagnostics.kktComplementarity<=c.numerics.kktTolerance);
        if ~strcmp(name{1},'central')
            assert(~s.usedCentralFallback && abs(s.relativeDualityGap)<=c.numerics.dual.relativeGapTolerance);
        end
    end
    s=r.sampledReceiver; f=r.positionFeedback;
    assert(s.sampledReceiverImplemented && ~s.truthUsedByEstimator && ~s.steadyStateDllImplemented);
    assert(isequal(size(s.rangeM),[r.num.K,r.num.M]) && ~s.fullyMatched);
    assert(~r.prsComparison.fullyMatched && ~r.prsComparison.observationMatched);
    assert(all(isnan(s.positionErrorM(~s.validPosition))));
    assert(height(readtable(fullfile(folder,'sampled_ranges.csv')))==r.num.K*r.num.M);
    assert(height(readtable(fullfile(folder,'sampled_positions.csv')))==r.num.M);
    assert(height(readtable(fullfile(folder,'sampled_packet_audit.csv')))==r.num.K);
    assert(f.firstFixUnitLambda && ~f.truthGeometryUsed && ~f.usedOracleFallback && f.knownChannelCsiAssumed);
    assert(all(f.history(1).lambda==1,'all'));
    for j=2:f.iterations
        previous=f.history(j-1).acquisition;
        estimate=struct('bsPosM',r.scene.bsPosM,'pUserPosM',previous.estimatedPositionM);
        spec=c.models.geometry; spec.source='estimated_sampled_ranges';
        expected=msnoma.geometry.weights(estimate,true(r.num.K,r.num.M),spec);
        assert(isequal(expected.lambda,f.history(j).lambda));
        assert(strcmp(f.history(j).geometrySource,'estimated_sampled_ranges'));
    end
    if strcmp(f.status,'acquisition_or_position_failure')
        a=f.history(end).acquisition;
        assert(~all(a.validPosition) || ~all(a.validRanges,'all'));
        assert(~f.converged);
    end
    if b==5
        assert(all(s.validPosition) && all(s.validRanges,'all'));
        assert(f.iterations>=2 && f.converged, ...
            'Interior complete sampled feedback must reach estimated geometry and convergence.');
    end
    expectedFigures=13+double(b==3);
    assert(numel(dir(fullfile(folder,'figure_*.png')))==expectedFigures);
    if b==3
        step=r.stepSensitivity; assert(step.verifiedCount==numel(step.trials));
        assert(step.objectiveRelativeSpread<1e-6 && ~step.authorStepsKnown);
        steps=readtable(fullfile(folder,'step_sensitivity.csv'));
        assert(height(steps)==3 && all(steps.lastAcceptedInnerStep>0));
        assert(all(steps.lastAcceptedOuterStep>0));
    end
    actualErrors=s.positionErrorM(s.validPosition);
    last=f.history(end).acquisition;
    metrics{b}=struct('experimentName',c.experiment.name,'domainM',c.assumed.domainM, ...
        'effectiveBandwidthHz',r.num.BHz,'budgetPerBSW',c.paper.positioningBudgetW, ...
        'codeMode',s.packet.codeMode,'includeCommunication',s.options.includeCommunication, ...
        'noiseMultiplier',s.options.noiseMultiplier,'whitening',s.options.whitening, ...
        'activeDurationS',s.packet.activeDurationS,'sampleRateHz',s.packet.sampleRateHz, ...
        'validRanges',nnz(s.validRanges),'totalRanges',numel(s.validRanges), ...
        'validPositions',nnz(s.validPosition),'totalPositions',numel(s.validPosition), ...
        'conditionalMeanActualPositionErrorM',mean(actualErrors), ...
        'conditionalMaxActualPositionErrorM',max(actualErrors), ...
        'feedbackStatus',f.status,'feedbackIterations',f.iterations, ...
        'feedbackConverged',f.converged,'finalFeedbackValidPositions',nnz(last.validPosition), ...
        'finalFeedbackConditionalMeanErrorM',mean(last.positionErrorM,'omitnan'), ...
        'stepSensitivityStatus',r.stepSensitivity.status,'figures',expectedFigures);
    runs{b}=folder; fprintf('Verified stage7 main %d: %s\n',b,folder);
end
hashes=msnoma.report.sourceHashes(root);
for b=1:numel(runs)
    manifest=jsondecode(fileread(fullfile(runs{b},'manifest.json')));
    assert(numel(manifest.sourceHashes)==numel(hashes));
    for i=1:numel(hashes)
        assert(strcmp(hashes(i).file,manifest.sourceHashes(i).file));
        assert(strcmp(hashes(i).sha256,manifest.sourceHashes(i).sha256));
    end
end
counts=[height(basic),height(precision),height(dual),height(projected), ...
    height(comparison),height(bias),height(sampled)];
verification=struct('date','2026-11-01','timezone','Asia/Shanghai', ...
    'basicTests',counts(1),'precisionTests',counts(2),'dualTests',counts(3), ...
    'projectedTests',counts(4),'prsComparisonTests',counts(5),'biasAndFilterTests',counts(6), ...
    'sampledTests',counts(7),'totalTests',sum(counts),'codeAnalyzerDiagnostics',0, ...
    'sourceFiles',sourceCount,'sourceHashesMatchedAllRuns',true,'runs',{runs},'scenarios',{metrics});
save(fullfile(root,'docs','STAGE7_VERIFICATION.mat'),'verification','basic','precision','dual', ...
    'projected','comparison','bias','sampled','hashes');
fid=fopen(fullfile(root,'docs','STAGE7_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true)); disp(verification);
end
