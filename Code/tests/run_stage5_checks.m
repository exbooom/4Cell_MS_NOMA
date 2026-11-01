function verification = run_stage5_checks()
%RUN_STAGE5_CHECKS Preserve optimizer certificates and audit PRS comparability.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
basic=run_tests(); precision=run_precision_tests();
dual=run_dual_tests(); projected=run_subgradient_tests(); comparison=run_prs_comparison_tests();
files=dir(fullfile(root,'**','*.m')); sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep])
        continue;
    end
    messages=checkcode(path,'-id');
    if ~isempty(messages), disp(path); disp(messages); end
    assert(isempty(messages),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,3); metrics=cell(1,3);
for b=1:3
    cfg=msnoma.config.base; cfg.io.showFigures=false;
    cfg.models.dllBias.enabled=false;
    cfg.models.prs.comparison.receiverAuditEnabled=false;
    if b>1, cfg.paper.requestedBandwidthHz=20e6; end
    if b==3
        cfg.paper.positioningBudgetW=1e-3;
        cfg.numerics.powerScaleW=cfg.paper.positioningBudgetW/cfg.paper.M;
        cfg.experiment.name='stage5_binding_budget_and_service_audit';
    end
    runDir=main(cfg); data=load(fullfile(runDir,'run.mat'),'run'); r=data.run;
    central=r.algorithms(strcmp({r.algorithms.id},'central')).solution;
    assert(strcmp(central.status,'converged_verified'));
    for id={'central','dual','cp4a_consistent_candidate','cp4a_subgradient_candidate'}
        a=r.algorithms(strcmp({r.algorithms.id},id{1})); s=a.solution;
        assert(strcmp(s.status,'converged_verified'),'%s: %s',id{1},s.status);
        assert(all(a.evaluation.valid) && a.evaluation.communicationQosSatisfied);
        assert(abs(s.objectiveM2/central.objectiveM2-1)<1e-6);
        assert(abs(a.evaluation.objectiveAllBsM2/s.objectiveM2-1)<1e-12);
        assert(s.diagnostics.maxNormalizedViolation<=cfg.numerics.feasibilityTolerance);
        assert(s.diagnostics.maxRelativeHearabilityViolation<=cfg.models.hearability.comparisonTolerance);
        assert(s.diagnostics.kktStationarity<=cfg.numerics.kktTolerance);
        assert(s.diagnostics.kktComplementarity<=cfg.numerics.kktTolerance);
        if ~strcmp(id{1},'central')
            assert(~s.usedCentralFallback);
            assert(abs(s.relativeDualityGap)<=cfg.numerics.dual.relativeGapTolerance);
        end
        assert(height(readtable(fullfile(runDir,[id{1} '_dll_validation.csv'])))==12);
    end
    assert(all(r.integralKernel.converged) && all(r.prsBenchmark.budget.feasible));
    c=r.prsComparison;
    assert(strcmp(c.status,'computed_partial_match_service_and_observation_audit'));
    assert(c.positioningBroadcastEnergyMatched && c.geometryConventionMatched);
    assert(~c.observationMatched && ~c.fullyMatched && ~c.detectionCoverageEvaluated);
    assert(all(c.msTdoaCommon.valid) && all(c.msTdoaHearable.valid));
    assert(numel(c.policies)==2);
    for j=1:2
        p=c.policies(j); s=p.service;
        assert(all(p.tdoa.valid) && all(isfinite(p.tdoa.positionStdM)));
        assert(all(isfinite(s.meanActiveBer(:))) && all(isfinite(s.worstActiveBer(:))));
        assert(all(s.operatingPowerFeasible) && ~any(s.attemptedRateContractMet));
        assert(~any(s.serviceContractMet) && ~any(s.communicationEnergyMatched));
        assert(~any(s.energyAndServiceMatched));
        assert(all(s.worstActiveBer(:)>=s.meanActiveBer(:)-1e-15));
        if j==1
            fraction=(1-r.prsBenchmark.budget.jointDutyCycle)*ones(r.num.K,1);
            assert(isequal(p.crb.rangeVarianceM2,r.prsBenchmark.crb.rangeVarianceM2));
        else
            fraction=1-r.prsBenchmark.budget.dutyCycle;
            assert(all(p.crb.rangeVarianceM2(:)>=c.policies(1).crb.rangeVarianceM2(:)));
            assert(all(p.tdoa.positionStdM>=c.policies(1).tdoa.positionStdM*(1-1e-12)));
        end
        assert(max(abs(s.communicationOpportunityFraction-fraction))<1e-14);
        assert(max(abs(s.requiredPowerForEqualCommunicationEnergyW.*fraction-cfg.paper.communicationTotalPowerW))<1e-12);
        spec=cfg.models.prs.comparison.tdoa; spec.referencePolicy='fixed'; spec.referenceBs=1;
        fixed=msnoma.geometry.tdoaAccuracy(r.scene,p.crb.rangeVarianceM2,spec);
        assert(max(abs(fixed.positionStdM./p.tdoa.positionStdM-1))<1e-11);
        assert(height(readtable(fullfile(runDir,['prs_' p.id '_service.csv'])))==r.num.K);
        assert(height(readtable(fullfile(runDir,['prs_' p.id '_communication.csv'])))==r.num.K*r.num.N);
        assert(height(readtable(fullfile(runDir,['prs_' p.id '_tdoa.csv'])))==r.num.M);
        assert(height(readtable(fullfile(runDir,['prs_' p.id '_ranging.csv'])))==r.num.K*r.num.M);
    end
    assert(numel(dir(fullfile(runDir,'figure_*.png')))==9);
    tab=readtable(fullfile(runDir,'prs_comparison_summary.csv'));
    assert(height(tab)==2 && ~any(tab.fullyMatched) && ~any(tab.communicationServiceMet));
    runs{b}=runDir;
    policies=cell(1,2);
    for j=1:2
        p=c.policies(j); s=p.service;
        policies{j}=struct('id',p.id,'meanTdoaBoundM',mean(p.tdoa.positionStdM), ...
            'minimumCommunicationFraction',min(s.communicationOpportunityFraction), ...
            'maximumCommunicationFraction',max(s.communicationOpportunityFraction), ...
            'maxRequiredEqualEnergyPowerW',max(s.requiredPowerForEqualCommunicationEnergyW), ...
            'energyCompensationPeakFeasible',all(s.energyCompensationPeakFeasible), ...
            'worstActiveBer',max(s.worstActiveBer(:)), ...
            'activeBerQosMet',all(s.activeSlotBerQosMet), ...
            'serviceContractMet',all(s.serviceContractMet));
    end
    metrics{b}=struct('effectiveBandwidthHz',r.num.BHz,'budgetPerBSW',cfg.paper.positioningBudgetW, ...
        'jointPrsDutyCycle',r.prsBenchmark.budget.jointDutyCycle, ...
        'meanMsTdoaVarianceProxyStdM',mean(c.msTdoaCommon.positionStdM), ...
        'fullyMatched',false,'policies',{policies});
    fprintf('Verified stage5 main %d: %s\n',b,runDir);
end
hashes=msnoma.report.sourceHashes(root);
for b=1:3
    manifest=jsondecode(fileread(fullfile(runs{b},'manifest.json')));
    assert(numel(manifest.sourceHashes)==numel(hashes));
    for i=1:numel(hashes)
        assert(strcmp(manifest.sourceHashes(i).file,hashes(i).file));
        assert(strcmp(manifest.sourceHashes(i).sha256,hashes(i).sha256));
    end
end
verification=struct('date','2026-11-01','timezone','Asia/Shanghai', ...
    'basicTests',height(basic),'precisionTests',height(precision), ...
    'dualTests',height(dual),'projectedTests',height(projected),'prsComparisonTests',height(comparison), ...
    'totalTests',height(basic)+height(precision)+height(dual)+height(projected)+height(comparison), ...
    'codeAnalyzerDiagnostics',0,'sourceFiles',sourceCount,'sourceHashesMatchedAllRuns',true, ...
    'run50',runs{1},'run20',runs{2},'budgetStressRun',runs{3},'scenarios',{metrics});
save(fullfile(root,'docs','STAGE5_VERIFICATION.mat'), ...
    'verification','basic','precision','dual','projected','comparison','hashes');
fid=fopen(fullfile(root,'docs','STAGE5_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true)); disp(verification);
end
