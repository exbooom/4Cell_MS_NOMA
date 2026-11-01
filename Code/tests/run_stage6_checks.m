function verification = run_stage6_checks()
%RUN_STAGE6_CHECKS Controlled DLL differences, receiver PSD integration, main.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
basic=run_tests(); precision=run_precision_tests(); dual=run_dual_tests();
projected=run_subgradient_tests(); comparison=run_prs_comparison_tests(); bias=run_bias_tests();
files=dir(fullfile(root,'**','*.m')); sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep]), continue; end
    diagnostics=checkcode(path,'-id');
    if ~isempty(diagnostics), disp(path); disp(diagnostics); end
    assert(isempty(diagnostics),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,3); metrics=cell(1,3);
for b=1:3
    c=msnoma.config.base; c.io.showFigures=false;
    if b>1, c.paper.requestedBandwidthHz=20e6; end
    if b==3
        c.paper.positioningBudgetW=1e-3;
        c.numerics.powerScaleW=c.paper.positioningBudgetW/c.paper.M;
        c.experiment.name='stage6_binding_budget_fixed_power_bias_study';
    end
    folder=main(c); loaded=load(fullfile(folder,'run.mat'),'run'); r=loaded.run;
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
    study=r.dllBias; v=study.result; bundle=study.bundle; L=numel(v.userIds);
    assert(~study.optimizerChanged && ~study.sampledReceiverImplemented);
    assert(bundle.allConverged && numel(bundle.variants)==8 && nnz([bundle.variants.reused])==1);
    assert(all(isfinite(v.delayVarianceS2) & v.delayVarianceS2>0,'all'));
    assert(v.maximumRelativeClosureResidual<1e-12);
    assert(isequal(v.delayVarianceS2(:,:,1),central.evaluation.integral.delayVarianceS2));
    assert(max(abs((v.noiseDelayVarianceS2+v.communicationDelayVarianceS2+v.positioningDelayVarianceS2) ...
        ./v.delayVarianceS2-1),[],'all')<1e-12);
    assert(numel(v.pairedContrasts)==12 && size(v.chainVarianceS2,3)==7);
    assert(max(abs(sum(v.chainRelativeToBaseline,3)-v.totalRelativeBias),[],'all')<1e-12);
    for i=1:8
        assert(height(readtable(fullfile(folder,'dll_bias_kernels.csv')))==8*L);
        assert(all(bundle.variants(i).kernel.converged));
    end
    assert(height(readtable(fullfile(folder,'dll_bias_variants.csv')))==8*r.num.K*L);
    assert(height(readtable(fullfile(folder,'dll_bias_paired_contrasts.csv')))==12*r.num.K*L);
    assert(height(readtable(fullfile(folder,'dll_bias_chain.csv')))==6*r.num.K*L);
    assert(height(readtable(fullfile(folder,'dll_bias_positions.csv')))==8*L);
    assert(height(readtable(fullfile(folder,'dll_bias_effect_summary.csv')))==9);
    prs=r.prsComparison; a=prs.receiverAudit; p=a.pointService; f=a.filteredService;
    assert(prs.positioningBroadcastEnergyMatched && prs.geometryConventionMatched);
    assert(~prs.observationMatched && ~prs.fullyMatched && ~prs.detectionCoverageEvaluated);
    assert(~a.sampledReceiverImplemented && ~a.finiteObservationMatched && ~a.prsRangingRecomputed);
    assert(all(a.filteredUnitPsdWPerHzPerW>0,'all'));
    assert(all(isfinite(f.worstActiveBer),'all') && all(isfinite(f.meanActiveBer),'all'));
    assert(isequal(p.communicationOpportunityFraction,f.communicationOpportunityFraction));
    assert(isequal(p.activeCommunicationPowerW,f.activeCommunicationPowerW));
    assert(isequal(p.worstActiveBer,prs.policies(2).service.worstActiveBer));
    assert(~any(p.serviceContractMet) && ~any(f.serviceContractMet));
    if b==1
        assert(max(p.worstActiveBer(:))==0 && max(f.worstActiveBer(:))>0);
    end
    assert(height(readtable(fullfile(folder,'prs_receiver_filter_links.csv')))==r.num.K^2*r.num.N);
    assert(height(readtable(fullfile(folder,'prs_receiver_filter_service.csv')))==r.num.K);
    assert(numel(dir(fullfile(folder,'figure_*.png')))==11);
    effects=reshape(v.mainEffectRelativeToBaseline,[],3);
    chain=reshape(v.chainRelativeToBaseline,[],6);
    metrics{b}=struct('effectiveBandwidthHz',r.num.BHz,'budgetPerBSW',c.paper.positioningBudgetW, ...
        'selectedUsers',v.userIds,'maximumAbsSingleFactorEffects',max(abs(effects),[],1), ...
        'meanSingleFactorEffects',mean(effects,1),'meanOrderedChainContributions',mean(chain,1), ...
        'maximumAbsTotalRelativeVarianceBias',max(abs(v.totalRelativeBias(:))), ...
        'maximumRelativeClosureResidual',v.maximumRelativeClosureResidual, ...
        'maximumQuadratureChange',max(arrayfun(@(q) max(q.kernel.relativeChange),bundle.variants)), ...
        'pointWorstActiveBer',max(p.worstActiveBer(:)),'filteredWorstActiveBer',max(f.worstActiveBer(:)), ...
        'pointActiveBerQosMet',all(p.activeSlotBerQosMet),'filteredActiveBerQosMet',all(f.activeSlotBerQosMet), ...
        'fullyMatched',false,'sampledReceiverImplemented',false);
    runs{b}=folder; fprintf('Verified stage6 main %d: %s\n',b,folder);
end
hashes=msnoma.report.sourceHashes(root);
for b=1:3
    manifest=jsondecode(fileread(fullfile(runs{b},'manifest.json')));
    assert(numel(manifest.sourceHashes)==numel(hashes));
    for i=1:numel(hashes)
        assert(strcmp(hashes(i).file,manifest.sourceHashes(i).file));
        assert(strcmp(hashes(i).sha256,manifest.sourceHashes(i).sha256));
    end
end
verification=struct('date','2026-11-01','timezone','Asia/Shanghai', ...
    'basicTests',height(basic),'precisionTests',height(precision),'dualTests',height(dual), ...
    'projectedTests',height(projected),'prsComparisonTests',height(comparison),'biasAndFilterTests',height(bias), ...
    'totalTests',height(basic)+height(precision)+height(dual)+height(projected)+height(comparison)+height(bias), ...
    'codeAnalyzerDiagnostics',0,'sourceFiles',sourceCount,'sourceHashesMatchedAllRuns',true, ...
    'run50',runs{1},'run20',runs{2},'budgetStressRun',runs{3},'scenarios',{metrics});
save(fullfile(root,'docs','STAGE6_VERIFICATION.mat'),'verification','basic','precision','dual', ...
    'projected','comparison','bias','hashes');
fid=fopen(fullfile(root,'docs','STAGE6_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true)); disp(verification);
end
