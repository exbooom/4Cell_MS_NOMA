function verification = run_stage4_checks()
%RUN_STAGE4_CHECKS First-order hierarchy, full metrics and failure semantics.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
basic=run_tests(); precision=run_precision_tests();
dual=run_dual_tests(); projected=run_subgradient_tests();
files=dir(fullfile(root,'**','*.m')); sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep])
        continue;
    end
    assert(isempty(checkcode(path,'-id')),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,3); metrics=cell(1,3);
for b=1:3
    cfg=msnoma.config.base; cfg.io.showFigures=false;
    cfg.models.dllBias.enabled=false;
    cfg.models.prs.comparison.enabled=false;
    if b>1, cfg.paper.requestedBandwidthHz=20e6; end
    if b==3
        cfg.paper.positioningBudgetW=1e-3;
        cfg.numerics.powerScaleW=cfg.paper.positioningBudgetW/cfg.paper.M;
        cfg.experiment.name='projected_budget_binding_stress';
    end
    runDir=main(cfg); loaded=load(fullfile(runDir,'run.mat'),'run'); r=loaded.run;
    central=r.algorithms(strcmp({r.algorithms.id},'central')).solution;
    assert(strcmp(central.status,'converged_verified'));
    for name={'dual','cp4a_consistent_candidate','cp4a_subgradient_candidate'}
        a=r.algorithms(strcmp({r.algorithms.id},name{1})); s=a.solution;
        assert(strcmp(s.status,'converged_verified'),'%s failed: %s',name{1},s.status);
        assert(~s.usedCentralFallback && all(a.evaluation.valid));
        assert(a.evaluation.communicationQosSatisfied);
        assert(abs(s.objectiveM2/central.objectiveM2-1)<1e-6);
        assert(abs(s.relativeDualityGap)<=cfg.numerics.dual.relativeGapTolerance);
        assert(s.dualBoundM2<=central.objectiveM2*(1+1e-7));
        assert(abs(a.evaluation.objectiveAllBsM2/s.objectiveM2-1)<1e-12);
        assert(s.diagnostics.maxNormalizedViolation<=cfg.numerics.feasibilityTolerance);
        assert(s.diagnostics.maxRelativeHearabilityViolation<=cfg.models.hearability.comparisonTolerance);
        assert(s.diagnostics.kktStationarity<=cfg.numerics.kktTolerance);
        assert(s.diagnostics.kktComplementarity<=cfg.numerics.kktTolerance);
        assert(height(readtable(fullfile(runDir,[name{1} '_constraints.csv'])))==size(r.problem.Aineq,1));
        assert(height(readtable(fullfile(runDir,[name{1} '_dll_validation.csv'])))==12);
    end
    s=r.algorithms(strcmp({r.algorithms.id},'cp4a_subgradient_candidate')).solution;
    assert(~s.usedPriceOptimizer && ~s.literalPaperAlgorithm);
    assert(s.innerCalls==numel(s.innerTrace));
    assert(all(cellfun(@(item) item.state.converged,s.innerTrace)));
    assert(s.innerIterations==sum(cellfun(@(item) item.state.iterations,s.innerTrace)));
    assert(height(readtable(fullfile(runDir,'subgradient_inner_calls.csv')))==s.innerCalls);
    assert(height(readtable(fullfile(runDir,'subgradient_outer_prices.csv')))==size(s.history,1));
    assert(height(readtable(fullfile(runDir,'optimization_comparison.csv')))==6);
    assert(numel(dir(fullfile(runDir,'figure_*.png')))==8);
    assert(all(r.integralKernel.converged) && all(r.prsBenchmark.budget.feasible));
    assert(all(isfinite(r.prsBenchmark.position.positionStdM)));
    if b==3
        assert(s.solverOutput.iterations>3 && s.innerCalls>3);
        assert(any(s.multipliers.ineqlin(r.problem.rowMap.budget)>1e-5));
    end
    runs{b}=runDir;
    metrics{b}=struct('effectiveBandwidthHz',r.num.BHz, ...
        'budgetPerBSW',cfg.paper.positioningBudgetW, ...
        'relativeToCentral',s.objectiveM2/central.objectiveM2-1, ...
        'relativeDualityGap',s.relativeDualityGap, ...
        'maxNormalizedViolation',s.diagnostics.maxNormalizedViolation, ...
        'maxRelativeHearabilityViolation',s.diagnostics.maxRelativeHearabilityViolation, ...
        'kktStationarity',s.diagnostics.kktStationarity, ...
        'kktComplementarity',s.diagnostics.kktComplementarity, ...
        'outerIterations',s.solverOutput.iterations,'innerCalls',s.innerCalls, ...
        'innerIterations',s.innerIterations,'elapsedS',s.elapsedS);
    fprintf('Verified projected-price stage %d: %s\n',b,runDir);
end
verification=struct('date','2026-11-01','timezone','Asia/Shanghai', ...
    'basicTests',height(basic),'precisionTests',height(precision), ...
    'dualTests',height(dual),'projectedTests',height(projected), ...
    'totalTests',height(basic)+height(precision)+height(dual)+height(projected), ...
    'codeAnalyzerDiagnostics',0,'sourceFiles',sourceCount, ...
    'run50',runs{1},'run20',runs{2},'budgetStressRun',runs{3},'scenarios',{metrics});
save(fullfile(root,'docs','STAGE4_VERIFICATION.mat'), ...
    'verification','basic','precision','dual','projected');
fid=fopen(fullfile(root,'docs','STAGE4_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true)); disp(verification);
end
