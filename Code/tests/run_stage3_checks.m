function verification = run_stage3_checks()
%RUN_STAGE3_CHECKS Certify both dual solvers with bandwidth and budget stress.
codeRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(codeRoot);
basic=run_tests(); precision=run_precision_tests(); dual=run_dual_tests();
files=dir(fullfile(codeRoot,'**','*.m'));
sourceCount=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep])
        continue;
    end
    assert(isempty(checkcode(path,'-id')),'Code Analyzer diagnostics in %s',path);
    sourceCount=sourceCount+1;
end
runs=cell(1,3); allMetrics=cell(1,3);
for b=1:3
    c=msnoma.config.base; c.io.showFigures=false;
    c.models.dllBias.enabled=false;
    c.models.prs.comparison.enabled=false;
    c.experiment.algorithmIds={'budget_equal','qos_scaled_equal','central', ...
        'dual','cp4a_consistent_candidate'};
    if b>1, c.paper.requestedBandwidthHz=20e6; end
    if b==3
        c.paper.positioningBudgetW=1e-3;
        c.numerics.powerScaleW=c.paper.positioningBudgetW/c.paper.M;
        c.experiment.name='budget_binding_stress';
    end
    runDir=main(c);
    loaded=load(fullfile(runDir,'run.mat'),'run'); r=loaded.run;
    central=r.algorithms(strcmp({r.algorithms.id},'central')).solution;
    assert(strcmp(central.status,'converged_verified'));
    m=repmat(struct('algorithm','','relativeToCentral',0,'relativeDualityGap',0, ...
        'maxNormalizedViolation',0,'kktStationarity',0,'kktComplementarity',0, ...
        'objectiveM2',0,'elapsedS',0),1,2);
    ids={'dual','cp4a_consistent_candidate'};
    for j=1:2
        a=r.algorithms(strcmp({r.algorithms.id},ids{j})); s=a.solution;
        assert(strcmp(s.status,'converged_verified'));
        assert(~s.usedCentralFallback && all(a.evaluation.valid));
        assert(a.evaluation.communicationQosSatisfied);
        assert(abs(s.objectiveM2/central.objectiveM2-1)<1e-6);
        assert(abs(s.relativeDualityGap)<c.numerics.dual.relativeGapTolerance);
        assert(s.dualBoundM2<=central.objectiveM2*(1+1e-7));
        assert(abs(a.evaluation.objectiveAllBsM2/s.objectiveM2-1)<1e-12);
        d=msnoma.optim.diagnostics(r.problem,s);
        assert(d.maxNormalizedViolation<=c.numerics.feasibilityTolerance);
        assert(d.kktStationarity<=c.numerics.kktTolerance);
        assert(d.kktComplementarity<=c.numerics.kktTolerance);
        m(j)=struct('algorithm',ids{j},'relativeToCentral',s.objectiveM2/central.objectiveM2-1, ...
            'relativeDualityGap',s.relativeDualityGap, ...
            'maxNormalizedViolation',d.maxNormalizedViolation, ...
            'kktStationarity',d.kktStationarity,'kktComplementarity',d.kktComplementarity, ...
            'objectiveM2',s.objectiveM2,'elapsedS',s.elapsedS);
        constraints=readtable(fullfile(runDir,[ids{j} '_constraints.csv']));
        assert(height(constraints)==size(r.problem.Aineq,1));
        assert(height(readtable(fullfile(runDir,[ids{j} '_dll_validation.csv'])))==12);
    end
    candidate=r.algorithms(strcmp({r.algorithms.id},ids{2})).solution;
    assert(~candidate.literalPaperAlgorithm && ~candidate.groupAudit.solverUsesPrintedGroups);
    if b==3
        assert(candidate.innerEvaluations>3 && size(candidate.history,1)>3);
        assert(any(candidate.multipliers.ineqlin(r.problem.rowMap.budget)>1e-5));
        assert(max(sum(candidate.powerW,2))<=c.paper.positioningBudgetW*(1+1e-7));
    end
    assert(all(r.integralKernel.converged));
    assert(all(isfinite(r.prsBenchmark.position.positionStdM)));
    assert(all(r.prsBenchmark.budget.feasible));
    assert(numel(dir(fullfile(runDir,'figure_*.png')))==7);
    comparison=readtable(fullfile(runDir,'optimization_comparison.csv'));
    assert(height(comparison)==5);
    assert(height(readtable(fullfile(runDir,'paper_carrier_group_audit.csv')))==r.num.M);
    runs{b}=runDir;
    allMetrics{b}=struct('effectiveBandwidthHz',r.num.BHz, ...
        'budgetPerBSW',c.paper.positioningBudgetW,'metrics',m, ...
        'candidateOuterIterations',candidate.solverOutput.iterations, ...
        'candidateInnerEvaluations',candidate.innerEvaluations, ...
        'printedCarrierOverrun',candidate.groupAudit.outOfBandCount);
    fprintf('Verified optimization stage %d: %s\n',b,runDir);
end
verification=struct('date','2026-11-01','timezone','Asia/Shanghai', ...
    'basicTests',height(basic),'precisionTests',height(precision),'dualTests',height(dual), ...
    'totalTests',height(basic)+height(precision)+height(dual), ...
    'codeAnalyzerDiagnostics',0,'sourceFiles',sourceCount, ...
    'run50',runs{1},'run20',runs{2},'budgetStressRun',runs{3}, ...
    'scenarios',{allMetrics});
save(fullfile(codeRoot,'docs','STAGE3_VERIFICATION.mat'),'verification','basic','precision','dual');
fid=fopen(fullfile(codeRoot,'docs','STAGE3_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true));
disp(verification);
end
