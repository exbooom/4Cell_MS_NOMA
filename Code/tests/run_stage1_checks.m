function verification = run_stage1_checks()
%RUN_STAGE1_CHECKS Model tests, two real-size runs, replay/floor and export QA.
codeRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(codeRoot);
testReport=run_tests();
files=dir(fullfile(codeRoot,'**','*.m'));
diagnostics=struct('file',{},'messages',{});
for i=1:numel(files)
    f=fullfile(files(i).folder,files(i).name);
    if contains(f,[filesep 'MS_NOMA' filesep]) || contains(f,[filesep '__MACOSX' filesep])
        continue;
    end
    messages=checkcode(f,'-id');
    if ~isempty(messages), diagnostics(end+1)=struct('file',f,'messages',messages); end %#ok<AGROW>
end
assert(isempty(diagnostics),'New MATLAB files must pass Code Analyzer.');
fprintf('Code Analyzer: zero diagnostics.\n');
bandwidths=[50e6,20e6];
runs=cell(1,2);
for b=1:2
    cfg=msnoma.config.base;
    cfg.models.dllBias.enabled=false;
    cfg.experiment.algorithmIds={'budget_equal','qos_scaled_equal','central'};
    cfg.paper.requestedBandwidthHz=bandwidths(b);
    cfg.io.showFigures=false;
    cfg.models.rangingIntegral.enabled=false;
    cfg.models.prs.enabled=false;
    runDir=main(cfg);
    saved=load(fullfile(runDir,'run.mat'),'run');
    r=saved.run;
    central=r.algorithms(strcmp({r.algorithms.id},'central'));
    assert(strcmp(central.solution.status,'converged_verified'));
    assert(central.evaluation.communicationQosSatisfied && all(central.evaluation.valid));
    assert(central.solution.diagnostics.maxNormalizedViolation<cfg.numerics.feasibilityTolerance);
    assert(abs(central.evaluation.objectiveAllBsM2/central.solution.objectiveM2-1)<1e-12);
    assert(nnz(central.evaluation.hearability.mask)==r.num.K*r.num.M);
    assert(min(central.solution.powerW(:))>100*cfg.numerics.minimumPowerW);
    manifest=jsondecode(fileread(fullfile(runDir,'manifest.json')));
    assert(numel(manifest.sourceHashes)>=30 && manifest.numerology.N==r.num.N);
    constraints=readtable(fullfile(runDir,'central_constraints.csv'));
    assert(height(constraints)==r.num.K*r.num.N+r.num.K+r.num.M*r.num.K*(r.num.K-1));
    assert(numel(dir(fullfile(runDir,'figure_*.png')))==4);
    runs{b}=runDir;
    fprintf('Verified full-size %.2f MHz run: %s\n',r.num.BHz/1e6,runDir);
end
% Independent rerun at a 100x smaller lower bound, using exactly the same scene.
cfg=r.cfg;
cfg.numerics.minimumPowerW=cfg.numerics.minimumPowerW/100;
replay=msnoma.experiment.runSingle(cfg,{'central'},r.scene);
reference=r.algorithms(strcmp({r.algorithms.id},'central')).solution;
candidate=replay.algorithms(1).solution;
assert(strcmp(candidate.status,'converged_verified'));
floorRelativeObjectiveChange=abs(candidate.objectiveM2/reference.objectiveM2-1);
assert(floorRelativeObjectiveChange<1e-5);
assert(isequal(r.scene,replay.scene));
verification=struct('date','2026-10-31','mathematicalTests',height(testReport), ...
    'codeAnalyzerDiagnostics',0,'run50',runs{1},'run20',runs{2}, ...
    'floorRelativeObjectiveChange',floorRelativeObjectiveChange);
save(fullfile(codeRoot,'docs','STAGE1_VERIFICATION.mat'),'verification','testReport');
fid=fopen(fullfile(codeRoot,'docs','STAGE1_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true));
disp(verification);
end
