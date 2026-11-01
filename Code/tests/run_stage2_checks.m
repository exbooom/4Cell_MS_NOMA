function verification = run_stage2_checks()
%RUN_STAGE2_CHECKS End-to-end DLL and PRS validation with two bandwidths.
codeRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(codeRoot);
basic=run_tests();
precision=run_precision_tests();
files=dir(fullfile(codeRoot,'**','*.m'));
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    if contains(path,[filesep 'MS_NOMA' filesep]) || contains(path,[filesep '__MACOSX' filesep])
        continue;
    end
    messages=checkcode(path,'-id');
    assert(isempty(messages),'Code Analyzer diagnostics in %s',path);
end
runs=cell(1,2);
metrics=repmat(struct('effectiveBandwidthHz',0,'maxQuadratureChange',0, ...
    'maxTwoTermRelativeVarianceError',0,'maxThreeTermRelativeVarianceError',0, ...
    'maxFrontendRelativeVarianceDifference',0,'prsJointDutyCycle',0, ...
    'prsMeanPositionBoundM',0),1,2);
bandwidths=[50e6,20e6];
for b=1:2
    c=msnoma.config.base;
    c.models.dllBias.enabled=false;
    c.models.prs.comparison.enabled=false;
    c.experiment.algorithmIds={'budget_equal','qos_scaled_equal','central'};
    c.paper.requestedBandwidthHz=bandwidths(b);
    c.io.showFigures=false;
    runDir=main(c);
    saved=load(fullfile(runDir,'run.mat'),'run');
    r=saved.run;
    central=r.algorithms(strcmp({r.algorithms.id},'central'));
    assert(strcmp(central.solution.status,'converged_verified'));
    assert(central.evaluation.communicationQosSatisfied && all(central.evaluation.valid));
    assert(all(r.integralKernel.converged));
    exact=central.evaluation.integral;
    assert(all(isfinite(exact.delayVarianceS2(:))) && all(exact.delayVarianceS2(:)>0));
    reconstructed=exact.noiseDelayVarianceS2+exact.communicationDelayVarianceS2 ...
        +exact.positioningDelayVarianceS2;
    assert(max(abs(reconstructed(:)./exact.delayVarianceS2(:)-1))<1e-12);
    p=r.prsBenchmark;
    assert(strcmp(p.status,'computed_assumed_AWGN_broadcast_bound'));
    assert(~p.detectionCoverageEvaluated && ~p.budget.communicationQosMatched);
    assert(all(p.budget.feasible) && p.budget.jointDutyCycle<=1);
    assert(max(abs(p.budget.averagePowerW-sum(central.solution.powerW,2)))<1e-14);
    assert(max(abs(p.budget.instantaneousPowerW.*p.budget.observationTimeS ...
        -p.budget.energyPerMeasurementJ))<1e-14);
    assert(numel(dir(fullfile(runDir,'figure_*.png')))==6);
    assert(height(readtable(fullfile(runDir,'central_dll_validation.csv')))==12);
    assert(height(readtable(fullfile(runDir,'prs_broadcast_budget.csv')))==r.num.K);
    % Isolate frontend interpretation on identical powers/channels.
    o=c.models.rangingIntegral;
    o.userIds=r.integralKernel.userIds;
    o.frontendMode='local_symmetric';
    localKernel=msnoma.ranging.integralKernels(r.ch,r.num,c.paper.dll,o);
    local=msnoma.ranging.varianceIntegral(central.solution.powerW,r.ch,localKernel,c.assumed.noisePsdWPerHz);
    difference=exact.delayVarianceS2./local.delayVarianceS2-1;
    [bs,user]=ndgrid(1:r.num.K,o.userIds);
    t=table(bs(:),user(:),exact.delayVarianceS2(:),local.delayVarianceS2(:),difference(:), ...
        'VariableNames',{'baseStation','positioningUser','globalDelayVarianceS2', ...
        'localSymmetricDelayVarianceS2','relativeGlobalMinusLocal'});
    writetable(t,fullfile(runDir,'dll_frontend_sensitivity.csv'));
    save(fullfile(runDir,'dll_frontend_sensitivity.mat'),'localKernel','local');
    metrics(b)=struct('effectiveBandwidthHz',r.num.BHz, ...
        'maxQuadratureChange',max(r.integralKernel.relativeChange), ...
        'maxTwoTermRelativeVarianceError',max(abs(exact.twoTermRelativeVarianceError(:))), ...
        'maxThreeTermRelativeVarianceError',max(abs(exact.threeTermRelativeVarianceError(:))), ...
        'maxFrontendRelativeVarianceDifference',max(abs(difference(:))), ...
        'prsJointDutyCycle',p.budget.jointDutyCycle, ...
        'prsMeanPositionBoundM',mean(p.position.positionStdM(p.position.valid)));
    runs{b}=runDir;
    fprintf('Verified precision stage %.2f MHz: %s\n',r.num.BHz/1e6,runDir);
end
verification=struct('date','2026-11-01','basicTests',height(basic), ...
    'precisionTests',height(precision),'totalTests',height(basic)+height(precision), ...
    'codeAnalyzerDiagnostics',0,'run50',runs{1},'run20',runs{2},'metrics',metrics);
save(fullfile(codeRoot,'docs','STAGE2_VERIFICATION.mat'),'verification','basic','precision');
fid=fopen(fullfile(codeRoot,'docs','STAGE2_VERIFICATION.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,'PrettyPrint',true));
disp(verification);
end
