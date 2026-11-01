function runDir = saveRun(run,outputRoot)
%SAVERUN Preserve inputs, raw results, tabular metrics and provenance.
if nargin<2, outputRoot=run.cfg.io.resultsRoot; end
if ~exist(outputRoot,'dir'), mkdir(outputRoot); end
runId = char(java.util.UUID.randomUUID());
stamp=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss'));
runDir = fullfile(outputRoot,['run_' stamp '_' runId(1:8)]);
mkdir(runDir);
save(fullfile(runDir,'run.mat'),'run','-v7.3');
writetable(run.summary,fullfile(runDir,'summary.csv'));
scene=run.scene;
coordinates=table((1:run.num.M)',scene.pUserPosM(:,1),scene.pUserPosM(:,2), ...
    'VariableNames',{'positioningUser','xM','yM'});
writetable(coordinates,fullfile(runDir,'positioning_coordinates.csv'));
writeJson(fullfile(runDir,'config.json'),run.cfg);
for i=1:numel(run.algorithms)
    a=run.algorithms(i);
    [k,m]=ndgrid(1:run.num.K,1:run.num.M);
    powers=table(k(:),m(:),a.solution.powerW(:), ...
        'VariableNames',{'baseStation','positioningUser','powerW'});
    writetable(powers,fullfile(runDir,[a.id '_powers.csv']));
    msnoma.report.exportConstraints(runDir,a,run.problem);
    if ~isempty(a.evaluation)
        users=table((1:run.num.M)',a.evaluation.valid, ...
            a.evaluation.positionStdM,a.evaluation.hearability.count, ...
            a.evaluation.failureReason,'VariableNames', ...
            {'positioningUser','valid','positionStdM','hearableBs','failureReason'});
        writetable(users,fullfile(runDir,[a.id '_users.csv']));
        [bs,user]=ndgrid(1:run.num.K,1:run.num.M);
        e=a.evaluation;
        ranging=table(bs(:),user(:),e.delayVarianceS2(:),e.rangeVarianceM2(:), ...
            e.threeTerm.rangeVarianceM2(:),e.hearability.receivedPowerW(:), ...
            e.hearability.mask(:),'VariableNames',{'baseStation','positioningUser', ...
            'twoTermDelayVarianceS2','twoTermRangeVarianceM2', ...
            'threeTermRangeVarianceM2','receivedPowerW','relativeHearable'});
        writetable(ranging,fullfile(runDir,[a.id '_ranging.csv']));
        if ~isempty(a.solution.history)
            h=array2table(a.solution.history,'VariableNames',a.solution.historyColumns);
            writetable(h,fullfile(runDir,[a.id '_solver_history.csv']));
        end
    end
end
msnoma.report.exportPrecision(run,runDir);
msnoma.report.exportOptimization(run,runDir);
msnoma.report.exportSubgradient(run,runDir);
msnoma.report.exportPrsComparison(run,runDir);
msnoma.report.exportDllBias(run,runDir);
msnoma.report.exportReceiverAudit(run,runDir);
msnoma.report.exportSampled(run,runDir);
manifest = struct;
manifest.schemaVersion = run.schemaVersion;
manifest.timezone = 'Asia/Shanghai';
manifest.paperDoi = run.cfg.paper.doi;
manifest.label = run.label;
manifest.profile = run.cfg.profile;
manifest.requestedBandwidthHz = run.num.requestedBandwidthHz;
manifest.effectiveBandwidthHz = run.num.effectiveBandwidthHz;
manifest.numerology = struct('K',run.num.K,'M',run.num.M,'N',run.num.N,'G',run.num.G);
manifest.provenance = run.cfg.provenance;
manifest.feasibilityStatus = run.feasibility.status;
manifest.algorithms = {run.algorithms.id};
manifest.optimizationLabel=run.optimizationLabel;
manifest.carrierGroupAudit=msnoma.optim.carrierGroupAudit(run.num);
manifest.matlabVersion = version;
manifest.toolboxVersions = ver;
manifest.summary = table2struct(run.summary);
manifest.prsStatus=run.prsBenchmark.status;
manifest.prsComparisonStatus=run.prsComparison.status;
manifest.dllBiasStatus=run.dllBias.status;
manifest.sampledReceiver=struct('status',run.sampledReceiver.status, ...
    'implemented',true,'steadyStateDll',false,'prsObservationMatched',false);
if isfield(run.sampledReceiver,'validCoverage')
    s=run.sampledReceiver;
    manifest.sampledReceiver.userIds=s.userIds;
    manifest.sampledReceiver.validCoverage=s.validCoverage;
    manifest.sampledReceiver.sampleRateHz=s.packet.sampleRateHz;
    manifest.sampledReceiver.activeDurationS=s.packet.activeDurationS;
    manifest.sampledReceiver.codeMode=s.packet.codeMode;
    manifest.sampledReceiver.codeReusedAcrossDifferentCarriers=s.packet.codeReusedAcrossDifferentCarriers;
    manifest.sampledReceiver.truthUsedByEstimator=false;
    manifest.sampledReceiver.estimatedPCancellation=s.options.interferenceCancellation.enabled;
    manifest.sampledReceiver.allPositioningCodesKnown=s.allPositioningCodesKnown;
    manifest.sampledReceiver.communicationDataUsedByEstimator=false;
    manifest.sampledReceiver.cancellationUncertaintyPropagated=false;
    manifest.sampledReceiver.maximumRelativeResidualEnergy= ...
        s.options.interferenceCancellation.maximumRelativeResidualEnergy;
    manifest.sampledReceiver.cancellationGateUsesTruth=false;
    manifest.sampledReceiver.cancellationGateIsAccuracyGuarantee=false;
end
if isfield(run,'interferenceAblation')
    manifest.interferenceAblationStatus=run.interferenceAblation.status;
end
manifest.positionFeedback=struct('status',run.positionFeedback.status);
if isfield(run.positionFeedback,'iterations')
    manifest.positionFeedback.iterations=run.positionFeedback.iterations;
    manifest.positionFeedback.converged=run.positionFeedback.converged;
    manifest.positionFeedback.firstFixUnitLambda=true;
    manifest.positionFeedback.usedOracleFallback=false;
end
manifest.stepSensitivity=struct('status',run.stepSensitivity.status,'authorStepsKnown',false);
if isfield(run.dllBias,'result')
    b=run.dllBias;
    manifest.dllBias=struct('referenceAlgorithm',b.referenceAlgorithm,'userIds',b.bundle.userIds, ...
        'variantCount',numel(b.bundle.variants),'quadratureAllConverged',b.bundle.allConverged, ...
        'maximumRelativeClosureResidual',b.result.maximumRelativeClosureResidual, ...
        'sampledReceiverImplemented',false,'optimizerChanged',false);
end
if isfield(run.prsComparison,'policies') && ~isempty(run.prsComparison.policies)
    c=run.prsComparison;
    manifest.prsMatching=struct('positioningBroadcastEnergyMatched',c.positioningBroadcastEnergyMatched, ...
        'geometryConventionMatched',c.geometryConventionMatched,'observationMatched',c.observationMatched, ...
        'detectionCoverageEvaluated',c.detectionCoverageEvaluated,'fullyMatched',c.fullyMatched, ...
        'communicationServiceMet',arrayfun(@(p) all(p.service.serviceContractMet),c.policies), ...
        'communicationEnergyMatched',arrayfun(@(p) all(p.service.communicationEnergyMatched),c.policies));
    manifest.prsMatching.communicationReceiverModel=c.kernel.communicationReceiverModel;
    manifest.prsMatching.gridCaveat=c.kernel.gridCaveat;
    if isfield(c,'receiverAudit') && isfield(c.receiverAudit,'filteredService')
        a=c.receiverAudit;
        manifest.receiverAudit=struct('status',a.status,'assumptions',a.assumptions, ...
            'pointWorstActiveBer',max(a.pointService.worstActiveBer(:)), ...
            'filteredWorstActiveBer',max(a.filteredService.worstActiveBer(:)), ...
            'sampledReceiverImplemented',false,'finiteObservationMatched',false);
    end
end
if ~isempty(run.integralKernel)
    manifest.dllIntegral=struct('frontendMode',run.integralKernel.frontendMode, ...
        'discriminator',run.integralKernel.discriminator,'userIds',run.integralKernel.userIds, ...
        'allConverged',all(run.integralKernel.converged), ...
        'maxRelativeChange',max(run.integralKernel.relativeChange));
end
codeRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
manifest.sourceHashes = msnoma.report.sourceHashes(codeRoot);
writeJson(fullfile(runDir,'manifest.json'),manifest);
end

function writeJson(path,value)
fid=fopen(path,'w','n','UTF-8');
assert(fid>0);
cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
