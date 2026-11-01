function runDir = main(cfg)
%MAIN Four-cell theoretical reference and optional sampled/feedback experiments.
% Run from Code: main; or pass an edited cfg made by msnoma.config.base.
if nargin<1, cfg=msnoma.config.base('smoke_assumed'); end
validation=msnoma.config.validate(cfg);
assert(validation.canRun,validation.message);
fprintf('MS-NOMA 4-cell ASSUMED profile; analytical central reference, not CP4A.\n');
if any(strcmp(cfg.experiment.algorithmIds,'cp4a_consistent_candidate'))
    fprintf('Hierarchical candidate: corrected full coupling, fixed oracle; not literal paper Algorithm 1.\n');
end
if any(strcmp(cfg.experiment.algorithmIds,'cp4a_subgradient_candidate'))
    fprintf('Projected-price experiment: explicit step rules, full coupling, no fmincon fallback.\n');
end
run=msnoma.experiment.runSingle(cfg,cfg.experiment.algorithmIds,[]);
disp(run.summary);
fprintf('Feasibility: %s; effective bandwidth: %.3f MHz; elapsed: %.2f s\n', ...
    run.feasibility.status,run.num.effectiveBandwidthHz/1e6,run.elapsedS);
if ~isempty(run.integralKernel)
    fprintf('DLL integral: %s, users %s, max quadrature change %.3g\n', ...
        run.integralKernel.frontendMode,mat2str(run.integralKernel.userIds), ...
        max(run.integralKernel.relativeChange));
end
fprintf('PRS benchmark: %s\n',run.prsBenchmark.status);
fprintf('PRS comparison: %s\n',run.prsComparison.status);
fprintf('Sampled receiver: %s; estimated-geometry feedback: %s\n', ...
    run.sampledReceiver.status,run.positionFeedback.status);
if isfield(run.sampledReceiver,'validCoverage')
    a=run.sampledReceiver;
    fprintf('  Finite packet %.3f ms, %.3f Msps; valid fixes %d/%d; not steady-state DLL.\n', ...
        1e3*a.packet.activeDurationS,1e-6*a.packet.sampleRateHz,nnz(a.validPosition),numel(a.userIds));
end
if isfield(run.prsComparison,'policies')
    for i=1:numel(run.prsComparison.policies)
        p=run.prsComparison.policies(i);
        fprintf('  %s: service=%d, C energy=%d, equal-energy peak feasible=%d; observation matched=0\n', ...
            p.id,all(p.service.serviceContractMet),all(p.service.communicationEnergyMatched), ...
            all(p.service.energyCompensationPeakFeasible));
    end
end
runDir=msnoma.report.saveRun(run,cfg.io.resultsRoot);
fprintf('DLL bias study: %s\n',run.dllBias.status);
if isfield(run.dllBias,'result')
    fprintf('  8 controlled variants; max relative chain closure residual %.3g\n', ...
        run.dllBias.result.maximumRelativeClosureResidual);
end
if isfield(run.prsComparison,'receiverAudit') && isfield(run.prsComparison.receiverAudit,'filteredService')
    a=run.prsComparison.receiverAudit;
    fprintf('Analytical C receiver audit: point BER %.4g; average-timing filtered BER %.4g; sampled receiver=0\n', ...
        max(a.pointService.worstActiveBer(:)),max(a.filteredService.worstActiveBer(:)));
end
if cfg.io.showFigures || cfg.io.saveFigures
    figures=msnoma.report.plotRun(run,cfg.io);
    figures=[figures,msnoma.report.plotPrecision(run,cfg.io)];
    figures=[figures,msnoma.report.plotOptimization(run,cfg.io)];
    figures=[figures,msnoma.report.plotSubgradient(run,cfg.io)];
    figures=[figures,msnoma.report.plotPrsComparison(run,cfg.io)];
    figures=[figures,msnoma.report.plotDllBias(run,cfg.io)];
    figures=[figures,msnoma.report.plotReceiverAudit(run,cfg.io)];
    figures=[figures,msnoma.report.plotSampled(run,cfg.io)];
    if cfg.io.saveFigures
        for i=1:numel(figures)
            exportgraphics(figures(i),fullfile(runDir,sprintf('figure_%02d.png',i)),'Resolution',160);
        end
    end
    if ~cfg.io.showFigures, close(figures); end
end
fprintf('Run directory: %s\n',runDir);
end
