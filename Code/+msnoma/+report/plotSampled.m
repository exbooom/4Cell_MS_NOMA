function figures = plotSampled(run,io)
%PLOTSAMPLED Actual finite-packet errors and estimated-geometry iterations.
figures=gobjects(0); visibility='off'; if io.showFigures, visibility='on'; end
if isfield(run.sampledReceiver,'rangeM')
    s=run.sampledReceiver; ids=s.userIds; good=s.validPosition;
    f=figure('Visible',visibility,'Position',[80,80,1200,760]);
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile; plot(run.scene.bsPosM(:,1),run.scene.bsPosM(:,2),'ks','DisplayName','BS'); hold on;
    truePosition=run.scene.pUserPosM(ids,:);
    plot(truePosition(:,1),truePosition(:,2),'bo','DisplayName','Truth for evaluation');
    plot(s.estimatedPositionM(good,1),s.estimatedPositionM(good,2),'r+', ...
        'DisplayName','Estimated from ranges');
    axis equal; grid on; legend('Location','best'); xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Actual fixes %d/%d; known receiver clock',nnz(good),numel(good)));
    nexttile; plot(s.rangeErrorM(:),'.-'); grid on;
    xlabel('Link (BS varies first)'); ylabel('Estimated range - true range (m)');
    title('All other P components and optional C traffic present');
    nexttile; plot(ids(good),s.positionErrorM(good),'o-'); grid on;
    xlabel('Positioning user'); ylabel('Actual position error (m)');
    title('Valid estimates only; failures remain missing');
    nexttile;
    metric=cellfun(@(x) x.normalizedPeak,s.receivers);
    threshold=cellfun(@(x) x.threshold,s.receivers);
    semilogy(metric(:),'.-'); hold on; semilogy(threshold(:),'--'); grid on;
    xlabel('Link (BS varies first)'); ylabel('Coarse peak / assumed correlation variance');
    title('Threshold is not calibrated for code/C interference');
    legend('Observed','Assumed-Gaussian threshold','Location','best');
    sgtitle(sprintf('Finite sampled MS acquisition: %.3f ms active; %.2f Msps; not a steady-state DLL', ...
        1e3*s.packet.activeDurationS,1e-6*s.packet.sampleRateHz));
    figures(end+1)=f;
end
feedback=run.positionFeedback;
if isfield(feedback,'history') && ~isempty(feedback.history)
    h=feedback.history; L=numel(h); meanError=nan(1,L); coverage=zeros(1,L);
    for i=1:L
        if isempty(h(i).acquisition), continue; end
        a=h(i).acquisition; coverage(i)=mean(a.validPosition);
        if any(a.validPosition), meanError(i)=mean(a.positionErrorM(a.validPosition)); end
    end
    f=figure('Visible',visibility,'Position',[100,100,1150,720]);
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile; plot(1:L,meanError,'o-'); grid on;
    xlabel('Fix / power solve'); ylabel('Mean actual position error (m)'); title('Truth used only for evaluation');
    nexttile; plot(1:L,coverage,'o-'); ylim([0,1.05]); grid on;
    xlabel('Fix / power solve'); ylabel('Valid position fraction'); title('No truth or oracle fallback');
    nexttile; plot(1:L,[h.positionChangeM],'o-'); grid on;
    xlabel('Fix / power solve'); ylabel('Max estimated position change (m)');
    title('First solve uses lambda = 1');
    nexttile; plot(1:L,[h.relativePowerChange],'o-'); grid on;
    xlabel('Fix / power solve'); ylabel('Relative Frobenius power change');
    title('Geometry is frozen within each certified power solve');
    sgtitle(sprintf('Estimated-geometry experiment: %s; fixed CSI; common noise = %d', ...
        feedback.status,feedback.commonRandomNumbers),'Interpreter','none');
    figures(end+1)=f;
end
if isfield(run.stepSensitivity,'trials') && ~isempty(run.stepSensitivity.trials)
    trials=run.stepSensitivity.trials;
    factor=[trials.multiplier]; elapsed=arrayfun(@(x) x.solution.elapsedS,trials);
    iterations=arrayfun(@(x) x.solution.innerIterations,trials);
    gaps=arrayfun(@(x) abs(x.solution.relativeDualityGap),trials);
    f=figure('Visible',visibility,'Position',[120,120,1100,380]);
    tiledlayout(1,3,'TileSpacing','compact');
    nexttile; semilogx(factor,iterations,'o-'); grid on;
    xlabel('Initial step multiplier'); ylabel('Total inner iterations');
    nexttile; semilogx(factor,elapsed,'o-'); grid on;
    xlabel('Initial step multiplier'); ylabel('Solve elapsed (s)');
    nexttile; loglog(factor,max(gaps,1e-16),'o-'); grid on;
    xlabel('Initial step multiplier'); ylabel('Absolute relative duality gap (floor 1e-16)');
    sgtitle('Engineering initial steps; accepted steps change under backtracking');
    figures(end+1)=f;
end
figures=[figures,msnoma.report.plotCancellation(run,io)];
end
