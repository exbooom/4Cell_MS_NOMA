function figures = plotCancellation(run,io)
%PLOTCANCELLATION P fit validity/residuals and controlled transmitter masks.
figures=gobjects(0); visibility='off'; if io.showFigures, visibility='on'; end
s=run.sampledReceiver;
if isfield(s,'allPositioningCodesKnown') && s.allPositioningCodesKnown
    f=figure('Visible',visibility,'Position',[120,120,1150,700]);
    fits=s.estimatedPChannels; valid=cellfun(@(x) x.valid,fits);
    eligible=cellfun(@(x) x.cancellationEligible,fits);
    residual=cellfun(@(x) x.relativeResidualEnergy,fits);
    iterations=cellfun(@(x) x.iterations,fits);
    applied=cellfun(@(x) x.pCancellationApplied,s.receivers);
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile; stem(s.userIds,valid,'DisplayName','Local fit converged'); hold on;
    plot(s.userIds,eligible,'rx','DisplayName','All cancellation gates passed');
    ylim([-.1,1.1]); grid on; legend('Location','best');
    xlabel('Positioning user'); ylabel('Joint P-channel fit valid');
    nexttile; semilogy(s.userIds,max(residual,realmin),'o-'); grid on;
    yline(s.options.interferenceCancellation.maximumRelativeResidualEnergy,'r--','Residual limit');
    xlabel('Positioning user'); ylabel('Weighted residual / received energy');
    nexttile; plot(s.userIds,iterations,'o-'); grid on;
    xlabel('Positioning user'); ylabel('Delay refinement iterations');
    nexttile; imagesc(s.userIds,1:run.num.K,applied); clim([0,1]); colorbar;
    xlabel('Positioning user'); ylabel('BS'); title('Estimated P cancellation applied');
    sgtitle({'Known all-P code/power bank: estimated delays and complex gains', ...
        'No true channels or C data in fit; detection false alarms and cancellation uncertainty uncalibrated'});
    figures(end+1)=f;
end
if isfield(run,'interferenceAblation') && isfield(run.interferenceAblation,'rows')
    a=run.interferenceAblation.rows; modes=unique(a.transmitterMask,'stable');
    f=figure('Visible',visibility,'Position',[140,140,1150,550]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;
    for k=1:run.num.K
        selected=a.baseStation==k;
        plot(find(selected),a.rangeErrorForEvaluationM(selected),'o-','DisplayName',sprintf('BS %d',k)); hold on;
    end
    grid on; legend('Location','best'); xlabel('Ablation row within exported table');
    ylabel('Estimated range - true range (m)');
    nexttile; count=zeros(numel(modes),1);
    for j=1:numel(modes), count(j)=nnz(a.valid(a.transmitterMask==modes(j))); end
    bar(count); xticks(1:numel(modes)); xticklabels(strrep(cellstr(modes),'_',' '));
    xtickangle(35); ylabel('Valid reported ranges (not a correct-lock guarantee)'); grid on;
    sgtitle('Simulation transmitter masks; original single-component receiver; no truth cancellation');
    figures(end+1)=f;
end
end
