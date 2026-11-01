function figures = plotPrsComparison(run,options)
%PLOTPRSCOMPARISON A partial comparison and explicit communication costs.
figures=gobjects(0); c=run.prsComparison;
if ~isfield(c,'policies') || isempty(c.policies), return; end
visible='off'; if options.showFigures, visible='on'; end
figures=figure('Name','PRS TDoA and service audit','Visible',visible, ...
    'Color','w','Position',[100,70,1250,850]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
semilogy(1:run.num.M,c.msTdoaCommon.positionStdM,'s-'); hold on;
for i=1:2, semilogy(1:run.num.M,c.policies(i).tdoa.positionStdM,'o-'); end
grid on; xlabel('Positioning user'); ylabel('Local TDoA RMS position (m)');
legend('MS steady-state variance proxy','PRS globally muted CRB','PRS other-BS C / PSD approx.', ...
    'Location','best'); title('Common oracle mask; correlated TDoA covariance');
nexttile;
fraction=[c.policies(1).service.communicationOpportunityFraction,c.policies(2).service.communicationOpportunityFraction];
bar(1:run.num.K,fraction); hold on; yline(run.cfg.models.prs.comparison.minimumAttemptedRateRatio,'r--');
ylim([0,1.15]); grid on; xlabel('Base station'); ylabel('C symbol opportunity fraction');
legend('Global muting','Other-BS C','Service requirement','Location','best');
title({'Same QPSK symbol clock';'Muted slots lose data opportunities'});
nexttile;
worst=[max(c.policies(1).service.worstActiveBer,[],2),max(c.policies(2).service.worstActiveBer,[],2)];
% A display floor is used only on the logarithmic chart, never in exported data.
semilogy(1:run.num.K,max(worst(:,1),1e-16),'o-','LineWidth',1.1); hold on;
semilogy(1:run.num.K,max(worst(:,2),1e-16),'s-','LineWidth',1.1);
yline(run.cfg.paper.berThreshold,'r--'); grid on; set(gca,'YMinorGrid','off');
ylim([1e-17,max(.02,1.5*max(worst(:)))]);
xlabel('Base station'); ylabel('Worst active-slot BER (display floor 1e-16)');
legend('Global muting','Other-BS C','BER limit','Location','best');
title({'Carrier-PSD approximation; grid nulls possible';'Silent slots are excluded from BER'});
nexttile;
power=[c.policies(1).service.requiredPowerForEqualCommunicationEnergyW, ...
    c.policies(2).service.requiredPowerForEqualCommunicationEnergyW];
bar(1:run.num.K,power); hold on;
yline(run.cfg.models.prs.comparison.communicationPeakLimitW,'r--'); grid on;
xlabel('Base station'); ylabel('Required active C power (W)');
legend('Global muting','Other-BS C','C peak limit','Location','best');
title({'Equal C energy compensation audit';'Actual active C power stays original'});
sgtitle('PRS partial-match audit: finite observation differs; service and energy flags exported','FontSize',12);
end
