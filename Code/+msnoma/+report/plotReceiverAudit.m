function figures = plotReceiverAudit(run,options)
%PLOTRECEIVERAUDIT Carrier samples versus stationary average-timing filtering.
figures=gobjects(0); c=run.prsComparison;
if ~isfield(c,'receiverAudit') || ~isfield(c.receiverAudit,'filteredService'), return; end
a=c.receiverAudit; visible='off'; if options.showFigures, visible='on'; end
figures=figure('Name','PRS C receiver filtering audit','Visible',visible, ...
    'Color','w','Position',[70,100,1400,500]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
difference=abs(a.filteredUnitPsdWPerHzPerW-a.pointUnitPsdWPerHzPerW);
for j=1:run.num.K, difference(j,:,j)=0; end
[~,index]=max(difference(:)); [bs,~,tx]=ind2sub(size(difference),index);
nexttile;
p=reshape(a.pointUnitPsdWPerHzPerW(bs,:,tx),[],1);
f=reshape(a.filteredUnitPsdWPerHzPerW(bs,:,tx),[],1);
semilogy(1:run.num.N,max(p,realmin)); hold on;
semilogy(1:run.num.N,max(f,realmin)); grid on;
xlabel('C carrier / receiver index (users differ)'); ylabel('Received interference PSD per TX watt (1/Hz)');
legend('Carrier PSD sample','Rectangular-filter integral','Location','best');
title({sprintf('PRS BS %d to C users served by BS %d',tx,bs),'Same physical channels and PRS power'});
nexttile;
p=a.pointService; f=a.filteredService;
semilogy(1:run.num.K,max(max(p.worstActiveBer,[],2),1e-16),'o-'); hold on;
semilogy(1:run.num.K,max(max(f.worstActiveBer,[],2),1e-16),'s-');
yline(run.cfg.paper.berThreshold,'r--'); grid on;
ylim([1e-17,.6]); xticks(1:run.num.K);
xlabel('Base station'); ylabel('Worst active-slot BER (display floor 1e-16)');
legend('Point PSD model','Average-timing filter model','BER limit','Location','best');
title({'Distinct receiver assumptions','No silent-slot BER in either model'});
nexttile;
reference=sum(p.referenceExpectedCorrectBitRatePerCarrierBps,2);
bar(1:run.num.K,[sum(p.expectedCorrectBitRatePerCarrierBps,2)./reference, ...
    sum(f.expectedCorrectBitRatePerCarrierBps,2)./reference]); hold on; yline(1,'r--'); grid on;
xlabel('Base station'); ylabel('Expected correct bit rate / MS reference');
legend('Point PSD model','Average-timing filter model','MS reference','Location','best');
title({'Same communication opportunity fraction','Expected correct bits are not packet goodput'});
sgtitle('Analytical PRS receiver audit: stationary independent symbols, average timing; sampled receiver not implemented','FontSize',11);
end
