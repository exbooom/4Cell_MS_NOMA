function figures = plotDllBias(run,options)
%PLOTDLLBIAS Independent interventions and a documented attribution path.
figures=gobjects(0); b=run.dllBias;
if ~isfield(b,'result'), return; end
r=b.result; visible='off'; if options.showFigures, visible='on'; end
figures=figure('Name','DLL controlled variance decomposition','Visible',visible, ...
    'Color','w','Position',[80,60,1250,850]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
value=reshape(r.mainEffectRelativeToBaseline,[],3);
plot(100*value,'LineWidth',1.1); grid on; yline(0,'k:');
xlabel('Selected link (BS varies first)'); ylabel('Predicted variance difference (%)');
legend('Sine linearized','Finite C PSD to spectral lines','Global to local frontend','Location','best');
title({'One factor changed from the common baseline';'Other channels, powers and loop settings fixed'});
nexttile;
value=reshape(r.chainRelativeToBaseline,[],6);
bar(100*mean(value,1)); grid on; yline(0,'k:');
xticks(1:6); xticklabels({'Sine','C linewidth','Frontend','Other asympt.','N to N+1','Drop P'});
xtickangle(25); ylabel('Mean step / baseline variance (%)');
title({'Ordered path to paper two-term variance';'Signed contributions sum to the final relative error'});
nexttile;
values=reshape(r.chainVarianceS2,[],7);
semilogy(run.num.c0Mps*sqrt(values(:,[1,5,7])),'LineWidth',1.1); grid on;
xlabel('Selected link (BS varies first)'); ylabel('Predicted range standard deviation (m)');
legend('Exact sine / finite PSD','Appendix B / exact finite N','Paper two-term / nominal B','Location','best');
title('Variance models evaluated at identical reference power');
nexttile;
changes=arrayfun(@(v) max(v.kernel.relativeChange),b.bundle.variants);
semilogy(1:8,max(changes,1e-16),'o-'); hold on;
yline(b.bundle.variants(1).kernel.quadratureOptions.relativeTolerance,'r--'); grid on;
xticks(1:8); xlabel('Variant index (see CSV for flags)'); ylabel('Last quadrature relative change');
title({'Numerical refinement checks';'Grid convergence is separate from model agreement'});
sgtitle(['DLL variance audit: ' strrep(b.referenceAlgorithm,'_',' ') ...
    '; oracle geometry, selected users; no sampled-delay bias'],'FontSize',12);
end
