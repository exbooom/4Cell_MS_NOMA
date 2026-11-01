function figures = plotOptimization(run,options)
%PLOTOPTIMIZATION Keep infeasible dual iterates separate from certificates.
ids={run.algorithms.id};
d=find(strcmp(ids,'dual'),1);
h=find(strcmp(ids,'cp4a_consistent_candidate'),1);
figures=gobjects(0);
if isempty(d) || isempty(h), return; end
a=run.algorithms(d).solution; b=run.algorithms(h).solution;
if isempty(a.history) || isempty(b.history), return; end
visible='off'; if options.showFigures, visible='on'; end
figures=figure('Name','Full-coupling dual and hierarchical candidate', ...
    'Visible',visible,'Color','w','Position',[180,100,1180,780]);
tiledlayout(2,2);
nexttile;
plot(a.history(:,1),a.history(:,3),'o-'); hold on;
plot(a.history(:,1),a.history(:,4),'s-');
xlabel('Dual constraint-generation pass'); ylabel('Scaled objective'); grid on;
legend('Recovered power objective (may be infeasible)','Dual lower bound','Location','best');
nexttile;
semilogy(a.history(:,1),max(1e-14,a.history(:,5)),'o-'); hold on;
yline(run.cfg.numerics.feasibilityTolerance,'k--','Feasibility tolerance');
xlabel('Dual constraint-generation pass'); ylabel('Maximum normalized violation'); grid on;
nexttile;
semilogy(b.history(:,1),max(1e-14,abs(b.history(:,4))),'o-'); hold on;
semilogy(b.history(:,1),max(1e-14,b.history(:,3)),'s-');
xlabel('Outer BS-price iteration'); ylabel('Residual / relative gap'); grid on;
legend('Absolute relative gap','Maximum violation','Location','best');
nexttile;
values=[sum(a.powerW,2),sum(b.powerW,2),run.problem.budgetW];
bar(1:run.num.K,values);
xlabel('Base station'); ylabel('Total positioning power (W)'); grid on;
legend('Full dual','Hierarchical candidate','BS budget','Location','best');
sgtitle('Full coupling and fixed oracle geometry: corrected candidate, not literal paper CP4A');
end
