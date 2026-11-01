function figures = plotSubgradient(run,options)
%PLOTSUBGRADIENT Separate accepted outer iterates and all inner trial solves.
index=find(strcmp({run.algorithms.id},'cp4a_subgradient_candidate'),1);
figures=gobjects(0);
if isempty(index), return; end
s=run.algorithms(index).solution;
if ~isfield(s,'innerTrace') || isempty(s.history), return; end
visible='off'; if options.showFigures, visible='on'; end
figures=figure('Name','Projected-price hierarchical experiment', ...
    'Visible',visible,'Color','w','Position',[200,100,1180,780]);
tiledlayout(2,2);
h=s.history;
nexttile;
semilogy(h(:,1),max(1e-14,abs(h(:,4))),'o-'); hold on;
semilogy(h(:,1),max(1e-14,h(:,3)),'s-');
xlabel('Accepted outer update'); ylabel('Absolute gap / maximum violation'); grid on;
legend('Relative duality gap','Normalized constraint violation','Location','best');
nexttile;
plot(h(:,1),s.outerPriceHistory,'o-'); grid on;
xlabel('Accepted outer update'); ylabel('Scaled BS budget price');
legend(compose('BS %d',1:run.num.K),'Location','best');
nexttile;
steps=cellfun(@(x) x.state.iterations,s.innerTrace);
semilogy(1:numel(steps),max(1,steps),'o-'); grid on;
xlabel('Inner call (includes rejected outer trials)'); ylabel('Inner updates (minimum display 1)');
nexttile;
bar(1:run.num.K,[sum(s.powerW,2),run.problem.budgetW]);
xlabel('Base station'); ylabel('Total positioning power (W)'); grid on;
legend('Projected candidate','BS budget','Location','best');
sgtitle(['Projected full-coupling candidate | ' strrep(s.stepRule,'_',' ') ...
    ' | ' strrep(s.status,'_',' ')]);
end
