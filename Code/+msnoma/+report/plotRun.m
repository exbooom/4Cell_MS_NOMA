function figures = plotRun(run,options)
%PLOTRUN Present stored data only; invalid users are explicitly marked.
if nargin<2, options=run.cfg.io; end
visible = 'off';
if options.showFigures, visible='on'; end
figures = gobjects(0);
count = numel(run.algorithms);
if count==0, return; end
labels = strrep({run.algorithms.id},'_',' ');
labels(strcmp({run.algorithms.id},'cp4a_consistent_candidate'))={'Hierarchical candidate'};
labels(strcmp({run.algorithms.id},'cp4a_subgradient_candidate'))={'Projected candidate'};
width=1240; height=440;
if count>3, height=780; end
figures(1) = figure('Name','4Cell topology and position bounds','Visible',visible, ...
    'Color','w','Position',[80,100,width,height]);
allBounds=[];
for i=1:count
    e=run.algorithms(i).evaluation;
    if ~isempty(e), allBounds=[allBounds;e.positionStdM(e.valid)]; end %#ok<AGROW>
end
tiledlayout(ceil(count/3),min(count,3));
for i=1:count
    nexttile; hold on;
    scene=run.scene;
    plot(scene.bsPosM(:,1),scene.bsPosM(:,2),'ks','MarkerFaceColor','k','MarkerSize',8);
    e=run.algorithms(i).evaluation;
    if ~isempty(e)
        good=e.valid;
        scatter(scene.pUserPosM(good,1),scene.pUserPosM(good,2),55,e.positionStdM(good),'filled');
        plot(scene.pUserPosM(~good,1),scene.pUserPosM(~good,2),'rx','LineWidth',2,'MarkerSize',10);
        colorbar;
        if ~isempty(allBounds) && max(allBounds)>min(allBounds)
            clim([min(allBounds),max(allBounds)]);
        end
        title(sprintf('%s | coverage %.0f%%',labels{i},100*e.coverage));
    else
        title([labels{i} ': no solver result']);
    end
    axis equal; xlim(scene.domainM(1,:)); ylim(scene.domainM(2,:)); grid on;
    xlabel('x (m)'); ylabel('y (m)');
end
sgtitle('Assumed 4-cell scenario: shared color scale in meters; red x = unavailable');
figures(2) = figure('Name','Power allocation and precision','Visible',visible, ...
    'Color','w','Position',[100,100,980,700]);
tiledlayout(2,1);
nexttile; hold on;
for i=1:count
    plot(1:run.num.M,sum(run.algorithms(i).solution.powerW,1),'LineWidth',1.2);
end
ylabel('Sum of BS powers per user (W)'); xlabel('Positioning user');
grid on; legend(labels,'Location','best');
nexttile; hold on;
for i=1:count
    a=run.algorithms(i);
    values=nan(run.num.M,1);
    if ~isempty(a.evaluation), values=a.evaluation.positionStdM; end
    plot(1:run.num.M,values,'o-','LineWidth',1.0);
end
ylabel('Position bound on valid users (m)'); xlabel('Positioning user');
grid on; legend(labels,'Location','best');
sgtitle('Power and position bounds on each algorithm''s valid users');
figures(3) = figure('Name','Communication QoS and coverage CDF','Visible',visible, ...
    'Color','w','Position',[120,100,1080,480]);
tiledlayout(1,2);
nexttile; hold on;
for i=1:count
    a=run.algorithms(i);
    if isempty(a.evaluation), continue; end
    sortedBer=sort(a.evaluation.ber(:));
    semilogy(linspace(0,1,numel(sortedBer)),max(realmin,sortedBer),'LineWidth',1.2);
end
yline(run.cfg.paper.berThreshold,'k--','QoS threshold');
set(gca,'YScale','log'); ylim([1e-8,0.2]);
xlabel('Communication-link quantile'); ylabel('BER'); grid on;
legend(labels,'Location','best');
nexttile; hold on;
for i=1:count
    a=run.algorithms(i);
    if isempty(a.evaluation), continue; end
    [x,F]=msnoma.report.unconditionalCdf(a.evaluation.positionStdM,a.evaluation.valid);
    stairs([0;x],[0;F],'LineWidth',1.2);
end
xlabel('Position bound (m)'); ylabel('Unconditional CDF (end = coverage)');
ylim([0,1.02]); grid on; legend(labels,'Location','best');
sgtitle('ASSUMED profile | relative hearability only');
central=find(strcmp({run.algorithms.id},'central'),1);
if ~isempty(central) && ~isempty(run.algorithms(central).solution.history)
    a=run.algorithms(central);
    figures(4)=figure('Name','Central reference solution','Visible',visible, ...
        'Color','w','Position',[140,100,980,700]);
    tiledlayout(2,1);
    nexttile;
    plot(1:run.num.M,a.solution.powerW','o-');
    xlabel('Positioning user'); ylabel('Per-BS power (W)'); grid on;
    legend(compose('BS %d',1:run.num.K),'Location','best');
    nexttile;
    h=a.solution.history;
    semilogy(h(:,1),h(:,2),'LineWidth',1.2);
    xlabel('fmincon iteration'); ylabel('Scaled mean position variance'); grid on;
    sgtitle('Centralized convex model: this is not the CP4A implementation');
end
end
