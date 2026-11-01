function figures = plotPrecision(run,options)
%PLOTPRECISION Show approximation error and a separately labelled PRS bound.
visible='off';
if options.showFigures, visible='on'; end
figures=gobjects(0);
reference=find(strcmp({run.algorithms.id},'central'),1);
if ~isempty(reference)
    e=run.algorithms(reference).evaluation;
    if ~isempty(e) && isfield(e,'integral')
        v=e.integral;
        figures(end+1)=figure('Name','DLL finite-band integral validation', ...
            'Visible',visible,'Color','w','Position',[150,100,1100,500]);
        tiledlayout(1,2);
        nexttile;
        approximate=e.twoTerm.rangeVarianceM2(:,v.userIds);
        values=[v.rangeStdM(:),sqrt(approximate(:))];
        semilogy(values(:,1),'o-'); hold on; semilogy(values(:,2),'s-');
        xlabel('Selected link (BS varies first)'); ylabel('Range standard deviation (m)');
        legend('Finite-PSD integral','Two-term optimizer model','Location','best'); grid on;
        nexttile;
        plot(100*v.twoTermRelativeVarianceError(:),'o-'); hold on;
        plot(100*v.threeTermRelativeVarianceError(:),'s-');
        yline(0,'k:'); xlabel('Selected link (BS varies first)'); ylabel('Relative variance error (%)');
        legend('Two-term vs integral','Three-term vs integral','Location','best'); grid on;
        sgtitle(['DLL validation: ' strrep(v.frontendMode,'_',' ') ', ' strrep(v.discriminator,'_',' ')]);
    end
end
p=run.prsBenchmark;
if isfield(p,'position') && ~isempty(p.position)
    figures(end+1)=figure('Name','Explicit PRS observation benchmark', ...
        'Visible',visible,'Color','w','Position',[170,100,1100,500]);
    tiledlayout(1,2);
    nexttile;
    semilogy(1:run.num.M,p.position.positionStdM,'o-','LineWidth',1.1); hold on;
    if ~isempty(reference) && ~isempty(run.algorithms(reference).evaluation)
        semilogy(1:run.num.M,run.algorithms(reference).evaluation.positionStdM,'s-','LineWidth',1.1);
        legend('PRS finite-observation AWGN CRB','MS steady-state two-term DLL','Location','best');
    end
    grid on; xlabel('Positioning user'); ylabel('Position bound (m)');
    nexttile;
    yyaxis left; bar(1:run.num.K,p.budget.averagePowerW); ylabel('Matched average PRS power (W)');
    yyaxis right; plot(1:run.num.K,1e3*p.budget.observationTimeS,'ko-'); ylabel('PRS observation (ms)');
    xlabel('Base station'); grid on;
    sgtitle('Assumed broadcast PRS: matched positioning power; communication QoS and estimators differ');
end
end
