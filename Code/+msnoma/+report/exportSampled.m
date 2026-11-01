function exportSampled(run,folder)
%EXPORTSAMPLED Finite acquisition, estimated geometry history and step trials.
if isfield(run,'interferenceAblation') && isfield(run.interferenceAblation,'rows')
    writetable(run.interferenceAblation.rows,fullfile(folder,'interference_ablation.csv'));
end
if isfield(run.sampledReceiver,'rangeM')
    s=run.sampledReceiver;
    [ranges,positions]=msnoma.report.sampledTables(s);
    writetable(ranges,fullfile(folder,'sampled_ranges.csv'));
    writetable(positions,fullfile(folder,'sampled_positions.csv'));
    a=s.txAudit; K=size(s.powerW,1);
    packet=table((1:K)',repmat(s.packet.sampleRateHz,K,1), ...
        repmat(s.packet.activeDurationS,K,1),repmat(s.packet.recordDurationS,K,1), ...
        a.nominalPositioningPowerW,a.finitePacketPositioningMeanPowerW, ...
        a.finitePacketCombinedMeanPowerW,'VariableNames', ...
        {'baseStation','sampleRateHz','activeDurationS','recordDurationS', ...
         'nominalPositioningPowerW','finitePositioningPowerW','finiteCombinedPowerW'});
    writetable(packet,fullfile(folder,'sampled_packet_audit.csv'));
    msnoma.report.exportPChannels(s,folder,'sampled');
end
f=run.positionFeedback;
if isfield(f,'history') && ~isempty(f.history)
    L=numel(f.history); rows=cell(L,1); rangeRows=cell(L,1); positionRows=cell(L,1);
    geometryRows=cell(L,1);
    for i=1:L
        h=f.history(i); s=h.solution; a=h.acquisition;
        validFixes=0; meanError=NaN; maximumError=NaN;
        if ~isempty(a)
            msnoma.report.exportPChannels(a,folder,sprintf('feedback_%02d',i));
            validFixes=nnz(a.validPosition);
            values=a.positionErrorM(a.validPosition);
            if ~isempty(values), meanError=mean(values); maximumError=max(values); end
            [r,p]=msnoma.report.sampledTables(a);
            r.iteration=repmat(i,height(r),1); p.iteration=repmat(i,height(p),1);
            rangeRows{i}=r; positionRows{i}=p;
        end
        [k,m]=ndgrid(1:run.num.K,1:run.num.M);
        geometryRows{i}=table(repmat(i,numel(k),1),k(:),m(:),h.lambda(:),s.powerW(:), ...
            repmat(string(h.geometrySource),numel(k),1),'VariableNames', ...
            {'iteration','baseStation','positioningUser','lambda','powerW','geometrySource'});
        rows{i}=table(i,string(h.geometrySource),string(s.status),s.objectiveM2, ...
            validFixes,meanError,maximumError,h.positionChangeM,h.relativePowerChange, ...
            'VariableNames',{'iteration','geometrySource','powerStatus','currentGeometryObjectiveM2', ...
            'validFixes','meanActualPositionErrorM','maxActualPositionErrorM', ...
            'maximumPositionChangeM','relativePowerChange'});
    end
    writetable(vertcat(rows{:}),fullfile(folder,'feedback_iterations.csv'));
    writetable(vertcat(geometryRows{:}),fullfile(folder,'feedback_geometry_powers.csv'));
    nonempty=~cellfun(@isempty,rangeRows);
    if any(nonempty)
        writetable(vertcat(rangeRows{nonempty}),fullfile(folder,'feedback_ranges.csv'));
        writetable(vertcat(positionRows{nonempty}),fullfile(folder,'feedback_positions.csv'));
    end
end
step=run.stepSensitivity;
if isfield(step,'trials') && ~isempty(step.trials)
    rows=cell(numel(step.trials),1);
    for i=1:numel(step.trials)
        t=step.trials(i); s=t.solution; d=s.diagnostics;
        acceptedInner=NaN; acceptedOuter=NaN;
        if isfield(s,'innerTrace')
            % A warm-started final inner solve may accept no new step.
            % Report the last step actually accepted across the whole solve.
            for j=numel(s.innerTrace):-1:1
                history=s.innerTrace{j}.state.history;
                last=find(history(:,7)>0,1,'last');
                if ~isempty(last), acceptedInner=history(last,7); break; end
            end
        end
        if ~isempty(s.history)
            last=find(s.history(:,7)>0,1,'last');
            if ~isempty(last), acceptedOuter=s.history(last,7); end
        end
        rows{i}=table(t.multiplier,t.innerInitialStep,t.outerInitialStep,acceptedInner,acceptedOuter, ...
            string(s.status),s.objectiveM2,s.relativeDualityGap,d.maxNormalizedViolation, ...
            d.maxRelativeHearabilityViolation,d.kktStationarity,d.kktComplementarity,s.elapsedS, ...
            s.innerIterations,s.innerCalls,'VariableNames',{'initialStepMultiplier','innerInitialStep', ...
            'outerInitialStep','lastAcceptedInnerStep','lastAcceptedOuterStep','status','objectiveM2', ...
            'relativeDualityGap','maxNormalizedViolation','maxRelativeHearabilityViolation', ...
            'kktStationarity','kktComplementarity','elapsedS','innerIterations','innerCalls'});
    end
    writetable(vertcat(rows{:}),fullfile(folder,'step_sensitivity.csv'));
end
end
