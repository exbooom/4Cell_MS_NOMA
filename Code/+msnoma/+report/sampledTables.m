function [ranges,positions] = sampledTables(s)
%SAMPLEDTABLES Keep measured results, truth-only diagnostics and failures apart.
[k,u]=ndgrid(1:size(s.rangeM,1),1:numel(s.userIds));
m=reshape(s.userIds(u(:)),[],1); receiver=s.receivers(:);
metric=cellfun(@(x) x.normalizedPeak,receiver);
threshold=cellfun(@(x) x.threshold,receiver);
status=string(cellfun(@(x) x.status,receiver,'UniformOutput',false));
bound=cellfun(@(x) x.noiseOnlyRangeCrbM2,receiver);
equivalentBound=cellfun(@(x) x.gaussianEquivalentRangeCrbM2,receiver);
ranges=table(k(:),m,s.detected(:),s.validRanges(:),s.delayS(:),s.rangeM(:), ...
    s.rangeErrorM(:),metric,threshold,bound,equivalentBound,status,'VariableNames', ...
    {'baseStation','positioningUser','detected','valid','estimatedDelayS','estimatedRangeM', ...
     'rangeErrorForEvaluationM','noiseNormalizedPeak','detectionThreshold', ...
     'noiseOnlyRangeCrbM2','gaussianEquivalentRangeCrbM2','status'});
if isfield(s.receivers{1},'pCancellationApplied')
    ranges.pCancellationApplied=cellfun(@(x) x.pCancellationApplied,receiver);
    ranges.pChannelFitStatus=string(cellfun(@(x) x.pChannelFitStatus,receiver,'UniformOutput',false));
end
residual=cellfun(@(x) x.residualRmsM,s.fixes);
used=cellfun(@(x) numel(x.usedBs),s.fixes);
status=string(cellfun(@(x) x.status,s.fixes,'UniformOutput',false));
positions=table(s.userIds(:),s.validPosition,s.estimatedPositionM(:,1),s.estimatedPositionM(:,2), ...
    s.positionErrorM,residual,used,status,'VariableNames', ...
    {'positioningUser','valid','estimatedXM','estimatedYM','positionErrorForEvaluationM', ...
     'rangeResidualRmsM','usedBs','status'});
end
