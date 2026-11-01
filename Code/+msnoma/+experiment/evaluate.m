function evaluation = evaluate(powerW,context)
%EVALUATE Recompute physical metrics from powers, independently of solver state.
num = context.num;
cfg = context.cfg;
ch = context.ch;
evaluation.interferenceWPerHz = msnoma.comm.interference(powerW,context.A,num);
evaluation.ber = msnoma.comm.ber(ch.gCServing,num.PcW,num.TcS, ...
    cfg.assumed.noisePsdWPerHz,evaluation.interferenceWPerHz,context.berModel);
evaluation.communicationQosSatisfied = all(evaluation.ber(:)<=cfg.paper.berThreshold*(1+1e-7));
evaluation.twoTerm = msnoma.ranging.varianceTwoTerm(powerW,context.rangeFactor);
evaluation.threeTerm = msnoma.ranging.varianceThreeTerm(powerW,ch,num,context.rangeFactor);
evaluation.hearability = msnoma.geometry.hearability(powerW,ch,cfg.models.hearability);
geom = msnoma.geometry.weights(context.scene,evaluation.hearability.mask,cfg.models.geometry);
position = msnoma.geometry.positionAccuracy(evaluation.twoTerm.rangeVarianceM2,geom);
position.valid = position.valid & evaluation.hearability.enoughBs;
position.positionStdM(~position.valid) = NaN;
position.positionVarianceM2(~position.valid) = NaN;
evaluation.position = position;
evaluation.geometry = geom;
evaluation.coverage = mean(position.valid);
evaluation.valid = position.valid;
evaluation.positionStdM = position.positionStdM;
evaluation.delayVarianceS2 = evaluation.twoTerm.delayVarianceS2;
evaluation.rangeVarianceM2 = evaluation.twoTerm.rangeVarianceM2;
evaluation.objectiveAllBsM2 = mean(sum(context.geom.lambda.^2 ...
    .*evaluation.rangeVarianceM2,1));
evaluation.failureReason = repmat({'valid'},num.M,1);
evaluation.failureReason(~geom.valid) = {'geometry_rank_or_range_failure'};
evaluation.failureReason(~evaluation.hearability.enoughBs) = {'insufficient_relative_hearable_BSs'};
if isfield(context,'integralKernel')
    exact=msnoma.ranging.varianceIntegral(powerW,ch,context.integralKernel,cfg.assumed.noisePsdWPerHz);
    ids=exact.userIds;
    exact.twoTermRelativeVarianceError=evaluation.twoTerm.delayVarianceS2(:,ids)./exact.delayVarianceS2-1;
    exact.threeTermRelativeVarianceError=evaluation.threeTerm.delayVarianceS2(:,ids)./exact.delayVarianceS2-1;
    exact.positionStdM=sqrt(sum(geom.lambda(:,ids).^2.*exact.rangeVarianceM2,1))';
    exact.positionStdM(~position.valid(ids))=NaN;
    evaluation.integral=exact;
end
end
