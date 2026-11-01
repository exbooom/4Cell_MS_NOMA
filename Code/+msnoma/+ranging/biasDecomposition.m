function out = biasDecomposition(powerW,ch,num,dll,N0,bundle,factor,geom)
%BIASDECOMPOSITION Fixed-power interventions, paired contrasts and closure.
ids=bundle.userIds; count=numel(bundle.variants);
assert(all(powerW(:,ids)>0,'all') && all(isfinite(powerW(:))));
shape=[num.K,numel(ids),count]; values=zeros(shape);
noise=zeros(shape); communication=zeros(shape); positioning=zeros(shape);
for i=1:count
    v=msnoma.ranging.varianceIntegral(powerW,ch,bundle.variants(i).kernel,N0);
    values(:,:,i)=v.delayVarianceS2; noise(:,:,i)=v.noiseDelayVarianceS2;
    communication(:,:,i)=v.communicationDelayVarianceS2;
    positioning(:,:,i)=v.positioningDelayVarianceS2;
end
flags=vertcat(bundle.variants.flags);
baseline=values(:,:,bundle.baselineIndex);
out.userIds=ids; out.variantIds={bundle.variants.id}; out.flags=flags;
out.delayVarianceS2=values; out.noiseDelayVarianceS2=noise;
out.communicationDelayVarianceS2=communication; out.positioningDelayVarianceS2=positioning;
out.rangeVarianceM2=num.c0Mps^2*values;
out.relativeToBaseline=values./baseline-1;
out.positionStdM=nan(numel(ids),count);
for i=1:count
    p=sqrt(sum(geom.lambda(:,ids).^2.*out.rangeVarianceM2(:,:,i),1))';
    p(~geom.valid(ids))=NaN; out.positionStdM(:,i)=p;
end
out.pairedContrasts=struct('factor',{},'fromIndex',{},'toIndex',{}, ...
    'absoluteDifferenceS2',{},'relativeToFrom',{},'logVarianceRatio',{});
q=0;
for axis=1:3
    origins=find(flags(:,axis)==0)';
    for source=origins
        targetFlags=flags(source,:); targetFlags(axis)=1;
        target=find(all(flags==targetFlags,2)); q=q+1;
        out.pairedContrasts(q)=struct('factor',bundle.factorOrder{axis}, ...
            'fromIndex',source,'toIndex',target, ...
            'absoluteDifferenceS2',values(:,:,target)-values(:,:,source), ...
            'relativeToFrom',values(:,:,target)./values(:,:,source)-1, ...
            'logVarianceRatio',log(values(:,:,target)./values(:,:,source)));
    end
end
out.mainEffectRelativeToBaseline=zeros(num.K,numel(ids),3);
for axis=1:3
    one=zeros(1,3); one(axis)=1; i=all(flags==one,2);
    out.mainEffectRelativeToBaseline(:,:,axis)=values(:,:,i)./baseline-1;
end
% A named telescoping path keeps other Appendix B approximations separate.
pathFlags=[0,0,0;1,0,0;1,1,0;1,1,1]; path=zeros(1,4);
for i=1:4, path(i)=find(all(flags==pathFlags(i,:),2)); end
analytic=msnoma.ranging.asymptoticVariance(powerW,ch,num,N0,dll);
nominal=msnoma.ranging.varianceThreeTerm(powerW,ch,num,factor);
two=msnoma.ranging.varianceTwoTerm(powerW,factor);
out.chainVarianceS2=cat(3,values(:,:,path),analytic.delayVarianceS2(:,ids), ...
    nominal.delayVarianceS2(:,ids),two.delayVarianceS2(:,ids));
out.chainStateIds=[out.variantIds(path),{'AppendixB_asymptotic_finite_N', ...
    'paper_Eq9_nominal_B','paper_Eq10_drop_positioning_interference'}];
out.chainStepIds=[bundle.factorOrder,{'remaining_A0_noise_positioning_asymptotics', ...
    'finite_N_to_nominal_B_Nplus1','drop_positioning_interference'}];
out.chainDifferenceS2=diff(out.chainVarianceS2,1,3);
out.chainRelativeToBaseline=out.chainDifferenceS2./baseline;
out.chainRelativeToPrevious=out.chainDifferenceS2./out.chainVarianceS2(:,:,1:end-1);
out.totalRelativeBias=out.chainVarianceS2(:,:,end)./baseline-1;
out.closureResidualS2=sum(out.chainDifferenceS2,3) ...
    -(out.chainVarianceS2(:,:,end)-out.chainVarianceS2(:,:,1));
out.maximumRelativeClosureResidual=max(abs(out.closureResidualS2./baseline),[],'all');
out.allConverged=bundle.allConverged;
out.frontendIsInterpretationSensitivity=true;
out.biasDefinition='difference_in_predicted_variance_not_empirical_mean_delay_bias';
out.attributionCaveat='paired_effects_depend_on_context_chain_attribution_depends_on_order';
out.positionConvention='fixed_oracle_no_clock_OLS_selected_users_not_detection_coverage';
out.model='fixed_power_factorial_variance_approximation_audit_not_sampled_receiver';
end
