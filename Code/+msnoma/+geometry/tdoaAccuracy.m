function out = tdoaAccuracy(scene,rangeVarianceM2,spec,bsMask)
%TDOAACCURACY Local GLS bound, correlated reference differences [28] Eq20-23.
% Independent TOA errors, synchronized BS clocks, common receiver clock removed.
% The TOA covariance is frozen at the true location; no covariance derivatives.
K=size(scene.bsPosM,1); M=size(scene.pUserPosM,1);
assert(isequal(size(rangeVarianceM2),[K,M]));
if nargin<4, bsMask=true(K,M); end
assert(islogical(bsMask) && isequal(size(bsMask),[K,M]));
assert(any(strcmp(spec.referencePolicy,{'minimum_variance','fixed'})));
validateattributes(spec.rankTolerance,{'numeric'},{'scalar','positive','<',1});
if strcmp(spec.referencePolicy,'fixed')
    validateattributes(spec.referenceBs,{'numeric'},{'scalar','integer','>=',1,'<=',K});
end
out.positionVarianceM2=nan(M,1); out.positionStdM=nan(M,1);
out.covarianceM2=cell(M,1); out.rangeDifferenceCovarianceM2=cell(M,1);
out.differenceMatrix=cell(M,1); out.H=cell(M,1); out.bsIds=cell(M,1);
out.referenceBs=nan(M,1); out.rank=zeros(M,1);
out.conditionNumber=inf(M,1); out.valid=false(M,1);
out.usableBsMask=bsMask & isfinite(rangeVarianceM2) & rangeVarianceM2>0;
out.failureReason=repmat({'insufficient_usable_BSs'},M,1);
for m=1:M
    ids=find(out.usableBsMask(:,m)); out.bsIds{m}=ids;
    if strcmp(spec.referencePolicy,'fixed') && ~ismember(spec.referenceBs,ids)
        out.failureReason{m}='fixed_reference_unavailable'; continue;
    end
    if numel(ids)<3, continue; end
    v=rangeVarianceM2(ids,m);
    if strcmp(spec.referencePolicy,'fixed')
        r=find(ids==spec.referenceBs);
    else
        [~,r]=min(v);
    end
    out.referenceBs(m)=ids(r);
    offset=scene.pUserPosM(m,:)-scene.bsPosM(ids,:);
    distance=sqrt(sum(offset.^2,2));
    if any(distance==0)
        out.failureReason{m}='user_at_BS_undefined_LOS'; continue;
    end
    other=setdiff(1:numel(ids),r,'stable');
    D=zeros(numel(other),numel(ids));
    D(:,r)=-1; D(sub2ind(size(D),(1:numel(other))',other'))=1;
    C=(D.*v')*D'; H=D*(offset./distance);
    out.differenceMatrix{m}=D; out.H{m}=H;
    out.rangeDifferenceCovarianceM2{m}=C;
    [L,flag]=chol(C,'lower');
    if flag~=0, out.failureReason{m}='difference_covariance_not_positive_definite'; continue; end
    W=L\H; s=svd(W,0);
    out.rank(m)=sum(s>spec.rankTolerance*max(s));
    if out.rank(m)<2, out.failureReason{m}='tdoa_geometry_rank_deficient'; continue; end
    [~,R]=qr(W,0); invR=R\eye(2); covariance=invR*invR';
    if any(~isfinite(covariance(:)))
        out.failureReason{m}='nonfinite_position_covariance'; continue;
    end
    out.covarianceM2{m}=covariance;
    out.positionVarianceM2(m)=trace(covariance);
    out.positionStdM(m)=sqrt(trace(covariance));
    out.conditionNumber(m)=s(1)/s(end); out.valid(m)=true;
    out.failureReason{m}='valid_local_frozen_covariance_bound';
end
out.model='reference28_Eq20_23_correlated_TDoA_local_GLS_bound';
out.clockModel='synchronized_BS_common_receiver_offset_eliminated';
out.geometrySource='oracle'; out.covarianceSpatialDerivativeIncluded=false;
out.detectionCoverageEvaluated=false;
end
