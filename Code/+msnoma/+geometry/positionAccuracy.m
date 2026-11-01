function out = positionAccuracy(rangeVarianceM2,geom)
%POSITIONACCURACY Eq1/50; variance is not squared a second time.
M = size(rangeVarianceM2,2);
out.positionVarianceM2 = nan(M,1);
out.positionStdM = nan(M,1);
out.covarianceM2 = cell(M,1);
out.valid = geom.valid;
for m = 1:M
    if ~geom.valid(m), continue; end
    ids = geom.bsIds{m};
    variance = rangeVarianceM2(ids,m);
    if any(~isfinite(variance)), out.valid(m)=false; continue; end
    H = geom.H{m};
    covariance = (H.*variance')*H';
    out.covarianceM2{m} = covariance;
    out.positionVarianceM2(m) = trace(covariance);
    out.positionStdM(m) = sqrt(out.positionVarianceM2(m));
end
end
