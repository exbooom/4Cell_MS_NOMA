function geom = weights(scene,bsMask,spec)
%WEIGHTS Appendix A OLS propagation, using QR and explicit rank checks.
K = size(scene.bsPosM,1);
M = size(scene.pUserPosM,1);
assert(isequal(size(bsMask),[K,M]));
assert(spec.dimension==2,'Only 2D known-height geometry is implemented.');
geom.lambda = zeros(K,M);
geom.H = cell(M,1);
geom.bsIds = cell(M,1);
geom.rank = zeros(M,1);
geom.conditionNumber = inf(M,1);
geom.valid = false(M,1);
for m = 1:M
    ids = find(bsMask(:,m));
    geom.bsIds{m} = ids;
    if numel(ids)<2, continue; end
    offsets = scene.pUserPosM(m,:)-scene.bsPosM(ids,:);
    distance = sqrt(sum(offsets.^2,2));
    if any(distance==0), continue; end
    G = offsets./distance;
    s = svd(G,0);
    rankG = sum(s>spec.rankTolerance*max(s));
    geom.rank(m) = rankG;
    if rankG<2, continue; end
    [Q,R] = qr(G,0);
    H = R\Q';
    geom.H{m} = H;
    geom.lambda(ids,m) = sqrt(sum(H.^2,1))';
    geom.conditionNumber(m) = s(1)/s(end);
    geom.valid(m) = true;
end
geom.source = spec.source;
geom.clockModel = spec.clockModel;
geom.model = 'paper_AppendixA_OLS_2D_known_height_candidate';
end
