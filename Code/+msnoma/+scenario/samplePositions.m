function positions = samplePositions(count,bsPosM,domainM,minDistanceM,stream,servingBs)
%SAMPLEPOSITIONS Rejection sampling with an explicit local random stream.
if nargin<6, servingBs = []; end
positions = zeros(count,2);
accepted = 0;
attempts = 0;
while accepted<count
    attempts = attempts+1;
    assert(attempts<=10000,'msnoma:scenario:Sampling','User sampling failed.');
    batch = max(64,2*(count-accepted));
    x = domainM(:,1)'+rand(stream,batch,2).*(domainM(:,2)-domainM(:,1))';
    d2 = sum((reshape(x,[],1,2)-reshape(bsPosM,1,[],2)).^2,3);
    [distance2,nearest] = min(d2,[],2);
    keep = distance2>=minDistanceM^2;
    if ~isempty(servingBs), keep = keep & nearest==servingBs; end
    x = x(keep,:);
    take = min(size(x,1),count-accepted);
    positions(accepted+(1:take),:) = x(1:take,:);
    accepted = accepted+take;
end
end
