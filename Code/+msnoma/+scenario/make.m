function scene = make(cfg,num,stream)
%MAKE Separate positioning users from per-carrier communication receivers.
if nargin<3, stream = RandStream('mt19937ar','Seed',cfg.assumed.seed); end
scene.bsPosM = cfg.paper.bsPosM;
scene.pUserPosM = msnoma.scenario.samplePositions(num.M,scene.bsPosM, ...
    cfg.assumed.domainM,cfg.assumed.minimumDistanceM,stream);
scene.cUserPosM = zeros(num.K,num.N,2);
for k = 1:num.K
    positions = msnoma.scenario.samplePositions(num.N,scene.bsPosM, ...
        cfg.assumed.domainM,cfg.assumed.minimumDistanceM,stream,k);
    scene.cUserPosM(k,:,:) = reshape(positions,1,num.N,2);
end
scene.seed = cfg.assumed.seed;
scene.randomGenerator = 'mt19937ar';
scene.initialStreamSeed = stream.Seed;
scene.domainM = cfg.assumed.domainM;
scene.provenance = cfg.provenance.users;
end
