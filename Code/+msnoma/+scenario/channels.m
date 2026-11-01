function ch = channels(scene,num,prop)
%CHANNELS Logical channel tensors retain different receiver identities.
[gp,dp] = msnoma.scenario.freeSpaceGain(scene.bsPosM,scene.pUserPosM,prop);
ch.gP = gp';                                  % K-by-M
ch.positioningDistanceM = dp';
rx = reshape(scene.cUserPosM,num.K*num.N,2);
[gc,dc] = msnoma.scenario.freeSpaceGain(scene.bsPosM,rx,prop);
byReceiver = reshape(gc,num.K,num.N,num.K);
ch.gCServing = zeros(num.K,num.N);
for k = 1:num.K
    ch.gCServing(k,:) = byReceiver(k,:,k);
end
ch.gPC = repmat(byReceiver,1,1,1,num.M);         % K-by-N-by-K-by-M
ch.gCP = repmat(reshape(gp,num.M,num.K,1),1,1,num.N);
ch.communicationDistanceM = reshape(dc,num.K,num.N,num.K);
ch.model = 'assumed_2D_Friis_power_gain_frequency_flat_no_fading';
end
