function [gain,distanceM] = freeSpaceGain(txPosM,rxPosM,prop)
%FREESPACEGAIN Receiver-by-transmitter Friis POWER gains, not amplitudes.
assert(size(txPosM,2)==size(rxPosM,2),'Coordinate dimensions differ.');
distanceM = sqrt(sum((reshape(rxPosM,size(rxPosM,1),1,[]) ...
    -reshape(txPosM,1,size(txPosM,1),[])).^2,3));
assert(all(distanceM(:)>=prop.minimumDistanceM*(1-1e-12)), ...
    'msnoma:channel:MinimumDistance','A receiver violates the minimum distance.');
gain = prop.txGainLinear*prop.rxGainLinear ...
    *(prop.speedMps/(4*pi*prop.carrierHz))^2./distanceM.^2;
end
