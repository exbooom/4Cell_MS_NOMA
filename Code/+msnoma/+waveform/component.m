function signal = component(packet,k,m)
%COMPONENT One unit-power rectangular positioning component before frontend.
chips=reshape(packet.chips(k,m,:),1,[]);
pulse=repelem(chips,packet.samplesPerChip);
t=(0:packet.activeSamples-1)/packet.sampleRateHz;
signal=pulse.*exp(1i*2*pi*packet.positioningCarrierHz(m)*t);
end
