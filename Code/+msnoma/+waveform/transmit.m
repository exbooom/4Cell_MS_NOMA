function tx = transmit(packet,powerW,includeCommunication)
%TRANSMIT Add independent per-component powers before the common channel.
[K,M]=size(powerW);
assert(isequal([K,M],[size(packet.chips,1),size(packet.chips,2)]));
assert(all(isfinite(powerW) & powerW>=0,'all'));
signal=complex(zeros(K,packet.activeSamples));
for k=1:K
    for m=1:M
        unit=msnoma.waveform.component(packet,k,m);
        signal(k,:)=signal(k,:)+sqrt(powerW(k,m))*unit;
    end
end
positioningMeanPowerW=mean(abs(signal).^2,2);
if includeCommunication, signal=signal+packet.communication; end
tx=struct('spectrum',fft(signal,packet.fftLength,2), ...
    'nominalPositioningPowerW',sum(powerW,2), ...
    'finitePacketPositioningMeanPowerW',positioningMeanPowerW, ...
    'finitePacketCombinedMeanPowerW',mean(abs(signal).^2,2), ...
    'includeCommunication',logical(includeCommunication));
end
