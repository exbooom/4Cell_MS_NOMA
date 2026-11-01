function bank = positioningBank(packet,powerW)
%POSITIONINGBANK Known P templates only; communication data are not read.
[K,M]=size(powerW);
assert(isequal([K,M],[size(packet.chips,1),size(packet.chips,2)]));
assert(all(isfinite(powerW) & powerW>=0,'all'));
bank=complex(zeros(K,packet.fftLength));
for k=1:K
    x=complex(zeros(1,packet.activeSamples));
    for m=1:M
        x=x+sqrt(powerW(k,m))*msnoma.waveform.component(packet,k,m);
    end
    bank(k,:)=fft(x,packet.fftLength).*packet.frontendMask;
end
end
