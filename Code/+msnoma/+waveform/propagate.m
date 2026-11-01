function received = propagate(tx,packet,gain,delayS,phaseRad,noisePsd,seed)
%PROPAGATE Physical simulation boundary: flat channels, delays, complex AWGN.
% Returned receiver spectrum contains no delay labels or true coordinates.
K=size(tx.spectrum,1);
assert(numel(gain)==K && numel(delayS)==K && numel(phaseRad)==K);
assert(all(gain>=0) && all(delayS>=0) && noisePsd>=0);
Y=complex(zeros(1,packet.fftLength));
for k=1:K
    Y=Y+sqrt(gain(k))*exp(1i*phaseRad(k))*tx.spectrum(k,:) ...
        .*exp(-1i*2*pi*packet.frequencyHz*delayS(k));
end
if noisePsd>0
    stream=RandStream('mt19937ar','Seed',seed);
    w=sqrt(noisePsd*packet.sampleRateHz/2)*(randn(stream,1,packet.fftLength) ...
        +1i*randn(stream,1,packet.fftLength));
    Y=Y+fft(w);
end
received=Y.*packet.frontendMask;
end
