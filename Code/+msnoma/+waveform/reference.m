function [spectrum,signal] = reference(packet,k,m)
%REFERENCE Known unit-amplitude positioning template, without channel truth.
signal=msnoma.waveform.component(packet,k,m);
spectrum=fft(signal,packet.fftLength).*packet.frontendMask;
end
