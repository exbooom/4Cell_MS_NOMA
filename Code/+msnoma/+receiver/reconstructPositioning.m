function prediction = reconstructPositioning(bank,fit,packet)
%RECONSTRUCTPOSITIONING Rebuild known P waveforms using estimated channels.
assert(fit.valid,'Only a valid estimated fit may be used for cancellation.');
assert(size(bank,1)==numel(fit.delayS));
prediction=complex(zeros(1,packet.fftLength));
for k=1:size(bank,1)
    prediction=prediction+fit.complexGain(k)*bank(k,:) ...
        .*exp(-1i*2*pi*packet.frequencyHz*fit.delayS(k));
end
end
