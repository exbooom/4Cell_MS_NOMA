function out = asymptoticVariance(powerW,ch,num,noisePsdWPerHz,dll)
%ASYMPTOTICVARIANCE Appendix B Eq56-61 with the exact finite N carrier sum.
assert(isequal(size(powerW),size(ch.gP)) && all(isfinite(powerW(:)) & powerW(:)>=0));
validateattributes(noisePsdWPerHz,{'numeric'},{'scalar','positive','finite'});
a=dll.loopBandwidthHz*(1-.5*dll.loopBandwidthHz*dll.coherentTimeS);
assert(a>0);
received=ch.gP.*powerW;
weightedGain=sum(reshape(sum(ch.gCP,2),num.M,num.N) ...
    .*sin(pi*(1:num.N)/num.G).^2,2)';
other=max(0,sum(received,1)-received);
out.noiseDelayVarianceS2=a*num.TpS*noisePsdWPerHz./(2*num.BfeHz*received);
out.communicationDelayVarianceS2=a*num.TpS*num.PcW*weightedGain./(num.BfeHz^2*received);
out.positioningDelayVarianceS2=a*num.TpS*other./(2*num.BfeHz^2*received);
out.delayVarianceS2=out.noiseDelayVarianceS2+out.communicationDelayVarianceS2+out.positioningDelayVarianceS2;
out.delayVarianceS2(received==0)=Inf;
out.rangeVarianceM2=num.c0Mps^2*out.delayVarianceS2;
out.model='AppendixB_asymptotic_A0_noise_P_interference_exact_finite_N_C_sum';
out.nominalBCommunicationMultiplier=(num.N+1)/num.N;
end
