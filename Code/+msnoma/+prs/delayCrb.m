function out = delayCrb(allocation,cn0Hz,symbolCount,c0Mps)
%DELAYCRB Reference [28] Eq5, information summed over known repeated symbols.
% SNR=(C/N0)/B and sum(weightsSquared)=Nc. Phase is assumed known.
assert(all(cn0Hz(:)>=0 & isfinite(cn0Hz(:))));
assert(all(symbolCount(:)>=1 & symbolCount(:)==fix(symbolCount(:))));
assert(all(allocation.weightsSquared>=0) && ...
    abs(sum(allocation.weightsSquared)/allocation.numCarriers-1)<1e-12);
moment=sum(allocation.weightsSquared.*allocation.activeK.^2);
assert(moment>0);
snr=cn0Hz/allocation.occupiedBandwidthHz;
information=8*pi^2*snr*moment/allocation.symbolPeriodS^2.*symbolCount;
out.delayVarianceS2=1./information;
out.rangeVarianceM2=c0Mps^2*out.delayVarianceS2;
out.rangeStdM=sqrt(out.rangeVarianceM2);
out.fisherInformationPerS2=information;
out.frequencyMoment=moment;
out.snrLinear=snr;
out.phaseModel='reference_known_phase_no_nuisance_phase_projection';
out.model='PeralRosado2012_reference28_Eq5_with_explicit_information_accumulation';
end
