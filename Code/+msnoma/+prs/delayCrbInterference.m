function out = delayCrbInterference(allocation,gP,budget,noisePsdWPerHz,kernel,c0Mps)
%DELAYCRBINTERFERENCE Known-phase information with frequency-dependent PSD.
assert(all(isfinite(gP(:)) & gP(:)>=0) && all(budget.feasible));
validateattributes(noisePsdWPerHz,{'numeric'},{'scalar','positive','finite'});
I=kernel.communicationPsdWPerHz;
assert(size(I,1)==size(gP,1) && size(I,2)==size(gP,2));
assert(all(isfinite(I(:)) & I(:)>=0));
assert(size(I,3)==numel(allocation.activeK));
f=allocation.activeK*allocation.deltaFHz;
moment=reshape(allocation.weightsSquared/allocation.numCarriers.*f.^2,1,1,[]);
weightedMoment=sum(moment./(noisePsdWPerHz+I),3);
information=8*pi^2*allocation.symbolPeriodS* ...
    (budget.instantaneousPowerW.*gP).*budget.symbolCount.*weightedMoment;
out.fisherInformationPerS2=information;
out.delayVarianceS2=1./information;
out.rangeVarianceM2=c0Mps^2*out.delayVarianceS2;
out.rangeStdM=sqrt(out.rangeVarianceM2);
out.phaseModel='reference_known_phase_no_nuisance_phase_projection';
out.model='known_signal_diagonal_frequency_noise_PSD_approximation';
out.fullColoredNoiseLikelihood=false;
out.communicationInterferenceIncluded=true;
end
