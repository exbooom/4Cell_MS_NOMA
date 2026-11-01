function out = receiverAudit(ch,num,allocation,budget,N0,berModel,berThreshold,spec,pointKernel,referenceBer)
%RECEIVERAUDIT C carrier sampling versus average-timing rectangular filtering.
% This is an analytical PSD audit, not generated samples or a timing receiver.
shape=zeros(num.N,1);
for first=1:spec.chunkSize:allocation.numPilots
    ids=first:min(first+spec.chunkSize-1,allocation.numPilots);
    frequency=pointKernel.communicationFrequencyHz-pointKernel.pilotFrequencyHz(ids)';
    overlap=msnoma.receiver.rectangularPsdOverlap(frequency,allocation.symbolPeriodS,num.TcS);
    shape=shape+overlap*(allocation.weightsSquared(ids)/allocation.numCarriers);
end
filtered=pointKernel;
filtered.prsPsdPerWAtCommunication=ch.gPC(:,:,:,1).*reshape(shape,1,num.N,1);
out.pointService=msnoma.prs.serviceAudit(ch,num,budget,pointKernel,N0,berModel, ...
    berThreshold,spec,'other_bs_communication',referenceBer);
out.filteredService=msnoma.prs.serviceAudit(ch,num,budget,filtered,N0,berModel, ...
    berThreshold,spec,'other_bs_communication',referenceBer);
out.filteredService.receiverCaveat='stationary_independent_PRS_symbols_average_timing_rectangular_C_filter';
out.pointUnitPsdWPerHzPerW=pointKernel.prsPsdPerWAtCommunication;
out.filteredUnitPsdWPerHzPerW=filtered.prsPsdPerWAtCommunication;
out.status='computed_analytical_average_timing_filter_audit';
out.assumptions='independent_stationary_symbol_PSD_average_timing_not_coherent_aligned_OFDM';
out.sampledReceiverImplemented=false; out.finiteObservationMatched=false;
out.prsRangingRecomputed=false; out.cToCInterferenceIncluded=false;
out.comparisonCaveat='filtering_changes_the_C_receiver_model_not_a_universal_BER_correction';
end
