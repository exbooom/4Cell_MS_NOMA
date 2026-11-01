function out = interferenceKernels(ch,num,allocation,spec)
%INTERFERENCEKERNELS Finite OFDM PSD in both directions, sequential PRS slots.
% C-to-C interference is omitted, consistently with the main MS model.
% Diagonal tone noise is an explicit approximation, not a colored likelihood.
validateattributes(spec.chunkSize,{'numeric'},{'scalar','integer','positive'});
assert(isequal(size(ch.gCP),[num.M,num.K,num.N]));
assert(size(ch.gPC,1)==num.K && size(ch.gPC,2)==num.N && size(ch.gPC,3)==num.K);
assert(all(isfinite(ch.gCP(:)) & ch.gCP(:)>=0));
assert(all(isfinite(ch.gPC(:)) & ch.gPC(:)>=0));
fC=(1:num.N)'*num.deltaFcHz;
fBase=allocation.activeK*allocation.deltaFHz;
fPilot=(num.N+1)*num.deltaFcHz/2+fBase;
L=numel(fPilot); out.communicationPsdWPerHz=zeros(num.K,num.M,L);
shape=zeros(num.N,1);
for first=1:spec.chunkSize:L
    ids=first:min(first+spec.chunkSize-1,L);
    cShape=num.TcS*msnoma.model.sincpi((fPilot(ids)-fC')*num.TcS).^2;
    pShape=allocation.symbolPeriodS* ...
        msnoma.model.sincpi((fC-fPilot(ids)')*allocation.symbolPeriodS).^2;
    shape=shape+pShape*(allocation.weightsSquared(ids)/allocation.numCarriers);
    for k=1:num.K
        other=setdiff(1:num.K,k);
        gain=reshape(sum(ch.gCP(:,other,:),2),num.M,num.N);
        psd=num.PcW*(cShape*gain');
        out.communicationPsdWPerHz(k,:,ids)=reshape(psd',1,num.M,numel(ids));
    end
end
out.prsPsdPerWAtCommunication=ch.gPC(:,:,:,1).*reshape(shape,1,num.N,1);
out.activeCommunicationPowerPerBSW=num.PcW*num.N;
out.pilotFrequencyHz=fPilot; out.basebandFrequencyHz=fBase;
out.communicationFrequencyHz=fC;
out.model='finite_rectangular_OFDM_PSD_independent_Gaussian_tone_noise_approximation';
out.prsInterferencePolicy='sequential_PRS_no_simultaneous_other_PRS';
out.communicationInterferencePolicy='other_BS_C_active_at_original_power_own_BS_C_muted';
out.channelAssumption='existing_frequency_flat_gPC_gain_shared_across_positioning_users';
out.communicationReceiverModel='PSD_sample_at_C_carrier_not_matched_filter_band_integral';
out.gridCaveat='aligned_rectangular_symbol_grids_can_have_exact_PSD_sample_nulls';
out.cToCInterferenceIncluded=false; out.toneNoiseCorrelationsIncluded=false;
end
