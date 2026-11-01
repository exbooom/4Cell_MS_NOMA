function out = serviceAudit(ch,num,budget,kernel,noisePsdWPerHz,berModel,berThreshold,spec,policy,referenceBer)
%SERVICEAUDIT Explicit time, energy, peak and BER contracts for a PRS schedule.
% Active C power stays at the original value. Required energy compensation
% is audited separately, never silently capped or used as a feasible signal.
assert(any(strcmp(policy,{'global_muting','other_bs_communication'})));
assert(all(budget.feasible));
assert(isequal(size(referenceBer),[num.K,num.N]));
assert(strcmp(berModel.id,'qpsk_bit_energy'),'Service rate currently requires the explicit QPSK model.');
validateattributes(spec.minimumAttemptedRateRatio,{'numeric'},{'scalar','positive','<=',1});
validateattributes(spec.communicationPeakLimitW,{'numeric'},{'scalar','positive','finite'});
idle=max(0,1-budget.jointDutyCycle);
if strcmp(policy,'global_muting')
    fraction=repmat(idle,num.K,1);
else
    fraction=1-budget.dutyCycle;
end
PcTotal=num.N*num.PcW;
noiseBer=msnoma.comm.ber(ch.gCServing,num.PcW,num.TcS,noisePsdWPerHz,zeros(num.K,num.N),berModel);
weightedBer=idle*noiseBer;
worstActiveBer=nan(num.K,num.N);
if idle>0, worstActiveBer=noiseBer; end
slotBer=nan(num.K,num.N,num.K);
if strcmp(policy,'other_bs_communication')
    assert(~isempty(kernel));
    for j=1:num.K
        active=setdiff(1:num.K,j);
        interference=kernel.prsPsdPerWAtCommunication(:,:,j)*budget.instantaneousPowerW(j);
        value=msnoma.comm.ber(ch.gCServing,num.PcW,num.TcS,noisePsdWPerHz,interference,berModel);
        slotBer(active,:,j)=value(active,:);
        weightedBer(active,:)=weightedBer(active,:)+budget.dutyCycle(j)*value(active,:);
        if budget.dutyCycle(j)>0
            for k=active
                if all(isnan(worstActiveBer(k,:)))
                    worstActiveBer(k,:)=value(k,:);
                else
                    worstActiveBer(k,:)=max(worstActiveBer(k,:),value(k,:));
                end
            end
        end
    end
end
meanActiveBer=weightedBer./fraction;
meanActiveBer(fraction==0,:)=NaN;
requiredPower=PcTotal./fraction;
operating=PcTotal<=spec.communicationPeakLimitW*(1+1e-12) & all(budget.feasible);
activeQos=all(worstActiveBer<=berThreshold*(1+1e-7),2) & fraction>0;
attemptsMet=fraction>=spec.minimumAttemptedRateRatio*(1-1e-12);
energyMet=abs(fraction-1)<=1e-12;
compensationPossible=requiredPower<=spec.communicationPeakLimitW*(1+1e-12) & fraction>0;
out.policy=policy; out.communicationOpportunityFraction=fraction;
out.activeCommunicationPowerW=PcTotal*ones(num.K,1);
out.averageCommunicationPowerW=PcTotal*fraction;
out.requiredPowerForEqualCommunicationEnergyW=requiredPower;
out.energyCompensationPeakFeasible=compensationPossible;
out.operatingPowerFeasible=repmat(operating,num.K,1);
out.communicationEnergyMatched=energyMet;
out.attemptedRateContractMet=attemptsMet;
out.activeSlotBerQosMet=activeQos;
out.serviceContractMet=operating & activeQos & attemptsMet;
out.energyAndServiceMatched=out.serviceContractMet & energyMet;
out.minimumAttemptedRateRatio=spec.minimumAttemptedRateRatio;
out.communicationPeakLimitW=spec.communicationPeakLimitW;
out.meanActiveBer=meanActiveBer; out.worstActiveBer=worstActiveBer;
out.berDuringPrsSlot=slotBer; out.berNoPrs=noiseBer;
out.referenceBer=referenceBer;
bitRate=2/num.TcS;
out.attemptedBitRatePerCarrierBps=bitRate*fraction;
out.referenceAttemptedBitRatePerCarrierBps=bitRate;
out.expectedCorrectBitRatePerCarrierBps=bitRate*(fraction-weightedBer);
out.referenceExpectedCorrectBitRatePerCarrierBps=bitRate*(1-referenceBer);
out.failureReason=cell(num.K,1);
for k=1:num.K
    reasons={};
    if ~operating, reasons{end+1}='active_communication_peak_exceeded'; end %#ok<AGROW>
    if fraction(k)==0, reasons{end+1}='no_communication_opportunity'; end %#ok<AGROW>
    if ~activeQos(k), reasons{end+1}='active_slot_BER_QoS_failed'; end %#ok<AGROW>
    if ~attemptsMet(k), reasons{end+1}='attempted_rate_contract_failed'; end %#ok<AGROW>
    if isempty(reasons), reasons={'service_contract_met'}; end
    out.failureReason{k}=strjoin(reasons,';');
end
out.energyPolicy='original_active_C_power_energy_loss_explicit_no_silent_compensation';
out.cToCInterferenceIncluded=false;
out.rateDefinition='QPSK_symbol_opportunities_and_expected_correct_bits_not_packet_goodput';
out.berDefinition='maximum_over_active_slot_states_silent_slots_not_zero_BER';
out.receiverCaveat='C_interference_uses_carrier_PSD_samples_grid_nulls_are_not_zero_physical_overlap';
end
