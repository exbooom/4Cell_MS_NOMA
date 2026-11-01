function budget = powerBudget(averagePowerW,measurementPeriodS,peakLimitW,allocation,spec)
%POWERBUDGET Broadcast energy per BS, counted once across all P-users.
validateattributes(averagePowerW,{'numeric'},{'vector','nonnegative','finite'});
validateattributes(measurementPeriodS,{'numeric'},{'scalar','positive','finite'});
validateattributes(peakLimitW,{'numeric'},{'scalar','positive','finite'});
averagePowerW=averagePowerW(:);
energy=averagePowerW*measurementPeriodS;
Ts=allocation.symbolPeriodS;
switch spec.symbolPolicy
    case 'fit_peak_limit'
        assert(strcmp(allocation.mode,'generalized_ofdm'), ...
            'Repeated generic symbols do not represent a standard LTE occasion.');
        count=max(1,ceil(energy/(peakLimitW*Ts)));
    case 'single_reference_symbol'
        count=ones(size(energy));
    otherwise
        error('msnoma:prs:SymbolPolicy','Unsupported PRS observation schedule.');
end
duration=count*Ts;
power=energy./duration;
feasible=duration<=measurementPeriodS*(1+1e-12) & power<=peakLimitW*(1+1e-12);
jointDuty=sum(duration)/measurementPeriodS;
feasible=feasible & jointDuty<=1+1e-12;
budget.averagePowerW=averagePowerW;
budget.energyPerMeasurementJ=energy;
budget.symbolCount=count;
budget.observationTimeS=duration;
budget.instantaneousPowerW=power;
budget.dutyCycle=duration/measurementPeriodS;
budget.measurementPeriodS=measurementPeriodS;
budget.peakLimitW=peakLimitW;
budget.feasible=feasible;
budget.slotStartS=[0;cumsum(duration(1:end-1))];
budget.jointDutyCycle=jointDuty;
budget.energyCount='one_broadcast_energy_per_BS_not_multiplied_by_user_count';
budget.interference='assumed_sequential_BS_slots_all_other_transmissions_muted';
budget.communicationQosMatched=false;
budget.schedule='sequential_BS_pilot_only_symbols_without_CP_or_LTE_occasion_timing';
end
