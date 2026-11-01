function allocation = allocation(spec,bandwidthHz)
%ALLOCATION Explicit generalized comb, or the LTE one-symbol reference [28].
df=spec.deltaFHz;
validateattributes(df,{'numeric'},{'scalar','positive','finite'});
switch spec.mode
    case 'generalized_ofdm'
        Nc=2*floor(bandwidthHz/(2*df));
        assert(Nc>=spec.combSpacing*2,'Requested PRS bandwidth is too small.');
        allK=[(-Nc/2:-1),(1:Nc/2)]';
        activeK=allK(mod(allK,spec.combSpacing)==spec.combOffset);
        symbol=NaN;
        isStandard=false;
    case 'lte_single_symbol_reference'
        assert(~isempty(which('ltePRSIndices')),'LTE Toolbox is required for LTE allocation.');
        assert(any(spec.ndlrb==[6,15,25,50,75,100]),'Invalid LTE PRS RB count.');
        assert(df==15e3,'LTE reference requires 15 kHz subcarrier spacing.');
        Nc=12*spec.ndlrb;
        assert(Nc*df<=bandwidthHz*(1+1e-12),'LTE occupied bandwidth exceeds requested bandwidth.');
        enb=struct('NDLRB',spec.ndlrb,'NPRSRB',spec.ndlrb,'CellRefP',1, ...
            'NCellID',spec.cellId,'NSubframe',spec.subframe,'PRSPeriod','On', ...
            'CyclicPrefix','Normal','DuplexMode','FDD');
        indices=double(ltePRSIndices(enb,'sub'));
        assert(~isempty(indices),'No PRS in the requested reference subframe.');
        symbol=min(indices(:,2));
        carriers=indices(indices(:,2)==symbol,1);
        activeK=carriers-(Nc/2+1)+(carriers>Nc/2);
        isStandard=true;
    otherwise
        error('msnoma:prs:Mode','Unsupported PRS allocation mode.');
end
allocation.mode=spec.mode;
allocation.activeK=activeK(:);
allocation.weightsSquared=ones(numel(activeK),1)*Nc/numel(activeK);
allocation.numCarriers=Nc;
allocation.numPilots=numel(activeK);
allocation.deltaFHz=df;
allocation.symbolPeriodS=1/df;
allocation.occupiedBandwidthHz=Nc*df;
allocation.requestedBandwidthHz=bandwidthHz;
allocation.isStandardLtePattern=isStandard;
allocation.referenceSymbol=symbol;
allocation.model='reference28_Eq5_Eq16_pilot_only_AWGN_no_CP';
end
