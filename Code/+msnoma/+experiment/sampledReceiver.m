function study = sampledReceiver(scene,ch,num,powerW,opts,noisePsd,packet)
%SAMPLEDRECEIVER Waveforms -> channels -> receiver -> ranges -> actual fix.
% Truth is used only for propagation and separately labelled error metrics.
if nargin<7, packet=msnoma.waveform.makePacket(num,opts); end
ids=opts.userIds; if isempty(ids), ids=1:num.M; end
L=numel(ids); tx=msnoma.waveform.transmit(packet,powerW,opts.includeCommunication);
rangeM=nan(num.K,L); delayS=rangeM; detected=false(num.K,L); valid=detected;
receivers=cell(num.K,L); fixes=cell(L,1); estimated=nan(L,2);
positionValid=false(L,1); noise=noisePsd*opts.noiseMultiplier;
channelFits=cell(L,1); bank=[];
if opts.interferenceCancellation.enabled
    bank=msnoma.waveform.positioningBank(packet,powerW);
end
for u=1:L
    m=ids(u); trueDelay=ch.positioningDistanceM(:,m)/num.c0Mps;
    phase=-2*pi*opts.carrierHz*trueDelay;
    Y=msnoma.waveform.propagate(tx,packet,ch.gP(:,m),trueDelay,phase,noise, ...
        mod(opts.noiseSeed+104729*m,2^32-1));
    cPsd=zeros(size(Y)); pPsd=cPsd;
    if opts.includeCommunication
        cPsd=sum(ch.gP(:,m))*packet.communicationPsdWPerHz;
    end
    if strcmp(opts.whitening,'known_average_CP_psd')
        pPsd=(ch.gP(:,m)'*powerW)*packet.positioningPsdPerW;
    end
    cancelled=Y; fit=[]; cancellationApplied=false;
    cWeight=ones(size(Y)); cDetectorPsd=noise;
    if opts.interferenceCancellation.enabled
        totalC=noise+cPsd;
        if ~strcmp(opts.whitening,'none'), cDetectorPsd=max(totalC(packet.frontendMask)); end
        if cDetectorPsd>0 && ~strcmp(opts.whitening,'none')
            floorC=max(noise,cDetectorPsd*opts.whiteningFloorFraction);
            cWeight=sqrt(cDetectorPsd./max(totalC,floorC));
        end
        fitOpts=opts; fitOpts.speedMps=num.c0Mps;
        fit=msnoma.receiver.estimatePChannels(Y.*cWeight,bank.*cWeight,packet,fitOpts,cDetectorPsd);
        if fit.valid && fit.cancellationEligible
            prediction=msnoma.receiver.reconstructPositioning(bank,fit,packet);
            cancelled=Y-prediction; cancellationApplied=true;
        end
        channelFits{u}=fit;
    end
    for k=1:num.K
        reference=msnoma.waveform.reference(packet,k,m);
        weight=ones(size(Y)); detectorPsd=noise;
        whitened=~strcmp(opts.whitening,'none');
        if whitened
            totalPsd=noise+cPsd;
            if strcmp(opts.whitening,'known_average_CP_psd')
                otherPsd=max(0,pPsd-ch.gP(k,m)*powerW(k,m)*packet.positioningPsdPerW(m,:));
                totalPsd=totalPsd+otherPsd;
            end
            detectorPsd=max(totalPsd(packet.frontendMask));
            if detectorPsd>0
                floorPsd=max(noise,detectorPsd*opts.whiteningFloorFraction);
                weight=sqrt(detectorPsd./max(totalPsd,floorPsd));
            else
                whitened=false;
            end
        end
        observation=Y;
        if cancellationApplied
            desired=fit.complexGain(k)*sqrt(powerW(k,m))*reference ...
                .*exp(-1i*2*pi*packet.frequencyHz*fit.delayS(k));
            observation=cancelled+desired;
            % After successful cancellation, use the residual C+noise proxy.
            % This changes whitening and detection; it is explicitly logged.
            weight=cWeight; detectorPsd=cDetectorPsd; whitened=detectorPsd>0;
        end
        receiverOpts=opts; receiverOpts.speedMps=num.c0Mps;
        r=msnoma.receiver.estimateDelay(observation.*weight,reference.*weight,packet,receiverOpts,detectorPsd);
        r.whitening=opts.whitening; r.whitened=whitened;
        r.pCancellationApplied=cancellationApplied;
        r.pChannelFitStatus='disabled';
        if ~isempty(fit), r.pChannelFitStatus=fit.status; end
        r.gaussianEquivalentRangeCrbM2=r.noiseOnlyRangeCrbM2;
        if whitened
            r.noiseOnlyRangeCrbM2=NaN;
            r.detector='average_PSD_Gaussian_equivalent_not_interference_calibrated';
        end
        if cancellationApplied
            r.noiseOnlyRangeCrbM2=NaN;
            r.whitening='known_average_C_psd_after_estimated_P_cancellation';
            if strcmp(opts.whitening,'none'), r.whitening='none_after_estimated_P_cancellation'; end
            r.detector='Gaussian_proxy_after_data_dependent_P_cancellation_not_false_alarm_calibrated';
        end
        receivers{k,u}=r; rangeM(k,u)=r.rangeM; delayS(k,u)=r.delayS;
        detected(k,u)=r.detected; valid(k,u)=r.valid;
    end
    fixes{u}=msnoma.geometry.estimateRanges(scene.bsPosM,rangeM(:,u),valid(:,u),opts.position);
    estimated(u,:)=fixes{u}.positionM; positionValid(u)=fixes{u}.valid;
end
metadata=rmfield(packet,{'frequencyHz','frontendMask','communication', ...
    'communicationPsdWPerHz','positioningPsdPerW'});
rangeError=rangeM-ch.positioningDistanceM(:,ids); rangeError(~valid)=NaN;
positionError=sqrt(sum((estimated-scene.pUserPosM(ids,:)).^2,2));
positionError(~positionValid)=NaN;
study=struct('status','completed','userIds',ids,'referenceAlgorithm',opts.referenceAlgorithm, ...
    'packet',metadata,'options',opts,'powerW',powerW,'txAudit',rmfield(tx,'spectrum'), ...
    'receivers',{receivers},'estimatedPChannels',{channelFits}, ...
    'allPositioningCodesKnown',opts.interferenceCancellation.enabled, ...
    'communicationDataUsedByEstimator',false, ...
    'estimatedCancellationUncertaintyPropagated',false, ...
    'rangeM',rangeM,'delayS',delayS,'detected',detected, ...
    'validRanges',valid,'fixes',{fixes},'estimatedPositionM',estimated,'validPosition',positionValid, ...
    'rangeErrorM',rangeError,'positionErrorM',positionError, ...
    'validCoverage',mean(positionValid),'truthUsedByEstimator',false, ...
    'sampledReceiverImplemented',true,'steadyStateDllImplemented',false, ...
    'clockModel','synchronized_BS_known_receiver_clock_zero', ...
    'channelModel','known_flat_CSI_for_simulation_unknown_complex_gain_in_delay_estimator', ...
    'observationMatchedToPrs',false,'fullyMatched',false);
end
