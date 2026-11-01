function audit = interferenceAblation(ch,num,powerW,opts,noisePsd,userIds)
%INTERFERENCEABLATION Controlled transmitter masks; truth only in propagation.
% This is a simulation diagnosis, never a selectable oracle receiver.
packet=msnoma.waveform.makePacket(num,opts);
allTx=msnoma.waveform.transmit(packet,powerW,opts.includeCommunication);
labels={'desired_only','same_carrier_P','same_BS_P','all_P','all_P_C','all_P_C_AWGN'};
rows=cell(numel(userIds)*num.K*numel(labels),1); row=0;
for m=userIds(:)'
    delay=ch.positioningDistanceM(:,m)/num.c0Mps;
    phase=-2*pi*opts.carrierHz*delay;
    for k=1:num.K
        S=msnoma.waveform.reference(packet,k,m);
        for v=1:numel(labels)
            mask=zeros(num.K,num.M); includeC=false; N0=0;
            switch v
                case 1, mask(k,m)=1;
                case 2, mask(:,m)=1;
                case 3, mask(k,:)=1;
                otherwise, mask(:)=1;
            end
            if v>=5, includeC=opts.includeCommunication; end
            if v==6, N0=noisePsd*opts.noiseMultiplier; end
            if v>=5
                tx=allTx;
            else
                tx=msnoma.waveform.transmit(packet,powerW.*mask,includeC);
            end
            Y=msnoma.waveform.propagate(tx,packet,ch.gP(:,m),delay,phase,N0, ...
                mod(opts.noiseSeed+104729*m,2^32-1));
            psd=N0*ones(size(Y));
            if includeC, psd=psd+sum(ch.gP(:,m))*packet.communicationPsdWPerHz; end
            if strcmp(opts.whitening,'known_average_CP_psd')
                other=powerW.*mask; other(k,m)=0;
                psd=psd+(ch.gP(:,m)'*other)*packet.positioningPsdPerW;
            end
            detector=N0; weight=ones(size(Y));
            if ~strcmp(opts.whitening,'none')
                detector=max(psd(packet.frontendMask));
                if detector>0
                    weight=sqrt(detector./max(psd,max(N0,detector*opts.whiteningFloorFraction)));
                end
            end
            o=opts; o.speedMps=num.c0Mps;
            r=msnoma.receiver.estimateDelay(Y.*weight,S.*weight,packet,o,detector);
            receivedPower=ch.gP(:,m).*powerW(:,m);
            nearFar=max(receivedPower)/max(receivedPower(k),realmin);
            row=row+1;
            rows{row}=table(m,k,string(labels{v}),includeC,N0,r.valid,r.rangeM, ...
                r.rangeM-ch.positioningDistanceM(k,m),r.normalizedPeak,r.threshold, ...
                nearFar,string(r.status),'VariableNames',{'positioningUser','baseStation', ...
                'transmitterMask','communicationPresent','noisePsdWPerHz','valid', ...
                'estimatedRangeM','rangeErrorForEvaluationM','normalizedPeak','threshold', ...
                'sameCarrierNearFarRatio','status'});
        end
    end
end
audit=struct('status','completed_controlled_transmitter_ablation', ...
    'userIds',userIds,'rows',vertcat(rows{:}),'truthUsedForCancellation',false, ...
    'simulationTransmitterMasks',true,'receiverUsesTrueDelay',false, ...
    'notAnOperationalReceiver',true,'options',opts);
end
