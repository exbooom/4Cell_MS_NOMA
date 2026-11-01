function qos = qosLimit(gServing,PcW,TcS,N0WPerHz,threshold,model)
%QOSLIMIT Exact inverse of the selected forward BER; negative b is retained.
if strcmp(model.id,'qpsk_bit_energy')
    assert(threshold>0 && threshold<0.5);
    qos.limitWPerHz = (gServing*PcW*TcS/2)/erfcinv(2*threshold)^2-N0WPerHz;
else
    assert(threshold>0 && threshold<model.Gamma);
    qos.limitWPerHz = model.gamma*gServing*PcW*TcS/erfcinv(threshold/model.Gamma)-2*N0WPerHz;
end
qos.cleanBer = msnoma.comm.ber(gServing,PcW,TcS,N0WPerHz,zeros(size(gServing)),model);
qos.cleanFeasible = all(qos.limitWPerHz(:)>=0);
qos.threshold = threshold;
qos.model = model;
end
