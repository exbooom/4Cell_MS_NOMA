function value = ber(gServing,PcW,TcS,N0WPerHz,IWPerHz,model)
%BER Forward communication model; use qosLimit with the same model.
assert(all(IWPerHz(:)>=0) && all(N0WPerHz(:)>0));
if strcmp(model.id,'qpsk_bit_energy')
    value = 0.5*erfc(sqrt((gServing*PcW*TcS/2)./(N0WPerHz+IWPerHz)));
else
    value = model.Gamma*erfc(model.gamma*gServing*PcW*TcS./(IWPerHz+2*N0WPerHz));
end
end
