function factor = factorTwoTerm(ch,num,noisePsdWPerHz,dll)
%FACTORTWOTERM Eq10 factor v_tau*p; noise and communication terms retained.
a = dll.loopBandwidthHz*(1-0.5*dll.loopBandwidthHz*dll.coherentTimeS);
assert(a>0,'Invalid DLL parameters.');
gCeq = msnoma.ranging.equivalentCommunicationGain(ch,num);
noiseTerm = noisePsdWPerHz/(num.BfeHz*num.TpS);
commTerm = num.BHz*num.G*num.PcW*sum(gCeq,2)'/num.BfeHz^2;
factor.noiseRtauS2W = a*num.TpS^2./(2*ch.gP)*noiseTerm;
factor.communicationRtauS2W = a*num.TpS^2./(2*ch.gP).*commTerm;
factor.RtauS2W = factor.noiseRtauS2W+factor.communicationRtauS2W;
factor.RrangeM2W = num.c0Mps^2*factor.RtauS2W;
factor.equivalentCommunicationGain = gCeq;
factor.aPerS = a;
factor.c0Mps = num.c0Mps;
factor.model = 'paper_Eq10_two_term_nominal_largeN_form_using_effective_B';
end
