function model = makeBerModel(spec)
%MAKEBERMODEL Matching forward and inverse models, with explicit conventions.
assert(any(strcmp(spec.id,{'qpsk_bit_energy','paper_literal'})),'Unknown BER model.');
model = spec;
if strcmp(model.id,'qpsk_bit_energy')
    model.convention = 'Eb=received_Pc*Tc/2; denominator=N0+I; Gaussian interference candidate';
else
    validateattributes(spec.Gamma,{'numeric'},{'scalar','positive'});
    validateattributes(spec.gamma,{'numeric'},{'scalar','positive'});
    model.convention = 'printed Eq2: Gamma*erfc(gamma*received_Pc*Tc/(I+2*N0))';
end
end
