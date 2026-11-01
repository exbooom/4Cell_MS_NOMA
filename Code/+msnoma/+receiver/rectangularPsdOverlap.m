function value = rectangularPsdOverlap(deltaFrequencyHz,transmitSymbolS,receiveSymbolS)
%RECTANGULARPSDOVERLAP Integrated PSD through a unit-area rectangular filter.
% Integral Ttx*sinc^2((f-ftx)*Ttx)*Trx*sinc^2((f-frx)*Trx) df.
% Independent stationary symbols / average timing, not coherent aligned OFDM.
validateattributes(deltaFrequencyHz,{'numeric'},{'real','finite'});
validateattributes(transmitSymbolS,{'numeric'},{'scalar','positive','finite'});
validateattributes(receiveSymbolS,{'numeric'},{'scalar','positive','finite'});
a=min(transmitSymbolS,receiveSymbolS); r=a/max(transmitSymbolS,receiveSymbolS);
z=2*pi*a*deltaFrequencyHz; value=zeros(size(z)); small=abs(z)<.25;
zs=z(small); series=zeros(size(zs));
for q=0:7
    moment=1/(2*q+1)-(1+r)/(2*q+2)+r/(2*q+3);
    series=series+(-1)^q*zs.^(2*q)/factorial(2*q)*moment;
end
value(small)=2*a*series;
zl=z(~small);
value(~small)=2*a./zl.^2.*((1-r)*2*sin(zl/2).^2+2*r*(1-sin(zl)./zl));
assert(all(isfinite(value(:)) & value(:)>=0));
end
