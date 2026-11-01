function out = varianceIntegral(powerW,ch,kernel,N0WPerHz)
%VARIANCEINTEGRAL Evaluate cached Eq6 kernels for every BS of selected users.
validateattributes(N0WPerHz,{'numeric'},{'scalar','positive','finite'});
assert(isequal(size(powerW),size(ch.gP)) && all(powerW(:)>=0));
ids=kernel.userIds;
received=ch.gP(:,ids).*powerW(:,ids);
others=max(0,sum(received,1)-received);
denominator=4*pi^2*received.*kernel.A0PerS.^2;
noiseNumerator=N0WPerHz*kernel.noiseIntegral;
communicationNumerator=kernel.PcW*kernel.communicationIntegralS;
positioningNumerator=others.*kernel.positioningIntegralS;
out.userIds=ids;
out.noiseDelayVarianceS2=kernel.aPerS*noiseNumerator./denominator;
out.communicationDelayVarianceS2=kernel.aPerS*communicationNumerator./denominator;
out.positioningDelayVarianceS2=kernel.aPerS*positioningNumerator./denominator;
out.delayVarianceS2=out.noiseDelayVarianceS2+out.communicationDelayVarianceS2 ...
    +out.positioningDelayVarianceS2;
out.delayVarianceS2(received==0)=Inf;
out.rangeVarianceM2=out.delayVarianceS2*kernel.c0Mps^2;
out.rangeStdM=sqrt(out.rangeVarianceM2);
out.converged=kernel.converged;
out.frontendMode=kernel.frontendMode;
out.discriminator=kernel.discriminator;
out.model=kernel.model;
end
