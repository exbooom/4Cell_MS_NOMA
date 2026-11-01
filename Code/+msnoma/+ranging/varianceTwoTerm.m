function out = varianceTwoTerm(powerW,factor)
%VARIANCETWOTERM Delay in s^2, distance in m^2; zero power is infinite error.
assert(isequal(size(powerW),size(factor.RtauS2W)) && all(powerW(:)>=0));
out.delayVarianceS2 = factor.RtauS2W./powerW;
out.rangeVarianceM2 = factor.RrangeM2W./powerW;
out.rangeStdM = sqrt(out.rangeVarianceM2);
out.noiseDelayVarianceS2 = factor.noiseRtauS2W./powerW;
out.communicationDelayVarianceS2 = factor.communicationRtauS2W./powerW;
out.model = factor.model;
end
