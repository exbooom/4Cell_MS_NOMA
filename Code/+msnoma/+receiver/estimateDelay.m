function result = estimateDelay(receivedSpectrum,templateSpectrum,packet,opts,noisePsd)
%ESTIMATEDELAY Finite matched-filter acquisition with unknown complex gain.
% Input interface intentionally excludes the true delay, location and gain.
L=packet.fftLength; Fs=packet.sampleRateHz;
assert(numel(receivedSpectrum)==L && numel(templateSpectrum)==L && noisePsd>=0);
maximumDelayS=opts.maximumRangeM/opts.speedMps;
maxLag=ceil(maximumDelayS*Fs);
assert(maxLag+packet.activeSamples<L/2,'Increase packet guard / FFT length.');
cross=receivedSpectrum.*conj(templateSpectrum);
q=opts.coarseOversampling;
% Interpolate the correlation before selecting a narrow fractional peak.
if q==1
    corr=ifft(cross);
else
    padded=[cross(1:L/2),zeros(1,L*(q-1)),cross(L/2+1:end)];
    corr=ifft(padded)*q;
end
energy=sum(abs(templateSpectrum).^2)/L;
assert(energy>0);
searchPoints=floor(maximumDelayS*Fs*q)+1;
score=abs(corr(1:searchPoints)).^2;
[peak,index]=max(score); coarse=(index-1)/q;
noiseVariance=noisePsd*Fs;
threshold=log(searchPoints/opts.falseAlarmProbability);
if noiseVariance>0
    metric=peak/(noiseVariance*energy);
else
    metric=Inf*double(peak>0);
    if peak==0, metric=0; end
end
detected=metric>=threshold;
result=struct('status','not_detected','valid',false,'detected',detected, ...
    'delayS',NaN,'rangeM',NaN,'coarseLagSamples',coarse,'fractionalLagSamples',NaN, ...
    'normalizedPeak',metric,'threshold',threshold,'atSearchBoundary',false, ...
    'complexAmplitude',complex(NaN,NaN),'noiseOnlyRangeCrbM2',NaN, ...
    'detector','discrete_lag_AWGN_union_bound_before_fractional_refinement', ...
    'coarseOversampling',q,'interferenceFalseAlarmCalibrated',false);
if ~detected, return; end
active=packet.frontendMask & abs(cross)>0;
f=packet.frequencyHz(active); z=cross(active)/L;
lo=max(0,coarse-1/q); hi=min(maximumDelayS*Fs,coarse+1/q);
if hi<=lo, result.status='empty_refinement_interval'; return; end
objective=@(lag) -abs(sum(z.*exp(1i*2*pi*f*(lag/Fs)))).^2;
options=optimset('TolX',opts.delayToleranceSamples,'Display','off');
[lag,negativePeak,exitflag]=fminbnd(objective,lo,hi,options);
if exitflag<=0 || ~isfinite(negativePeak)
    result.status='delay_refinement_failed'; return;
end
amplitude=sum(z.*exp(1i*2*pi*f*(lag/Fs)))/energy;
result.fractionalLagSamples=lag; result.delayS=lag/Fs;
result.rangeM=result.delayS*opts.speedMps; result.complexAmplitude=amplitude;
result.atSearchBoundary=lag<=10*opts.delayToleranceSamples ...
    || maximumDelayS*Fs-lag<=10*opts.delayToleranceSamples;
powerSpectrum=abs(templateSpectrum).^2;
fMean=sum(packet.frequencyHz.*powerSpectrum)/sum(powerSpectrum);
information=2*abs(amplitude)^2*(2*pi)^2 ...
    *sum((packet.frequencyHz-fMean).^2.*powerSpectrum)/L;
if information>0
    result.noiseOnlyRangeCrbM2=noiseVariance/information*opts.speedMps^2;
end
if result.atSearchBoundary
    result.status='search_boundary_peak';
else
    result.status='estimated'; result.valid=true;
end
end
