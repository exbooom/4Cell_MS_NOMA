function report = run_cancellation_tests()
%RUN_CANCELLATION_TESTS Independent channel fits and no-truth cancellation.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
tests={@testBank,@testCommunicationNotRead,@testIntegerCommonChannels, ...
    @testFractionalNearFar,@testDescent,@testNoSignal,@testRankFailure, ...
    @testTruthFree,@testSmallChain,@testCancellationDisabled,@testAblation, ...
    @testExport,@testFitFailureFallback,@testPrecisionGate,@testPrecisionScaling};
elapsed=zeros(numel(tests),1); names=cellfun(@func2str,tests,'UniformOutput',false);
for i=1:numel(tests)
    started=tic; tests{i}(); elapsed(i)=toc(started);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
end

function [c,n,p,bank,o]=setup()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.sampledReceiver.codeChips=511;
c.models.sampledReceiver.codeMode='seeded_independent_bpsk';
n=msnoma.model.numerology(c); o=c.models.sampledReceiver; o.speedMps=n.c0Mps;
p=msnoma.waveform.makePacket(n,o);
bank=msnoma.waveform.positioningBank(p,ones(n.K,n.M)*0.003);
end

function testBank()
[~,n,p,bank]=setup(); power=ones(n.K,n.M)*0.003;
tx=msnoma.waveform.transmit(p,power,false);
assert(norm(bank-tx.spectrum.*p.frontendMask,'fro')/norm(bank,'fro')<1e-13);
end

function testCommunicationNotRead()
[~,n,p,bank]=setup(); p.communication(:)=NaN; p.communicationSymbols(:)=NaN;
actual=msnoma.waveform.positioningBank(p,ones(n.K,n.M)*0.003);
assert(isequal(actual,bank));
end

function testIntegerCommonChannels()
[~,n,p,bank,o]=setup(); lags=[1;2;3;4]; amplitude=[1;.4i;-.2;.1*exp(.7i)];
Y=complex(zeros(1,p.fftLength));
for k=1:n.K
    x=zeros(1,p.activeSamples);
    for m=1:n.M, x=x+sqrt(.003)*msnoma.waveform.component(p,k,m); end
    shifted=[zeros(1,lags(k)),x,zeros(1,p.fftLength-numel(x)-lags(k))];
    Y=Y+amplitude(k)*fft(shifted).*p.frontendMask;
end
fit=msnoma.receiver.estimatePChannels(Y,bank,p,o,0);
assert(fit.valid && max(abs(fit.delayS*p.sampleRateHz-lags))<1e-3);
assert(max(abs(fit.complexGain-amplitude))<1e-3);
end

function [fit,Y,bank,p,o]=fractional()
[~,~,p,bank,o]=setup(); lags=[.73;1.61;2.23;3.41];
amplitude=[1;.01*exp(.3i);-.03i;-.005];
Y=sum(amplitude.*bank.*exp(-1i*2*pi*p.frequencyHz.*(lags/p.sampleRateHz)),1);
fit=msnoma.receiver.estimatePChannels(Y,bank,p,o,0);
assert(fit.valid && max(abs(fit.delayS*p.sampleRateHz-lags))<2e-3);
assert(norm(fit.complexGain-amplitude)/norm(amplitude)<2e-4);
end

function testFractionalNearFar()
[fit,Y,bank,p]=fractional(); prediction=msnoma.receiver.reconstructPositioning(bank,fit,p);
assert(norm(prediction-Y)/norm(Y)<1e-5);
end

function testDescent()
[fit,~,~,~,~]=fractional(); assert(all(diff(fit.history(:,2))<=0));
end

function testNoSignal()
[~,~,p,bank,o]=setup(); fit=msnoma.receiver.estimatePChannels(zeros(1,p.fftLength),bank,p,o,0);
assert(~fit.valid && strcmp(fit.status,'signal_absent'));
assert(all(isnan(fit.delayS)));
end

function testRankFailure()
[~,~,p,bank,o]=setup(); bank(2,:)=bank(1,:);
fit=msnoma.receiver.estimatePChannels(bank(1,:),bank,p,o,0);
assert(~fit.valid);
end

function testTruthFree()
for name={'msnoma.receiver.estimatePChannels','msnoma.receiver.reconstructPositioning','msnoma.waveform.positioningBank'}
    body=fileread(which(name{1})); first=strtok(body,newline);
    assert(~contains(first,{'trueDelay','scene','ch,','gainTruth'}));
    assert(~contains(body,{'pUserPosM','positioningDistanceM','packet.communication'}));
end
end

function testSmallChain()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false; c.models.dllBias.enabled=false; c.models.prs.enabled=false;
c.models.sampledReceiver.enabled=true; c.models.sampledReceiver.codeChips=511;
c.models.sampledReceiver.includeCommunication=false; c.models.sampledReceiver.noiseMultiplier=0;
c.models.sampledReceiver.interferenceCancellation.enabled=true;
c.experiment.algorithmIds={'central'};
r=msnoma.experiment.runSingle(c); s=r.sampledReceiver;
assert(all(s.validPosition) && all(s.validRanges,'all'));
assert(max(s.positionErrorM)<0.01 && ~s.communicationDataUsedByEstimator);
assert(all(cellfun(@(x) x.valid,s.estimatedPChannels)));
assert(all(cellfun(@(x) x.pCancellationApplied,s.receivers),'all'));
assert(all(cellfun(@(x) isnan(x.noiseOnlyRangeCrbM2),s.receivers),'all'));
end

function testCancellationDisabled()
c=msnoma.config.base; assert(~c.models.sampledReceiver.interferenceCancellation.enabled);
end

function r=smallRun()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false; c.models.dllBias.enabled=false; c.models.prs.enabled=false;
c.models.sampledReceiver.enabled=true; c.models.sampledReceiver.codeChips=511;
c.models.sampledReceiver.includeCommunication=false; c.models.sampledReceiver.noiseMultiplier=0;
c.models.sampledReceiver.interferenceCancellation.enabled=true;
c.experiment.interferenceAblation.enabled=true; c.experiment.interferenceAblation.userIds=1;
c.experiment.algorithmIds={'central'}; r=msnoma.experiment.runSingle(c);
end

function testAblation()
r=smallRun(); a=r.interferenceAblation;
assert(height(a.rows)==24 && a.notAnOperationalReceiver && ~a.truthUsedForCancellation);
desired=a.rows(a.rows.transmitterMask=="desired_only",:);
assert(all(desired.valid) && max(abs(desired.rangeErrorForEvaluationM))<.01);
assert(~any(a.rows.communicationPresent) && all(a.rows.noisePsdWPerHz==0));
end

function testExport()
r=smallRun(); folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
msnoma.report.exportSampled(r,folder);
t=readtable(fullfile(folder,'sampled_estimated_P_channels.csv'));
assert(height(t)==12 && all(t.fitValid));
t=readtable(fullfile(folder,'sampled_ranges.csv'));
assert(all(t.pCancellationApplied));
assert(height(readtable(fullfile(folder,'interference_ablation.csv')))==24);
end

function testFitFailureFallback()
[c,n,p,~,o]=setup(); scene=msnoma.scenario.make(c,n);
ch=msnoma.scenario.channels(scene,n,c.assumed.propagation);
o.userIds=1; o.maximumRangeM=.01; o.includeCommunication=false; o.noiseMultiplier=0;
o.interferenceCancellation.enabled=true;
s=msnoma.experiment.sampledReceiver(scene,ch,n,ones(n.K,n.M)*1e-3,o,0,p);
assert(~s.estimatedPChannels{1}.valid);
assert(~any(cellfun(@(x) x.pCancellationApplied,s.receivers),'all'));
o.interferenceCancellation.enabled=false;
b=msnoma.experiment.sampledReceiver(scene,ch,n,ones(n.K,n.M)*1e-3,o,0,p);
assert(isequaln(s.rangeM,b.rangeM) && isequal(s.validRanges,b.validRanges));
end

function testPrecisionGate()
[fit,Y,bank,p,o]=fractional(); assert(fit.cancellationEligible);
large=msnoma.receiver.estimatePChannels(Y,bank,p,o,1e3);
assert(large.valid && ~large.cancellationEligible);
assert(strcmp(large.eligibilityStatus,'Gaussian_precision_gate_failed'));
assert(isequal(large.delayS,fit.delayS));
end

function testPrecisionScaling()
[~,Y,bank,p,o]=fractional();
a=msnoma.receiver.estimatePChannels(Y,bank,p,o,1e-12);
b=msnoma.receiver.estimatePChannels(Y,bank,p,o,4e-12);
assert(max(abs(b.gaussianDelayStdSamples./a.gaussianDelayStdSamples-2))<1e-12);
assert(max(abs(b.gaussianRelativeGainStd./a.gaussianRelativeGainStd-2))<1e-12);
end
