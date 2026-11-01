function report = run_sampled_tests()
%RUN_SAMPLED_TESTS Finite receiver, independent ranges, feedback interfaces.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
tests={@testPacket,@testCommunicationPsd,@testCommunicationNulls,@testPropagation,@testRandomIsolation,@testAuthorCodes,@testIntegerDelay, ...
    @testFractionalDelay,@testNoSignal,@testNoisyDelay,@testRangeFix, ...
    @testMissingAndRank,@testTruthFreeInterfaces,@testSuppliedPacket,@testSmallFeedback,@testFeedbackFailure,@testStepTrials};
names=cellfun(@func2str,tests,'UniformOutput',false); elapsed=zeros(numel(tests),1);
for i=1:numel(tests)
    t=tic; tests{i}(); elapsed(i)=toc(t);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
end

function testPropagation()
[~,n,p]=setup();
x=msnoma.waveform.component(p,1,1); y=msnoma.waveform.component(p,2,1);
tx=struct('spectrum',complex(zeros(n.K,p.fftLength)));
tx.spectrum(1,:)=fft(x,p.fftLength); tx.spectrum(2,:)=fft(y,p.fftLength);
g=[1;0.25;0;0]; delay=[2;5;0;0]/p.sampleRateHz; phase=[0;pi/3;0;0];
Y=msnoma.waveform.propagate(tx,p,g,delay,phase,0,1);
expected=[zeros(1,2),x,zeros(1,p.fftLength-numel(x)-2)] ...
    +0.5*exp(1i*pi/3)*[zeros(1,5),y,zeros(1,p.fftLength-numel(y)-5)];
actual=ifft(Y); expected=ifft(fft(expected).*p.frontendMask);
assert(norm(actual-expected)/norm(expected)<1e-12);
end

function testCommunicationPsd()
[~,n,p]=setup();
f=linspace(-p.sampleRateHz/2,p.sampleRateHz/2,301);
expected=n.PcW*n.TcS*sum(msnoma.model.sincpi(f(:)*n.TcS+(n.N+1)/2-(1:n.N)).^2,2);
actual=msnoma.waveform.communicationPsd(f,n);
assert(max(abs(actual(:)-expected))/max(expected)<1e-12);
f=((1:n.N)-(n.N+1)/2)*n.deltaFcHz;
assert(max(abs(msnoma.waveform.communicationPsd(f,n)/(n.PcW*n.TcS)-1))<1e-12);
end

function [c,n,p] = setup()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.sampledReceiver.codeChips=511;
c.models.sampledReceiver.codeMode='seeded_independent_bpsk';
n=msnoma.model.numerology(c); p=msnoma.waveform.makePacket(n,c.models.sampledReceiver);
end

function testPacket()
[~,n,p]=setup();
assert(p.sampleRateHz==4*n.BHz && p.samplesPerChip==16);
assert(p.actualCodeChips==ceil(511/n.G)*n.G && p.activeSamples==p.actualCodeChips*16);
assert(abs(p.activeDurationS-p.actualCodeChips*n.TpS)<1e-16);
assert(all(abs(p.chips)==1,'all') && all(abs(abs(p.communicationSymbols)-1)<1e-14,'all'));
x=msnoma.waveform.component(p,2,3); assert(all(abs(abs(x)-1)<1e-13));
[S,y]=msnoma.waveform.reference(p,2,3);
assert(isequal(x,y) && norm(S(~p.frontendMask))==0);
assert(p.fftLength>2*p.activeSamples);
end

function testCommunicationNulls()
[~,n,p]=setup();
assert(mod(p.actualCodeChips,n.G)==0);
t=(0:p.activeSamples-1)/p.sampleRateHz;
f=(n.N+3-(n.N+1)/2)*n.deltaFcHz;
projection=sum(p.communication(1,:).*exp(-1i*2*pi*f*t));
assert(abs(projection)/sum(abs(p.communication(1,:)))<1e-12);
end

function testRandomIsolation()
[c,n,p]=setup(); before=rng;
q=msnoma.waveform.makePacket(n,c.models.sampledReceiver); after=rng;
assert(isequal(before,after) && isequal(q,p));
tx=msnoma.waveform.transmit(p,zeros(n.K,n.M),false);
Y=msnoma.waveform.propagate(tx,p,ones(n.K,1),zeros(n.K,1),zeros(n.K,1),1e-8,1);
Z=msnoma.waveform.propagate(tx,p,ones(n.K,1),zeros(n.K,1),zeros(n.K,1),1e-8,1);
assert(isequal(Y,Z) && isequal(before,rng));
end

function testAuthorCodes()
c=msnoma.config.base; c.paper.M=9; c.models.sampledReceiver.codeChips=31;
n=msnoma.model.numerology(c); p=msnoma.waveform.makePacket(n,c.models.sampledReceiver);
assert(p.codeReusedAcrossDifferentCarriers && all(abs(p.chips)==1,'all'));
for m=1:n.M, assert(numel(unique(p.prn(:,m)))==n.K); end
bank=load(p.codeSourcePath,'L1CPcodes');
assert(isequal(reshape(p.chips(2,3,:),1,[]),bank.L1CPcodes(p.prn(2,3),1:p.actualCodeChips)));
end

function testIntegerDelay()
[c,n,p]=setup(); [S,x]=msnoma.waveform.reference(p,1,2);
lag=3; y=[zeros(1,lag),x,zeros(1,p.fftLength-numel(x)-lag)];
Y=fft(y).*p.frontendMask;
o=c.models.sampledReceiver; o.speedMps=n.c0Mps;
r=msnoma.receiver.estimateDelay(Y,S,p,o,0);
assert(r.valid && abs(r.fractionalLagSamples-lag)<5e-5);
end

function testFractionalDelay()
[c,n,p]=setup(); S=msnoma.waveform.reference(p,1,2);
lag=2.37; amplitude=0.3*exp(1i*1.2);
Y=amplitude*S.*exp(-1i*2*pi*p.frequencyHz*(lag/p.sampleRateHz));
o=c.models.sampledReceiver; o.speedMps=n.c0Mps;
r=msnoma.receiver.estimateDelay(Y,S,p,o,0);
assert(r.valid && abs(r.fractionalLagSamples-lag)<5e-5);
assert(abs(r.complexAmplitude/amplitude-1)<1e-4);
end

function testNoSignal()
[c,n,p]=setup(); S=msnoma.waveform.reference(p,1,2);
o=c.models.sampledReceiver; o.speedMps=n.c0Mps;
r=msnoma.receiver.estimateDelay(zeros(size(S)),S,p,o,0);
assert(~r.valid && ~r.detected && isnan(r.rangeM) && r.normalizedPeak==0);
end

function testNoisyDelay()
[c,n,p]=setup(); S=msnoma.waveform.reference(p,1,2); lag=2.4;
stream=RandStream('mt19937ar','Seed',800);
noise=sqrt(1e-3/2)*(randn(stream,1,p.fftLength)+1i*randn(stream,1,p.fftLength));
Y=(S.*exp(-1i*2*pi*p.frequencyHz*(lag/p.sampleRateHz))+fft(noise)).*p.frontendMask;
o=c.models.sampledReceiver; o.speedMps=n.c0Mps;
r=msnoma.receiver.estimateDelay(Y,S,p,o,1e-3/p.sampleRateHz);
assert(r.valid && abs(r.fractionalLagSamples-lag)<0.01 && r.noiseOnlyRangeCrbM2>0);
end

function testRangeFix()
c=msnoma.config.base; b=c.paper.bsPosM; o=c.models.sampledReceiver.position;
for position=[31,47;100,100;178,8]'
    r=sqrt(sum((position'-b).^2,2));
    fix=msnoma.geometry.estimateRanges(b,r,true(4,1),o);
    assert(fix.valid && norm(fix.positionM-position')<1e-7 && fix.residualRmsM<1e-7);
end
end

function testMissingAndRank()
c=msnoma.config.base; o=c.models.sampledReceiver.position;
a=msnoma.geometry.estimateRanges(c.paper.bsPosM,ones(4,1),[1;1;0;0],o);
assert(strcmp(a.status,'insufficient_ranges') && all(isnan(a.positionM)));
a=msnoma.geometry.estimateRanges([0,0;10,0;20,0],ones(3,1),true(3,1),o);
assert(strcmp(a.status,'anchor_rank_deficient'));
outside=[250,250]; ranges=sqrt(sum((outside-c.paper.bsPosM).^2,2));
a=msnoma.geometry.estimateRanges(c.paper.bsPosM,ranges,true(4,1),o);
assert(strcmp(a.status,'outside_known_domain') && norm(a.positionM-outside)<1e-7);
end

function testTruthFreeInterfaces()
for name={'msnoma.receiver.estimateDelay','msnoma.geometry.estimateRanges'}
    text=fileread(which(name{1})); first=strtok(text,newline);
    assert(~contains(first,{'trueDelay','scene','ch,','positionTruth'}));
    assert(~contains(text,'pUserPosM'));
end
end

function testSuppliedPacket()
[c,n,p]=setup(); scene=msnoma.scenario.make(c,n);
ch=msnoma.scenario.channels(scene,n,c.assumed.propagation);
o=c.models.sampledReceiver; o.userIds=1; o.codeChips=1023;
o.includeCommunication=false; o.noiseMultiplier=0;
s=msnoma.experiment.sampledReceiver(scene,ch,n,ones(n.K,n.M)*1e-3,o,0,p);
assert(s.packet.requestedCodeChips==p.requestedCodeChips);
assert(isequal(s.packet.chips,p.chips));
assert(s.options.codeChips~=s.packet.requestedCodeChips);
end

function r = smallRun(maximumRangeM)
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false; c.models.dllBias.enabled=false; c.models.prs.enabled=false;
c.models.sampledReceiver.enabled=true; c.models.sampledReceiver.codeChips=2047;
c.models.sampledReceiver.includeCommunication=false; c.models.sampledReceiver.noiseMultiplier=0;
if nargin>0, c.models.sampledReceiver.maximumRangeM=maximumRangeM; end
c.experiment.positionFeedback.enabled=true; c.experiment.positionFeedback.solverId='central';
c.experiment.positionFeedback.maxFixes=2;
c.experiment.algorithmIds={'central'};
r=msnoma.experiment.runSingle(c);
end

function testFeedbackFailure()
r=smallRun(0.1); f=r.positionFeedback;
assert(strcmp(f.status,'acquisition_or_position_failure') && f.iterations==1);
assert(~f.usedOracleFallback && ~f.truthGeometryUsed && ~f.converged);
assert(~all(f.history(1).acquisition.validPosition));
end

function testSmallFeedback()
r=smallRun(); f=r.positionFeedback;
assert(r.sampledReceiver.sampledReceiverImplemented && ~r.sampledReceiver.truthUsedByEstimator);
assert(all(r.sampledReceiver.validPosition));
assert(f.firstFixUnitLambda && ~f.usedOracleFallback && ~f.truthGeometryUsed);
assert(all(f.history(1).lambda==1,'all') && numel(f.history)==2);
assert(strcmp(f.history(2).geometrySource,'estimated_sampled_ranges'));
assert(all(isfinite(f.history(2).lambda),'all'));
end

function testStepTrials()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false; c.models.dllBias.enabled=false; c.models.prs.enabled=false;
c.experiment.algorithmIds={'central'}; c.experiment.stepSensitivity.enabled=true;
c.experiment.stepSensitivity.multipliers=[0.1,1,10];
r=msnoma.experiment.runSingle(c); s=r.stepSensitivity;
assert(strcmp(s.status,'completed') && numel(s.trials)==3 && ~s.authorStepsKnown);
assert(s.verifiedCount==3 && s.objectiveRelativeSpread<1e-6);
end
