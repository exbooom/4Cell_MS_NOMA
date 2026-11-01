function report = run_robust_cancellation_tests()
%RUN_ROBUST_CANCELLATION_TESTS Data-only rejection and exact fallback.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
tests={@testStrongCommunication,@testThresholdSensitivity,@testLegacyGate, ...
    @testInvalidDiagnostics,@testScaleInvariance};
names=cellfun(@func2str,tests,'UniformOutput',false); elapsed=zeros(numel(tests),1);
for i=1:numel(tests)
    started=tic; tests{i}(); elapsed(i)=toc(started);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
end

function [s,b,c]=paired(limit,power)
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.sampledReceiver.codeChips=511;
o=c.models.sampledReceiver; o.userIds=1:3;
o.interferenceCancellation.enabled=true;
o.interferenceCancellation.maximumRelativeResidualEnergy=limit;
n=msnoma.model.numerology(c); scene=msnoma.scenario.make(c,n);
ch=msnoma.scenario.channels(scene,n,c.assumed.propagation);
p=msnoma.waveform.makePacket(n,o); powerW=ones(n.K,n.M)*power;
s=msnoma.experiment.sampledReceiver(scene,ch,n,powerW,o,c.assumed.noisePsdWPerHz,p);
o.interferenceCancellation.enabled=false;
b=msnoma.experiment.sampledReceiver(scene,ch,n,powerW,o,c.assumed.noisePsdWPerHz,p);
end

function testStrongCommunication()
[s,b]=paired(.5,1e-7);
assert(~any(cellfun(@(x) x.pCancellationApplied,s.receivers),'all'));
assert(isequaln(s.rangeM,b.rangeM) && isequal(s.validRanges,b.validRanges));
assert(isequaln(s.estimatedPositionM,b.estimatedPositionM));
assert(~s.truthUsedByEstimator && ~s.communicationDataUsedByEstimator);
end

function testThresholdSensitivity()
[s,b]=paired(1e-12,1e-3);
assert(any(cellfun(@(x) x.valid,s.estimatedPChannels)));
assert(~any(cellfun(@(x) x.residualGatePassed,s.estimatedPChannels)));
assert(isequaln(s.rangeM,b.rangeM));
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
msnoma.report.exportPChannels(s,folder,'sampled');
t=readtable(fullfile(folder,'sampled_estimated_P_channels.csv'));
assert(all(t.maximumRelativeResidualEnergy==1e-12) && ~any(t.residualGatePassed));
end

function testLegacyGate()
[s,~,c]=paired(1,1e-3); o=c.models.sampledReceiver.interferenceCancellation;
o.maximumRelativeResidualEnergy=1;
for u=1:numel(s.estimatedPChannels)
    f=s.estimatedPChannels{u}; gate=msnoma.receiver.cancellationGate(f,o);
    expected=f.valid && all(f.gaussianDelayStdSamples<=o.maximumDelayStdSamples) ...
        && all(f.gaussianRelativeGainStd<=o.maximumRelativeGainStd);
    assert(gate.eligible==expected);
end
end

function testInvalidDiagnostics()
c=msnoma.config.base; o=c.models.sampledReceiver.interferenceCancellation;
f=struct('valid',true,'gaussianDelayStdSamples',0,'gaussianRelativeGainStd',0, ...
    'relativeResidualEnergy',NaN);
g=msnoma.receiver.cancellationGate(f,o); assert(~g.eligible);
f.relativeResidualEnergy=.1; f.gaussianRelativeGainStd=NaN;
g=msnoma.receiver.cancellationGate(f,o); assert(~g.eligible);
f.valid=false; g=msnoma.receiver.cancellationGate(f,o);
assert(~g.eligible && strcmp(g.status,'fit_invalid'));
end

function testScaleInvariance()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
o=c.models.sampledReceiver; o.codeChips=511; o.codeMode='seeded_independent_bpsk';
n=msnoma.model.numerology(c); o.speedMps=n.c0Mps;
p=msnoma.waveform.makePacket(n,o);
bank=msnoma.waveform.positioningBank(p,ones(n.K,n.M)*.003);
lag=[.7;1.6;2.3;3.4]; gain=[1;.4i;-.2;.1];
Y=sum(gain.*bank.*exp(-1i*2*pi*p.frequencyHz.*(lag/p.sampleRateHz)),1);
a=msnoma.receiver.estimatePChannels(Y,bank,p,o,1e-12);
b=msnoma.receiver.estimatePChannels(7*Y,7*bank,p,o,49e-12);
assert(a.valid && b.valid && a.cancellationEligible==b.cancellationEligible);
assert(max(abs(a.delayS-b.delayS))*p.sampleRateHz<1e-10);
assert(max(abs(a.gaussianDelayStdSamples-b.gaussianDelayStdSamples))<1e-12);
assert(abs(a.relativeResidualEnergy-b.relativeResidualEnergy)<1e-12);
end
