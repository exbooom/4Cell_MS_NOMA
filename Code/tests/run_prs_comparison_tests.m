function report = run_prs_comparison_tests()
%RUN_PRS_COMPARISON_TESTS Independent covariance, service and PRS PSD checks.
tests={@testTdoaCovariance,@testReferenceInvariance,@testClockSchur, ...
    @testTdoaScaling,@testTdoaFailures,@testPsdDirections,@testPsdChunking, ...
    @testZeroInterference,@testInterferenceMonotonicity,@testServiceGlobal, ...
    @testServiceOtherBs,@testServiceStrictEnergy,@testServiceZeroOpportunity, ...
    @testServicePeakFailure,@testCarrierGridNullAndOffset,@testActiveBerFailure, ...
    @testComparisonIntegration,@testReferenceUnavailable};
names=cellfun(@func2str,tests,'UniformOutput',false); elapsed=zeros(numel(tests),1);
for i=1:numel(tests)
    started=tic; tests{i}(); elapsed(i)=toc(started);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
fprintf('All %d PRS comparison checks passed.\n',height(report));
end

function [scene,v,spec]=geometryFixture()
scene=struct('bsPosM',[0,0;200,0;200,200;0,200],'pUserPosM',[80,110]);
v=[1;4;9;16];
spec=struct('referencePolicy','fixed','referenceBs',1,'rankTolerance',1e-10);
end

function testTdoaCovariance()
[scene,v,spec]=geometryFixture(); a=msnoma.geometry.tdoaAccuracy(scene,v,spec);
expected=diag(v(2:end))+v(1)*ones(3);
assert(a.valid && isequal(a.rangeDifferenceCovarianceM2{1},expected));
H=a.H{1}; expectedPosition=(H'*(expected\H))\eye(2);
assert(norm(a.covarianceM2{1}-expectedPosition,'fro')<1e-12);
end

function testReferenceInvariance()
[scene,v,spec]=geometryFixture(); a=msnoma.geometry.tdoaAccuracy(scene,v,spec);
for r=2:4
    spec.referenceBs=r; b=msnoma.geometry.tdoaAccuracy(scene,v,spec);
    assert(norm(a.covarianceM2{1}-b.covarianceM2{1},'fro')<1e-12);
end
spec.referencePolicy='minimum_variance'; v=[4;1;9;16];
a=msnoma.geometry.tdoaAccuracy(scene,v,spec); assert(a.referenceBs==2);
end

function testClockSchur()
[scene,v,spec]=geometryFixture(); a=msnoma.geometry.tdoaAccuracy(scene,v,spec);
offset=scene.pUserPosM-scene.bsPosM; G=offset./sqrt(sum(offset.^2,2));
J=[G,ones(4,1)]'*(diag(1./v)*[G,ones(4,1)]);
S=J(1:2,1:2)-J(1:2,3)/J(3,3)*J(3,1:2);
assert(norm(a.covarianceM2{1}-(S\eye(2)),'fro')<1e-12);
end

function testTdoaScaling()
[scene,v,spec]=geometryFixture(); a=msnoma.geometry.tdoaAccuracy(scene,v,spec);
b=msnoma.geometry.tdoaAccuracy(scene,7*v,spec);
assert(abs(b.positionVarianceM2/a.positionVarianceM2-7)<1e-12);
end

function testTdoaFailures()
[scene,v,spec]=geometryFixture(); mask=true(4,1); mask(1)=false;
a=msnoma.geometry.tdoaAccuracy(scene,v,spec,mask);
assert(~a.valid && strcmp(a.failureReason{1},'fixed_reference_unavailable'));
spec.referencePolicy='minimum_variance'; mask(2)=false;
a=msnoma.geometry.tdoaAccuracy(scene,v,spec,mask);
assert(~a.valid && isnan(a.positionStdM));
scene.bsPosM=[0,0;1,0;2,0;3,0]; scene.pUserPosM=[5,0];
a=msnoma.geometry.tdoaAccuracy(scene,v,spec);
assert(~a.valid && a.rank==0 && strcmp(a.failureReason{1},'tdoa_geometry_rank_deficient'));
[scene,v,spec]=geometryFixture(); v(2)=Inf;
a=msnoma.geometry.tdoaAccuracy(scene,v,spec); assert(a.valid && numel(a.bsIds{1})==3);
v(3)=NaN; a=msnoma.geometry.tdoaAccuracy(scene,v,spec); assert(~a.valid);
scene.pUserPosM=[0,0]; v=ones(4,1);
a=msnoma.geometry.tdoaAccuracy(scene,v,spec);
assert(~a.valid && strcmp(a.failureReason{1},'user_at_BS_undefined_LOS'));
end

function [cfg,num,scene,ch,allocation,budget,kernel,model]=fixture()
cfg=msnoma.config.base; cfg.paper.M=3; cfg.paper.requestedBandwidthHz=1.2e6;
num=msnoma.model.numerology(cfg); scene=msnoma.scenario.make(cfg,num);
ch=msnoma.scenario.channels(scene,num,cfg.assumed.propagation);
allocation=msnoma.prs.allocation(cfg.models.prs,num.BHz);
budget=msnoma.prs.powerBudget([.04;.03;.02;.01],cfg.models.prs.measurementPeriodS, ...
    cfg.models.prs.peakLimitW,allocation,cfg.models.prs);
kernel=msnoma.prs.interferenceKernels(ch,num,allocation,cfg.models.prs.comparison);
model=msnoma.comm.makeBerModel(cfg.models.ber);
end

function testPsdDirections()
[c,n,~,ch,a,~,kernel]=fixture(); k=2; m=3; ell=3;
expected=0; fp=(n.N+1)*n.deltaFcHz/2+a.activeK(ell)*a.deltaFHz;
for j=setdiff(1:n.K,k)
    for carrier=1:n.N
        x=(fp-carrier*n.deltaFcHz)*n.TcS;
        if x==0, s=1; else, s=sin(pi*x)/(pi*x); end
        expected=expected+n.PcW*n.TcS*ch.gCP(m,j,carrier)*s^2;
    end
end
assert(abs(kernel.communicationPsdWPerHz(k,m,ell)/expected-1)<1e-12);
carrier=5; j=3; expected=0;
for ell=1:a.numPilots
    fp=(n.N+1)*n.deltaFcHz/2+a.activeK(ell)*a.deltaFHz;
    x=(carrier*n.deltaFcHz-fp)*a.symbolPeriodS;
    if x==0, s=1; else, s=sin(pi*x)/(pi*x); end
    expected=expected+a.symbolPeriodS*a.weightsSquared(ell)/a.numCarriers*s^2;
end
expected=expected*ch.gPC(k,carrier,j,1);
assert(abs(kernel.prsPsdPerWAtCommunication(k,carrier,j)/expected-1)<1e-12);
ch.gCP(:,k,:)=ch.gCP(:,k,:)*1e6;
b=msnoma.prs.interferenceKernels(ch,n,a,c.models.prs.comparison);
assert(isequal(b.communicationPsdWPerHz(k,:,:),kernel.communicationPsdWPerHz(k,:,:)));
end

function testPsdChunking()
[c,n,~,ch,a,~,k]=fixture(); c.models.prs.comparison.chunkSize=1;
b=msnoma.prs.interferenceKernels(ch,n,a,c.models.prs.comparison);
assert(max(abs(b.communicationPsdWPerHz(:)./k.communicationPsdWPerHz(:)-1))<1e-12);
assert(max(abs(b.prsPsdPerWAtCommunication(:)./k.prsPsdPerWAtCommunication(:)-1))<1e-12);
end

function testZeroInterference()
[c,n,~,ch,a,b,k]=fixture(); k.communicationPsdWPerHz(:)=0;
v=msnoma.prs.delayCrbInterference(a,ch.gP,b,c.assumed.noisePsdWPerHz,k,n.c0Mps);
r=msnoma.prs.delayCrb(a,ch.gP.*b.instantaneousPowerW/c.assumed.noisePsdWPerHz,b.symbolCount,n.c0Mps);
assert(max(abs(v.delayVarianceS2(:)./r.delayVarianceS2(:)-1))<1e-12);
end

function testInterferenceMonotonicity()
[c,n,~,ch,a,b,k]=fixture(); N0=c.assumed.noisePsdWPerHz;
v=msnoma.prs.delayCrbInterference(a,ch.gP,b,N0,k,n.c0Mps);
k.communicationPsdWPerHz=3*k.communicationPsdWPerHz;
w=msnoma.prs.delayCrbInterference(a,ch.gP,b,N0,k,n.c0Mps);
assert(all(w.delayVarianceS2(:)>=v.delayVarianceS2(:)));
k.communicationPsdWPerHz(:)=N0;
u=msnoma.prs.delayCrbInterference(a,ch.gP,b,N0,k,n.c0Mps);
k.communicationPsdWPerHz(:)=0;
z=msnoma.prs.delayCrbInterference(a,ch.gP,b,N0,k,n.c0Mps);
assert(max(abs(u.delayVarianceS2(:)./z.delayVarianceS2(:)-2))<1e-12);
end

function testServiceGlobal()
[c,n,~,ch,~,b,k,model]=fixture(); ref=ones(n.K,n.N)*.005;
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    c.models.prs.comparison,'global_muting',ref);
L=1-b.jointDutyCycle;
assert(max(abs(s.communicationOpportunityFraction-L))<1e-15);
expected=msnoma.comm.ber(ch.gCServing,n.PcW,n.TcS,c.assumed.noisePsdWPerHz,zeros(n.K,n.N),model);
assert(max(abs(s.meanActiveBer(:)-expected(:)))<1e-15);
assert(all(s.activeSlotBerQosMet) && ~any(s.attemptedRateContractMet) && ~any(s.serviceContractMet));
assert(max(abs(s.expectedCorrectBitRatePerCarrierBps(:)-2/n.TcS*L*(1-expected(:))))<1e-9);
assert(max(abs(s.requiredPowerForEqualCommunicationEnergyW-n.N*n.PcW/L))<1e-12);
end

function testServiceOtherBs()
[c,n,~,ch,~,b,k,model]=fixture(); ref=zeros(n.K,n.N);
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    c.models.prs.comparison,'other_bs_communication',ref);
assert(max(abs(s.communicationOpportunityFraction-(1-b.dutyCycle)))<1e-15);
bs=2; carrier=5; expected=(1-b.jointDutyCycle)*s.berNoPrs(bs,carrier);
for j=setdiff(1:n.K,bs)
    noise=c.assumed.noisePsdWPerHz+k.prsPsdPerWAtCommunication(bs,carrier,j)*b.instantaneousPowerW(j);
    ber=.5*erfc(sqrt(ch.gCServing(bs,carrier)*n.PcW*n.TcS/2/noise));
    assert(abs(s.berDuringPrsSlot(bs,carrier,j)-ber)<1e-15);
    expected=expected+b.dutyCycle(j)*ber;
end
assert(isnan(s.berDuringPrsSlot(bs,carrier,bs)));
assert(abs(s.meanActiveBer(bs,carrier)-expected/(1-b.dutyCycle(bs)))<1e-15);
assert(abs(s.expectedCorrectBitRatePerCarrierBps(bs,carrier)-2/n.TcS*(1-b.dutyCycle(bs)-expected))<1e-9);
assert(all(s.worstActiveBer(:)>=s.meanActiveBer(:)-1e-15));
end

function testServiceStrictEnergy()
[c,n,~,ch,~,b,k,model]=fixture(); b.dutyCycle(:)=.05; b.jointDutyCycle=.2;
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    c.models.prs.comparison,'global_muting',zeros(n.K,n.N));
assert(~any(s.energyCompensationPeakFeasible) && ~any(s.communicationEnergyMatched));
assert(all(abs(s.requiredPowerForEqualCommunicationEnergyW-1.25)<1e-12));
spec=c.models.prs.comparison; spec.minimumAttemptedRateRatio=.8;
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    spec,'global_muting',zeros(n.K,n.N));
assert(all(s.serviceContractMet) && ~any(s.energyAndServiceMatched));
end

function testServiceZeroOpportunity()
[c,n,~,ch,~,b,k,model]=fixture(); b.dutyCycle(:)=.25; b.jointDutyCycle=1;
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    c.models.prs.comparison,'global_muting',zeros(n.K,n.N));
assert(all(isnan(s.meanActiveBer(:))) && all(isnan(s.worstActiveBer(:))));
assert(all(isinf(s.requiredPowerForEqualCommunicationEnergyW)) && ~any(s.serviceContractMet));
assert(all(s.expectedCorrectBitRatePerCarrierBps(:)==0));
end

function testServicePeakFailure()
[c,n,~,ch,~,b,k,model]=fixture(); spec=c.models.prs.comparison; spec.communicationPeakLimitW=.5;
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    spec,'other_bs_communication',zeros(n.K,n.N));
assert(~any(s.operatingPowerFeasible) && ~any(s.serviceContractMet));
assert(all(s.activeCommunicationPowerW==1)); % No silent cap to the requested peak.
end

function testCarrierGridNullAndOffset()
[c,n,~,~,~,~,~,~]=fixture(); n.N=4; n.PcW=1/n.N;
ch.gCP=ones(n.M,n.K,n.N); ch.gPC=ones(n.K,n.N,n.K,n.M);
spec=c.models.prs; spec.combSpacing=2;
a=msnoma.prs.allocation(spec,150e3);
k=msnoma.prs.interferenceKernels(ch,n,a,c.models.prs.comparison);
% C grid 30/60/90/120 kHz; PRS grid 15/45/105/135 kHz.
% All differences are nonzero integer PRS symbol frequencies: point nulls.
assert(max(k.prsPsdPerWAtCommunication(:))<1e-30);
spec.deltaFHz=17e3; a=msnoma.prs.allocation(spec,150e3);
q=msnoma.prs.interferenceKernels(ch,n,a,c.models.prs.comparison);
assert(max(q.prsPsdPerWAtCommunication(:))>1e-8);
assert(strcmp(k.communicationReceiverModel,'PSD_sample_at_C_carrier_not_matched_filter_band_integral'));
end

function testActiveBerFailure()
[c,n,~,ch,~,b,k,model]=fixture();
% Synthetic nonzero received PSD isolates service state weighting from grid nulls.
k.prsPsdPerWAtCommunication=repmat(ch.gCServing*n.PcW*n.TcS*100,1,1,n.K) ...
    ./reshape(b.instantaneousPowerW,1,1,n.K);
s=msnoma.prs.serviceAudit(ch,n,b,k,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    c.models.prs.comparison,'other_bs_communication',zeros(n.K,n.N));
assert(~any(s.activeSlotBerQosMet) && ~any(s.serviceContractMet));
assert(all(s.worstActiveBer(:)>.4));
assert(all(s.worstActiveBer(:)>s.meanActiveBer(:)));
assert(all(contains(s.failureReason,'active_slot_BER_QoS_failed')));
end

function testComparisonIntegration()
[c,~,~,~,~,~,~,~]=fixture(); c.models.rangingIntegral.enabled=false;
r=msnoma.experiment.runSingle(c,{'central'});
a=r.prsComparison;
assert(strcmp(r.algorithms(1).solution.status,'converged_verified'));
assert(strcmp(r.prsBenchmark.status,'computed_assumed_AWGN_broadcast_bound'));
assert(numel(a.policies)==2 && all(a.msTdoaCommon.valid));
assert(a.positioningBroadcastEnergyMatched && a.geometryConventionMatched);
assert(~a.observationMatched && ~a.fullyMatched && ~a.detectionCoverageEvaluated);
v1=a.policies(1).crb.rangeVarianceM2; v2=a.policies(2).crb.rangeVarianceM2;
assert(isequal(v1,r.prsBenchmark.crb.rangeVarianceM2) && all(v2(:)>=v1(:)));
assert(~any(a.policies(1).service.serviceContractMet));
directory=tempname; mkdir(directory);
cleanup=onCleanup(@() rmdir(directory,'s'));
msnoma.report.exportPrsComparison(r,directory);
tab=readtable(fullfile(directory,'prs_comparison_summary.csv'));
assert(height(tab)==2 && ~any(tab.communicationServiceMet) && ~any(tab.fullyMatched));
assert(numel(dir(fullfile(directory,'*.csv')))==9);
end

function testReferenceUnavailable()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false;
r=msnoma.experiment.runSingle(c,{'budget_equal'});
assert(strcmp(r.prsComparison.status,'reference_solution_unavailable'));
assert(~isfield(r.prsComparison,'policies'));
c.models.prs.measurementPeriodS=1e-6;
r=msnoma.experiment.runSingle(c,{'central'});
assert(strcmp(r.prsComparison.status,'prs_observation_schedule_infeasible'));
assert(isempty(r.prsComparison.policies));
end
