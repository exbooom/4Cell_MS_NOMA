function report = run_precision_tests()
%RUN_PRECISION_TESTS Independent DLL quadrature and PRS energy/CRB checks.
tests={@testAdaptiveDll,@testDllScaling,@testFrontendCoordinates, ...
    @testSmallSpacingLimit,@testPrsNormalization,@testPrsPhysicalEnergy, ...
    @testPrsSchedule,@testLteReference,@testBroadcastEnergy};
names=cellfun(@func2str,tests,'UniformOutput',false);
elapsed=zeros(numel(tests),1);
for i=1:numel(tests)
    start=tic; tests{i}(); elapsed(i)=toc(start);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
fprintf('All %d precision and PRS checks passed.\n',height(report));
end

function [c,n,ch,o]=context()
c=msnoma.config.base;
c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
n=msnoma.model.numerology(c);
s=msnoma.scenario.make(c,n);
ch=msnoma.scenario.channels(s,n,c.assumed.propagation);
o=c.models.rangingIntegral;
o.userIds=1:3; o.relativeTolerance=1e-7; o.maxRefinements=4;
end

function testAdaptiveDll()
[c,n,ch,o]=context();
k=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
m=1; D=c.paper.dll.earlyLateSpacingChips;
b=k.localBoundsHz(1,:)*n.TpS;
gain=reshape(sum(ch.gCP(m,:,:),2),1,n.N);
q=@(x) msnoma.model.sincpi(x).^2;
u=@(x) sin(pi*D*x);
psd=@(x) reshape(msnoma.model.sincpi(n.G*x(:)+m*n.G-(1:n.N)).^2*gain',size(x));
v=[integral(@(x) x.*q(x).*u(x),b(1),b(2),'RelTol',1e-11,'AbsTol',1e-15)/n.TpS, ...
    integral(@(x) q(x).*u(x).^2,b(1),b(2),'RelTol',1e-11,'AbsTol',1e-15), ...
    n.TcS*integral(@(x) psd(x).*q(x).*u(x).^2,b(1),b(2),'RelTol',1e-10,'AbsTol',1e-27), ...
    n.TpS*integral(@(x) q(x).^2.*u(x).^2,b(1),b(2),'RelTol',1e-11,'AbsTol',1e-20)];
actual=[k.A0PerS(1),k.noiseIntegral(1),k.communicationIntegralS(1),k.positioningIntegralS(1)];
assert(max(abs(actual./v-1))<3e-7 && all(k.converged));
end

function testDllScaling()
[c,n,ch,o]=context(); o.userIds=2;
k=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
P=0.002*ones(n.K,n.M); N0=c.assumed.noisePsdWPerHz;
a=msnoma.ranging.varianceIntegral(P,ch,k,N0);
k2=k; k2.PcW=2*k.PcW;
b=msnoma.ranging.varianceIntegral(2*P,ch,k2,2*N0);
assert(max(abs(a.delayVarianceS2(:)./b.delayVarianceS2(:)-1))<1e-12);
P(2:end,:)=0;
a=msnoma.ranging.varianceIntegral(P,ch,k,N0);
b=msnoma.ranging.varianceIntegral(2*P,ch,k,N0);
assert(abs(a.delayVarianceS2(1)/b.delayVarianceS2(1)-2)<1e-12);
assert(all(isinf(a.delayVarianceS2(2:end))));
end

function testFrontendCoordinates()
[c,n,ch,o]=context(); ch.gCP(:)=1;
globalK=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
o.frontendMode='local_symmetric';
localK=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
assert(max(abs(globalK.localBoundsHz(2,:)-localK.localBoundsHz(2,:)))<1e-6);
assert(abs(globalK.A0PerS(2)/localK.A0PerS(2)-1)<1e-12);
assert(abs(globalK.communicationIntegralS(1)/globalK.communicationIntegralS(3)-1)<1e-10);
assert(globalK.localBoundsHz(1,1)~=localK.localBoundsHz(1,1));
end

function testSmallSpacingLimit()
[c,n,ch,o]=context(); o.userIds=2;
c.paper.dll.earlyLateSpacingChips=1e-4;
a=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
o.discriminator='small_spacing';
b=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
P=ones(n.K,n.M)*0.003;
va=msnoma.ranging.varianceIntegral(P,ch,a,c.assumed.noisePsdWPerHz);
vb=msnoma.ranging.varianceIntegral(P,ch,b,c.assumed.noisePsdWPerHz);
assert(max(abs(va.delayVarianceS2./vb.delayVarianceS2-1))<1e-5);
end

function testPrsNormalization()
c=msnoma.config.base;
a=msnoma.prs.allocation(c.models.prs,49.77e6);
assert(a.numCarriers==3318 && ~a.isStandardLtePattern);
assert(abs(sum(a.weightsSquared)/a.numCarriers-1)<1e-12);
assert(~any(a.activeK==0));
assert(all(mod(a.activeK,c.models.prs.combSpacing)==c.models.prs.combOffset));
x=msnoma.prs.delayCrb(a,1e8,1,299792458);
y=msnoma.prs.delayCrb(a,2e8,3,299792458);
assert(abs(x.delayVarianceS2/y.delayVarianceS2-6)<1e-12);
end

function testPrsPhysicalEnergy()
c=msnoma.config.base; a=msnoma.prs.allocation(c.models.prs,20.16e6);
cn0=1e8; observations=7;
v=msnoma.prs.delayCrb(a,cn0,observations,299792458);
frequencyHz=a.activeK*a.deltaFHz;
meanSquareFrequency=sum(a.weightsSquared.*frequencyHz.^2)/sum(a.weightsSquared);
energyOverN0=cn0*a.symbolPeriodS*observations;
independent=1/(8*pi^2*energyOverN0*meanSquareFrequency);
assert(abs(v.delayVarianceS2/independent-1)<1e-12);
end

function testPrsSchedule()
c=msnoma.config.base; a=msnoma.prs.allocation(c.models.prs,20.16e6);
b=msnoma.prs.powerBudget([0.04;0.08],0.02,1.1,a,c.models.prs);
assert(all(b.feasible) && all(b.instantaneousPowerW<=1.1));
assert(max(abs(b.instantaneousPowerW.*b.observationTimeS-b.energyPerMeasurementJ))<1e-15);
assert(max(abs(b.energyPerMeasurementJ-[0.04;0.08]*0.02))<1e-15);
spec=c.models.prs; spec.symbolPolicy='single_reference_symbol';
b=msnoma.prs.powerBudget(0.08,0.02,1.1,a,spec);
assert(~b.feasible && b.instantaneousPowerW>1.1);
b=msnoma.prs.powerBudget(0.3*ones(4,1),0.02,1.1,a,c.models.prs);
assert(all(b.dutyCycle<1) && b.jointDutyCycle>1 && ~any(b.feasible));
end

function testLteReference()
c=msnoma.config.base; p=c.models.prs;
p.mode='lte_single_symbol_reference';
a=msnoma.prs.allocation(p,1.4e6);
assert(a.isStandardLtePattern && a.numCarriers==72 && a.numPilots==12);
assert(a.occupiedBandwidthHz==1.08e6 && abs(sum(a.weightsSquared)-72)<1e-12);
end

function testBroadcastEnergy()
[c,n,ch]=context(); s=msnoma.scenario.make(c,n);
p=c.models.prs;
b=msnoma.prs.positionBenchmark(s,ch,n,0.08*ones(n.K,1), ...
    c.assumed.noisePsdWPerHz,p,c.models.geometry);
assert(strcmp(b.status,'computed_assumed_AWGN_broadcast_bound'));
assert(numel(b.budget.energyPerMeasurementJ)==n.K);
assert(all(b.position.valid) && ~b.detectionCoverageEvaluated);
assert(abs(sum(b.budget.energyPerMeasurementJ)-n.K*0.08*p.measurementPeriodS)<1e-15);
end
