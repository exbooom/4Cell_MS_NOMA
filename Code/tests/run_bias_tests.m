function report = run_bias_tests()
%RUN_BIAS_TESTS Controlled DLL approximations and independent filter integrals.
tests={@testFinitePsd,@testSpectralLineFormula,@testFiniteNAsymptotics, ...
    @testFamilyAndCache,@testPairedContrasts,@testChainClosure,@testComponentIsolation, ...
    @testSmallSpacingLimit,@testFrontendCenter,@testZeroCommunication, ...
    @testOverlapFrequencyIntegral,@testOverlapSpecialValues,@testOverlapScaleSymmetry, ...
    @testOverlapNarrowFilter,@testReceiverAudit,@testBiasIntegration,@testBiasFailures};
names=cellfun(@func2str,tests,'UniformOutput',false); elapsed=zeros(numel(tests),1);
for i=1:numel(tests)
    started=tic; tests{i}(); elapsed(i)=toc(started);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
fprintf('All %d DLL bias and receiver-filter checks passed.\n',height(report));
end

function [c,n,s,ch,o,P,f,g]=fixture()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.dllBias.enabled=false; n=msnoma.model.numerology(c);
s=msnoma.scenario.make(c,n); ch=msnoma.scenario.channels(s,n,c.assumed.propagation);
ch.gCP=ch.gCP.*reshape(1+.2*sin(1:n.N),1,1,n.N);
o=c.models.rangingIntegral; o.userIds=1:3; o.relativeTolerance=1e-7; o.maxRefinements=4;
P=reshape(linspace(.001,.003,n.K*n.M),n.K,n.M);
f=msnoma.ranging.factorTwoTerm(ch,n,c.assumed.noisePsdWPerHz,c.paper.dll);
g=msnoma.geometry.weights(s,true(n.K,n.M),c.models.geometry);
end

function testFinitePsd()
[~,n,~,ch,~,~,~,~]=fixture(); freq=[-14000,1234,17000]; m=2;
actual=msnoma.ranging.communicationPsd(freq,m,ch,n,1); expected=zeros(size(freq));
for i=1:numel(freq)
    for j=1:n.K
        for carrier=1:n.N
            x=(freq(i)+m*n.deltaFpHz-carrier*n.deltaFcHz)*n.TcS;
            if x==0, v=1; else, v=sin(pi*x)/(pi*x); end
            expected(i)=expected(i)+n.TcS*ch.gCP(m,j,carrier)*v^2;
        end
    end
end
assert(max(abs(actual./expected-1))<1e-12);
end

function testSpectralLineFormula()
[c,n,~,ch,o,~,~,~]=fixture(); o.discriminator='small_spacing'; o.communicationModel='spectral_lines';
k=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
gain=reshape(sum(ch.gCP,2),n.M,n.N);
expected=n.TpS*c.paper.dll.earlyLateSpacingChips^2*sum(gain.*sin(pi*(1:n.N)/n.G).^2,2)';
assert(max(abs(k.communicationIntegralS./expected-1))<1e-12);
end

function testFiniteNAsymptotics()
[c,n,~,ch,~,P,f,~]=fixture();
a=msnoma.ranging.asymptoticVariance(P,ch,n,c.assumed.noisePsdWPerHz,c.paper.dll);
b=msnoma.ranging.varianceThreeTerm(P,ch,n,f);
assert(max(abs(a.noiseDelayVarianceS2./b.noiseDelayVarianceS2-1),[],'all')<1e-12);
assert(max(abs(a.positioningDelayVarianceS2./b.positioningDelayVarianceS2-1),[],'all')<1e-12);
assert(max(abs(b.communicationDelayVarianceS2./a.communicationDelayVarianceS2-(n.N+1)/n.N),[],'all')<1e-12);
end

function testFamilyAndCache()
[c,n,~,ch,o,~,~,~]=fixture(); k=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
b=msnoma.ranging.biasKernels(ch,n,c.paper.dll,o,k);
assert(numel(b.variants)==8 && b.allConverged && b.variants(1).reused);
assert(nnz([b.variants.reused])==1 && size(unique(vertcat(b.variants.flags),'rows'),1)==8);
assert(isequaln(k,b.variants(1).kernel));
end

function [r,b,c,n,ch,P,f,g]=result()
[c,n,~,ch,o,P,f,g]=fixture(); b=msnoma.ranging.biasKernels(ch,n,c.paper.dll,o);
r=msnoma.ranging.biasDecomposition(P,ch,n,c.paper.dll,c.assumed.noisePsdWPerHz,b,f,g);
end

function testPairedContrasts()
[r,b]=result(); assert(numel(r.pairedContrasts)==12);
for q=1:12
    p=r.pairedContrasts(q); difference=b.variants(p.toIndex).flags-b.variants(p.fromIndex).flags;
    assert(nnz(difference)==1 && all(difference>=0));
    ratio=r.delayVarianceS2(:,:,p.toIndex)./r.delayVarianceS2(:,:,p.fromIndex);
    assert(max(abs(p.relativeToFrom-(ratio-1)),[],'all')<1e-14);
    assert(max(abs(exp(p.logVarianceRatio)./ratio-1),[],'all')<1e-14);
end
end

function testChainClosure()
[r,~,~,~,~,P,f]=result(); expected=msnoma.ranging.varianceTwoTerm(P,f);
assert(r.maximumRelativeClosureResidual<1e-12);
assert(max(abs(sum(r.chainRelativeToBaseline,3)-r.totalRelativeBias),[],'all')<1e-12);
assert(isequal(r.chainVarianceS2(:,:,end),expected.delayVarianceS2));
assert(all(r.chainDifferenceS2(:,:,end)<0,'all'));
end

function testComponentIsolation()
[r,b]=result(); flags=vertcat(b.variants.flags);
target=find(all(flags==[0,1,0],2));
assert(max(abs(r.noiseDelayVarianceS2(:,:,target)./r.noiseDelayVarianceS2(:,:,1)-1),[],'all')<5e-7);
assert(max(abs(r.positioningDelayVarianceS2(:,:,target)./r.positioningDelayVarianceS2(:,:,1)-1),[],'all')<5e-7);
end

function testSmallSpacingLimit()
[c,n,~,ch,o,P,~,~]=fixture(); c.paper.dll.earlyLateSpacingChips=1e-5;
a=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o); o.discriminator='small_spacing';
b=msnoma.ranging.integralKernels(ch,n,c.paper.dll,o);
x=msnoma.ranging.varianceIntegral(P,ch,a,c.assumed.noisePsdWPerHz);
y=msnoma.ranging.varianceIntegral(P,ch,b,c.assumed.noisePsdWPerHz);
assert(max(abs(x.delayVarianceS2./y.delayVarianceS2-1),[],'all')<1e-7);
end

function testFrontendCenter()
[r]=result(); assert(max(abs(r.mainEffectRelativeToBaseline(:,2,3)))<1e-12);
end

function testZeroCommunication()
[c,n,~,ch,o,P,f,g]=fixture(); ch.gCP(:)=0;
b=msnoma.ranging.biasKernels(ch,n,c.paper.dll,o);
r=msnoma.ranging.biasDecomposition(P,ch,n,c.paper.dll,c.assumed.noisePsdWPerHz,b,f,g);
assert(all(r.communicationDelayVarianceS2==0,'all'));
assert(max(abs(r.mainEffectRelativeToBaseline(:,:,2)),[],'all')<1e-12);
end

function testOverlapFrequencyIntegral()
T=1/30000; offsets=[0,.5,3]/T;
actual=msnoma.receiver.rectangularPsdOverlap(offsets,2*T,T);
for i=1:numel(offsets)
    d=offsets(i);
    fun=@(f) 2*T*msnoma.model.sincpi((f-d)*2*T).^2.*T.*msnoma.model.sincpi(f*T).^2;
    expected=integral(fun,-150/T,150/T,'Waypoints',(-150:150)/T,'RelTol',1e-9,'AbsTol',1e-13);
    assert(abs(actual(i)/expected-1)<2e-6);
end
end

function testOverlapSpecialValues()
T=1/30000;
assert(abs(msnoma.receiver.rectangularPsdOverlap(0,T,T)/(2*T/3)-1)<1e-14);
assert(abs(msnoma.receiver.rectangularPsdOverlap(0,2*T,T)/(T*5/6)-1)<1e-14);
assert(abs(msnoma.receiver.rectangularPsdOverlap(1/T,T,T)/(T/pi^2)-1)<1e-14);
end

function testOverlapScaleSymmetry()
f=[-60000,-.001,0,.001,60000]; T=1/30000;
a=msnoma.receiver.rectangularPsdOverlap(f,T,2*T);
b=msnoma.receiver.rectangularPsdOverlap(-f,2*T,T);
c=msnoma.receiver.rectangularPsdOverlap(f/7,7*T,14*T);
assert(all(a>0) && max(abs(a./b-1))<1e-12 && max(abs(c./a-7))<1e-12);
end

function testOverlapNarrowFilter()
T=1/30000; f=.35/T;
a=msnoma.receiver.rectangularPsdOverlap(f,T,1e8*T);
expected=T*(sin(pi*f*T)/(pi*f*T))^2;
assert(abs(a/expected-1)<1e-7);
end

function testReceiverAudit()
[c,n,~,ch,~,~,~,~]=fixture(); a=msnoma.prs.allocation(c.models.prs,n.BHz);
b=msnoma.prs.powerBudget(.02*ones(n.K,1),.02,1.1,a,c.models.prs);
k=msnoma.prs.interferenceKernels(ch,n,a,c.models.prs.comparison);
model=msnoma.comm.makeBerModel(c.models.ber);
r=msnoma.prs.receiverAudit(ch,n,a,b,c.assumed.noisePsdWPerHz,model,c.paper.berThreshold, ...
    c.models.prs.comparison,k,zeros(n.K,n.N));
assert(all(r.filteredUnitPsdWPerHzPerW>0,'all'));
assert(~r.sampledReceiverImplemented && ~r.finiteObservationMatched && ~r.prsRangingRecomputed);
assert(isequal(r.pointService.communicationOpportunityFraction,r.filteredService.communicationOpportunityFraction));
end

function testBiasIntegration()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false;
r=msnoma.experiment.runSingle(c,{'central'});
assert(strcmp(r.dllBias.status,'computed_fixed_reference_power_variance_bias_decomposition'));
assert(r.dllBias.result.allConverged && ~r.dllBias.optimizerChanged);
assert(numel(r.dllBias.bundle.variants)==8 && all(r.algorithms(1).evaluation.valid));
directory=tempname; mkdir(directory); cleanup=onCleanup(@() rmdir(directory,'s'));
msnoma.report.exportDllBias(r,directory); msnoma.report.exportReceiverAudit(r,directory);
assert(numel(dir(fullfile(directory,'*.csv')))==8);
assert(height(readtable(fullfile(directory,'dll_bias_variants.csv')))==4*3*8);
assert(height(readtable(fullfile(directory,'dll_bias_paired_contrasts.csv')))==4*3*12);
assert(height(readtable(fullfile(directory,'dll_bias_chain.csv')))==4*3*6);
assert(height(readtable(fullfile(directory,'dll_bias_effect_summary.csv')))==9);
tab=readtable(fullfile(directory,'prs_receiver_filter_links.csv'));
assert(height(tab)==r.num.K^2*r.num.N);
assert(all(isnan(tab.filteredSlotBer(~tab.communicationActive))));
assert(all(isfinite(tab.filteredSlotBer(logical(tab.communicationActive)))));
end

function testBiasFailures()
c=msnoma.config.base; c.paper.M=3; c.paper.requestedBandwidthHz=1.2e6;
c.models.rangingIntegral.enabled=false;
r=msnoma.experiment.runSingle(c,{'budget_equal'});
assert(strcmp(r.dllBias.status,'reference_solution_unavailable'));
failed=false;
try
    msnoma.receiver.rectangularPsdOverlap(0,-1,1);
catch
    failed=true;
end
assert(failed);
c.models.dllBias.enabled=false; c.models.prs.comparison.receiverAuditEnabled=false;
r=msnoma.experiment.runSingle(c,{'central'});
assert(strcmp(r.dllBias.status,'disabled') && strcmp(r.prsComparison.receiverAudit.status,'disabled'));
end
