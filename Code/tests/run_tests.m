function report = run_tests()
%RUN_TESTS Mathematical invariants and independent solver checks.
codeRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(codeRoot);
tests={@testNumerology,@testRandomIsolation,@testUserAssociation,@testFriis,@testChannels, ...
    @testInterference,@testBerInverse,@testNegativeQos,@testRanging, ...
    @testGeometry,@testRankDeficiency,@testHearability,@testDerivative, ...
    @testAnalyticOptimum,@testImpossibleHearability,@testBaselines,@testCdf};
names=cellfun(@func2str,tests,'UniformOutput',false);
elapsed=zeros(numel(tests),1);
for i=1:numel(tests)
    start=tic;
    tests{i}();
    elapsed(i)=toc(start);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
fprintf('All %d mathematical checks passed.\n',height(report));
end

function [cfg,num,scene,ch,kernels] = smallContext()
cfg=msnoma.config.base;
cfg.paper.M=3;
cfg.paper.requestedBandwidthHz=1.2e6;
num=msnoma.model.numerology(cfg);
scene=msnoma.scenario.make(cfg,num);
ch=msnoma.scenario.channels(scene,num,cfg.assumed.propagation);
kernels=msnoma.model.spectralKernels(num);
end

function testNumerology()
c=msnoma.config.base;
n=msnoma.model.numerology(c);
assert(n.G==79 && n.N==1658 && n.effectiveBandwidthHz==49.77e6);
c.paper.requestedBandwidthHz=20e6;
n=msnoma.model.numerology(c);
assert(n.G==32 && n.N==671 && n.effectiveBandwidthHz==20.16e6);
assert(abs(n.TcS/n.TpS-n.G)<1e-12);
end

function testRandomIsolation()
[cfg,num,s]=smallContext();
before=rng;
t=msnoma.scenario.make(cfg,num);
after=rng;
assert(isequal(before,after) && isequal(s,t));
end

function testUserAssociation()
[~,num,s]=smallContext();
for k=1:num.K
    rx=reshape(s.cUserPosM(k,:,:),num.N,2);
    d=sum((reshape(rx,num.N,1,2)-reshape(s.bsPosM,1,num.K,2)).^2,3);
    [~,nearest]=min(d,[],2);
    assert(all(nearest==k));
end
end

function testFriis()
c=msnoma.config.base;
g=msnoma.scenario.freeSpaceGain([0,0],[10,0;20,0],c.assumed.propagation);
assert(abs(g(1)/g(2)-4)<1e-12);
lambda=c.assumed.propagation.speedMps/c.paper.carrierHz;
assert(abs(g(1)/(lambda/(4*pi*10))^2-1)<1e-12);
end

function testInterference()
[~,n,~,ch,w]=smallContext();
ch.gPC=reshape(1:numel(ch.gPC),size(ch.gPC))*1e-10;
A=msnoma.comm.interferenceMatrix(ch,n,w);
P=reshape(linspace(1e-4,0.01,n.K*n.M),n.K,n.M);
expected=zeros(n.K,n.N);
for k=1:n.K
    for carrier=1:n.N
        for j=1:n.K
            for m=1:n.M
                expected(k,carrier)=expected(k,carrier)+ ...
                    ch.gPC(k,carrier,j,m)*P(j,m)*n.TpS ...
                    *msnoma.model.sincpi(m-carrier/n.G)^2;
            end
        end
    end
end
actual=msnoma.comm.interference(P,A,n);
assert(norm(actual(:)-expected(:))/norm(expected(:))<1e-13);
end

function testChannels()
[c,n,s,ch]=smallContext();
k=2; carrier=7; j=4; m=3;
rx=reshape(s.cUserPosM(k,carrier,:),1,2);
distance=norm(rx-s.bsPosM(j,:));
expected=(c.assumed.propagation.speedMps/(4*pi*c.paper.carrierHz*distance))^2;
assert(abs(ch.gPC(k,carrier,j,m)/expected-1)<1e-12);
distance=norm(s.pUserPosM(m,:)-s.bsPosM(j,:));
expected=(c.assumed.propagation.speedMps/(4*pi*c.paper.carrierHz*distance))^2;
assert(abs(ch.gCP(m,j,carrier)/expected-1)<1e-12);
assert(ch.gP(j,m)==ch.gCP(m,j,n.N));
replay=msnoma.scenario.fromSnapshot(s,n);
assert(isequal(replay,s));
end

function testBerInverse()
c=msnoma.config.base;
for id={'qpsk_bit_energy','paper_literal'}
    spec=c.models.ber; spec.id=id{1}; model=msnoma.comm.makeBerModel(spec);
    g=[1,0.7;0.3,0.1];
    q=msnoma.comm.qosLimit(g,0.01,1/30000,1e-12,0.01,model);
    b=q.limitWPerHz;
    center=msnoma.comm.ber(g,0.01,1/30000,1e-12,b,model);
    low=msnoma.comm.ber(g,0.01,1/30000,1e-12,0.8*b,model);
    high=msnoma.comm.ber(g,0.01,1/30000,1e-12,1.2*b,model);
    assert(max(abs(center(:)-0.01))<1e-12 && all(low(:)<0.01) && all(high(:)>0.01));
end
end

function testNegativeQos()
c=msnoma.config.base; model=msnoma.comm.makeBerModel(c.models.ber);
q=msnoma.comm.qosLimit(1e-15,1e-3,1e-5,1,0.01,model);
assert(~q.cleanFeasible && q.limitWPerHz<0 && q.cleanBer>0.01);
p=toyProblem(); p.bqos=-1;
f=msnoma.optim.checkFeasibility(p,c.numerics);
assert(~f.feasible && strcmp(f.status,'noise_alone_violates_qos'));
end

function testRanging()
[c,n,~,ch]=smallContext();
f=msnoma.ranging.factorTwoTerm(ch,n,c.assumed.noisePsdWPerHz,c.paper.dll);
P=ones(n.K,n.M)*0.003;
a=msnoma.ranging.varianceTwoTerm(P,f);
b=msnoma.ranging.varianceTwoTerm(2*P,f);
assert(max(abs(a.delayVarianceS2(:)./b.delayVarianceS2(:)-2))<1e-12);
assert(max(abs(a.rangeVarianceM2(:)./a.delayVarianceS2(:)/n.c0Mps^2-1))<1e-12);
t=msnoma.ranging.varianceThreeTerm(P,ch,n,f);
rx=ch.gP.*P;
expected=f.aPerS*n.TpS/(2*n.BfeHz^2)*(sum(rx(:,1))-rx(1,1))/rx(1,1);
assert(abs(t.positioningDelayVarianceS2(1,1)/expected-1)<1e-12);
assert(all(t.delayVarianceS2(:)>=a.delayVarianceS2(:)));
end

function testGeometry()
c=msnoma.config.base;
s=struct('bsPosM',c.paper.bsPosM,'pUserPosM',[100,100]);
g=msnoma.geometry.weights(s,true(4,1),c.models.geometry);
assert(g.valid && max(abs(g.lambda-0.5))<1e-12);
p=msnoma.geometry.positionAccuracy(4*ones(4,1),g);
assert(abs(p.positionStdM-2)<1e-12 && norm(p.covarianceM2{1}-2*eye(2))<1e-12);
end

function testRankDeficiency()
c=msnoma.config.base;
s=struct('bsPosM',[0,0;1,0;2,0],'pUserPosM',[3,0]);
g=msnoma.geometry.weights(s,true(3,1),c.models.geometry);
p=msnoma.geometry.positionAccuracy(ones(3,1),g);
assert(~g.valid && ~p.valid && isnan(p.positionStdM));
end

function testHearability()
c=msnoma.config.base;
ch=struct('gP',ones(4,2));
P=[1,0;0.2,0;0.11,0;0.05,0];
h=msnoma.geometry.hearability(P,ch,c.models.hearability);
assert(isequal(h.mask(:,1),logical([1;1;1;0])) && ~any(h.mask(:,2)));
assert(h.count(1)==3 && ~h.enoughBs(1));
end

function testDerivative()
p=toyProblem(); z=[1.1;2.2;3.3];
[~,gradient,H]=msnoma.optim.objective(z,p);
h=1e-5; fd=zeros(3,1); fH=zeros(3);
for i=1:3
    e=zeros(3,1); e(i)=h;
    fp=msnoma.optim.objective(z+e,p); fm=msnoma.optim.objective(z-e,p);
    fd(i)=(fp-fm)/(2*h);
    [~,gp]=msnoma.optim.objective(z+e,p); [~,gm]=msnoma.optim.objective(z-e,p);
    fH(:,i)=(gp-gm)/(2*h);
end
assert(norm(fd-gradient,inf)<1e-8 && norm(fH-H,inf)<1e-8);
end

function testAnalyticOptimum()
c=msnoma.config.base; p=toyProblem();
s=msnoma.optim.solveCentral(p,c.numerics);
expected=p.budgetW*sqrt(p.c)/sum(sqrt(p.c));
assert(strcmp(s.status,'converged_verified'));
assert(norm(s.powerW(:)-expected,inf)<1e-6);
assert(abs(s.objectiveM2/(sum(sqrt(p.c))^2/p.budgetW)-1)<1e-7);
end

function testImpossibleHearability()
c=msnoma.config.base; p=toyProblem(); p.q=1.1;
f=msnoma.optim.checkFeasibility(p,c.numerics);
assert(~f.feasible && strcmp(f.status,'all_pair_hearability_impossible'));
end

function testBaselines()
p=toyProblem(); p.Aqos=[1,0,0]; p.bqos=0.1;
a=msnoma.optim.equalPower(p,'budget_equal');
b=msnoma.optim.equalPower(p,'qos_scaled_equal');
assert(max(a.powerW)-min(a.powerW)<1e-12 && a.powerW(1)>0.1);
assert(max(abs(b.powerW-0.1))<1e-12);
end

function testCdf()
[x,F]=msnoma.report.unconditionalCdf([0.5;NaN;1;NaN],[true;false;true;false]);
assert(isequal(x,[0.5;1]) && isequal(F,[0.25;0.5]));
end

function p = toyProblem()
p.K=1; p.M=3; p.N=0;
p.c=[1;4;9]; p.Aqos=zeros(0,3); p.bqos=zeros(0,1);
p.Bbudget=ones(1,3); p.budgetW=0.7; p.Chear=zeros(0,3); p.q=0;
p.Aineq=p.Bbudget; p.bineq=0.7;
p.powerScaleW=0.1; p.objectiveScaleM2=sum(p.c/p.powerScaleW);
p.rowScale=0.7; p.scaledA=p.Aineq*p.powerScaleW/p.rowScale;
p.scaledB=1; p.scaledC=p.c/(p.powerScaleW*p.objectiveScaleM2);
p.lbW=ones(3,1)*1e-10; p.ubW=ones(3,1)*0.7;
p.scaledLb=p.lbW/p.powerScaleW; p.scaledUb=p.ubW/p.powerScaleW;
end
