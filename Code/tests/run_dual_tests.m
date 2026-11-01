function report = run_dual_tests()
%RUN_DUAL_TESTS Independent dual algebra, coupling and hierarchical checks.
tests={@testScalarBox,@testDualDerivatives,@testWeakDuality, ...
    @testBudgetOptimum,@testCoupledOptimum,@testActiveSetCompleteness, ...
    @testDualInfeasible,@testCp4aMode,@testCarrierAudit,@testPriceCoupling};
elapsed=zeros(numel(tests),1);
names=cellfun(@func2str,tests,'UniformOutput',false);
for i=1:numel(tests)
    started=tic; tests{i}(); elapsed(i)=toc(started);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
fprintf('All %d dual and hierarchical checks passed.\n',height(report));
end

function p=makeProblem(K,c,A,b,B,budget,C)
p.K=K; p.M=numel(c)/K; p.N=size(A,1)/K; p.G=1; p.q=0.1;
p.c=c(:); p.Aqos=A; p.bqos=b(:); p.Bbudget=B;
p.budgetW=budget(:); p.Chear=C;
p.Aineq=[A;B;C]; p.bineq=[b(:);budget(:);zeros(size(C,1),1)];
p.powerScaleW=0.1; p.objectiveScaleM2=sum(c/p.powerScaleW);
p.lbW=ones(numel(c),1)*1e-10; p.ubW=repmat(budget(:),p.M,1);
p.rowScale=max(max(abs(p.Aineq*p.powerScaleW),[],2),abs(p.bineq));
p.rowScale(p.rowScale==0)=1;
p.scaledA=p.Aineq*p.powerScaleW./p.rowScale;
p.scaledB=p.bineq./p.rowScale;
p.scaledC=p.c/(p.powerScaleW*p.objectiveScaleM2);
p.scaledLb=p.lbW/p.powerScaleW; p.scaledUb=p.ubW/p.powerScaleW;
p.rowMap.qos=(1:size(A,1))';
p.rowMap.budget=size(A,1)+(1:size(B,1))';
p.rowMap.hearability=size(A,1)+size(B,1)+(1:size(C,1))';
end

function p=budgetProblem()
p=makeProblem(1,[1;4;9],zeros(0,3),[],ones(1,3),0.7,zeros(0,3));
end

function p=coupledProblem()
A=[1,.8,.1,.2;.2,1,.3,.7]; b=[.5;.7];
B=[1,0,1,0;0,1,0,1];
C=[-1,.2,0,0;.1,-1,0,0;0,0,-1,.3;0,0,.2,-1];
p=makeProblem(2,[1;4;9;16],A,b,B,[.7;.9],C);
end

function testScalarBox()
p=budgetProblem(); p.scaledC=[1;0;4];
p.scaledA=[-1,2,0]; p.scaledB=0;
[value,~,z]=msnoma.optim.dualValue(1,p);
assert(z(1)==p.scaledUb(1) && z(2)==p.scaledLb(2) && z(3)==p.scaledUb(3));
independent=sum(p.scaledC./z+(p.scaledA'*1).*z);
assert(abs(value+independent)<1e-12);
end

function testDualDerivatives()
p=coupledProblem(); l=[.2;.3;2;2;.01;.02;.01;.02];
[~,g,~,H]=msnoma.optim.dualValue(l,p);
h=1e-6; gn=zeros(size(g)); Hn=zeros(size(H));
for j=1:numel(l)
    a=l; b=l; a(j)=a(j)+h; b(j)=b(j)-h;
    [fa,ga]=msnoma.optim.dualValue(a,p); [fb,gb]=msnoma.optim.dualValue(b,p);
    gn(j)=(fa-fb)/(2*h); Hn(:,j)=(ga-gb)/(2*h);
end
assert(norm(gn-g,inf)<1e-8 && norm(Hn-H,inf)<1e-7);
assert(min(eig(H))>-1e-10);
end

function testWeakDuality()
p=coupledProblem(); c=msnoma.config.base;
feas=msnoma.optim.checkFeasibility(p,c.numerics);
assert(feas.feasible);
primal=sum(p.scaledC./feas.scaledPoint);
for scale=[0,.1,1,10]
    v=msnoma.optim.dualValue(scale*ones(size(p.scaledB)),p);
    assert(-v<=primal+1e-12);
end
end

function testBudgetOptimum()
p=budgetProblem(); c=msnoma.config.base;
expected=p.budgetW*sqrt(p.c)/sum(sqrt(p.c));
a=msnoma.optim.solveDual(p,c.numerics);
b=msnoma.optim.solveCp4a(p,c.numerics);
assert(strcmp(a.status,'converged_verified') && strcmp(b.status,'converged_verified'));
assert(norm(a.powerW(:)-expected,inf)<1e-6 && norm(b.powerW(:)-expected,inf)<1e-6);
assert(b.innerEvaluations>1 && b.multipliers.ineqlin(p.rowMap.budget)>0);
assert(~b.literalPaperAlgorithm && ~a.usedCentralFallback && ~b.usedCentralFallback);
end

function testCoupledOptimum()
p=coupledProblem(); c=msnoma.config.base;
a=msnoma.optim.solveCentral(p,c.numerics);
b=msnoma.optim.solveDual(p,c.numerics);
d=msnoma.optim.solveCp4a(p,c.numerics);
assert(strcmp(a.status,'converged_verified') && strcmp(b.status,'converged_verified'));
assert(strcmp(d.status,'converged_verified'));
assert(abs(b.objectiveM2/a.objectiveM2-1)<1e-6 && abs(d.objectiveM2/a.objectiveM2-1)<1e-6);
assert(b.dualBoundM2<=a.objectiveM2*(1+1e-8));
end

function testActiveSetCompleteness()
p=coupledProblem(); c=msnoma.config.base; c.numerics.dual.rowsPerPass=1;
s=msnoma.optim.solveDual(p,c.numerics);
assert(strcmp(s.status,'converged_verified') && s.dual.passes>1);
assert(numel(s.multipliers.ineqlin)==size(p.scaledA,1));
assert(max(p.scaledA*(s.powerW(:)/p.powerScaleW)-p.scaledB)<1e-7);
end

function testDualInfeasible()
p=coupledProblem(); p.bqos(1)=-1; c=msnoma.config.base;
a=msnoma.optim.solveDual(p,c.numerics); b=msnoma.optim.solveCp4a(p,c.numerics);
assert(strcmp(a.status,'noise_alone_violates_qos') && all(isnan(a.powerW(:))));
assert(strcmp(b.status,'noise_alone_violates_qos') && all(isnan(b.powerW(:))));
end

function testCp4aMode()
p=budgetProblem(); c=msnoma.config.base; c.numerics.cp4a.mode='paper_literal';
failed=false;
try
    msnoma.optim.solveCp4a(p,c.numerics);
catch exception
    failed=strcmp(exception.identifier,'msnoma:optim:Cp4aMode');
end
assert(failed);
end

function testCarrierAudit()
c=msnoma.config.base; n=msnoma.model.numerology(c);
a=msnoma.optim.carrierGroupAudit(n);
assert(a.printedCarrierCount==3140 && a.actualCarrierCount==1658 && a.outOfBandCount==1482);
assert(~a.isExactPartition && ~a.solverUsesPrintedGroups);
end

function testPriceCoupling()
p=coupledProblem(); l=zeros(size(p.scaledB));
l(p.rowMap.qos(2))=.5;
price=p.scaledA'*l;
assert(all(price>0)); % Incoming off-BS communication derivative is included.
[~,~,z]=msnoma.optim.dualValue(l,p);
stationary=max(p.scaledLb,min(p.scaledUb,sqrt(p.scaledC./price)));
assert(norm(z-stationary,inf)<1e-12);
l(p.rowMap.hearability(1))=.2;
price=p.scaledA'*l;
[~,~,z]=msnoma.optim.dualValue(l,p);
assert(price(1)<0 && z(1)==p.scaledUb(1));
end
