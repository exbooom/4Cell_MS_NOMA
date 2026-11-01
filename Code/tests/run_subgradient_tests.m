function report = run_subgradient_tests()
%RUN_SUBGRADIENT_TESTS Price update algebra and independent certificates.
tests={@testProjectedUpdate,@testPlainGradient,@testHierarchicalBudget, ...
    @testCoupledAndReset,@testDiminishingReplay,@testIterationFailure, ...
    @testInvalidRule,@testNoOptimizerCalls,@testRelativeHearabilityGate};
elapsed=zeros(numel(tests),1);
names=cellfun(@func2str,tests,'UniformOutput',false);
for i=1:numel(tests)
    started=tic; tests{i}(); elapsed(i)=toc(started);
    fprintf('PASS %s (%.3f s)\n',names{i},elapsed(i));
end
report=table(names(:),elapsed,'VariableNames',{'test','elapsedS'});
fprintf('All %d projected-price checks passed.\n',height(report));
end

function p=problem(coupled)
if coupled
    p.K=2; p.M=2; p.G=1; p.N=1; p.c=[1;4;9;16];
    p.Aqos=[1,.8,.1,.2;.2,1,.3,.7]; p.bqos=[.5;.7];
    p.Bbudget=[1,0,1,0;0,1,0,1]; p.budgetW=[.7;.9];
    p.Chear=[-1,.2,0,0;.1,-1,0,0;0,0,-1,.3;0,0,.2,-1];
else
    p.K=1; p.M=3; p.N=0; p.G=1; p.c=[1;4;9];
    p.Aqos=zeros(0,3); p.bqos=zeros(0,1);
    p.Bbudget=ones(1,3); p.budgetW=.7; p.Chear=zeros(0,3);
end
p.q=.1; p.powerScaleW=.1; p.objectiveScaleM2=sum(p.c/p.powerScaleW);
p.Aineq=[p.Aqos;p.Bbudget;p.Chear];
p.bineq=[p.bqos;p.budgetW;zeros(size(p.Chear,1),1)];
p.rowScale=max(max(abs(p.Aineq*p.powerScaleW),[],2),abs(p.bineq));
p.rowScale(p.rowScale==0)=1;
p.scaledA=p.Aineq*p.powerScaleW./p.rowScale;
p.scaledB=p.bineq./p.rowScale;
p.scaledC=p.c/(p.powerScaleW*p.objectiveScaleM2);
p.lbW=ones(size(p.c))*1e-10; p.ubW=repmat(p.budgetW,p.M,1);
p.scaledLb=p.lbW/p.powerScaleW; p.scaledUb=p.ubW/p.powerScaleW;
p.rowMap.qos=(1:size(p.Aqos,1))';
p.rowMap.budget=size(p.Aqos,1)+(1:p.K)';
p.rowMap.hearability=size(p.Aqos,1)+p.K+(1:size(p.Chear,1))';
end

function testProjectedUpdate()
p=problem(true); c=msnoma.config.base;
c.numerics.subgradient.stepRule='diminishing';
c.numerics.subgradient.maxInnerIterations=1;
s=msnoma.optim.projectedDual(p,c.numerics);
scale=1./sqrt(sum(p.scaledA.^2,2));
residual=p.scaledA*p.scaledUb-p.scaledB;
expected=scale.*max(0,c.numerics.subgradient.diminishingInitialStep*scale.*residual);
assert(norm(s.lambda-expected,inf)<1e-12 && s.iterations==1);
end

function testPlainGradient()
p=problem(false); c=msnoma.config.base;
c.numerics.subgradient.accelerate=false;
s=msnoma.optim.projectedDual(p,c.numerics);
a=msnoma.optim.dualSolution(p,c.numerics,s,'projected_test');
assert(s.converged && strcmp(a.status,'converged_verified'));
expected=p.budgetW*sqrt(p.c)/sum(sqrt(p.c));
assert(norm(a.powerW(:)-expected,inf)<1e-5);
assert(~s.output.accelerated);
end

function testHierarchicalBudget()
p=problem(false); c=msnoma.config.base;
s=msnoma.optim.solveSubgradient(p,c.numerics);
expected=p.budgetW*sqrt(p.c)/sum(sqrt(p.c));
assert(strcmp(s.status,'converged_verified') && norm(s.powerW(:)-expected,inf)<1e-5);
assert(s.innerCalls>1 && s.outerPriceHistory(end)>0);
assert(~s.usedCentralFallback && ~s.usedPriceOptimizer && ~s.literalPaperAlgorithm);
assert(size(s.outerPriceHistory,1)==size(s.history,1));
end

function testCoupledAndReset()
p=problem(true); c=msnoma.config.base;
a=msnoma.optim.solveCentral(p,c.numerics);
b=msnoma.optim.solveSubgradient(p,c.numerics);
c.numerics.subgradient.warmStartInner=false;
d=msnoma.optim.solveSubgradient(p,c.numerics);
assert(strcmp(b.status,'converged_verified') && strcmp(d.status,'converged_verified'));
assert(abs(b.objectiveM2/a.objectiveM2-1)<1e-6 && abs(d.objectiveM2/a.objectiveM2-1)<1e-6);
assert(b.diagnostics.maxNormalizedViolation<=c.numerics.feasibilityTolerance);
assert(~d.stepConfiguration.warmStartInner);
end

function testDiminishingReplay()
p=problem(true); c=msnoma.config.base;
c.numerics.subgradient.stepRule='diminishing';
c.numerics.subgradient.maxInnerIterations=20;
a=msnoma.optim.projectedDual(p,c.numerics);
b=msnoma.optim.projectedDual(p,c.numerics);
assert(isequal(a.lambda,b.lambda) && isequal(a.history,b.history));
assert(~a.converged && ~a.output.accelerated && a.iterations==20);
end

function testIterationFailure()
p=problem(true); c=msnoma.config.base;
c.numerics.subgradient.maxInnerIterations=0;
a=msnoma.optim.solveSubgradient(p,c.numerics);
assert(strcmp(a.status,'subgradient_inner_iteration_limit') && ~a.usedCentralFallback);
p.bqos(1)=-1;
b=msnoma.optim.solveSubgradient(p,c.numerics);
assert(strcmp(b.status,'noise_alone_violates_qos') && all(isnan(b.powerW(:))));
end

function testInvalidRule()
p=problem(false); c=msnoma.config.base;
c.numerics.subgradient.stepRule='unknown'; rejected=false;
try
    msnoma.optim.projectedDual(p,c.numerics);
catch exception
    rejected=strcmp(exception.identifier,'msnoma:optim:StepRule');
end
assert(rejected);
end

function testNoOptimizerCalls()
for name={'msnoma.optim.projectedDual','msnoma.optim.solveSubgradient'}
    source=fileread(which(name{1}));
    forbidden='\<(fmincon|solveCentral|solveDual|solveCp4a|dualActiveSet)\s*\(';
    assert(isempty(regexp(source,forbidden,'once')));
end
end

function testRelativeHearabilityGate()
p=problem(true); c=msnoma.config.base;
p.hearabilityComparisonTolerance=1e-7;
power=[2e-7*(1-1e-6);1e-6;1e-6;1e-6];
candidate=struct('powerW',reshape(power,p.K,p.M),'multipliers',[]);
d=msnoma.optim.diagnostics(p,candidate);
assert(d.maxNormalizedViolation<1e-7 && d.maxRelativeHearabilityViolation>1e-7);
% Isolate the relative-hearability gate with other optimality checks disabled.
state=struct('z',power/p.powerScaleW,'lambda',zeros(size(p.scaledB)), ...
    'negativeDual',msnoma.optim.dualValue(zeros(size(p.scaledB)),p), ...
    'exitflag',1,'output',struct,'history',[],'historyColumns',{{}});
c.numerics.kktTolerance=Inf; c.numerics.dual.relativeGapTolerance=Inf;
s=msnoma.optim.dualSolution(p,c.numerics,state,'isolated_gate_test');
assert(~strcmp(s.status,'converged_verified'));
end
