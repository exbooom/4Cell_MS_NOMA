function state = dualActiveSet(problem,opts,variableRows,fixed,warm)
%DUALACTIVESET Constraint generation in dual space; no primal solve/fallback.
% Every omitted dual coordinate is zero. Each value is a valid dual bound.
if nargin<3, variableRows=(1:size(problem.scaledA,1))'; end
if nargin<4, fixed=zeros(size(problem.scaledB)); end
if nargin<5, warm=zeros(size(problem.scaledB)); end
variableRows=variableRows(:);
assert(all(fixed(variableRows)==0));
active=intersect(variableRows,find(warm>0));
active=union(active,intersect(variableRows,problem.rowMap.budget));
lambda=fixed;
lambda(active)=warm(active);
history=zeros(0,6);
totalIterations=0;
flag=0;
output=struct;
for pass=1:opts.dual.maxPasses
    if ~isempty(active)
        options=optimoptions('fmincon','Algorithm','interior-point', ...
            'Display',opts.solverDisplay,'SpecifyObjectiveGradient',true, ...
            'HessianFcn',@(x,~) dualHessian(x), ...
            'OptimalityTolerance',opts.dual.optimalityTolerance, ...
            'StepTolerance',1e-14,'ConstraintTolerance',1e-12, ...
            'MaxIterations',opts.dual.maxIterations, ...
            'MaxFunctionEvaluations',20000);
        [x,~,flag,output]=fmincon(@dualObjective,max(lambda(active),1e-5), ...
            [],[],[],[],zeros(numel(active),1),[],[],options);
        lambda(active)=x;
        totalIterations=totalIterations+output.iterations;
    end
    [value,gradient,z]=msnoma.optim.dualValue(lambda(variableRows),problem,variableRows,fixed);
    residual=-gradient;
    projected=lambda(variableRows)-max(0,lambda(variableRows)-gradient);
    f=sum(problem.scaledC./z);
    gap=f+value;
    history(end+1,:)=[pass,numel(active),f,-value,max([0;residual]),norm(projected,inf)]; %#ok<AGROW>
    violating=variableRows(residual>opts.dual.constraintTolerance);
    missing=setdiff(violating,active);
    if isempty(missing)
        break;
    end
    % Add strongest violations in each physical family, retaining every
    % previously priced row. Screening never drops a row from final checks.
    additions=[];
    groups={problem.rowMap.qos,problem.rowMap.budget,problem.rowMap.hearability};
    for group=1:numel(groups)
        ids=intersect(missing,groups{group});
        if isempty(ids), continue; end
        violations=problem.scaledA(ids,:)*z-problem.scaledB(ids);
        [~,order]=sort(violations,'descend');
        additions=[additions;ids(order(1:min(numel(ids),opts.dual.rowsPerPass)))]; %#ok<AGROW>
    end
    active=union(active,additions);
end
state.lambda=lambda;
state.z=z;
state.negativeDual=value;
state.scaledObjective=f;
state.scaledGap=gap;
state.projectedGradient=norm(projected,inf);
state.innerComplementarity=norm(lambda(variableRows).*residual,inf);
state.maxVariableViolation=max([0;residual]);
state.activeRows=active;
state.passes=pass;
state.iterations=totalIterations;
state.exitflag=flag;
state.output=output;
state.history=history;
state.historyColumns={'pass','pricedRows','scaledPrimalObjective', ...
    'scaledDualBound','maxVariableViolation','dualProjectedGradient'};

    function [v,g]=dualObjective(x)
        [v,g]=msnoma.optim.dualValue(x,problem,active,fixed);
    end
    function H=dualHessian(x)
        [~,~,~,H]=msnoma.optim.dualValue(x,problem,active,fixed);
    end
end
