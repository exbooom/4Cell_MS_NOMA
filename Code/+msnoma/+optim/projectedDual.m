function state = projectedDual(problem,opts,variableRows,fixed,warm)
%PROJECTEDDUAL First-order nonnegative dual prices; no optimization solver.
% Default accelerated projected gradient with a checked quadratic line bound.
% Diminishing steps are an explicit experiment, not a convergence certificate.
if nargin<3, variableRows=(1:size(problem.scaledA,1))'; end
if nargin<4, fixed=zeros(size(problem.scaledB)); end
if nargin<5, warm=zeros(size(problem.scaledB)); end
s=opts.subgradient;
assert(any(strcmp(s.stepRule,{'backtracking','diminishing'})), ...
    'msnoma:optim:StepRule','Unsupported projected-price step rule.');
assert(s.diminishingExponent>0.5 && s.diminishingExponent<=1);
validateattributes(s.maxInnerIterations,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(s.recordEvery,{'numeric'},{'scalar','integer','positive'});
assert(s.initialStep>0 && s.maximumStep>=s.initialStep && s.maxBacktracks>=1);
assert(s.backtrackFactor>0 && s.backtrackFactor<1 && s.stepGrowth>=1);
variableRows=variableRows(:);
assert(all(fixed(variableRows)==0));
A=problem.scaledA(variableRows,:); b=problem.scaledB(variableRows);
hearLocal=find(ismember(variableRows,problem.rowMap.hearability));
positiveHear=max(A(hearLocal,:),0);
hearTolerance=s.innerRelativeHearabilityTolerance;
if isfield(problem,'hearabilityComparisonTolerance')
    hearTolerance=min(hearTolerance,0.1*problem.hearabilityComparisonTolerance);
end
fixedPrice=problem.scaledA'*fixed;
fixedConstant=problem.scaledB'*fixed;
% Diagonal change of dual coordinates; all matrix columns remain coupled.
scale=1./max(sqrt(sum(A.^2,2)),1e-12);
Q=A.*scale;
rhs=b.*scale;
x=warm(variableRows)./scale;
y=x; momentum=1; step=s.initialStep;
history=zeros(0,8); converged=false; evaluations=0; acceptedStep=0;
reason='iteration_limit';
for iteration=0:s.maxInnerIterations
    [value,g,z]=evaluate(x);
    lambdaVariable=scale.*x;
    residual=A*z-b;
    violation=max([0;residual]);
    complementarity=norm(lambdaVariable.*residual,inf);
    projected=norm(x-max(0,x-g),inf);
    relativeHear=max([0;residual(hearLocal)./max(positiveHear*z,realmin)]);
    if mod(iteration,s.recordEvery)==0 || iteration==s.maxInnerIterations
        history(end+1,:)=[iteration,-value,sum(problem.scaledC./z), ...
            violation,complementarity,projected,acceptedStep,evaluations]; %#ok<AGROW>
    end
    relativeComplementarity=sum(abs(lambdaVariable.*residual))/max(sum(problem.scaledC./z),realmin);
    if violation<=s.innerFeasibilityTolerance && relativeComplementarity<=s.innerRelativeGapTolerance ...
            && relativeHear<=hearTolerance
        converged=true; reason='inner_feasibility_and_complementarity'; break;
    end
    if iteration==s.maxInnerIterations, break; end
    if strcmp(s.stepRule,'diminishing')
        acceptedStep=s.diminishingInitialStep/(iteration+1)^s.diminishingExponent;
        x=max(0,x-acceptedStep*g);
        continue;
    end
    [~,gY,zY]=evaluate(y);
    step=min(s.maximumStep,step*s.stepGrowth);
    accepted=false;
    for backtrack=1:s.maxBacktracks
        candidate=max(0,y-step*gY);
        [~,~,zNext,wNext]=evaluate(candidate);
        delta=candidate-y;
        % Stable exact dual Bregman remainder. Subtracting two O(1) dual
        % values loses the O(1e-18) improvement required by tight residuals.
        remainder=-sum((zNext-zY).*(wNext-problem.scaledC./(zNext.*zY)));
        bound=sum(delta.^2)/(2*step);
        if remainder<=bound*(1+1e-10)+realmin
            accepted=true; break;
        end
        step=step*s.backtrackFactor;
    end
    if ~accepted, reason='line_search_failed'; break; end
    acceptedStep=step;
    nextMomentum=(1+sqrt(1+4*momentum^2))/2;
    nextY=candidate+(momentum-1)/nextMomentum*(candidate-x);
    % Gradient-based restart is stable when dual-value differences round off.
    if ~s.accelerate || (y-candidate)'*(candidate-x)>0
        nextY=candidate; nextMomentum=1;
    end
    x=candidate; y=nextY; momentum=nextMomentum;
end
[value,g,z]=evaluate(x);
lambda=fixed; lambda(variableRows)=scale.*x;
residual=A*z-b;
if isempty(history) || history(end,1)~=iteration
    history(end+1,:)=[iteration,-value,sum(problem.scaledC./z), ...
        max([0;residual]),norm(lambda(variableRows).*residual,inf), ...
        norm(x-max(0,x-g),inf),acceptedStep,evaluations];
end
state.lambda=lambda; state.z=z; state.negativeDual=value;
state.scaledObjective=sum(problem.scaledC./z);
state.scaledGap=state.scaledObjective+value;
state.projectedGradient=norm(x-max(0,x-g),inf);
state.innerComplementarity=norm(lambda(variableRows).*residual,inf);
state.innerRelativeComplementarity=sum(abs(lambda(variableRows).*residual))/max(state.scaledObjective,realmin);
state.maxVariableViolation=max([0;residual]);
state.maxRelativeHearabilityViolation=max([0;residual(hearLocal)./max(positiveHear*z,realmin)]);
state.iterations=iteration; state.evaluations=evaluations;
state.converged=converged; state.reason=reason;
state.exitflag=double(converged);
state.output=struct('iterations',iteration,'evaluations',evaluations, ...
    'message',reason,'stepRule',s.stepRule, ...
    'accelerated',s.accelerate && strcmp(s.stepRule,'backtracking'));
state.history=history;
state.historyColumns={'innerIteration','scaledDualBound','scaledRecoveredObjective', ...
    'maxVariableViolation','innerComplementarity','coordinateProjectedGradient', ...
    'acceptedStep','functionEvaluations'};
state.coordinateScale=scale;

    function [v,gradient,power,w]=evaluate(price)
        w=fixedPrice+Q'*price;
        power=problem.scaledUb;
        ids=w>0 & problem.scaledC>0;
        power(ids)=sqrt(problem.scaledC(ids)./w(ids));
        power=max(problem.scaledLb,min(problem.scaledUb,power));
        ids=problem.scaledC==0 & w>=0; power(ids)=problem.scaledLb(ids);
        v=fixedConstant+rhs'*price-sum(problem.scaledC./power+w.*power);
        gradient=rhs-Q*power;
        evaluations=evaluations+1;
    end
end
