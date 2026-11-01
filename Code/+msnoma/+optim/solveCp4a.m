function solution = solveCp4a(problem,opts,feasibility)
%SOLVECP4A Corrected hierarchical candidate, not literal paper Algorithm 1.
% Outer BS budget prices; inner full-coupling QoS/hearability dual. Both are
% minimized with analytic-gradient fmincon, NOT the printed subgradient rule.
% Oracle geometry remains fixed; this is not first-fix position feedback.
assert(strcmp(opts.cp4a.mode,'full_coupling_hierarchical_candidate'), ...
    'msnoma:optim:Cp4aMode','Only the explicitly corrected candidate is supported.');
if nargin<3, feasibility=msnoma.optim.checkFeasibility(problem,opts); end
solution=struct('algorithm','cp4a_consistent_candidate','status',feasibility.status, ...
    'powerW',nan(problem.K,problem.M),'objectiveM2',NaN,'exitflag',NaN, ...
    'multipliers',[],'history',[],'solverOutput',feasibility.message,'diagnostics',[]);
if ~feasibility.feasible, return; end
started=tic;
budgetRows=problem.rowMap.budget;
innerRows=setdiff((1:size(problem.scaledA,1))',budgetRows);
warm=zeros(size(problem.scaledB));
evaluations=0;
innerIterations=0;
history=zeros(0,7);
lastNu=[]; lastValue=[]; lastGradient=[]; lastState=[];
innerFailed=false;
options=optimoptions('fmincon','Algorithm','sqp','Display',opts.solverDisplay, ...
    'SpecifyObjectiveGradient',true,'OptimalityTolerance',opts.cp4a.outerTolerance, ...
    'ConstraintTolerance',1e-12,'StepTolerance',1e-12, ...
    'MaxIterations',opts.cp4a.maxOuterIterations,'MaxFunctionEvaluations',1000, ...
    'OutputFcn',@recordOuter);
try
    [nu,~,flag,output]=fmincon(@outerObjective,zeros(numel(budgetRows),1), ...
        [],[],[],[],zeros(numel(budgetRows),1),[],[],options);
catch exception
    if ~innerFailed, rethrow(exception); end
    solution.status='hierarchical_inner_dual_not_verified';
    solution.solverOutput=exception.message;
    solution.elapsedS=toc(started);
    return;
end
[~,~]=outerObjective(nu);
state=lastState;
state.exitflag=flag;
state.output=output;
state.history=history;
state.historyColumns={'outerIteration','scaledDualBound','maxPrimalViolation', ...
    'relativeGap','budgetProjectedGradient','innerProjectedGradient','innerEvaluations'};
solution=msnoma.optim.dualSolution(problem,opts,state,'cp4a_consistent_candidate');
solution.elapsedS=toc(started);
solution.mode=opts.cp4a.mode;
solution.literalPaperAlgorithm=false;
solution.innerMethod='complete_coupled_dual_constraint_generation_analytic_gradient_fmincon';
solution.outerMethod='BS_budget_price_value_function_SQP_analytic_envelope_gradient';
solution.geometryMode='fixed_oracle_no_first_fix_or_position_feedback';
solution.innerEvaluations=evaluations;
solution.innerIterations=innerIterations;
solution.usedCentralFallback=false;
if isfield(problem,'G')
    solution.groupAudit=msnoma.optim.carrierGroupAudit(problem);
else
    solution.groupAudit=struct('status','legacy_problem_missing_numerology_G');
end

    function [value,gradient]=outerObjective(price)
        if isequal(price,lastNu)
            value=lastValue; gradient=lastGradient; return;
        end
        fixed=zeros(size(problem.scaledB));
        fixed(budgetRows)=price;
        stateInner=msnoma.optim.dualActiveSet(problem,opts,innerRows,fixed,warm);
        if stateInner.maxVariableViolation>opts.cp4a.innerTolerance ...
                || stateInner.innerComplementarity>opts.cp4a.innerTolerance
            innerFailed=true;
            error('msnoma:optim:Cp4aInner','Inner dual failed feasibility/complementarity checks.');
        end
        warm=stateInner.lambda;
        warm(budgetRows)=0;
        evaluations=evaluations+1;
        innerIterations=innerIterations+stateInner.iterations;
        value=stateInner.negativeDual;
        gradient=problem.scaledB(budgetRows)-problem.scaledA(budgetRows,:)*stateInner.z;
        lastNu=price; lastValue=value; lastGradient=gradient; lastState=stateInner;
    end
    function stop=recordOuter(price,values,phase)
        stop=false;
        if strcmp(phase,'init') || strcmp(phase,'iter')
            [v,g]=outerObjective(price);
            violation=max([0;problem.scaledA*lastState.z-problem.scaledB]);
            relativeGap=(lastState.scaledObjective+v)/max(lastState.scaledObjective,realmin);
            history(end+1,:)=[values.iteration,-v,violation,relativeGap, ...
                norm(price-max(0,price-g),inf),lastState.projectedGradient,evaluations];
        end
    end
end
