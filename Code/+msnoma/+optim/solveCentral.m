function solution = solveCentral(problem,opts,feasibility)
%SOLVECENTRAL Full-coupling convex model solved numerically with fmincon.
% This is a centralized reference solver, NOT the paper's CP4A algorithm.
if nargin<3, feasibility = msnoma.optim.checkFeasibility(problem,opts); end
solution = struct('algorithm','central','status','not_started', ...
    'powerW',nan(problem.K,problem.M),'objectiveM2',NaN,'exitflag',NaN, ...
    'multipliers',[],'history',[],'solverOutput',[],'diagnostics',[]);
if ~feasibility.feasible
    solution.status = feasibility.status;
    solution.solverOutput = feasibility.message;
    return;
end
history = zeros(0,4);
options = optimoptions('fmincon','Algorithm','interior-point', ...
    'Display',opts.solverDisplay,'SpecifyObjectiveGradient',true, ...
    'HessianFcn',@(z,~) diag(2*problem.scaledC./z.^3), ...
    'ConstraintTolerance',min(1e-9,opts.feasibilityTolerance/10), ...
    'OptimalityTolerance',1e-8,'StepTolerance',1e-12, ...
    'MaxIterations',opts.maxIterations,'MaxFunctionEvaluations',10000, ...
    'OutputFcn',@recordIteration);
started = tic;
[z,~,flag,output,lambda] = fmincon(@(z) msnoma.optim.objective(z,problem), ...
    feasibility.scaledPoint,problem.scaledA,problem.scaledB,[],[], ...
    problem.scaledLb,problem.scaledUb,[],options);
solution.elapsedS = toc(started);
solution.powerW = reshape(z*problem.powerScaleW,problem.K,problem.M);
solution.objectiveM2 = sum(problem.c./solution.powerW(:));
solution.exitflag = flag;
solution.multipliers = lambda;
solution.solverOutput = output;
solution.history = history;
solution.historyColumns = {'iteration','scaledObjective','constraintViolation','firstOrderOptimality'};
solution.diagnostics = msnoma.optim.diagnostics(problem,solution);
d = solution.diagnostics;
if flag>0 && d.maxNormalizedViolation<=opts.feasibilityTolerance ...
        && d.kktStationarity<=opts.kktTolerance ...
        && d.kktComplementarity<=opts.kktTolerance
    solution.status = 'converged_verified';
elseif d.maxNormalizedViolation<=opts.feasibilityTolerance
    solution.status = 'feasible_not_verified_optimal';
else
    solution.status = 'solver_result_infeasible';
end

    function stop = recordIteration(~,values,state)
        stop = false;
        if strcmp(state,'init') || strcmp(state,'iter')
            firstOrder = NaN;
            if isfield(values,'firstorderopt'), firstOrder=values.firstorderopt; end
            history(end+1,:) = [values.iteration,values.fval,values.constrviolation,firstOrder];
        end
    end
end
