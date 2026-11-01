function solution = solveSubgradient(problem,opts,feasibility)
%SOLVESUBGRADIENT Hierarchical projected-price experiment, no fmincon calls.
% Full coupling, explicit step rules, oracle geometry; not literal Eq24/33.
if nargin<3, feasibility=msnoma.optim.checkFeasibility(problem,opts); end
solution=struct('algorithm','cp4a_subgradient_candidate','status',feasibility.status, ...
    'powerW',nan(problem.K,problem.M),'objectiveM2',NaN,'exitflag',NaN, ...
    'multipliers',[],'history',[],'solverOutput',feasibility.message,'diagnostics',[]);
if ~feasibility.feasible, return; end
started=tic; s=opts.subgradient;
budget=problem.rowMap.budget;
innerRows=setdiff((1:size(problem.scaledA,1))',budget);
nu=zeros(numel(budget),1); warm=zeros(size(problem.scaledB));
innerCalls=0; innerIterations=0; innerTrace=cell(0,1);
history=zeros(0,8); outerPrices=zeros(0,numel(budget));
reason='outer_iteration_limit'; acceptedStep=0;
step=s.outerInitialStep; lastNu=[]; lastState=[];
[state,value,gradient]=evaluate(nu);
for outer=0:s.maxOuterIterations
    state.history=history;
    state.historyColumns=columns();
    certificate=msnoma.optim.dualSolution(problem,opts,state,'cp4a_subgradient_candidate');
    violation=certificate.diagnostics.maxNormalizedViolation;
    projected=norm(nu-max(0,nu-gradient),inf);
    history(end+1,:)=[outer,-value,violation,certificate.relativeDualityGap, ...
        projected,state.innerComplementarity,acceptedStep,innerCalls]; %#ok<AGROW>
    outerPrices(end+1,:)=nu'; %#ok<AGROW>
    if strcmp(certificate.status,'converged_verified') && state.converged
        reason='primal_dual_certificate_passed'; break;
    end
    if ~state.converged, reason=['inner_' state.reason]; break; end
    if outer==s.maxOuterIterations, break; end
    if strcmp(s.stepRule,'diminishing')
        acceptedStep=s.outerInitialStep/(outer+1)^s.diminishingExponent;
        nu=max(0,nu-acceptedStep*gradient);
        [state,value,gradient]=evaluate(nu);
        continue;
    end
    step=min(s.maximumStep,step*s.stepGrowth);
    accepted=false;
    for backtrack=1:s.maxBacktracks
        nextNu=max(0,nu-step*gradient);
        delta=nextNu-nu;
        [nextState,nextValue,nextGradient]=evaluate(nextNu);
        if ~nextState.converged, reason=['inner_' nextState.reason]; break; end
        % Reduced dual uses approximate inner solves: final acceptance is
        % always by the full certificate, not a small outer price change.
        if nextValue<=value+1e-4*(gradient'*delta)+1e-13*max(1,abs(value))
            accepted=true; break;
        end
        step=step*s.backtrackFactor;
    end
    if ~accepted
        if nextState.converged, reason='outer_line_search_failed'; end
        break;
    end
    nu=nextNu; state=nextState; value=nextValue; gradient=nextGradient;
    acceptedStep=step;
end
state.history=history; state.historyColumns=columns();
state.output=struct('iterations',outer,'innerCalls',innerCalls, ...
    'innerIterations',innerIterations,'message',reason,'stepRule',s.stepRule);
state.exitflag=double(strcmp(reason,'primal_dual_certificate_passed'));
solution=msnoma.optim.dualSolution(problem,opts,state,'cp4a_subgradient_candidate');
if ~strcmp(reason,'primal_dual_certificate_passed')
    solution.status=['subgradient_' reason];
end
solution.elapsedS=toc(started);
solution.innerCalls=innerCalls; solution.innerIterations=innerIterations;
solution.innerTrace=innerTrace;
solution.outerPriceHistory=outerPrices;
solution.stepConfiguration=s;
solution.literalPaperAlgorithm=false;
solution.usedCentralFallback=false;
solution.usedPriceOptimizer=false;
solution.feasibilityPrecheck='shared_linprog_or_supplied_result';
solution.stepRule=s.stepRule;
solution.geometryMode='fixed_oracle_no_first_fix_or_position_feedback';
if isfield(problem,'geometrySource')
    solution.geometryMode=['fixed_' problem.geometrySource '_within_power_solve'];
end
solution.innerMethod='projected_full_coupling_prices_with_explicit_coordinate_scale';
solution.outerMethod='projected_BS_budget_prices';
solution.groupAudit=msnoma.optim.carrierGroupAudit(problem);

    function [st,v,g]=evaluate(price)
        if isequal(price,lastNu)
            st=lastState; v=st.negativeDual;
        else
            fixed=zeros(size(problem.scaledB)); fixed(budget)=price;
            if ~s.warmStartInner, warm=zeros(size(warm)); end
            st=msnoma.optim.projectedDual(problem,opts,innerRows,fixed,warm);
            warm=st.lambda; warm(budget)=0;
            innerCalls=innerCalls+1; innerIterations=innerIterations+st.iterations;
            innerTrace{end+1}=struct('budgetPrice',price,'state',st);
            lastNu=price; lastState=st;
            v=st.negativeDual;
        end
        g=problem.scaledB(budget)-problem.scaledA(budget,:)*st.z;
    end
    function names=columns()
        names={'outerIteration','scaledDualBound','maxPrimalViolation', ...
            'relativeGap','budgetProjectedGradient','innerComplementarity', ...
            'acceptedOuterStep','innerCalls'};
    end
end
