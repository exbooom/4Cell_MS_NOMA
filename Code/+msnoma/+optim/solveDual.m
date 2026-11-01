function solution = solveDual(problem,opts,feasibility)
%SOLVEDUAL Full-coupling box-constrained dual reference, independent of central.
if nargin<3, feasibility=msnoma.optim.checkFeasibility(problem,opts); end
solution=struct('algorithm','dual','status',feasibility.status, ...
    'powerW',nan(problem.K,problem.M),'objectiveM2',NaN,'exitflag',NaN, ...
    'multipliers',[],'history',[],'solverOutput',feasibility.message,'diagnostics',[]);
if ~feasibility.feasible, return; end
started=tic;
state=msnoma.optim.dualActiveSet(problem,opts);
solution=msnoma.optim.dualSolution(problem,opts,state,'dual');
solution.elapsedS=toc(started);
solution.model='full_coupling_Lagrange_dual_box_scalar_minimizer';
solution.usedCentralFallback=false;
end
