function solution = equalPower(problem,policy)
%EQUALPOWER Budget equal or globally QoS-scaled equal, without balancing gains.
powerW = repmat(problem.budgetW/problem.M,1,problem.M);
scale = 1;
switch char(policy)
    case 'budget_equal'
        % Intentionally retain possible QoS/hearability violations for auditing.
    case 'qos_scaled_equal'
        interference = problem.Aqos*powerW(:);
        active = interference>0;
        scale = max(0,min([1;problem.bqos(active)./interference(active)]));
        powerW = powerW*scale;
    otherwise
        error('msnoma:optim:Policy','Unknown equal-power policy.');
end
solution.algorithm = char(policy);
solution.status = 'baseline_evaluated_not_optimized';
solution.powerW = powerW;
solution.scale = scale;
solution.objectiveM2 = sum(problem.c./powerW(:));
solution.exitflag = NaN;
solution.multipliers = [];
solution.history = [];
solution.diagnostics = msnoma.optim.diagnostics(problem,solution);
end
