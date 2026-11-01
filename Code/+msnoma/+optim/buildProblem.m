function problem = buildProblem(A,qos,factor,geom,ch,cfg,num)
%BUILDPROBLEM Complete convex model: QoS, budgets, and all-pair hearability.
assert(all(geom.valid),'Optimization requires full-rank fixed geometry.');
nv = num.K*num.M;
B = kron(ones(1,num.M),eye(num.K));
count = num.M*num.K*(num.K-1);
columnIds = zeros(count,2);
coefficients = zeros(count,2);
hearMap = zeros(count,3);
r = 0;
for m = 1:num.M
    for k = 1:num.K
        for j = [1:k-1,k+1:num.K]
            r = r+1;
            columnIds(r,:) = [k+(m-1)*num.K,j+(m-1)*num.K];
            coefficients(r,:) = [-ch.gP(k,m),cfg.models.hearability.ratio*ch.gP(j,m)];
            hearMap(r,:) = [k,j,m];
        end
    end
end
rowIds = repmat((1:count)',1,2);
C = sparse(rowIds(:),columnIds(:),coefficients(:),count,nv);
problem.K = num.K;
problem.M = num.M;
problem.N = num.N;
problem.G = num.G;
problem.c = reshape(geom.lambda.^2.*factor.RrangeM2W/num.M,[],1);
problem.Aqos = A;
problem.bqos = qos.limitWPerHz(:);
problem.Bbudget = B;
problem.budgetW = repmat(cfg.paper.positioningBudgetW,num.K,1);
problem.Chear = C;
problem.q = cfg.models.hearability.ratio;
problem.hearabilityComparisonTolerance=cfg.models.hearability.comparisonTolerance;
problem.Aineq = [A;B;full(C)];
problem.bineq = [problem.bqos;problem.budgetW;zeros(count,1)];
problem.lbW = repmat(cfg.numerics.minimumPowerW,nv,1);
problem.ubW = repmat(problem.budgetW,num.M,1);
problem.powerScaleW = cfg.numerics.powerScaleW;
problem.objectiveScaleM2 = sum(problem.c/problem.powerScaleW);
problem.rowScale = max(max(abs(problem.Aineq*problem.powerScaleW),[],2),abs(problem.bineq));
problem.rowScale(problem.rowScale==0) = 1;
problem.scaledA = problem.Aineq*problem.powerScaleW./problem.rowScale;
problem.scaledB = problem.bineq./problem.rowScale;
problem.scaledC = problem.c/(problem.powerScaleW*problem.objectiveScaleM2);
problem.scaledLb = problem.lbW/problem.powerScaleW;
problem.scaledUb = problem.ubW/problem.powerScaleW;
problem.rowMap.qos = (1:num.K*num.N)';
problem.rowMap.budget = num.K*num.N+(1:num.K)';
problem.rowMap.hearability = num.K*num.N+num.K+(1:count)';
problem.rowMap.hearabilityIds = hearMap;
problem.variableMap = reshape(1:nv,num.K,num.M);
problem.geometrySource=geom.source;
problem.objective = ['mean_horizontal_position_variance_m2_fixed_' geom.source '_Eq10'];
problem.provenance = cfg.provenance;
end
