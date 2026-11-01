function exportConstraints(runDir,algorithm,problem)
%EXPORTCONSTRAINTS Signed SI residuals and normalized constraints, with IDs.
if isempty(algorithm.solution.diagnostics), return; end
rows=size(problem.Aineq,1);
kind=repmat({''},rows,1);
unit=repmat({''},rows,1);
k=zeros(rows,1); n=zeros(rows,1); m=zeros(rows,1); j=zeros(rows,1);
ids=problem.rowMap.qos;
kind(ids)={'qos'}; unit(ids)={'W/Hz'};
k(ids)=problem.rowMap.qosIds(:,1); n(ids)=problem.rowMap.qosIds(:,2);
ids=problem.rowMap.budget;
kind(ids)={'budget'}; unit(ids)={'W'}; k(ids)=(1:problem.K)';
ids=problem.rowMap.hearability;
kind(ids)={'relative_hearability'}; unit(ids)={'W'};
map=problem.rowMap.hearabilityIds;
k(ids)=map(:,1); j(ids)=map(:,2); m(ids)=map(:,3);
signedSi=problem.Aineq*algorithm.solution.powerW(:)-problem.bineq;
signedNormalized=algorithm.solution.diagnostics.scaledResidual;
t=table(kind,k,n,j,m,unit,signedSi,signedNormalized,'VariableNames', ...
    {'constraint','baseStation','communicationCarrier','otherBs','positioningUser', ...
    'siUnit','signedSiResidual','signedNormalizedResidual'});
writetable(t,fullfile(runDir,[algorithm.id '_constraints.csv']));
end
