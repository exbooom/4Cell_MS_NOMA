function exportOptimization(run,runDir)
%EXPORTOPTIMIZATION Independent comparison of primal and dual certificates.
count=numel(run.algorithms);
rows=repmat(struct('algorithm','','status','','objectiveM2',NaN, ...
    'dualBoundM2',NaN,'relativeDualityGap',NaN,'relativeToCentral',NaN, ...
    'maxNormalizedViolation',NaN,'maxRelativeHearabilityViolation',NaN, ...
    'kktStationarity',NaN,'kktComplementarity',NaN, ...
    'elapsedS',NaN,'usedCentralFallback',false,'usesFmincon',false,'stepRule','not_applicable'),count,1);
reference=find(strcmp({run.algorithms.id},'central'),1);
central=NaN;
if ~isempty(reference), central=run.algorithms(reference).solution.objectiveM2; end
for i=1:count
    a=run.algorithms(i); s=a.solution;
    rows(i).algorithm=a.id; rows(i).status=s.status;
    rows(i).usesFmincon=ismember(a.id,{'central','dual','cp4a_consistent_candidate'});
    if isfield(s,'stepRule'), rows(i).stepRule=s.stepRule; end
    rows(i).objectiveM2=s.objectiveM2;
    if isfinite(central), rows(i).relativeToCentral=s.objectiveM2/central-1; end
    if isfield(s,'dualBoundM2'), rows(i).dualBoundM2=s.dualBoundM2; end
    if isfield(s,'relativeDualityGap'), rows(i).relativeDualityGap=s.relativeDualityGap; end
    if isfield(s,'elapsedS'), rows(i).elapsedS=s.elapsedS; end
    if isfield(s,'usedCentralFallback'), rows(i).usedCentralFallback=s.usedCentralFallback; end
    if ~isempty(s.diagnostics)
        rows(i).maxNormalizedViolation=s.diagnostics.maxNormalizedViolation;
        rows(i).maxRelativeHearabilityViolation=s.diagnostics.maxRelativeHearabilityViolation;
        rows(i).kktStationarity=s.diagnostics.kktStationarity;
        rows(i).kktComplementarity=s.diagnostics.kktComplementarity;
    end
end
writetable(struct2table(rows),fullfile(runDir,'optimization_comparison.csv'));
audit=msnoma.optim.carrierGroupAudit(run.num);
t=table((1:run.num.M)',audit.printedFirst,audit.printedLast, ...
    audit.printedLast>run.num.N,'VariableNames', ...
    {'positioningUser','printedFirstCarrier','printedLastCarrier','extendsBeyondActualBand'});
writetable(t,fullfile(runDir,'paper_carrier_group_audit.csv'));
end
