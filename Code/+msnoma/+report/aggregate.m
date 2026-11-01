function summary = aggregate(run)
%AGGREGATE Coverage, conditional accuracy, full objective and solver quality.
count = numel(run.algorithms);
rows = repmat(struct('algorithm','','status','','nTotal',run.num.M,'nValid',0, ...
    'coverage',0,'conditionalMeanM',NaN,'conditionalMedianM',NaN, ...
    'objectiveMeanM2',NaN,'maxBer',NaN,'communicationQosSatisfied',false, ...
    'maxNormalizedViolation',NaN,'kktStationarity',NaN,'kktComplementarity',NaN),count,1);
for i = 1:count
    a = run.algorithms(i);
    rows(i).algorithm = a.id;
    rows(i).status = a.solution.status;
    rows(i).objectiveMeanM2 = a.solution.objectiveM2;
    if ~isempty(a.solution.diagnostics)
        rows(i).maxNormalizedViolation = a.solution.diagnostics.maxNormalizedViolation;
        rows(i).kktStationarity = a.solution.diagnostics.kktStationarity;
        rows(i).kktComplementarity = a.solution.diagnostics.kktComplementarity;
    end
    if isempty(a.evaluation), continue; end
    e = a.evaluation;
    rows(i).nValid = nnz(e.valid);
    rows(i).coverage = e.coverage;
    if any(e.valid)
        rows(i).conditionalMeanM = mean(e.positionStdM(e.valid));
        rows(i).conditionalMedianM = median(e.positionStdM(e.valid));
    end
    rows(i).maxBer = max(e.ber(:));
    rows(i).communicationQosSatisfied = e.communicationQosSatisfied;
end
summary = struct2table(rows);
end
