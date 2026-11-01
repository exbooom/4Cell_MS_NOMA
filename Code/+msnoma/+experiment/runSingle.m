function run = runSingle(cfg,algorithmIds,snapshot)
%RUNSINGLE Complete fixed-geometry single-scene path, shared by all algorithms.
if nargin<2, algorithmIds=cfg.experiment.algorithmIds; end
if nargin<3, snapshot=[]; end
msnoma.config.validate(cfg);
started = tic;
num = msnoma.model.numerology(cfg);
kernels = msnoma.model.spectralKernels(num);
if isempty(snapshot)
    scene = msnoma.scenario.make(cfg,num);
else
    scene = msnoma.scenario.fromSnapshot(snapshot,num);
end
ch = msnoma.scenario.channels(scene,num,cfg.assumed.propagation);
berModel = msnoma.comm.makeBerModel(cfg.models.ber);
[A,qosRowMap] = msnoma.comm.interferenceMatrix(ch,num,kernels);
qos = msnoma.comm.qosLimit(ch.gCServing,num.PcW,num.TcS, ...
    cfg.assumed.noisePsdWPerHz,cfg.paper.berThreshold,berModel);
rangeFactor = msnoma.ranging.factorTwoTerm(ch,num,cfg.assumed.noisePsdWPerHz,cfg.paper.dll);
geom = msnoma.geometry.weights(scene,true(num.K,num.M),cfg.models.geometry);
problem = msnoma.optim.buildProblem(A,qos,rangeFactor,geom,ch,cfg,num);
problem.rowMap.qosIds = qosRowMap;
feasibility = msnoma.optim.checkFeasibility(problem,cfg.numerics);
context = struct('cfg',cfg,'num',num,'scene',scene,'ch',ch,'A',A, ...
    'berModel',berModel,'qos',qos,'rangeFactor',rangeFactor,'geom',geom);
integralKernel=[];
if cfg.models.rangingIntegral.enabled
    integralOpts=cfg.models.rangingIntegral;
    if isempty(integralOpts.userIds), integralOpts.userIds=unique([1,ceil(num.M/2),num.M]); end
    integralKernel=msnoma.ranging.integralKernels(ch,num,cfg.paper.dll,integralOpts);
    context.integralKernel=integralKernel;
end
run.schemaVersion = cfg.schemaVersion;
run.label = 'ASSUMED research scenario; centralized reference is not CP4A';
run.optimizationLabel = 'dual, SQP hierarchy and projected-price candidates; fixed oracle, not literal paper CP4A';
run.cfg = cfg;
run.num = num;
run.scene = scene;
run.ch = ch;
run.kernels = kernels;
run.qos = qos;
run.rangeFactor = rangeFactor;
run.geometry = geom;
run.problem = problem;
run.feasibility = feasibility;
run.integralKernel=integralKernel;
run.algorithms = struct('id',{},'solution',{},'evaluation',{});
for i = 1:numel(algorithmIds)
    id = char(algorithmIds{i});
    if strcmp(id,'central')
        solution = msnoma.optim.solveCentral(problem,cfg.numerics,feasibility);
    elseif strcmp(id,'dual')
        solution = msnoma.optim.solveDual(problem,cfg.numerics,feasibility);
    elseif strcmp(id,'cp4a_consistent_candidate')
        solution = msnoma.optim.solveCp4a(problem,cfg.numerics,feasibility);
    elseif strcmp(id,'cp4a_subgradient_candidate')
        solution = msnoma.optim.solveSubgradient(problem,cfg.numerics,feasibility);
    else
        solution = msnoma.optim.equalPower(problem,id);
    end
    if all(isfinite(solution.powerW(:))) && all(solution.powerW(:)>=0)
        evaluation = msnoma.experiment.evaluate(solution.powerW,context);
    else
        evaluation = [];
    end
    run.algorithms(i) = struct('id',id,'solution',solution,'evaluation',evaluation);
end
run.prsBenchmark=struct('status','disabled');
run.prsComparison=struct('status','disabled');
if cfg.models.prs.enabled
    reference=find(strcmp({run.algorithms.id},cfg.models.prs.referenceAlgorithm),1);
    if ~isempty(reference) && ~isempty(run.algorithms(reference).evaluation) ...
            && all(isfinite(run.algorithms(reference).solution.powerW(:)))
        averagePowerW=sum(run.algorithms(reference).solution.powerW,2);
        run.prsBenchmark=msnoma.prs.positionBenchmark(scene,ch,num,averagePowerW, ...
            cfg.assumed.noisePsdWPerHz,cfg.models.prs,cfg.models.geometry);
        run.prsBenchmark.referenceAlgorithm=cfg.models.prs.referenceAlgorithm;
        if cfg.models.prs.comparison.enabled
            run.prsComparison=msnoma.prs.comparison(scene,ch,num,run.prsBenchmark, ...
                run.algorithms(reference),cfg,berModel);
        end
    else
        run.prsBenchmark.status='reference_solution_unavailable';
        if cfg.models.prs.comparison.enabled
            run.prsComparison.status='reference_solution_unavailable';
        end
    end
end
run.dllBias=msnoma.experiment.dllBias(run);
run.sampledReceiver=struct('status','disabled','sampledReceiverImplemented',true);
run.positionFeedback=struct('status','disabled');
if cfg.models.sampledReceiver.enabled || cfg.experiment.positionFeedback.enabled
    receiverOpts=cfg.models.sampledReceiver;
    packet=msnoma.waveform.makePacket(num,receiverOpts);
    if receiverOpts.enabled
        reference=find(strcmp({run.algorithms.id},receiverOpts.referenceAlgorithm),1);
        if ~isempty(reference) && ~isempty(run.algorithms(reference).evaluation)
            run.sampledReceiver=msnoma.experiment.sampledReceiver(scene,ch,num, ...
                run.algorithms(reference).solution.powerW,receiverOpts,cfg.assumed.noisePsdWPerHz,packet);
        else
            run.sampledReceiver.status='reference_solution_unavailable';
        end
    end
    run.positionFeedback=msnoma.experiment.positionFeedback(context,cfg,feasibility,packet);
end
run.stepSensitivity=msnoma.experiment.stepSensitivity(problem,cfg,feasibility);
run.interferenceAblation=struct('status','disabled');
if cfg.experiment.interferenceAblation.enabled
    reference=find(strcmp({run.algorithms.id},cfg.models.sampledReceiver.referenceAlgorithm),1);
    if ~isempty(reference) && ~isempty(run.algorithms(reference).evaluation)
        run.interferenceAblation=msnoma.experiment.interferenceAblation(ch,num, ...
            run.algorithms(reference).solution.powerW,cfg.models.sampledReceiver, ...
            cfg.assumed.noisePsdWPerHz,cfg.experiment.interferenceAblation.userIds);
    else
        run.interferenceAblation.status='reference_solution_unavailable';
    end
end
run.elapsedS = toc(started);
run.summary = msnoma.report.aggregate(run);
end
