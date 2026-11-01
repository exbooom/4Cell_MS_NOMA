function feedback = positionFeedback(context,cfg,feasibility,packet)
%POSITIONFEEDBACK Unit-lambda first fix, then estimated geometry, no truth fill.
o=cfg.experiment.positionFeedback; s=cfg.models.sampledReceiver;
feedback=struct('status','disabled','iterations',0,'history',[], ...
    'converged',false,'firstFixUnitLambda',true,'truthGeometryUsed',false, ...
    'knownChannelCsiAssumed',true,'usedOracleFallback',false, ...
    'commonRandomNumbers',o.commonRandomNumbers,'options',o);
if ~o.enabled, return; end
if nargin<4, packet=msnoma.waveform.makePacket(context.num,s); end
K=context.num.K; M=context.num.M;
geom=struct('lambda',ones(K,M),'valid',true(M,1),'source','first_fix_unit_lambda');
history=struct('iteration',{},'geometrySource',{},'lambda',{},'solution',{}, ...
    'acquisition',{},'positionChangeM',{},'relativePowerChange',{});
previousPosition=[]; previousPower=[];
s.userIds=1:M; s.referenceAlgorithm=['feedback_' o.solverId];
for iteration=1:o.maxFixes
    problem=msnoma.optim.buildProblem(context.A,context.qos,context.rangeFactor, ...
        geom,context.ch,cfg,context.num);
    switch o.solverId
        case 'central'
            solution=msnoma.optim.solveCentral(problem,cfg.numerics,feasibility);
        case 'dual'
            solution=msnoma.optim.solveDual(problem,cfg.numerics,feasibility);
        case 'cp4a_subgradient_candidate'
            solution=msnoma.optim.solveSubgradient(problem,cfg.numerics,feasibility);
        otherwise
            error('msnoma:feedback:Solver','Unsupported feedback solver.');
    end
    acquisition=[]; positionChange=NaN; powerChange=NaN;
    if ~strcmp(solution.status,'converged_verified')
        feedback.status=['power_' solution.status];
    else
        receiverOpts=s;
        if ~o.commonRandomNumbers, receiverOpts.noiseSeed=s.noiseSeed+iteration-1; end
        acquisition=msnoma.experiment.sampledReceiver(context.scene,context.ch,context.num, ...
            solution.powerW,receiverOpts,cfg.assumed.noisePsdWPerHz,packet);
        if ~all(acquisition.validPosition) || ~all(acquisition.validRanges,'all')
            feedback.status='acquisition_or_position_failure';
        elseif ~isempty(previousPosition)
            positionChange=max(sqrt(sum((acquisition.estimatedPositionM-previousPosition).^2,2)));
            powerChange=norm(solution.powerW-previousPower,'fro')/norm(previousPower,'fro');
            if iteration>=o.minimumFixes && positionChange<=o.positionToleranceM ...
                    && powerChange<=o.relativePowerTolerance
                feedback.converged=true; feedback.status='converged';
            end
        end
    end
    history(end+1)=struct('iteration',iteration,'geometrySource',geom.source, ...
        'lambda',geom.lambda,'solution',solution,'acquisition',acquisition, ...
        'positionChangeM',positionChange,'relativePowerChange',powerChange); %#ok<AGROW>
    if feedback.converged || ~strcmp(solution.status,'converged_verified') ...
            || isempty(acquisition) || ~all(acquisition.validPosition) || ~all(acquisition.validRanges,'all')
        break;
    end
    estimatedScene=struct('bsPosM',context.scene.bsPosM, ...
        'pUserPosM',acquisition.estimatedPositionM);
    spec=cfg.models.geometry; spec.source='estimated_sampled_ranges';
    geom=msnoma.geometry.weights(estimatedScene,true(K,M),spec);
    if ~all(geom.valid), feedback.status='estimated_geometry_rank_failure'; break; end
    previousPosition=acquisition.estimatedPositionM; previousPower=solution.powerW;
    feedback.status='fix_limit_reached';
end
feedback.history=history; feedback.iterations=numel(history);
end
