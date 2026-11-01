function fit = estimatePChannels(receivedSpectrum,bank,packet,opts,noisePsd)
%ESTIMATEPCHANNELS Joint P-bank delay/gain fit, without channel/position truth.
% All P codes and relative powers are assumed known; C data are not inputs.
K=size(bank,1); L=packet.fftLength;
assert(size(bank,2)==L && numel(receivedSpectrum)==L);
fit=struct('status','not_fitted','valid',false,'delayS',nan(K,1), ...
    'complexGain',complex(nan(K,1)), 'coarseOrder',[], ...
    'history',zeros(0,4),'iterations',0,'rankRatio',NaN, ...
    'truthUsed',false,'communicationDataUsed',false,'relativeResidualEnergy',NaN, ...
    'gaussianDelayStdSamples',nan(K,1),'gaussianRelativeGainStd',nan(K,1), ...
    'cancellationEligible',false,'eligibilityStatus','fit_invalid', ...
    'precisionGatePassed',false,'residualGatePassed',false, ...
    'maximumRelativeResidualEnergy',opts.interferenceCancellation.maximumRelativeResidualEnergy);
if norm(receivedSpectrum)==0, fit.status='signal_absent'; return; end
residual=receivedSpectrum; lag=nan(K,1); chosen=false(K,1);
o=opts; o.falseAlarmProbability=opts.falseAlarmProbability;
% Coarse SIC selects the observed component giving the largest reduction.
for stage=1:K
    best=-Inf; winner=0; result=[];
    for k=find(~chosen)'
        candidate=msnoma.receiver.estimateDelay(residual,bank(k,:),packet,o,0);
        if ~candidate.valid, continue; end
        reduction=abs(candidate.complexAmplitude)^2*sum(abs(bank(k,:)).^2);
        if reduction>best, best=reduction; winner=k; result=candidate; end
    end
    if winner==0, fit.status='coarse_acquisition_failed'; return; end
    lag(winner)=result.fractionalLagSamples; chosen(winner)=true;
    shifted=bank(winner,:).*exp(-1i*2*pi*packet.frequencyHz*result.delayS);
    residual=residual-result.complexAmplitude*shifted;
    fit.coarseOrder(end+1)=winner;
end
active=packet.frontendMask & any(abs(bank)>0,1);
f=packet.frequencyHz(active)'/packet.sampleRateHz;
B=bank(:,active).'; y=receivedSpectrum(active).';
yScale=norm(y); y=y/yScale;
columnScale=sqrt(sum(abs(B).^2,1));
if any(columnScale==0), fit.status='zero_template_energy'; return; end
B=B./columnScale;
maximum=opts.maximumRangeM/opts.speedMps*packet.sampleRateHz;
n=opts.interferenceCancellation;
[cost,g,r,Q,A,ratio]=evaluate(lag);
if ratio<n.rankTolerance, fit.status='template_rank_deficient'; return; end
converged=false; accepted=0;
for iteration=1:n.maxIterations
    D=(-1i*2*pi*f).*A.*g.';
    J=D-Q*(Q'*D); gradient=real(J'*r);
    H=real(J'*J);
    singular=svd(H);
    if singular(end)<=n.rankTolerance*singular(1)
        fit.status='delay_information_rank_deficient'; break;
    end
    delta=H\gradient;
    delta=delta*min(1,n.maximumStepSamples/max(norm(delta,inf),realmin));
    fit.history(end+1,:)=[iteration,cost,norm(delta,inf),accepted];
    if norm(delta,inf)<=n.delayToleranceSamples
        converged=true; fit.status='converged'; break;
    end
    accepted=0;
    for trial=1:n.maxBacktracks
        step=2^(1-trial); candidate=max(0,min(maximum,lag+step*delta));
        [nextCost,nextG,nextR,nextQ,nextA,nextRatio]=evaluate(candidate);
        if nextRatio>=n.rankTolerance && nextCost<cost
            accepted=step; break;
        end
    end
    if accepted==0, fit.status='line_search_failed'; break; end
    lag=candidate; cost=nextCost; g=nextG; r=nextR; Q=nextQ; A=nextA; ratio=nextRatio;
    fit.status='iteration_limit';
end
fit.iterations=iteration; fit.rankRatio=ratio;
fit.delayS=lag/packet.sampleRateHz;
fit.complexGain=g./columnScale.'*yScale;
fit.relativeResidualEnergy=cost;
fit.assumedNoisePsd=noisePsd;
boundary=any(lag<=10*opts.delayToleranceSamples | maximum-lag<=10*opts.delayToleranceSamples);
if boundary, fit.status='search_boundary_peak'; end
fit.valid=converged && ~boundary;
if fit.valid
    D=(-1i*2*pi*f).*A.*g.'; J=D-Q*(Q'*D);
    H=real(J'*J);
    variance=noisePsd*packet.sampleRateHz*L/yScale^2;
    fit.gaussianDelayStdSamples=sqrt(max(0,diag(H\eye(K)))*variance/2);
    [~,R]=qr(A,0); inverse=R\eye(K);
    fit.gaussianRelativeGainStd=sqrt(sum(abs(inverse).^2,2)*variance)./max(abs(g),realmin);
    gate=msnoma.receiver.cancellationGate(fit,n);
    fit.cancellationEligible=gate.eligible; fit.eligibilityStatus=gate.status;
    fit.precisionGatePassed=gate.precisionPassed;
    fit.residualGatePassed=gate.residualPassed;
end

    function [value,coefficient,errorVector,basis,templates,rankRatio]=evaluate(delays)
        templates=B.*exp(-1i*2*pi*f*delays.');
        [basis,R]=qr(templates,0);
        singularValues=svd(R); rankRatio=singularValues(end)/singularValues(1);
        coefficient=complex(zeros(K,1));
        if rankRatio>=n.rankTolerance, coefficient=R\(basis'*y); end
        errorVector=y-templates*coefficient; value=real(errorVector'*errorVector);
    end
end
