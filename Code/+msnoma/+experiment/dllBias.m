function out = dllBias(run)
%DLLBIAS Reference-power-only audit; oracle optimizer remains unchanged.
out=struct('status','disabled');
spec=run.cfg.models.dllBias;
if ~spec.enabled, return; end
out.referenceAlgorithm=spec.referenceAlgorithm;
index=find(strcmp({run.algorithms.id},spec.referenceAlgorithm),1);
if isempty(index) || isempty(run.algorithms(index).evaluation)
    out.status='reference_solution_unavailable'; return;
end
power=run.algorithms(index).solution.powerW;
if any(~isfinite(power(:)) | power(:)<=0)
    out.status='reference_power_outside_positive_domain'; return;
end
opts=run.cfg.models.rangingIntegral;
opts.userIds=spec.userIds;
if isempty(opts.userIds), opts.userIds=unique([1,ceil(run.num.M/2),run.num.M]); end
out.bundle=msnoma.ranging.biasKernels(run.ch,run.num,run.cfg.paper.dll,opts,run.integralKernel);
out.result=msnoma.ranging.biasDecomposition(power,run.ch,run.num,run.cfg.paper.dll, ...
    run.cfg.assumed.noisePsdWPerHz,out.bundle,run.rangeFactor,run.geometry);
out.status='computed_fixed_reference_power_variance_bias_decomposition';
out.optimizerChanged=false; out.sampledReceiverImplemented=false;
end
