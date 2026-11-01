function out = biasKernels(ch,num,dll,opts,existing)
%BIASKERNELS Full 2x2x2 controlled DLL family, independent of allocated power.
if nargin<5, existing=[]; end
frontends={'global_frontend','local_symmetric'};
discriminators={'exact_sine','small_spacing'};
communication={'finite_carriers','spectral_lines'};
out.variants=struct('id',{},'flags',{},'kernel',{},'reused',{});
index=0;
for f=1:2
    for d=1:2
        for c=1:2
            index=index+1; setting=opts;
            setting.frontendMode=frontends{f}; setting.discriminator=discriminators{d};
            setting.communicationModel=communication{c};
            reuse=~isempty(existing) && isfield(existing,'quadratureOptions') ...
                && isequaln(existing.quadratureOptions,setting);
            if reuse
                kernel=existing;
            else
                kernel=msnoma.ranging.integralKernels(ch,num,dll,setting);
            end
            out.variants(index)=struct('id',sprintf('%s__%s__%s',frontends{f}, ...
                discriminators{d},communication{c}),'flags',[d-1,c-1,f-1], ...
                'kernel',kernel,'reused',reuse);
        end
    end
end
out.userIds=opts.userIds(:)';
out.factorOrder={'sine_linearization','finite_C_PSD_to_spectral_lines','frontend_global_to_local'};
out.baselineIndex=1;
out.allConverged=all(arrayfun(@(v) all(v.kernel.converged),out.variants));
out.model='controlled_2x2x2_DLL_family_fixed_channels_power_and_loop';
end
