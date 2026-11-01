function kernel = integralKernels(ch,num,dll,opts)
%INTEGRALKERNELS Eq51-54 with exact sine and the finite communication PSD.
% Integrate x=f_local*Tp. RF carrier Hz never enters the discriminator.
ids=opts.userIds(:)';
if ~isfield(opts,'communicationModel'), opts.communicationModel='finite_carriers'; end
assert(any(strcmp(opts.communicationModel,{'finite_carriers','spectral_lines'})));
assert(~isempty(ids) && all(ids>=1 & ids<=num.M & ids==fix(ids)));
assert(numel(unique(ids))==numel(ids),'Duplicate selected positioning user.');
assert(any(strcmp(opts.frontendMode,{'local_symmetric','global_frontend'})));
assert(any(strcmp(opts.discriminator,{'exact_sine','small_spacing'})));
validateattributes(dll.earlyLateSpacingChips,{'numeric'},{'scalar','positive','finite'});
validateattributes(opts.initialSamplesPerCarrier,{'numeric'},{'scalar','integer','>=',4});
validateattributes(opts.maxRefinements,{'numeric'},{'scalar','integer','>=',1});
validateattributes(opts.relativeTolerance,{'numeric'},{'scalar','positive','<',1});
selected=numel(ids);
kernel.userIds=ids;
kernel.A0PerS=zeros(1,selected);
kernel.noiseIntegral=zeros(1,selected);
kernel.communicationIntegralS=zeros(1,selected);
kernel.positioningIntegralS=zeros(1,selected);
kernel.localBoundsHz=zeros(selected,2);
kernel.relativeChange=zeros(1,selected);
kernel.samplesPerCarrier=zeros(1,selected);
kernel.converged=false(1,selected);
kernel.refinementHistory=cell(1,selected);
for i=1:selected
    m=ids(i);
    if strcmp(opts.frontendMode,'local_symmetric')
        bounds=[-num.BfeHz/2,num.BfeHz/2];
    else
        center=(num.N+1)*num.deltaFcHz/2;
        bounds=center+[-num.BfeHz/2,num.BfeHz/2]-m*num.deltaFpHz;
    end
    points=opts.initialSamplesPerCarrier;
    previous=integrateOne(bounds*num.TpS,m,ch,num,dll,opts,points);
    history=[points,NaN,previous];
    for refinement=1:opts.maxRefinements
        points=points*2;
        current=integrateOne(bounds*num.TpS,m,ch,num,dll,opts,points);
        change=max(abs(current-previous)./max(abs(current),realmin));
        history(end+1,:)=[points,change,current]; %#ok<AGROW>
        if change<=opts.relativeTolerance
            kernel.converged(i)=true;
            break;
        end
        previous=current;
    end
    assert(isfinite(current(1)) && current(1)>0,'Nonpositive discriminator slope.');
    kernel.A0PerS(i)=current(1);
    kernel.noiseIntegral(i)=current(2);
    kernel.communicationIntegralS(i)=current(3);
    kernel.positioningIntegralS(i)=current(4);
    kernel.localBoundsHz(i,:)=bounds;
    kernel.relativeChange(i)=change;
    kernel.samplesPerCarrier(i)=points;
    kernel.refinementHistory{i}=history;
end
kernel.frontendMode=opts.frontendMode;
kernel.discriminator=opts.discriminator;
kernel.communicationModel=opts.communicationModel;
kernel.quadratureOptions=opts;
kernel.aPerS=dll.loopBandwidthHz*(1-0.5*dll.loopBandwidthHz*dll.coherentTimeS);
kernel.c0Mps=num.c0Mps;
kernel.PcW=num.PcW;
kernel.model='finite_PSD_finite_band_Eq6_Eq51_55';
if strcmp(opts.communicationModel,'spectral_lines')
    kernel.model='finite_C_tone_sum_zero_linewidth_limit_same_DLL_integrals';
end
kernel.historyColumns={'samplesPerCarrier','maxRelativeChange', ...
    'A0PerS','noiseIntegral','communicationIntegralS','positioningIntegralS'};
if opts.requireConvergence
    assert(all(kernel.converged),'msnoma:ranging:Convergence', ...
        'DLL quadrature failed refinement tolerance; inspect refinementHistory.');
end
end

function values=integrateOne(bounds,m,ch,num,dll,opts,samples)
intervals=ceil(diff(bounds)*num.G*samples);
x=linspace(bounds(1),bounds(2),intervals+1)';
q=msnoma.model.sincpi(x).^2;
if strcmp(opts.discriminator,'exact_sine')
    u=sin(pi*dll.earlyLateSpacingChips*x);
else
    u=pi*dll.earlyLateSpacingChips*x;
end
weighted=q.*u.^2;
if strcmp(opts.communicationModel,'finite_carriers')
    communication=msnoma.ranging.communicationPsd(x/num.TpS,m,ch,num,opts.chunkSize);
    communicationIntegral=trapz(x,communication.*weighted);
else
    % Zero-linewidth tones retain the finite list, gains and frontend bounds.
    toneX=(1:num.N)'/num.G-m;
    gain=reshape(sum(ch.gCP(m,:,:),2),num.N,1);
    if strcmp(opts.discriminator,'exact_sine')
        toneU=sin(pi*dll.earlyLateSpacingChips*toneX);
    else
        toneU=pi*dll.earlyLateSpacingChips*toneX;
    end
    boundaryTolerance=1e-12*max(1,max(abs(bounds)));
    inBand=double(toneX>bounds(1)+boundaryTolerance & toneX<bounds(2)-boundaryTolerance);
    onBoundary=abs(toneX-bounds(1))<=boundaryTolerance | abs(toneX-bounds(2))<=boundaryTolerance;
    inBand(onBoundary)=.5;
    communicationIntegral=num.TpS*sum(inBand.*gain.*msnoma.model.sincpi(toneX).^2.*toneU.^2);
end
values=[trapz(x,x.*q.*u)/num.TpS, ...
    trapz(x,weighted),communicationIntegral, ...
    num.TpS*trapz(x,q.*weighted)];
end
