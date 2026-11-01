function packet = makePacket(num,opts)
%MAKEPACKET Finite synchronous rectangular BPSK/QPSK complex-baseband packet.
% No CP, channel, user coordinates or receiver delay is used here.
validateattributes(opts.oversampling,{'numeric'},{'scalar','integer','>=',2});
validateattributes(opts.codeChips,{'numeric'},{'scalar','integer','>=',31});
requestedChips=opts.codeChips;
L=ceil(requestedChips/num.G)*num.G; Fs=opts.oversampling*num.BHz;
sps=(num.M+1)*opts.oversampling; active=L*sps;
guard=ceil(opts.maximumRangeM/num.c0Mps*Fs)+sps;
fftLength=2^nextpow2(2*(active+guard));
assert(fftLength<=opts.maximumFftLength,'msnoma:waveform:MemoryLimit', ...
    'Packet FFT length exceeds configured memory limit.');
f=(0:fftLength-1)*Fs/fftLength; f(f>=Fs/2)=f(f>=Fs/2)-Fs;
mask=abs(f)<=num.BfeHz/2;
assert(num.BfeHz<=Fs,'Sample rate must cover the full frontend bandwidth.');
stream=RandStream('mt19937ar','Seed',opts.sequenceSeed);
prn=zeros(num.K,num.M); sourcePath=''; reused=false;
if strcmp(opts.codeMode,'author_L1CP')
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    sourcePath=fullfile(root,'MS_NOMA','L1CPcodes.mat');
    data=load(sourcePath,'L1CPcodes'); bank=data.L1CPcodes;
    assert(size(bank,1)>=num.K && all(abs(bank(:))==1));
    chips=zeros(num.K,num.M,L);
    for m=1:num.M
        for k=1:num.K
            prn(k,m)=mod(k+(m-1)*num.K-1,size(bank,1))+1;
            chips(k,m,:)=bank(prn(k,m),mod(0:L-1,size(bank,2))+1);
        end
    end
    reused=num.K*num.M>size(bank,1);
elseif strcmp(opts.codeMode,'seeded_independent_bpsk')
    chips=2*randi(stream,[0,1],num.K,num.M,L)-1;
else
    error('msnoma:waveform:CodeMode','Unsupported sequence mode.');
end
Nsymbol=(num.N+1)*opts.oversampling;
symbolsCount=ceil(active/Nsymbol); center=(num.N+1)/2;
integerCenter=floor(center); bins=mod((1:num.N)-integerCenter,Nsymbol)+1;
t=(0:active-1)/Fs;
communication=complex(zeros(num.K,active));
symbols=complex(zeros(num.K,num.N,symbolsCount));
for k=1:num.K
    q=((2*randi(stream,[0,1],num.N,symbolsCount)-1) ...
        +1i*(2*randi(stream,[0,1],num.N,symbolsCount)-1))/sqrt(2);
    symbols(k,:,:)=q;
    grid=complex(zeros(Nsymbol,symbolsCount)); grid(bins,:)=sqrt(num.PcW)*q;
    signal=reshape(ifft(grid,[],1)*Nsymbol,1,[]);
    communication(k,:)=signal(1:active).*exp(-1i*2*pi*(center-integerCenter)*num.deltaFcHz*t);
end
packet=struct('sampleRateHz',Fs,'samplesPerChip',sps,'activeSamples',active, ...
    'fftLength',fftLength,'frequencyHz',f,'frontendMask',mask, ...
    'positioningCarrierHz',((1:num.M)*num.G-center)*num.deltaFcHz, ...
    'chips',chips,'communicationSymbols',symbols,'communication',communication, ...
    'prn',prn,'codeMode',opts.codeMode,'codeSourcePath',sourcePath, ...
    'requestedCodeChips',requestedChips,'actualCodeChips',L, ...
    'codeReusedAcrossDifferentCarriers',reused,'sequenceSeed',opts.sequenceSeed, ...
    'activeDurationS',active/Fs,'recordDurationS',fftLength/Fs, ...
    'frontendBandwidthHz',num.BfeHz, ...
    'communicationPsdWPerHz',msnoma.waveform.communicationPsd(f,num), ...
    'model','finite_sampled_rectangular_packet_no_CP_no_periodic_code_continuation', ...
    'frontend','global_ideal_rectangular_DFT_projection', ...
    'delayModel','zero_padded_finite_record_DFT_bandlimited_fractional_shift');
packet.positioningPsdPerW=zeros(num.M,fftLength);
for m=1:num.M
    packet.positioningPsdPerW(m,:)=num.TpS*msnoma.model.sincpi( ...
        (f-packet.positioningCarrierHz(m))*num.TpS).^2;
end
end
