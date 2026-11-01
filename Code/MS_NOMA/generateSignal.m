function signal = generateSignal(cfg)
validateattributes(cfg.M,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.G,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.numSymbols,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.oversampling,{'numeric'},{'scalar','integer','>=',2});
validateattributes(cfg.dfC,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.communicationPower,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.CPRdB,{'numeric'},{'scalar','real','finite'});
root = fileparts(mfilename('fullpath'));
data = load(fullfile(root,'L1CPcodes.mat'),'L1CPcodes');
assert(numel(cfg.PRN)==cfg.M && all(cfg.PRN>=1 & cfg.PRN<=size(data.L1CPcodes,1)) ...
    && all(cfg.PRN==fix(cfg.PRN)),'Invalid code selection.');

N = (cfg.M+1)*cfg.G-1;
NFFT = (N+1)*cfg.oversampling;
fs = NFFT*cfg.dfC;
Tc = 1/cfg.dfC;
Tp = Tc/cfg.G;
Pc = cfg.communicationPower/N;
Pp = 2*cfg.G*Pc/10^(cfg.CPRdB/10);   % CPR = 2G Pc/Pp
code = data.L1CPcodes(cfg.PRN,:);

% QPSK 通信符号按 Tc 更新。
symbols = ((2*randi([0,1],N,cfg.numSymbols)-1) ...
    +1i*(2*randi([0,1],N,cfg.numSymbols)-1))/sqrt(2);
centerTone = floor((N+1)/2);
tones = (1:N)-centerTone;
bins = mod(tones,NFFT)+1;
grid = complex(zeros(NFFT,cfg.numSymbols));
grid(bins,:) = sqrt(Pc)*symbols;
communication = reshape(ifft(grid,[],1)*NFFT,1,[]);

% 定位码片按 Tp 更新，并调制到间隔为 G*dfC 的子载波。
numChips = cfg.G*cfg.numSymbols;
chips = code(:,mod(0:numChips-1,size(code,2))+1);
samplesPerChip = NFFT/cfg.G;
numSamples = NFFT*cfg.numSymbols;
t = (0:numSamples-1)/fs;
positioning = complex(zeros(1,numSamples));
fp = ((1:cfg.M)*cfg.G-centerTone)*cfg.dfC;
for m = 1:cfg.M
    pulse = repelem(chips(m,:),samplesPerChip);
    positioning = positioning+sqrt(Pp)*pulse.*exp(1i*2*pi*fp(m)*t);
end
signal = struct('communication',communication,'positioning',positioning, ...
    'combined',communication+positioning,'sampleRate',fs,'N',N,'NFFT',NFFT, ...
    'Tc',Tc,'Tp',Tp,'Pc',Pc,'Pp',Pp*ones(1,cfg.M), ...
    'symbols',symbols,'chips',chips,'code',code,'positionFrequencies',fp);
end
