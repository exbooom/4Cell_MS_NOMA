cfg.M = 20;                        % 定位子载波数
cfg.G = 80;                        % 定位/通信子载波间隔比
cfg.dfC = 30e3;                    % 通信子载波间隔/Hz
cfg.communicationPower = 1;        % 通信总功率，归一化
cfg.CPRdB = 20;                    % CPR = 2G Pc/Pp，dB
cfg.PRN = 1:cfg.M;                 % 定位码编号
cfg.numSymbols = 128;              % 通信符号数
cfg.oversampling = 4;              % 采样倍数
cfg.EbN0dB = (0:2:20)';            % 通信 Eb/N0/dB
cfg.berCPRdB = [15,17,20,30,Inf];   % Inf：无定位信号

% 生成通信、定位和叠加信号。
signal = generateSignal(cfg);

% 比较通信符号与定位码片的时间尺度。
numShow = min(4,cfg.numSymbols);
fig = figure('Color','w','Position',[100,100,760,540]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
symbols = signal.symbols(1,[1:numShow,numShow]);
stairs((0:numShow)*signal.Tc*1e6,[real(symbols);imag(symbols)].','LineWidth',1.2);
grid on; ylim([-1.2,1.2]); xlim([0,numShow*signal.Tc*1e6]);
xlabel('Time (us)'); ylabel('QPSK symbol'); legend('I','Q');
nexttile;
numChips = numShow*cfg.G;
stairs((0:numChips)*signal.Tp*1e6, ...
    signal.chips(1,[1:numChips,numChips]),'LineWidth',1.1);
grid on; ylim([-1.2,1.2]); xlim([0,numShow*signal.Tc*1e6]);
xlabel('Time (us)'); ylabel('Positioning chip');
formatFigure(fig);

% 查看叠加前后的基带频谱。
nfft = 8*signal.NFFT;
[f,Pc] = calPSD(signal.communication,signal.sampleRate,nfft);
[~,Pp] = calPSD(signal.positioning,signal.sampleRate,nfft);
[~,Pt] = calPSD(signal.combined,signal.sampleRate,nfft);
fig = figure('Color','w','Position',[100,100,760,480]);
plot(f/1e6,10*log10(max([Pc,Pp,Pt],realmin)),'LineWidth',1.1);
grid on; xlabel('Baseband frequency (MHz)'); ylabel('PSD (dB/Hz)');
B = (signal.N+1)*cfg.dfC;
xlim([-B,B]/1e6);
legend('Communication','Positioning','Combined','Location','southwest');
formatFigure(fig);

% 原 L1C 码的相关特性。
[autoCorr,lags] = xcorr(signal.code(1,:),'coeff');
crossCorr = xcorr(signal.code(1,:),signal.code(min(2,cfg.M),:),'coeff');
fig = figure('Color','w','Position',[100,100,720,540]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(lags,abs(autoCorr));
grid on; xlabel('Code lag (chips)'); ylabel('|Autocorrelation|');
nexttile; plot(lags,abs(crossCorr));
grid on; xlabel('Code lag (chips)'); ylabel('|Cross-correlation|');
formatFigure(fig);

% 单小区下，定位分量对 QPSK 通信的等效干扰。
ber = zeros(numel(cfg.EbN0dB),numel(cfg.berCPRdB));
for p = 1:numel(cfg.berCPRdB)
    powerP = 2*cfg.G*signal.Pc/10^(cfg.berCPRdB(p)/10)*ones(1,cfg.M);
    ber(:,p) = mean(calBER(cfg,cfg.EbN0dB,powerP),2);
end
[berSubcarrier,interference] = calBER(cfg,5,signal.Pp);
clean = calBER(cfg,5,zeros(1,cfg.M));
fig = figure('Color','w','Position',[100,100,1100,480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
labels = cell(1,numel(cfg.berCPRdB));
for p = 1:numel(cfg.berCPRdB)
    semilogy(cfg.EbN0dB,ber(:,p),'LineWidth',1.2);
    if isinf(cfg.berCPRdB(p))
        labels{p} = 'No positioning';
    else
        labels{p} = sprintf('CPR = %g dB',cfg.berCPRdB(p));
    end
end
set(gca,'YScale','log'); ylim([1e-8,.2]); grid on;
xlabel('Communication E_b/N_0 (dB)'); ylabel('Average approximate BER');
legend(labels,'Location','southoutside','NumColumns',2);
nexttile;
plot(1:signal.N,berSubcarrier,'LineWidth',1.2); hold on;
plot(1:signal.N,clean,'--','LineWidth',1.2);
grid on; xlabel('C-Sub index'); ylabel('Approximate BER');
title(sprintf('E_b/N_0 = 5 dB, CPR = %g dB',cfg.CPRdB));
legend('With positioning','No positioning','Location','southoutside');
formatFigure(fig);
fprintf('Bandwidth: %.2f MHz; communication power: %.4f; positioning power: %.4f\n', ...
    B/1e6,mean(abs(signal.communication).^2),mean(abs(signal.positioning).^2));

function formatFigure(fig)
set(findall(fig,'Type','axes'),'Color','w','XColor','k','YColor','k', ...
    'GridColor',[.6,.6,.6],'FontSize',10);
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','legend'),'Color','w','TextColor','k');
end
