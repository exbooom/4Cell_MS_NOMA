function [ber,interference] = calBER(cfg,EbN0dB,positionPower)
N = (cfg.M+1)*cfg.G-1;
Tc = 1/cfg.dfC;
Tp = Tc/cfg.G;
Pc = cfg.communicationPower/N;
assert(numel(positionPower)==cfg.M && all(isfinite(positionPower)) ...
    && all(positionPower>=0),'Invalid positioning power.');

% sinc 平方谱叠加，计算各通信子载波上的等效干扰。
x = (1:cfg.M)-(1:N)'/cfg.G;
overlap = ones(size(x));
index = x~=0;
overlap(index) = sin(pi*x(index))./(pi*x(index));
interference = Tp*overlap.^2*positionPower(:);

% QPSK 按每比特能量计算 BER，定位干扰作高斯近似。
Eb = Pc*Tc/2;
N0 = Eb./10.^(EbN0dB(:)/10);
ber = 0.5*erfc(sqrt(Eb./(N0+interference.')));
end
