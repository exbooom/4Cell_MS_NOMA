function num = numerology(cfg)
%NUMEROLOGY Single source of carrier counts, bandwidth and SI units.
num.K = cfg.paper.K;
num.M = cfg.paper.M;
num.G = round(cfg.paper.requestedBandwidthHz/((num.M+1)*cfg.paper.deltaFcHz));
assert(num.G>=1,'Requested bandwidth is too small.');
num.N = (num.M+1)*num.G-1;
num.deltaFcHz = cfg.paper.deltaFcHz;
num.deltaFpHz = num.G*num.deltaFcHz;
num.TcS = 1/num.deltaFcHz;
num.TpS = 1/num.deltaFpHz;
num.requestedBandwidthHz = cfg.paper.requestedBandwidthHz;
num.effectiveBandwidthHz = (num.N+1)*num.deltaFcHz;
num.BHz = num.effectiveBandwidthHz;
num.BfeHz = cfg.paper.frontendMultiplier*num.BHz;
num.PcW = cfg.paper.communicationTotalPowerW/num.N;
num.c0Mps = cfg.assumed.propagation.speedMps;
num.bandwidthMode = cfg.models.bandwidth;
num.powerVectorOrder = 'P(:), variableId=k+(m-1)*K';
num.communicationRowOrder = 'I(:), rowId=k+(n-1)*K';
end
