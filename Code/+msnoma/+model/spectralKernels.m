function kernels = spectralKernels(num)
%SPECTRALKERNELS Position-to-communication overlap, N-by-M.
x = (1:num.M)-(1:num.N)'/num.G;
kernels.W = msnoma.model.sincpi(x).^2;
kernels.weightedOverlapS = num.TpS*kernels.W;
kernels.model = 'paper_Eq4_normalized_sinc_squared';
end
