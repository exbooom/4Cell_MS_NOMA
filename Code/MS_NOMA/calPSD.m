function [frequencyHz,psd]=calPSD(x,fs,nfft)
% 分段加窗后对功率谱取平均。
numBlocks=floor(numel(x)/nfft);
assert(numBlocks>=1,'Signal is too short for spectrum analysis.');
blocks=reshape(x(1:numBlocks*nfft),nfft,numBlocks);
window=0.5-0.5*cos(2*pi*(0:nfft-1)'/(nfft-1));
spectrum=fft(blocks.*window,nfft,1);
psd=fftshift(mean(abs(spectrum).^2,2)/(fs*sum(window.^2)));
frequencyHz=(-floor(nfft/2):ceil(nfft/2)-1)'*fs/nfft;
end
