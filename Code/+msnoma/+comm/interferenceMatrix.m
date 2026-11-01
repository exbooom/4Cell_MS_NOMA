function [A,rowMap] = interferenceMatrix(ch,num,kernels)
%INTERFERENCEMATRIX Eq3-4, full cross-BS/cross-subcarrier interference.
weighted = ch.gPC.*reshape(kernels.weightedOverlapS,1,num.N,1,num.M);
A = reshape(weighted,num.K*num.N,num.K*num.M);
[k,n] = ndgrid(1:num.K,1:num.N);
rowMap = [k(:),n(:)];
end
