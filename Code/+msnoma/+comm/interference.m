function IWPerHz = interference(powerW,A,num)
%INTERFERENCE Output is K-by-N PSD-like interference, not total watts.
assert(isequal(size(powerW),[num.K,num.M]) && all(powerW(:)>=0));
IWPerHz = reshape(A*powerW(:),num.K,num.N);
end
