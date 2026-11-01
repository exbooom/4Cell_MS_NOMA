function value = communicationPsd(localFrequencyHz,userId,ch,num,chunkSize)
%COMMUNICATIONPSD Exact finite-C-carrier PSD per unit per-carrier C power.
% Local frequency is relative to userId*deltaFpHz; gCP indexes P receivers.
if nargin<5, chunkSize=256; end
validateattributes(localFrequencyHz,{'numeric'},{'real','finite'});
validateattributes(userId,{'numeric'},{'scalar','integer','>=',1,'<=',num.M});
validateattributes(chunkSize,{'numeric'},{'scalar','integer','positive'});
gain=reshape(sum(ch.gCP(userId,:,:),2),num.N,1);
assert(all(isfinite(gain) & gain>=0));
frequency=localFrequencyHz(:); output=zeros(size(frequency));
for first=1:chunkSize:numel(frequency)
    ids=first:min(first+chunkSize-1,numel(frequency));
    x=(frequency(ids)+userId*num.deltaFpHz)*num.TcS-(1:num.N);
    output(ids)=num.TcS*(msnoma.model.sincpi(x).^2*gain);
end
value=reshape(output,size(localFrequencyHz));
end
