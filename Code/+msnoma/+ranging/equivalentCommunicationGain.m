function gCeq = equivalentCommunicationGain(ch,num)
%EQUIVALENTCOMMUNICATIONGAIN Eq62 finite-carrier weighted gain, M-by-K.
weights = reshape(sin(pi*(1:num.N)/num.G).^2,1,1,num.N);
gCeq = (2/num.N)*sum(ch.gCP.*weights,3);
end
