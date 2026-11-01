function psd = communicationPsd(frequencyHz,num)
%COMMUNICATIONPSD Exact finite equal-power rectangular C-carrier PSD.
% Positive-argument trigamma tails avoid an N-by-frequency dense matrix.
x=frequencyHz/num.deltaFcHz+(num.N+1)/2;
sine=sin(pi*(x-round(x))).^2/pi^2;
value=zeros(size(x)); inside=x>=1 & x<=num.N;
value(inside)=1-sine(inside).*(psi(1,x(inside))+psi(1,num.N+1-x(inside)));
left=x<1; z=1-x(left);
value(left)=msnoma.model.sincpi(z).^2 ...
    +sine(left).*(psi(1,z+1)-psi(1,num.N+z));
right=x>num.N; z=x(right)-num.N;
value(right)=msnoma.model.sincpi(z).^2 ...
    +sine(right).*(psi(1,z+1)-psi(1,num.N+z));
assert(all(value>=-1e-12,'all'));
psd=num.PcW*num.TcS*max(value,0);
end
