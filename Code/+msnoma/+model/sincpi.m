function y = sincpi(x)
%SINCPI Normalized sinc, including its removable singularity.
y = ones(size(x));
i = x~=0;
y(i) = sin(pi*x(i))./(pi*x(i));
end
