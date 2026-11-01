function [x,F] = unconditionalCdf(values,valid)
%UNCONDITIONALCDF End point is coverage, failed users never become zero error.
assert(numel(values)==numel(valid));
valid = valid(:) & isfinite(values(:));
x = sort(values(valid));
F = (1:numel(x))'/numel(values);
end
