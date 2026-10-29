function tf = isint(x)
%ISINT Logical mask for finite integers, including complex integers.
% Project implementation; distributed under the root BSD 3-Clause LICENSE.
% A complex entry is an integer when both components are integers.

narginchk(1, 1);
validateattributes(x, {'numeric'}, {}, mfilename, 'x', 1);
tf = ~isnan(x) & ~isinf(x) & (x == fix(x));
end
