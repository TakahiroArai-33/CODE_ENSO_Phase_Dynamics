function tf = isodd(x)
%ISODD Logical mask for odd integers in a real numeric array.
% Project implementation; distributed under the root BSD 3-Clause LICENSE.
% Both positive and negative odd integers are recognized.

narginchk(1, 1);
validateattributes(x, {'numeric'}, {'real'}, mfilename, 'x', 1);
tf = abs(rem(x, 2)) == 1;
end
