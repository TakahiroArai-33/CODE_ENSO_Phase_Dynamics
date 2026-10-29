function tf = iseven(x)
%ISEVEN Logical mask for even integers in a real numeric array.
% Project implementation; distributed under the root BSD 3-Clause LICENSE.
% Integer arithmetic preserves parity for int64 and uint64 inputs.

narginchk(1, 1);
validateattributes(x, {'numeric'}, {'real'}, mfilename, 'x', 1);
tf = rem(x, 2) == 0;
end
