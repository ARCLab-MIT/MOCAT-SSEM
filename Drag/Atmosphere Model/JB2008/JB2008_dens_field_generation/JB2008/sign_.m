% File retrieved on 11 June 2019 from:
% https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model
%
% Copyright (c) 2018, Meysam Mahooti
% See license in LICENSE_JB2008
%
%--------------------------------------------------------------------------
% sign: returns absolute value of a with sign of b
%--------------------------------------------------------------------------
function result = sign_(a, b)

if (b>=0)
    result = abs(a);
else
    result = -abs(a);
end
