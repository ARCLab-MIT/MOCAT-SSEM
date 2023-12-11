% File retrieved on 11 June 2019 from:
% https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model
%
%   - - - - - - -
%    i a u A n p
%   - - - - - - -
% 
%   Normalize angle into the range 0 <= a < 2pi.
% 
%   This function is part of the International Astronomical Union's
%   SOFA (Standards Of Fundamental Astronomy) software collection.
% 
%   Status:  vector/matrix support function.
% 
%   Given:
%      a        double     angle (radians)
% 
%   Returned (function value):
%               double     angle in range 0-2pi
% 
%   This revision:  2008 May 16
% 
%   SOFA release 2012-03-01
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = iauAnp(a)

global const

w = mod(a, const.D2PI);
if (w < 0)
    w = w + const.D2PI;
end

