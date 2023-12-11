%--------------------------------------------------------------------------
%
% Mjday: Modified Julian Date from calendar date and time
%
% Inputs:
%   year, month, day, hour, min, sec
%
% output:
%   Mjd		Modified Julian Date
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function Mjd = Mjday(year, month, day, hour, min, sec)

if (nargin < 4)
    hour = 0;
    min  = 0;
    sec  = 0;
end

y = year;
m = month;
b = 0;
c = 0;

if (m <= 2)
   y = y - 1;
   m = m + 12;
end

if (y < 0)
   c = -.75;
end

% check for valid calendar date
if (year < 1582)
   % null
elseif (year > 1582)
   a = fix(y / 100);
   b = 2 - a + floor(a / 4);
elseif (month < 10)
   % null
elseif (month > 10)
   a = fix(y / 100);
   b = 2 - a + floor(a / 4);
elseif (day <= 4)
   % null
elseif (day > 14)
   a = fix(y / 100);
   b = 2 - a + floor(a / 4);
else
    fprintf('\n\n  this is an invalid calendar date!!\n');
    return
end

jd = fix(365.25 * y + c) + fix(30.6001 * (m + 1));
jd = jd + day + b + 1720994.5;
jd = jd + (hour+min/60+sec/3600)/24;
Mjd = jd - 2400000.5;

