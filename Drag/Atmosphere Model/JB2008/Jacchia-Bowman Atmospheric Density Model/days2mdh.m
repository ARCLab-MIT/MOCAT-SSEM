%---------------------------------------------------------------------------
%
%  days2mdh.m
%
%  this function converts the day of the year, days, to the equivalent month
%  day, hour, minute and second.
%
%
%  inputs:
%    year        - year                           1900 .. 2100
%    days        - julian day of the year         0.0  .. 366.0
%
%  outputs:
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    minute      - minute                         0 .. 59
%    sec         - second                         0.0 .. 59.999
%
%---------------------------------------------------------------------------
function [mon, day, hr, minute, sec] = days2mdh (year, days)

% --------------- set up array of days in month  --------------
lmonth = zeros(12,1);
for i=1:12
    lmonth(i) = 31;
    if (i == 2)
        lmonth(i)= 28;
    end
    if (i == 4 || i == 6 || i == 9 || i == 11)
        lmonth(i)= 30;
    end
end

dayofyr= floor(days );

% ----------------- find month and day of month ---------------
if rem(year-1900,4) == 0
    lmonth(2)= 29;
end

i= 1;
inttemp= 0;
while (dayofyr > inttemp + lmonth(i) && i < 12)
    inttemp= inttemp + lmonth(i);
    i= i+1;
end

mon= i;
day= dayofyr - inttemp;

% ----------------- find hours minutes and seconds ------------
temp= (days - dayofyr )*24;
hr  = fix( temp );
temp= (temp-hr) * 60;
minute = fix( temp );
sec = (temp-minute) * 60;

