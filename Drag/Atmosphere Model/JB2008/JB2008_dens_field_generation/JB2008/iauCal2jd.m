% File retrieved on 11 June 2019 from:
% https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model
%
%  - - - - - - - - - -
%   i a u C a l 2 j d
%  - - - - - - - - - -
%
%  Gregorian Calendar to Julian Date.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     year, month, day in Gregorian calendar (Note 1)
%
%  Returned:
%     djm0      MJD zero-point: always 2400000.5
%     djm       Modified Julian Date for 0 hrs
%
%  Notes:
%  1) The algorithm used is valid from -4800 March 1, but this
%     implementation rejects dates before -4799 January 1.
%
%  2) The Julian Date is returned in two pieces, in the usual SOFA
%     manner, which is designed to preserve time resolution.  The
%     Julian Date is available as a single number by adding djm0 and
%     djm.
%
%  3) In early eras the conversion is from the "Proleptic Gregorian
%     Calendar";  no account is taken of the date(s) of adoption of
%     the Gregorian Calendar, nor is the AD/BC numbering convention
%     observed.
%
%  Reference:
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992),
%     Section 12.92 (p604).
%
%  This revision:  2009 October 19
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [djm0, djm] = iauCal2jd(year, month, day)

djm0 = 2400000.5;

b = 0;
c = 0;

if (month <= 2)
   year = year - 1;
   month = month + 12;
end

if (year < 0)
   c = -.75;
end

% check for valid calendar date
if (year < 1582)
   % null
elseif (year > 1582)
   a = fix(year / 100);
   b = 2 - a + floor(a / 4);
elseif (month < 10)
   % null
elseif (month > 10)
   a = fix(year / 100);
   b = 2 - a + floor(a / 4);
elseif (day <= 4)
   % null
elseif (day > 14)
   a = fix(year / 100);
   b = 2 - a + floor(a / 4);
else
    fprintf('\n\n  this is an invalid calendar date!!\n');
    return
end

jd = fix(365.25 * year + c) + fix(30.6001 * (month + 1));
jd = jd + day + b + 1720994.5;
djm = jd-djm0;

