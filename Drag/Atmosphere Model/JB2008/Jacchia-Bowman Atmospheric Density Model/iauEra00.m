%  - - - - - - - - -
%   i a u E r a 0 0
%  - - - - - - - - -
%
%  Earth rotation angle (IAU 2000 model).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     dj1,dj2       UT1 as a 2-part Julian Date (see note)
%
%  Returned (function value):
%                   Earth rotation angle (radians), range 0-2pi
%
%  Notes:
%  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
%     convenient way between the arguments dj1 and dj2.  For example,
%     JD(UT1)=2450123.7 could be expressed in any of these ways,
%     among others:
%
%             dj1            dj2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable.  The J2000 and MJD methods are good compromises
%     between resolution and convenience.  The date & time method is
%     best matched to the algorithm used:  maximum precision is
%     delivered when the dj1 argument is for 0hrs UT1 on the day in
%     question and the dj2 argument lies in the range 0 to 1, or vice
%     versa.
%
%  2) The algorithm is adapted from Expression 22 of Capitaine et al.
%     2000.  The time argument has been expressed in days directly,
%     and, to retain precision, integer contributions have been
%     eliminated.  The same formulation is given in IERS Conventions
%     (2003), Chap. 5, Eq. 14.
%
%  Called:
%     iauAnp       normalize angle into range 0 to 2pi
%
%  References:
%     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
%     Astrophys., 355, 398-405.
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2008 May 24
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = iauEra00(dj1, dj2)

global const

% Days since fundamental epoch.
if (dj1 < dj2)
    d1 = dj1;
    d2 = dj2;
else
    d1 = dj2;
    d2 = dj1;
end

t = d1 + (d2- const.DJ00);

% Fractional part of T (days).
f = mod(d1, 1) + mod(d2, 1);

% Earth rotation angle at this UT1.
theta = iauAnp(const.D2PI * (f + 0.7790572732640 + 0.00273781191135448 * t));

