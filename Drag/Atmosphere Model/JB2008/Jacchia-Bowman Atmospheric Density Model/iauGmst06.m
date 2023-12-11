%  - - - - - - - - - -
%   i a u G m s t 0 6
%  - - - - - - - - - -
%
%  Greenwich mean sidereal time (consistent with IAU 2006 precession).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     uta,utb        UT1 as a 2-part Julian Date (Notes 1,2)
%     tta,ttb        TT as a 2-part Julian Date (Notes 1,2)
%
%  Returned (function value):
%                    Greenwich mean sidereal time (radians)
%
%  Notes:
%  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
%     Julian Dates, apportioned in any convenient way between the
%     argument pairs.  For example, JD=2450123.7 could be expressed in
%     any of these ways, among others:
%
%            Part A        Part B
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable (in the case of UT;  the TT is not at all critical
%     in this respect).  The J2000 and MJD methods are good compromises
%     between resolution and convenience.  For UT, the date & time
%     method is best matched to the algorithm that is used by the Earth
%     rotation angle function, called internally:  maximum precision is
%     delivered when the uta argument is for 0hrs UT1 on the day in
%     question and the utb argument lies in the range 0 to 1, or vice
%     versa.
%
%  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
%     and TT to predict the effects of precession.  If UT1 is used for
%     both purposes, errors of order 100 microarcseconds result.
%
%  3) This GMST is compatible with the IAU 2006 precession and must not
%     be used with other precession models.
%
%  4) The result is returned in the range 0 to 2pi.
%
%  Called:
%     iauEra00     Earth rotation angle, IAU 2000
%     iauAnp       normalize angle into range 0 to 2pi
%
%  Reference:
%     Capitaine, N., Wallace, P.T. & Chapront, J., 2005,
%     Astron.Astrophys. 432, 355
%
%  This revision:  2008 May 24
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gmst = iauGmst06(uta, utb, tta, ttb)

global const

% TT Julian centuries since J2000.0.
t = ((tta - const.DJ00) + ttb) / const.DJC;

% Greenwich mean sidereal time, IAU 2006.
gmst = iauAnp(iauEra00(uta, utb) +...
             (    0.014506     +...
             (  4612.156534    +...
             (     1.3915817   +...
             (    -0.00000044  +...
             (    -0.000029956 +...
             (    -0.0000000368 )...
             * t) * t) * t) * t) * t) * const.DAS2R);

