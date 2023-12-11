%-Abstract
%
%   CSPICE_GFPOSC determines the time intervals for which a coordinate of
%   an observer-target position vector satisfies a numerical constraint.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      target   name of the target body.
%
%               [1,c1] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               Optionally, you may supply the integer ID code for the object
%               as an integer string. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               target body.
%
%               The target and observer define a position vector
%               that points from the observer to the target.
%
%      frame    name of the reference frame in which to perform state look-ups
%               and coordinate calculations.
%
%               [1,c2] = size(frame); char = class(frame)
%
%                  or
%
%               [1,1] = size(frame); cell = class(frame)
%
%               The SPICE frame subsystem must recognize the `frame' name.
%
%      abcorr   describes the aberration corrections to apply to the state
%               evaluations to account for one-way light time and stellar
%               aberration.
%
%               [1,c3] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               This routine accepts the same aberration corrections as does
%               the routine cspice_spkezr. See the header of cspice_spkezr for
%               a detailed description of the aberration correction options.
%               For convenience, the options are listed below:
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       "Reception" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     "Reception" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'CN'       "Reception" case: converged
%                             Newtonian light time correction.
%
%                  'CN+S'     "Reception" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%                  'XLT'      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'XCN'      "Transmission" case: converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    "Transmission" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               The `abcorr' string lacks sensitivity to case, and to
%               embedded, leading and trailing blanks.
%
%      obsrvr   name of the observing body.
%
%               [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Optionally, you may supply the ID code of the object as
%               an integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the
%               observer is the Earth.
%
%      crdsys   name of the coordinate system for which the coordinate of
%               interest is a member.
%
%               [1,c5] = size(crdsys); char = class(crdsys)
%
%                  or
%
%               [1,1] = size(crdsys); cell = class(crdsys)
%
%      coord    name of the coordinate of interest in `crdsys'.
%
%               [1,c6] = size(coord); char = class(coord)
%
%                  or
%
%               [1,1] = size(coord); cell = class(coord)
%
%               The supported coordinate systems and coordinate names are:
%
%                   crdsys              coord                Range
%                  ------------------  -------------------  -----------
%                   'RECTANGULAR'       'X'
%                                       'Y'
%                                       'Z'
%
%                   'LATITUDINAL'       'RADIUS'
%                                       'LONGITUDE'        (-Pi,Pi]
%                                       'LATITUDE'         [-Pi/2,Pi/2]
%
%                   'RA/DEC'            'RANGE'
%                                       'RIGHT ASCENSION'  [0,2Pi)
%                                       'DECLINATION'      [-Pi/2,Pi/2]
%
%                   'SPHERICAL'         'RADIUS'
%                                       'COLATITUDE'       [0,Pi]
%                                       'LONGITUDE'        (-Pi,Pi]
%
%                   'CYLINDRICAL'       'RADIUS'
%                                       'LONGITUDE'        [0,2Pi)
%                                       'Z'
%
%                   'GEODETIC'          'LONGITUDE'        (-Pi,Pi]
%                                       'LATITUDE'         [-Pi/2,Pi/2]
%                                       'ALTITUDE'
%
%                   'PLANETOGRAPHIC'    'LONGITUDE'        [0,2Pi)
%                                       'LATITUDE'         [-Pi/2,Pi/2]
%                                       'ALTITUDE'
%
%               Limit searches for coordinate events in the GEODETIC and
%               PLANETOGRAPHIC coordinate systems to `target' bodies with
%               axial symmetry in the equatorial plane, i.e. equality
%               of the body X and Y radii (oblate or prolate spheroids).
%
%      relate   describes the constraint relational on the selected coordinate
%               of the observer-target vector.
%
%               [1,c7] = size(relate); char = class(relate)
%
%                  or
%
%               [1,1] = size(relate); cell = class(relate)
%
%               The result window found by this routine indicates the time
%               intervals where the constraint is satisfied.
%
%               Supported values of `relate' and corresponding meanings are
%               shown below:
%
%                  '>'      Separation is greater than the reference
%                           value refval.
%
%                  '='      Separation is equal to the reference
%                           value refval.
%
%                  '<'      Separation is less than the reference
%                           value refval.
%
%                  'ABSMAX'  Separation is at an absolute maximum.
%
%                  'ABSMIN'  Separation is at an absolute  minimum.
%
%                  'LOCMAX'  Separation is at a local maximum.
%
%                  'LOCMIN'  Separation is at a local minimum.
%
%               The caller may indicate that the region of interest
%               is the set of time intervals where the quantity is
%               within a specified measure of an absolute extremum.
%               The argument `adjust' (described below) is used to
%               specify this measure.
%
%               Local extrema are considered to exist only in the
%               interiors of the intervals comprising the confinement
%               window: a local extremum cannot exist at a boundary
%               point of the confinement window.
%
%               The `relate' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%      refval   reference value used together with `relate' argument to
%               define an equality or inequality to satisfy by the selected
%               coordinate of the observer-target vector.
%
%               [1,1] = size(refval); double = class(refval)
%
%               See the discussion of `relate' above for further information.
%
%               The units of `refval' correspond to the type as defined
%               by `coord', radians for angular measures, kilometers for
%               distance measures.
%
%      adjust   value used to modify searches for absolute extrema.
%
%               [1,1] = size(adjust); double = class(adjust)
%
%               When `relate' is set to 'ABSMAX' or 'ABSMIN' and `adjust'
%               is set to a positive value, cspice_gfposc finds times
%               when the observer-target vector coordinate is within `adjust'
%               radians/kilometers of the specified extreme value.
%
%               For `relate' set to 'ABSMAX', the result window contains
%               time intervals when the observer-target vector coordinate has
%               values between ABSMAX - adjust and ABSMAX.
%
%               For `relate' set to `ABSMIN', the result window contains
%               time intervals when the observer-target vector coordinate has
%               values between ABSMIN and ABSMIN + adjust.
%
%               `adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     time step size to use in the search.
%
%               [1,1] = size(step); double = class(step)
%
%               `step' must be short enough to for a search using this step
%               size to locate the time intervals where coordinate
%               function of the observer-target vector is monotone increasing
%               or decreasing. However, `step' must not be *too* short, or
%               the search will take an unreasonable amount of time.
%
%               For coordinates other than LONGITUDE and RIGHT ASCENSION,
%               the step size must be shorter than the shortest interval,
%               within the confinement window, over which the coordinate
%               is monotone increasing or decreasing.
%
%               For LONGITUDE and RIGHT ASCENSION, the step size must
%               be shorter than the shortest interval, within the
%               confinement window, over which either the sin or cos
%               of the coordinate is monotone increasing or decreasing.
%
%               The choice of `step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               `step' has units of TDB seconds.
%
%      nintvls  value specifying the number of intervals in the internal
%               workspace array used by this routine.
%
%               [1,1] = size(nintvls); int32 = class(nintvls)
%
%               `nintvls' should be at least as large as the number of
%               intervals within the search region on which the specified
%               observer-target vector coordinate function is monotone
%               increasing or decreasing. It does no harm to pick a value of
%               `nintvls' larger than the minimum required to execute the
%               specified search, but if chosen too small, the search will
%               fail.
%
%      cnfine   a SPICE window that confines the time period over which the
%               specified search is conducted.
%
%               [2m,1] = size(cnfine); double = class(cnfine)
%
%               `cnfine' may consist of a single interval or a collection
%               of intervals.
%
%               In some cases the confinement window can be used to
%               greatly reduce the time period that must be searched
%               for the desired solution. See the -Particulars section
%               below for further discussion.
%
%               See the -Examples section below for a code example
%               that shows how to create a confinement window.
%
%               In some cases the observer's state may be computed at
%               times outside of `cnfine' by as much as 2 seconds. See
%               -Particulars for details.
%
%   the call:
%
%      [result] = cspice_gfposc( target, frame, abcorr,  obsrvr,           ...
%                                crdsys, coord, relate,  refval,           ...
%                                adjust, step,  nintvls, cnfine )
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window `cnfine', on which the specified
%               constraint is satisfied.
%
%               [2n,1] = size(result); double = class(result)
%
%               If the search is for local extrema, or for absolute
%               extrema with `adjust' set to zero, then normally each
%               interval of `result' will be a singleton: the left and
%               right endpoints of each interval will be identical.
%
%               If no times within the confinement window satisfy the
%               constraint, `result' will return with cardinality zero.
%
%-Parameters
%
%   All parameters described here are declared in the Mice include file
%   MiceGF.m. See that file for parameter values.
%
%   SPICE_GF_CNVTOL
%
%               is the convergence tolerance used for finding endpoints
%               of the intervals comprising the result window.
%               SPICE_GF_CNVTOL is used to determine when binary
%               searches for roots should terminate: when a root is
%               bracketed within an interval of length SPICE_GF_CNVTOL,
%               the root is considered to have been found.
%
%               The accuracy, as opposed to precision, of roots found
%               by this routine depends on the accuracy of the input
%               data. In most cases, the accuracy of solutions will be
%               inferior to their precision.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%   1) Find the time during 2007 for which the latitude of the
%      Earth-Sun vector in IAU_EARTH frame has the maximum value,
%      i.e. the latitude of the Tropic of Cancer.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfposc_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00009.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function gfposc_ex1()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfposc_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % The latitude varies relatively slowly, ~46 degrees during the
%         % year. The extrema occur approximately every six months.
%         % Search using a step size near half that value (180 days).
%         % For this example use ninety days (in units of seconds).
%         %
%         step   = 90.*cspice_spd;
%         adjust = 0.;
%         refval = 0.;
%
%         %
%         % Search for the date on which the `crdsys' system
%         % coordinate `coord' satisfies the `relate' constraint.
%         %
%         %
%         % Perform this search using the geometric position
%         % of the bodies; set the aberration correction to 'NONE'.
%         %
%         targ   = 'SUN';
%         obsrvr = 'EARTH';
%         frame  = 'IAU_EARTH';
%         abcorr = 'NONE';
%         relate = 'ABSMAX';
%         crdsys = 'LATITUDINAL';
%         coord  = 'LATITUDE';
%         nintvls= MAXWIN ;
%
%
%         result = cspice_gfposc( targ,   frame, abcorr,  obsrvr,          ...
%                                 crdsys, coord, relate,  refval,          ...
%                                 adjust, step,  nintvls, cnfine );
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         for i=1:numel(result)/2
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['From : ' output(1,:)] )
%               disp( ['To   : ' output(2,:)] )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Event time: 2007-JUN-21 17:54:13.172475 (TDB)
%
%
%   2) A minor modification of the program listed in Example 1; find
%      the time during 2007 for which the latitude of the Earth-Sun
%      vector in IAU_EARTH frame has the minimum value, i.e. the
%      latitude of the Tropic of Capricorn.
%
%
%      Use the meta-kernel from the first example.
%
%
%      Example code begins here.
%
%
%      function gfposc_ex2()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfposc_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % The latitude varies relatively slowly, ~46 degrees during the
%         % year. The extrema occur approximately every six months.
%         % Search using a step size near half that value (180 days).
%         % For this example use ninety days (in units of seconds).
%         %
%         step   = 90.*cspice_spd;
%         adjust = 0.;
%         refval = 0.;
%
%         %
%         % Search for the date on which the `crdsys' system
%         % coordinate `coord' satisfies the `relate' constraint.
%         %
%         %
%         % Perform this search using the geometric position
%         % of the bodies; set the aberration correction to 'NONE'.
%         %
%         targ   = 'SUN';
%         obsrvr = 'EARTH';
%         frame  = 'IAU_EARTH';
%         abcorr = 'NONE';
%         relate = 'ABSMIN';
%         crdsys = 'LATITUDINAL';
%         coord  = 'LATITUDE';
%         nintvls= MAXWIN ;
%
%
%         result = cspice_gfposc( targ,   frame, abcorr,  obsrvr,          ...
%                                 crdsys, coord, relate,  refval,          ...
%                                 adjust, step,  nintvls, cnfine );
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         for i=1:numel(result)/2
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['From : ' output(1,:)] )
%               disp( ['To   : ' output(2,:)] )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Event time: 2007-DEC-22 06:04:32.635539 (TDB)
%
%
%   3) Find the time during 2007 for which the Z component of the
%      Earth-Sun vector in IAU_EARTH frame has value 0, i.e. crosses
%      the equatorial plane (this also defines a zero latitude).
%      The search should return two times, one for an ascending
%      passage and one for descending.
%
%      Use the meta-kernel from the first example.
%
%
%      Example code begins here.
%
%
%      function gfposc_ex3()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfposc_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % The latitude varies relatively slowly, ~46 degrees during the
%         % year. The extrema occur approximately every six months.
%         % Search using a step size near half that value (180 days).
%         % For this example use ninety days (in units of seconds).
%         %
%         step   = 90.*cspice_spd;
%         adjust = 0.;
%         refval = 0.;
%
%         %
%         % Search for the date on which the `crdsys' system
%         % coordinate `coord' satisfies the `relate' constraint.
%         %
%         %
%         % Perform this search using the geometric position
%         % of the bodies; set the aberration correction to 'NONE'.
%         %
%         targ   = 'SUN';
%         obsrvr = 'EARTH';
%         frame  = 'IAU_EARTH';
%         abcorr = 'NONE';
%         relate = '=';
%         crdsys = 'RECTANGULAR';
%         coord  = 'Z';
%         nintvls= MAXWIN ;
%
%
%         result = cspice_gfposc( targ,   frame, abcorr,  obsrvr,          ...
%                                 crdsys, coord, relate,  refval,          ...
%                                 adjust, step,  nintvls, cnfine );
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         for i=1:numel(result)/2
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['From : ' output(1,:)] )
%               disp( ['To   : ' output(2,:)] )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Event time: 2007-MAR-21 00:01:25.500673 (TDB)
%      Event time: 2007-SEP-23 09:46:39.579484 (TDB)
%
%
%   4) Find the times between Jan 1, 2007 and Jan 1, 2008
%      corresponding to the apoapsis on the Moon's orbit around the
%      Earth (note, the cspice_gfdist routine can also perform this search).
%
%      Use the meta-kernel from the first example.
%
%
%      Example code begins here.
%
%
%      function gfposc_ex4()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfposc_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % This search requires a change in the step size since the
%         % Moon's orbit about the earth (earth-moon barycenter) has
%         % a twenty-eight day period. Use a step size something
%         % less than half that value. In this case, we use twelve
%         % days.
%         %
%         step   = 12.*cspice_spd;
%         adjust = 0.;
%         refval = 0.;
%
%         %
%         % Search for the date on which the `crdsys' system
%         % coordinate `coord' satisfies the `relate' constraint.
%         %
%         %
%         % Perform this search using the geometric position
%         % of the bodies; set the aberration correction to 'NONE'.
%         %
%         targ   = 'MOON';
%         obsrvr = 'EARTH';
%         frame  = 'J2000';
%         abcorr = 'NONE';
%         relate = 'LOCMAX';
%         crdsys = 'SPHERICAL';
%         coord  = 'RADIUS';
%         nintvls= MAXWIN ;
%
%
%         result = cspice_gfposc( targ,   frame, abcorr,  obsrvr,          ...
%                                 crdsys, coord, relate,  refval,          ...
%                                 adjust, step,  nintvls, cnfine );
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         for i=1:numel(result)/2
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['From : ' output(1,:)] )
%               disp( ['To   : ' output(2,:)] )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Event time: 2007-JAN-10 16:26:18.784521 (TDB)
%      Event time: 2007-FEB-07 12:39:35.055710 (TDB)
%      Event time: 2007-MAR-07 03:38:07.308330 (TDB)
%      Event time: 2007-APR-03 08:38:55.191516 (TDB)
%      Event time: 2007-APR-30 10:56:49.819340 (TDB)
%      Event time: 2007-MAY-27 22:03:28.834302 (TDB)
%      Event time: 2007-JUN-24 14:26:23.617432 (TDB)
%      Event time: 2007-JUL-22 08:43:50.113902 (TDB)
%      Event time: 2007-AUG-19 03:28:33.515939 (TDB)
%      Event time: 2007-SEP-15 21:07:13.940711 (TDB)
%      Event time: 2007-OCT-13 09:52:30.791223 (TDB)
%      Event time: 2007-NOV-09 12:32:50.039258 (TDB)
%      Event time: 2007-DEC-06 16:54:31.199770 (TDB)
%
%
%   5) Find times between Jan 1, 2007 and Jan 1, 2008 when the
%      latitude (elevation) of the observer-target vector between
%      DSS 17 and the Moon, as observed in the DSS 17 topocentric
%      (station) frame, exceeds 83 degrees.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfposc_ex5.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%            earthstns_itrf93_050714.bsp   SPK for DSN Station
%                                          Locations
%            earth_topo_050714.tf          Topocentric DSN stations
%                                          frame definitions
%            earth_latest_high_prec.bpc    High precision earth PCK
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00009.tpc',
%                                'naif0009.tls',
%                                'earthstns_itrf93_050714.bsp',
%                                'earth_topo_050714.tf',
%                                'earth_latest_high_prec.bpc'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function gfposc_ex5()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfposc_ex5.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % This search uses a step size of four hours since the
%         % time for all declination zero-to-max-to-zero passes
%         % within the search window exceeds eight hours.
%         %
%         % The example uses an 83 degree elevation because of its
%         % rare occurrence and short duration.
%         %
%         step   = (4./24.)*cspice_spd;
%         adjust = 0.;
%         refval = 83*cspice_rpd;
%
%         %
%         % Search for the date on which the `crdsys' system
%         % coordinate `coord' satisfies the `relate' constraint.
%         %
%         %
%         % Perform this search using the geometric position
%         % of the bodies; set the aberration correction to 'NONE'.
%         %
%         targ   = 'MOON';
%         obsrvr = 'DSS-17';
%         frame  = 'DSS-17_TOPO';
%         abcorr = 'NONE';
%         relate = '>';
%         crdsys = 'LATITUDINAL';
%         coord  = 'LATITUDE';
%         nintvls= MAXWIN ;
%
%
%         result = cspice_gfposc( targ,   frame, abcorr,  obsrvr,          ...
%                                 crdsys, coord, relate,  refval,          ...
%                                 adjust, step,  nintvls, cnfine );
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         for i=1:numel(result)/2
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['From : ' output(1,:)] )
%               disp( ['To   : ' output(2,:)] )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      From : 2007-FEB-26 03:18:48.229281 (TDB)
%      To   : 2007-FEB-26 03:31:29.734931 (TDB)
%      From : 2007-MAR-25 01:12:38.550572 (TDB)
%      To   : 2007-MAR-25 01:23:53.909469 (TDB)
%
%
%-Particulars
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when the selected coordinate of
%   the observer-target position vector satisfies a caller-specified
%   constraint. The resulting set of intervals is returned as a SPICE
%   window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%
%
%   The Search Process
%   ==================
%
%   Regardless of the type of constraint selected by the caller, this
%   routine starts the search for solutions by determining the time
%   periods, within the confinement window, over which the specified
%   coordinate function is monotone increasing and monotone
%   decreasing. Each of these time periods is represented by a SPICE
%   window. Having found these windows, all of the coordinate
%   function's local extrema within the confinement window are known.
%   Absolute extrema then can be found very easily.
%
%   Within any interval of these "monotone" windows, there will be at
%   most one solution of any equality constraint. Since the boundary
%   of the solution set for any inequality constraint is contained in
%   the union of
%
%   -  the set of points where an equality constraint is met
%
%   -  the boundary points of the confinement window
%
%   the solutions of both equality and inequality constraints can be
%   found easily once the monotone windows have been found.
%
%
%   Step Size
%   =========
%
%   The monotone windows (described above) are found using a two-step
%   search process. Each interval of the confinement window is
%   searched as follows: first, the input step size is used to
%   determine the time separation at which the sign of the rate of
%   change of coordinate will be sampled. Starting at
%   the left endpoint of an interval, samples will be taken at each
%   step. If a change of sign is found, a root has been bracketed; at
%   that point, the time at which the time derivative of the
%   coordinate is zero can be found by a refinement process, for
%   example, using a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the coordinate function is monotone:
%   the step size should be shorter than the shortest of these
%   intervals (within the confinement window).
%
%   The optimal step size is *not* necessarily related to the lengths
%   of the intervals comprising the result window. For example, if
%   the shortest monotone interval has length 10 days, and if the
%   shortest result window interval has length 5 minutes, a step size
%   of 9.9 days is still adequate to find all of the intervals in the
%   result window. In situations like this, the technique of using
%   monotone windows yields a dramatic efficiency improvement over a
%   state-based search that simply tests at each step whether the
%   specified constraint is satisfied. The latter type of search can
%   miss solution intervals if the step size is longer than the
%   shortest solution interval.
%
%   Having some knowledge of the relative geometry of the target and
%   observer can be a valuable aid in picking a reasonable step size.
%   In general, the user can compensate for lack of such knowledge by
%   picking a very short step size; the cost is increased computation
%   time.
%
%   Note that the step size is not related to the precision with which
%   the endpoints of the intervals of the result window are computed.
%   That precision level is controlled by the convergence tolerance.
%
%
%   Convergence Tolerance
%   =====================
%
%   As described above, the root-finding process used by this routine
%   involves first bracketing roots and then using a search process
%   to locate them. "Roots" are both times when local extrema are
%   attained and times when the coordinate function is equal to a
%   reference value. All endpoints of the intervals comprising the
%   result window are either endpoints of intervals of the
%   confinement window or roots.
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   "convergence tolerance." The default convergence tolerance
%   used by this routine is set by the parameter SPICE_GF_CNVTOL (defined
%   in MiceGF.m).
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't become the limiting factor in the accuracy of
%   solutions found by this routine. In general the accuracy of input
%   data will be the limiting factor.
%
%   The user may change the convergence tolerance from the default
%   SPICE_GF_CNVTOL value by calling the routine cspice_gfstol, e.g.
%
%      cspice_gfstol( tolerance value );
%
%   Call cspice_gfstol prior to calling this routine. All subsequent
%   searches will use the updated tolerance value.
%
%   Setting the tolerance tighter than SPICE_GF_CNVTOL is unlikely to be
%   useful, since the results are unlikely to be more accurate.
%   Making the tolerance looser will speed up searches somewhat,
%   since a few convergence steps will be omitted. However, in most
%   cases, the step size is likely to have a much greater effect
%   on processing time than would the convergence tolerance.
%
%
%   The Confinement Window
%   ======================
%
%   The simplest use of the confinement window is to specify a time
%   interval within which a solution is sought. However, the
%   confinement window can, in some cases, be used to make searches
%   more efficient. Sometimes it's possible to do an efficient search
%   to reduce the size of the time period over which a relatively
%   slow search of interest must be performed.
%
%   Practical use of the coordinate search capability would likely
%   consist of searches over multiple coordinate constraints to find
%   time intervals that satisfies the constraints. An
%   effective technique to accomplish such a search is
%   to use the result window from one search as the confinement window
%   of the next.
%
%   Certain types of searches require the state of the observer,
%   relative to the solar system barycenter, to be computed at times
%   slightly outside the confinement window `cnfine'. The time window
%   that is actually used is the result of "expanding" `cnfine' by a
%   specified amount "T": each time interval of `cnfine' is expanded by
%   shifting the interval's left endpoint to the left and the right
%   endpoint to the right by T seconds. Any overlapping intervals are
%   merged. (The input argument `cnfine' is not modified.)
%
%   The window expansions listed below are additive: if both
%   conditions apply, the window expansion amount is the sum of the
%   individual amounts.
%
%   -  If a search uses an equality constraint, the time window
%      over which the state of the observer is computed is expanded
%      by 1 second at both ends of all of the time intervals
%      comprising the window over which the search is conducted.
%
%   -  If a search uses stellar aberration corrections, the time
%      window over which the state of the observer is computed is
%      expanded as described above.
%
%   When light time corrections are used, expansion of the search
%   window also affects the set of times at which the light time-
%   corrected state of the target is computed.
%
%   In addition to the possible 2 second expansion of the search
%   window that occurs when both an equality constraint and stellar
%   aberration corrections are used, round-off error should be taken
%   into account when the need for data availability is analyzed.
%
%   Longitude and Right Ascension
%   =============================
%
%   The cyclic nature of the longitude and right ascension coordinates
%   produces branch cuts at +/- 180 degrees longitude and 0-360
%   right ascension. Round-off error may cause solutions near these
%   branches to cross the branch. Use of the Mice routine cspice_wncond
%   will contract solution windows by some epsilon, reducing the
%   measure of the windows and eliminating the branch crossing. A
%   one millisecond contraction will in most cases eliminate
%   numerical round-off caused branch crossings.
%
%-Exceptions
%
%   1)  In order for this routine to produce correct results,
%       the step size must be appropriate for the problem at hand.
%       Step sizes that are too large may cause this routine to miss
%       roots; step sizes that are too small may cause this routine
%       to run unacceptably slowly and in some cases, find spurious
%       roots.
%
%       This routine does not diagnose invalid step sizes, except
%       that if the step size is non-positive, an error is signaled
%       by a routine in the call tree of this routine.
%
%   2)  Due to numerical errors, in particular,
%
%          - truncation error in time values
%          - finite tolerance value
%          - errors in computed geometric quantities
%
%       it is *normal* for the condition of interest to not always be
%       satisfied near the endpoints of the intervals comprising the
%       `result' window. One technique to handle such a situation,
%       slightly contract `result' using the window routine cspice_wncond.
%
%   3)  If the output SPICE window `result' has insufficient capacity
%       to contain the number of intervals on which the specified
%       distance condition is met, an error is signaled
%       by a routine in the call tree of this routine.
%
%   4)  If an error (typically cell overflow) occurs during
%       window arithmetic, the error is signaled by a routine
%       in the call tree of this routine.
%
%   5)  If the relational operator `relate' is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If the aberration correction specifier contains an
%       unrecognized value, an error is signaled by a routine in the
%       call tree of this routine.
%
%   7)  If `adjust' is negative, an error is signaled by a routine in
%       the call tree of this routine.
%
%   8)  If either of the input body names do not map to NAIF ID
%       codes, an error is signaled by a routine in the call tree of
%       this routine.
%
%   9)  If required ephemerides or other kernel data are not
%       available, an error is signaled by a routine in the call tree
%       of this routine.
%
%   10) If the search uses GEODETIC or PLANETOGRAPHIC coordinates, and
%       the center body of the reference frame has unequal equatorial
%       radii, an error is signaled by a routine in the call tree of
%       this routine.
%
%   11) If any of the input arguments, `target', `frame', `abcorr',
%       `obsrvr', `crdsys', `coord', `relate', `refval', `adjust',
%       `step', `nintvls' or `cnfine', is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   12) If any of the input arguments, `target', `frame', `abcorr',
%       `obsrvr', `crdsys', `coord', `relate', `refval', `adjust',
%       `step', `nintvls' or `cnfine', is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPK and PCK kernels must be loaded by the calling
%   program before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: the calling application must load ephemeris data
%      for the targets, observer, and any intermediate objects in
%      a chain connecting the targets and observer that cover the
%      time period specified by the window `cnfine'. If aberration
%      corrections are used, the states of target and observer
%      relative to the solar system barycenter must be calculable
%      from the available ephemeris data. Typically ephemeris data
%      are made available by loading one or more SPK files using
%      cspice_furnsh.
%
%   -  If non-inertial reference frames are used, then PCK
%      files, frame kernels, C-kernels, and SCLK kernels may be
%      needed.
%
%   -  In some cases the observer's state may be computed at times
%      outside of `cnfine' by as much as 2 seconds; data required to
%      compute this state must be provided by loaded kernels. See
%      -Particulars for details.
%
%   Such kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   1)  The kernel files to be used by this routine must be loaded
%       (normally using the Mice routine cspice_furnsh) before this
%       routine is called.
%
%   2)  This routine has the side effect of re-initializing the
%       coordinate quantity utility package. Callers may
%       need to re-initialize the package after calling this routine.
%
%-Required_Reading
%
%   MICE.REQ
%   GF.REQ
%   SPK.REQ
%   CK.REQ
%   TIME.REQ
%   WINDOWS.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 03-NOV-2021 (EDW) (JDR)
%
%       Updated header to describe use of expanded confinement window.
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 13-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.2, 05-SEP-2012 (EDW)
%
%       Edit to comments to correct search description.
%
%       Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.1, 10-JUN-2009 (EDW)
%
%       Minor header edit "cosin" -> "cos".
%
%   -Mice Version 1.0.0, 15-APR-2009 (EDW)
%
%-Index_Entries
%
%   GF position coordinate search
%
%-&


function [result] = cspice_gfposc( target, frame, abcorr, obsrvr,  ...
                                   crdsys, coord, relate, refval,  ...
                                   adjust, step,  nintvls, cnfine )

   switch nargin

      case 12

         target  = zzmice_str(target);
         frame   = zzmice_str(frame);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         crdsys  = zzmice_str(crdsys);
         coord  = zzmice_str(coord);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfposc( `target`, `frame`, '  ...
                                      '`abcorr`, `obsrvr`, `crdsys`, '     ...
                                      '`coord`, `relate`, refval, '        ...
                                      'adjust, step, nintvls, cnfine )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfposc_c', target, frame, abcorr,  obsrvr,  ...
                                  crdsys, coord, relate,  refval,  ...
                                  adjust, step,  nintvls,          ...
                                  [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end
