%-Abstract
%
%   CSPICE_GFILUM determines the time intervals over which a specified
%   constraint on the observed phase, solar incidence, or emission angle
%   at a specified target body surface point is met.
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
%      method   name specifying the computation method to use.
%
%               [1,c1] = size(method); char = class(method)
%
%                 or
%
%               [1,1] = size(method); cell = class(method)
%
%               Parameters include, but are not limited to, the shape model
%               used to represent the surface of the target body.
%
%               The only choice currently supported is
%
%                  'Ellipsoid'        The illumination angle
%                                     computation uses a triaxial
%                                     ellipsoid to model the surface
%                                     of the target body. The
%                                     ellipsoid's radii must be
%                                     available in the kernel pool.
%
%               Neither case nor white space are significant in
%               `method'. For example, the string ' eLLipsoid ' is
%               valid.
%
%      angtyp   name specifying the illumination angle for which a search
%               is to be performed.
%
%               [1,c2] = size(angtyp); char = class(angtyp)
%
%                 or
%
%               [1,1] = size(angtyp); cell = class(angtyp)
%
%               The possible values of `angtyp' are:
%
%                  'PHASE'
%                  'INCIDENCE'
%                  'EMISSION'
%
%               See the -Particulars section below for a detailed
%               description of these angles.
%
%               Neither case nor white space are significant in
%               `angtyp'. For example, the string ' Incidence ' is
%               valid.
%
%      target   name of the target body.
%
%               [1,c3] = size(target); char = class(target)
%
%                 or
%
%               [1,1] = size(target); cell = class(target)
%
%               The point at which the illumination angles are defined is
%               located on the surface of this body.
%
%               Optionally, you may supply the integer ID code for
%               the object as an integer string. For example both
%               'MOON' and '301' are legitimate strings that indicate
%               the moon is the target body.
%
%      illmn    name of the illumination source.
%
%               [1,c4] = size(illmn); char = class(illmn)
%
%                 or
%
%               [1,1] = size(illmn); cell = class(illmn)
%
%               This source may be any ephemeris object. Case, blanks, and
%               numeric values are treated in the same way as for the
%               input `target'.
%
%      fixref   name of the body-fixed, body-centered reference frame
%               associated with the target body.
%
%               [1,c5] = size(fixref); char = class(fixref)
%
%                 or
%
%               [1,1] = size(fixref); cell = class(fixref)
%
%               The input surface point `spoint' is expressed relative to
%               this reference frame, and this frame is used to
%               define the orientation of the target body as a
%               function of time.
%
%               The string `fixref' is case-insensitive, and leading
%               and trailing blanks in `fixref' are not significant.
%
%      abcorr   describes the aberration corrections to apply to the state
%               evaluations to account for one-way light time and stellar
%               aberration.
%
%               [1,c6] = size(abcorr); char = class(abcorr)
%
%                 or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               Any 'reception' correction accepted by cspice_spkezr can be
%               used here. See the header of cspice_spkezr for a detailed
%               description of the aberration correction options.
%
%               For convenience, the options are listed below:
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       'Reception' case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     'Reception' case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'CN'       'Reception' case: converged
%                             Newtonian light time correction.
%
%                  'CN+S'     'Reception' case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               Case and blanks are not significant in the string
%               `abcorr'.
%
%      obsrvr   name of the observing body.
%
%               [1,c7] = size(obsrvr); char = class(obsrvr)
%
%                 or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate that the
%               observer is Earth.
%
%      spoint   a surface point on the target body, expressed in
%               Cartesian coordinates, relative to the body-fixed
%               target frame designated by `fixref'.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               `spoint' need not be visible from the observer's
%               location in order for the constraint specified by
%               `relate' and `refval' (see descriptions below) to be
%               satisfied.
%
%               The components of `spoint' have units of km.
%
%      relate   describes the constraint relational operator on a specified
%               illumination angle.
%
%               [1,c8] = size(relate); char = class(relate)
%
%                 or
%
%               [1,1] = size(relate); cell = class(relate)
%
%               The result window found by this routine indicates the time
%               intervals where the constraint is satisfied.
%
%               Supported values of `relate' and corresponding meanings
%               are shown below:
%
%                   '>'      The angle is greater than the reference
%                            value `refval'.
%
%                   '='      The angle is equal to the reference
%                            value `refval'.
%
%                   '<'      The angle is less than the reference
%                            value `refval'.
%
%                  'ABSMAX'  The angle is at an absolute maximum.
%
%                  'ABSMIN'  The angle is at an absolute  minimum.
%
%                  'LOCMAX'  The angle is at a local maximum.
%
%                  'LOCMIN'  The angle is at a local minimum.
%
%               The caller may indicate that the region of interest is
%               the set of time intervals where the angle is within a
%               specified separation from an absolute extremum. The
%               argument `adjust' (described below) is used to specify
%               this separation.
%
%               Local extrema are considered to exist only in the
%               interiors of the intervals comprising the confinement
%               window: a local extremum cannot exist at a boundary
%               point of the confinement window.
%
%               Case is not significant in the string `relate'.
%
%      refval   reference value used together with the argument
%               `relate' to define an equality or inequality to be
%               satisfied by the specified illumination angle.
%
%               [1,1] = size(refval); double = class(refval)
%
%               See the discussion of `relate' above for further information.
%
%               The units of `refval' are radians.
%
%      adjust   parameter used to modify searches for absolute extrema.
%
%               [1,1] = size(adjust); double = class(adjust)
%
%               When `relate' is set to 'ABSMAX' or 'ABSMIN' and `adjust'
%               is set to a positive value, cspice_gfilum will find times
%               when the observer-target distance is within `adjust' km of
%               the specified extreme value.
%
%               If `adjust' is non-zero and a search for an absolute
%               minimum `min' is performed, the result window contains
%               time intervals when the observer-target distance has
%               values between `min' and min+adjust.
%
%               If the search is for an absolute maximum `max', the
%               corresponding range is from max-adjust to `max'.
%
%               `adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     step size to use in the search.
%
%               [1,1] = size(step); double = class(step)
%
%               `step' must be short enough for a search using this step
%               size to locate the time intervals where the specified
%               illumination angle is monotone increasing or
%               decreasing. However, `step' must not be *too* short, or
%               the search will take an unreasonable amount of time.
%
%               The choice of `step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               `step' has units of seconds.
%
%      nintvls  a parameter specifying the number of intervals that
%               can be accommodated by each of the dynamically allocated
%               workspace windows used internally by this routine.
%
%               [1,1] = size(nintvls); int32 = class(nintvls)
%
%               In many cases, it's not necessary to compute an accurate
%               estimate of how many intervals are needed; rather, the
%               user can pick a size considerably larger than what's
%               really required.
%
%               However, since excessively large arrays can prevent
%               applications from compiling, linking, or running
%               properly, sometimes `nintvls' must be set according to
%               the actual workspace requirement. A rule of thumb for
%               the number of intervals needed is
%
%                  nintvls  =  2*n  +  ( m / step )
%
%               where
%
%                  n     is the number of intervals in the confinement
%                        window
%
%                  m     is the measure of the confinement window, in
%                        units of seconds
%
%                  step  is the search step size in seconds
%
%      cnfine   a SPICE window that confines the time period over
%               which the specified search is conducted.
%
%               [2r,1] = size(cnfine); double = class(cnfine)
%
%               `cnfine' may consist of a single interval or a collection of
%               intervals.
%
%               The endpoints of the time intervals comprising `cnfine'
%               are interpreted as seconds past J2000 TDB.
%
%               See the -Examples section below for a code example that
%               shows how to create a confinement window.
%
%               In some cases the observer's state may be computed at
%               times outside of `cnfine' by as much as 2 seconds. See
%               -Particulars for details.
%
%   the call:
%
%      [result] = cspice_gfilum( method, angtyp, target,  illmn,  fixref,  ...
%                                abcorr, obsrvr, spoint,  relate, refval,  ...
%                                adjust, step,   nintvls, cnfine )
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window `cnfine', on which the specified
%               constraint is satisfied.
%
%               [2s,1] = size(result); double = class(result)
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
%                is the convergence tolerance used for finding
%                endpoints of the intervals comprising the result
%                window.  SPICE_GF_CNVTOL is used to determine when
%                binary searches for roots should terminate: when a
%                root is bracketed within an interval of length
%                SPICE_GF_CNVTOL, the root is considered to have
%                been found.
%
%                The accuracy, as opposed to precision, of roots found
%                by this routine depends on the accuracy of the input
%                data. In most cases, the accuracy of solutions will be
%                inferior to their precision.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Determine time intervals over which the planned Mars Science
%      Laboratory (MSL) Gale Crater landing site satisfies certain
%      constraints on its illumination and visibility as seen from
%      the Mars Reconnaissance Orbiter (MRO) spacecraft. The
%      observation period will range from slightly before the planned
%      landing time to about 10 days later.
%
%      In this case we require the emission angle to be less than
%      30 degrees and the solar incidence angle to be less than
%      40 degrees.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfilum_ex1.tm
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
%            File name         Contents
%            ---------         --------
%            de421.bsp         Planetary ephemeris
%            pck00010.tpc      Planet orientation
%                              and radii
%            naif0012.tls      Leapseconds
%            mro_psp24.bsp     MRO ephemeris
%
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'de421.bsp',
%                             'pck00010.tpc',
%                             'naif0012.tls',
%                             'mro_psp24.bsp' )
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function gfilum_ex1()
%
%         %
%         % Output time format:
%         %
%         TIMFMT = 'YYYY MON DD HR:MN:SC.###### TDB::TDB';
%
%         %
%         % Meta-kernel name:
%         %
%         META = 'gfilum_ex1.tm';
%
%         %
%         % Maximum number of intervals in the windows
%         % used in this program:
%         %
%         MAXIVL = 1000;
%
%         %
%         % Local variables
%         %
%         r2d    = cspice_dpr();
%
%         %
%         % Initial values
%         %
%         % Mars planetodetic coordinates of landing site.
%         % Angular units are degrees; distance units are km.
%         %
%         gclat  =  -4.543182;
%         gclon  = 137.420000;
%         gcalt  =  -4.876405;
%
%         %
%         % Load kernels:
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert the landing site location from planetodetic
%         % to Cartesian coordinates for use with GFILUM.
%         %
%         radii = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%         re = radii(1);
%         rp = radii(3);
%
%         f  = ( re - rp ) / re;
%
%         gcpos = cspice_georec( gclon * cspice_rpd(),                     ...
%                                gclat * cspice_rpd(),                     ...
%                                gcalt, re, f);
%
%         %
%         % Set the search interval:
%         %
%         utcbeg = '2012 AUG 5 00:00:00 UTC';
%         et0    = cspice_str2et( utcbeg );
%
%         utcend = '2012 SEP 15 00:00:00 UTC';
%         et1    = cspice_str2et( utcend );
%
%         cnfine = cspice_wninsd( et0, et1 );
%
%
%         %
%         % Set observer, target, aberration correction, and the
%         % Mars body-fixed, body-centered reference frame. The
%         % lighting source is the sun.
%         %
%         % Aberration corrections are set for remote observations.
%         %
%         illmn  = 'sun';
%         obsrvr = 'mro';
%         target = 'mars';
%         abcorr = 'cn+s';
%         fixref = 'iau_mars';
%
%         %
%         % Initialize the adjustment value for absolute
%         % extremum searches. We're not performing
%         % such searches in this example, but this input
%         % to GFILUM must still be set.
%         %
%         adjust = 0.0;
%
%         %
%         % The computation uses an ellipsoidal model for the
%         % target body shape.
%         %
%         method = 'Ellipsoid';
%
%         %
%         % Set the reference value to use for the solar
%         % incidence angle search.
%         %
%         refval = 45.0 * cspice_rpd();
%
%         %
%         % Since the period of the solar incidence angle
%         % is about one Martian day, we can safely use 6 hours
%         % as the search step.
%         %
%         step   = 21600.0;
%
%         %
%         % Search over the confinement window for times
%         % when the solar incidence angle is less than
%         % the reference value.
%         %
%         [wnsolr] = cspice_gfilum( method, 'INCIDENCE', target,           ...
%                                   illmn,  fixref,      abcorr,           ...
%                                   obsrvr, gcpos,       '<',              ...
%                                   refval, adjust,      step,             ...
%                                   MAXIVL, cnfine );
%
%         %
%         % With the search on the incidence angle complete, perform
%         % a search on the emission angle.
%         %
%         % Set the reference value for the emission angle search.
%         %
%         refval = 80.0 * cspice_rpd();
%
%         %
%         % We'll use 15 minutes as the search step. This step
%         % is small enough to be suitable for Mars orbiters.
%         % Units are seconds.
%         %
%         step   = 900.0;
%
%         %
%         % Search over the previous result window for times when the
%         % emission angle is less than the reference value.
%         %
%         [result] = cspice_gfilum( method, 'EMISSION', target, illmn,     ...
%                                   fixref, abcorr,     obsrvr, gcpos,     ...
%                                   '<',    refval,     adjust, step,      ...
%                                   MAXIVL, wnsolr );
%
%         %
%         % Display the result window. Show the solar incidence
%         % and emission angles at the window's interval
%         % boundaries.
%         %
%         if ( cspice_wncard( result ) == 0 )
%
%            disp( '     Window is empty: condition is not met.' )
%
%         else
%
%            fprintf( '                                  '   )
%            fprintf( '       Solar Incidence   Emission\n'  )
%            fprintf( '                                  '   )
%            fprintf( '             (deg)         (deg)\n\n' )
%
%            for i=1:cspice_wncard( result )
%
%               [start, finish] = cspice_wnfetd( result, i );
%
%               %
%               % Compute the angles of interest at the boundary
%               % epochs.
%               %
%               timstr = cspice_timout( start, TIMFMT );
%               [trgepc, srfvec, phase, solar, emissn] =                   ...
%                                        cspice_ilumin( method, target,    ...
%                                                       start,  fixref,    ...
%                                                       abcorr, obsrvr,    ...
%                                                       gcpos );
%
%                  fprintf ( ' Start: %s %14.9f %14.9f\n', timstr,         ...
%                                                          solar*r2d,      ...
%                                                          emissn*r2d )
%
%
%               timstr = cspice_timout( finish, TIMFMT);
%               [trgepc, srfvec, phase, solar, emissn] =                   ...
%                                        cspice_ilumin( method, target,    ...
%                                                       finish, fixref,    ...
%                                                       abcorr, obsrvr,    ...
%                                                       gcpos );
%
%                  fprintf ( ' Start: %s %14.9f %14.9f\n\n', timstr,       ...
%                                                            solar*r2d,    ...
%                                                            emissn*r2d )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%                                               Solar Incidence   Emission
%                                                     (deg)         (deg)
%
%       Start: 2012 AUG 09 06:14:46.475539 TDB   41.793493032   80.000000000
%       Start: 2012 AUG 09 06:15:29.695045 TDB   41.954623385   80.000000002
%
%       Start: 2012 AUG 14 09:37:47.093234 TDB   42.772767813   80.000000007
%       Start: 2012 AUG 14 09:41:59.554719 TDB   43.729251675   79.999999998
%
%       Start: 2012 AUG 19 13:01:43.056249 TDB   44.000361046   80.000000017
%       Start: 2012 AUG 19 13:06:03.429007 TDB   44.999999999   75.754083310
%
%       Start: 2012 AUG 30 20:10:42.196910 TDB   42.214690783   79.999999993
%       Start: 2012 AUG 30 20:14:47.411493 TDB   43.170768309   79.999999996
%
%       Start: 2012 SEP 04 23:35:53.476437 TDB   43.804510481   79.999999983
%       Start: 2012 SEP 04 23:40:57.001978 TDB   45.000000001   77.221887661
%
%       Start: 2012 SEP 11 03:22:35.751759 TDB   41.115348965   80.000000009
%       Start: 2012 SEP 11 03:24:59.610628 TDB   41.684463728   79.999999996
%
%
%-Particulars
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when the specified illumination
%   angle satisfies a caller-specified constraint. The resulting set
%   of intervals is returned as a SPICE window.
%
%   The term "illumination angles" refers to the following set of
%   angles:
%
%
%      phase angle              Angle between the vectors from the
%                               surface point to the observer and
%                               from the surface point to the
%                               illumination source.
%
%      incidence angle          Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               illumination source. When the sun is
%                               the illumination source, this angle is
%                               commonly called the "solar incidence
%                               angle."
%
%      emission angle           Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               observer.
%
%   The diagram below illustrates the geometric relationships
%   defining these angles. The labels for the incidence, emission,
%   and phase angles are "inc.", "e.", and "phase".
%
%
%
%                                                    *
%                                            illumination source
%
%                  surface normal vector
%                            ._                 _.
%                            |\                 /|  illumination
%                              \    phase      /    source vector
%                               \   .    .    /
%                               .            .
%                                 \   ___   /
%                            .     \/     \/
%                                  _\ inc./
%                           .    /   \   /
%                           .   |  e. \ /
%       *             <--------------- *  surface point on
%    viewing            vector            target body
%    location           to viewing
%    (observer)         location
%
%
%
%   Note that if the target-observer vector, the target normal vector
%   at the surface point, and the target-illumination source vector
%   are coplanar, then phase is the sum of the incidence and emission
%   angles. This rarely occurs; usually
%
%      phase angle  <  incidence angle + emission angle
%
%   All of the above angles can be computed using light time
%   corrections, light time and stellar aberration corrections, or no
%   aberration corrections. In order to describe apparent geometry as
%   observed by a remote sensing instrument, both light time and
%   stellar aberration corrections should be used.
%
%   The way aberration corrections are applied by this routine
%   is described below.
%
%      Light time corrections
%      ======================
%
%         Observer-target surface point vector
%         ------------------------------------
%
%         Let `et' be the epoch at which an observation or remote
%         sensing measurement is made, and let et - lt (`lt' stands
%         for "light time") be the epoch at which the photons
%         received at `et' were emitted from the surface point `spoint'.
%         Note that the light time between the surface point and
%         observer will generally differ from the light time between
%         the target body's center and the observer.
%
%
%         Target body's orientation
%         -------------------------
%
%         Using the definitions of `et' and `lt' above, the target body's
%         orientation at et - lt is used. The surface normal is
%         dependent on the target body's orientation, so the body's
%         orientation model must be evaluated for the correct epoch.
%
%
%         Target body -- illumination source vector
%         -----------------------------------------
%
%         The surface features on the target body near `spoint' will
%         appear in a measurement made at `et' as they were at et-lt.
%         In particular, lighting on the target body is dependent on
%         the apparent location of the illumination source as seen
%         from the target body at et-lt. So, a second light time
%         correction is used to compute the position of the
%         illumination source relative to the surface point.
%
%
%      Stellar aberration corrections
%      ==============================
%
%      Stellar aberration corrections are applied only if
%      light time corrections are applied as well.
%
%         Observer-target surface point body vector
%         -----------------------------------------
%
%         When stellar aberration correction is performed, the
%         observer-to-surface point direction vector, which we'll
%         call SRFVEC, is adjusted so as to point to the apparent
%         position of `spoint': considering `spoint' to be an ephemeris
%         object, SRFVEC points from the observer's position at `et' to
%         the light time and stellar aberration
%         corrected position of `spoint'.
%
%         Target body-illumination source vector
%         --------------------------------------
%
%         The target body-illumination source vector is the apparent
%         position of the illumination source, corrected for light
%         time and stellar aberration, as seen from the surface point
%         `spoint' at time et-lt.
%
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
%   illumination angle is monotone increasing and monotone decreasing.
%   Each of these time periods is represented by a SPICE window.
%   Having found these windows, all of the illumination angle's local
%   extrema within the confinement window are known. Absolute extrema
%   then can be found very easily.
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
%   The monotone windows (described above) are found via a two-step
%   search process. Each interval of the confinement window is
%   searched as follows: first, the input step size is used to
%   determine the time separation at which the sign of the rate of
%   change of the illumination angle will be sampled. Starting at the
%   left endpoint of an interval, samples will be taken at each step.
%   If a change of sign is found, a root has been bracketed; at that
%   point, the time at which the rate of change of the selected
%   illumination angle is zero can be found by a refinement process,
%   for example, via binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the illumination angle is monotone:
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
%   Having some knowledge of the relative geometry of the target,
%   observer, and illumination source can be a valuable aid in
%   picking a reasonable step size. In general, the user can
%   compensate for lack of such knowledge by picking a very short
%   step size; the cost is increased computation time.
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
%   attained and times when the illumination angle is equal to a
%   reference value. All endpoints of the intervals comprising the
%   result window are either endpoints of intervals of the
%   confinement window or roots.
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   "convergence tolerance." The convergence tolerance used by this
%   routine is set via the parameter SPICE_GF_CNVTOL.
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't become the limiting factor in the accuracy of
%   solutions found by this routine. In general the accuracy of input
%   data will be the limiting factor.
%
%   The user may change the convergence tolerance from the default
%   SPICE_GF_CNVTOL value by calling the routine cspice_gfstol, e.g.
%
%      cspice_gfstol( tolerance value in seconds );
%
%   Call cspice_gfstol prior to calling this routine. All subsequent
%   searches will use the updated tolerance value.
%
%   Searches over time windows of long duration may require use of
%   larger tolerance values than the default: the tolerance must be
%   large enough so that it, when added to or subtracted from the
%   confinement window's lower and upper bounds, yields distinct time
%   values.
%
%   Setting the tolerance tighter than SPICE_GF_CNVTOL is unlikely to be
%   useful, since the results are unlikely to be more accurate.
%   Making the tolerance looser will speed up searches somewhat,
%   since a few convergence steps will be omitted.
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
%-Exceptions
%
%   1)  In order for this routine to produce correct results,
%       the step size must be appropriate for the problem at hand.
%       Step sizes that are too large may cause this routine to miss
%       roots; step sizes that are too small may cause this routine
%       to run unacceptably slowly and in some cases, find spurious
%       roots.
%
%       This routine does not diagnose invalid step sizes, except that
%       if the step size is non-positive, the error SPICE(INVALIDSTEP)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  Due to numerical errors, in particular,
%
%          - Truncation error in time values
%          - Finite tolerance value
%          - Errors in computed geometric quantities
%
%       it is *normal* for the condition of interest to not always be
%       satisfied near the endpoints of the intervals comprising the
%       result window.
%
%       The result window may need to be contracted slightly by the
%       caller to achieve desired results. The SPICE window routine
%       cspice_wncond can be used to contract the result window.
%
%   3)  If an error (typically cell overflow) occurs while performing
%       window arithmetic, the error is signaled by a routine
%       in the call tree of this routine.
%
%   4)  If the output SPICE window `result' has insufficient capacity to
%       hold the set of intervals on which the specified illumination
%       angle condition is met, an error is signaled by a routine in
%       the call tree of this routine.
%
%   5)  If the input target body-fixed frame `fixref' is not
%       recognized, an error is signaled by a routine in the call
%       tree of this routine. A frame name may fail to be recognized
%       because a required frame specification kernel has not been
%       loaded; another cause is a misspelling of the frame name.
%
%   6)  If the input frame `fixref' is not centered at the target body,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   7)  If the input argument `method' is not recognized, an error is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the illumination angle type `angtyp' is not recognized,
%       an error is signaled by a routine in the call tree
%       of this routine.
%
%   9)  If the relational operator `relate' is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   10) If the aberration correction specifier contains an
%       unrecognized value, an error is signaled by a routine in the
%       call tree of this routine.
%
%   11) If `adjust' is negative, an error is signaled by a routine in
%       the call tree of this routine.
%
%   12) If any of the input body names do not map to NAIF ID
%       codes, an error is signaled by a routine in the call tree of
%       this routine.
%
%   13) If the target coincides with the observer or the illumination
%       source, an error is signaled by a routine in the call tree
%       of this routine.
%
%   14) If required ephemerides or other kernel data are not
%       available, an error is signaled by a routine in the call tree
%       of this routine.
%
%   15) If any of the input arguments, `method', `angtyp', `target',
%       `illmn', `fixref', `abcorr', `obsrvr', `spoint', `relate',
%       `refval', `adjust', `step', `nintvls' or `cnfine', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   16) If any of the input arguments, `method', `angtyp', `target',
%       `illmn', `fixref', `abcorr', `obsrvr', `spoint', `relate',
%       `refval', `adjust', `step', `nintvls' or `cnfine', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   Appropriate kernels must be loaded by the calling program before
%   this routine is called.
%
%   The following data are required:
%
%   -  SPK data: ephemeris data for target, observer, and the
%      illumination source must be loaded. If aberration
%      corrections are used, the states of target, observer, and
%      the illumination source relative to the solar system
%      barycenter must be calculable from the available ephemeris
%      data. Typically ephemeris data are made available by loading
%      one or more SPK files via cspice_furnsh.
%
%   -  PCK data: if the target body shape is modeled as an
%      ellipsoid (currently no other shapes are supported),
%      triaxial radii for the target body must be loaded
%      into the kernel pool. Typically this is done by loading a
%      text PCK file via cspice_furnsh.
%
%   -  Further PCK data: rotation data for the target body must be
%      loaded. These may be provided in a text or binary PCK file.
%
%   -  Frame data: if a frame definition not built into SPICE
%      is required to convert the observer and target states to the
%      body-fixed frame of the target, that definition must be
%      available in the kernel pool. Typically the definition is
%      supplied by loading a frame kernel via cspice_furnsh.
%
%   -  In some cases the observer's state may be computed at times
%      outside of `cnfine' by as much as 2 seconds; data required to
%      compute this state must be provided by loaded kernels. See
%      -Particulars for details.
%
%   In all cases, kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   1)  The kernel files to be used by this routine must be loaded
%       (normally using the Mice routine cspice_furnsh) before this
%       routine is called.
%
%   2)  This routine has the side effect of re-initializing the
%       illumination angle utility package. Callers may
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
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Changed the input argument name "illum" to "illmn" for
%       consistency with other routines.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       edited -I/O section to comply with NAIF standard. Fixed minor typos in
%       header.
%
%       Updated header to describe use of expanded confinement window.
%
%       Edited the header to comply with NAIF standard.
%
%       Corrected error in header that listed 'SOLAR INCIDENCE' as an
%       allowed angle type rather than the correct value 'INCIDENCE'.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 11-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 07-NOV-2013 (EDW)
%
%-Index_Entries
%
%   solve for illumination_angle constraints
%   solve for phase_angle constraints
%   solve for solar_incidence_angle constraints
%   solve for incidence_angle constraints
%   solve for emission_angle constraints
%   search using illumination_angle constraints
%   search using lighting_angle constraints
%
%-&

function [result] = cspice_gfilum( method, angtyp,  target, illmn,  ...
                                   fixref, abcorr,  obsrvr, spoint, ...
                                   relate, refval,  adjust,         ...
                                   step,   nintvls, cnfine )

   switch nargin

      case 14

         method  = zzmice_str(method);
         angtyp  = zzmice_str(angtyp);
         target  = zzmice_str(target);
         illmn   = zzmice_str(illmn);
         fixref  = zzmice_str(fixref);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         spoint  = zzmice_dp(spoint);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfilum( `method`, `angtyp`, '...
                                    '`target`, `illmn`, `fixref`, '       ...
                                    '`abcorr`, `obsrvr`, spoint[3], '     ...
                                    '`relate`, refval, adjust, '          ...
                                    'step, nintvls, cnfine, )' ] )

   end

   %
   % Call the GF routine, add to `cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfilum_c', method,  angtyp, target, illmn,  ...
                                  fixref,  abcorr, obsrvr, spoint, ...
                                  relate,  refval, adjust, step,   ...
                                  nintvls, [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end
