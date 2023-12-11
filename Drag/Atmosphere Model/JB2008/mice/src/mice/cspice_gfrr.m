%-Abstract
%
%   CSPICE_GFRR determines the time intervals for which a specified constraint
%   on the observer-target range rate is met.
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
%               as an integer string. For example both 'MOON' and '301' are
%               legitimate strings that indicate the Moon is the target body.
%
%               Case and leading or trailing blanks are not significant
%               in the string `target'.
%
%               The target and observer define a position vector which
%               points from the observer to the target; the time derivative
%               length of this vector is the "range rate" that serves as
%               the subject of the search performed by this routine.
%
%      abcorr   describes the aberration corrections to apply to the state
%               evaluations to account for one-way light time and stellar
%               aberration.
%
%               [1,c2] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               This routine accepts the same aberration corrections as does
%               the Mice routine cspice_spkezr. See the header of
%               cspice_spkezr for a detailed description of the aberration
%               correction options. For convenience, the options are listed
%               below:
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
%               [1,c3] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example both 'MOON' and '301' are
%               legitimate strings that indicate the Moon is the observer.
%
%               Case and leading or trailing blanks are not significant
%               in the string `obsrvr'.
%
%      relate   the constraint relational operator on observer-target
%               distance.
%
%               [1,c4] = size(relate); char = class(relate)
%
%                  or
%
%               [1,1] = size(relate); cell = class(relate)
%
%               The result window found  by this routine indicates
%               the time intervals where the constraint is satisfied.
%
%               Supported values of `relate' and corresponding meanings are
%               shown below:
%
%                  '>'       Range rate is greater than the reference
%                            value `refval'.
%
%                  '='       Range rate is equal to the reference
%                            value `refval'.
%
%                  '<'       Range rate is less than the reference
%                            value `refval'.
%
%                  'ABSMAX'  Range rate is at an absolute maximum.
%
%                  'ABSMIN'  Range rate is at an absolute  minimum.
%
%                  'LOCMAX'  Range rate is at a local maximum.
%
%                  'LOCMIN'  Range rate is at a local minimum.
%
%               The caller may indicate that the region of interest
%               is the set of time intervals where the quantity is
%               within a specified distance of an absolute extremum.
%               The argument `adjust' (described below) is used to
%               specify this distance.
%
%               Local extrema are considered to exist only in the
%               interiors of the intervals comprising the confinement
%               window:  a local extremum cannot exist at a boundary
%               point of the confinement window.
%
%               The `relate' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      refval   reference value used together with `relate' argument to define
%               an equality or inequality to satisfy by the observer-target
%               distance.
%
%               [1,1] = size(refval); double = class(refval)
%
%               See the discussion of `relate' above for further information.
%
%               The units of `refval' are km.
%
%      adjust   value used to modify searches for absolute extrema.
%
%               [1,1] = size(adjust); double = class(adjust)
%
%               When `relate' is set to 'ABSMAX' or 'ABSMIN' and `adjust' is
%               set to a positive value, cspice_gfrr finds times when the
%               observer-target vector coordinate is within `adjust'
%               kilometers/second of the specified extreme value.
%
%               For `relate' set to 'ABSMAX', the result window contains
%               time intervals when the observer-target vector coordinate has
%               values between ABSMAX - adjust and ABSMAX.
%
%               For `relate' set to 'ABSMIN', the result window contains
%               time intervals when the phase angle has values between
%               ABSMIN and ABSMIN + adjust.
%
%               `adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     time step size to use in the search.
%
%               [1,1] = size(step); double = class(step)
%
%               `step' must be short enough for a search using this step size
%               to locate the time intervals where coordinate function of the
%               observer-target vector is monotone increasing or decreasing.
%               However, `step' must not be *too* short, or the search will
%               take an unreasonable amount of time.
%
%               The choice of `step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
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
%               `cnfine' may consist of a single interval or a collection of
%               intervals.
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
%      [result] = cspice_gfrr( target, abcorr, obsrvr,  relate, refval,    ...
%                              adjust, step,   nintvls, cnfine )
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
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Determine the time windows from January 1, 2007 UTC to
%      April 1, 2007 UTC for which the sun-moon range rate satisfies the
%      relation conditions with respect to a reference value of
%      0.3365 km/s radians (this range rate known to occur within the
%      search interval). Also determine the time windows corresponding
%      to the local maximum and minimum range rate, and the absolute
%      maximum and minimum range rate during the search interval.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfrr_ex1.tm
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
%      function gfrr_ex1()
%
%         MAXWIN  =  2000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###';
%
%         relate = { '=', '<', '>', 'LOCMIN', 'ABSMIN', 'LOCMAX', 'ABSMAX' };
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfrr_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2007 APR 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % Search using a step size of 1 day (in units of seconds).
%         % The reference value is .3365 km/s. We're not using the
%         % adjustment feature, so we set 'adjust' to zero.
%         %
%         step   = cspice_spd;
%         adjust = 0.D0;
%         refval = .3365D0;
%
%         target  = 'MOON';
%         abcorr  = 'NONE';
%         obsrvr  = 'SUN';
%         nintvls = MAXWIN;
%
%         for j=1:7
%
%            fprintf( 'Relation condition: %s\n',  char( relate(j) ) )
%
%            %
%            % Perform the search. The SPICE window 'result' contains
%            % the set of times when the condition is met.
%            %
%            result = cspice_gfrr( target,    abcorr,  obsrvr,             ...
%                                  relate(j), refval,  adjust,             ...
%                                  step,      nintvls, cnfine );
%
%            %
%            % List the beginning and ending times in each interval
%            % if 'result' contains data.
%            %
%            count = cspice_wncard( result );
%
%            if ( isequal(count,0) )
%
%               disp( 'Result window is empty.' )
%
%            else
%
%               for i= 1:count
%
%                  %
%                  % Fetch the endpoints of the Ith interval
%                  % of the result window.
%                  %
%                 [left, right] = cspice_wnfetd( result, i );
%
%                 timstr = cspice_timout( [left,right], TIMFMT );
%
%                 disp( ['Start time, drdt = ', timstr(1,:) ] )
%                 disp( ['Stop time,  drdt = ', timstr(2,:) ] )
%
%               end
%
%               disp( ' ' )
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
%      Relation condition: =
%      Start time, drdt = 2007-JAN-02 00:35:19.571
%      Stop time,  drdt = 2007-JAN-02 00:35:19.571
%      Start time, drdt = 2007-JAN-19 22:04:54.897
%      Stop time,  drdt = 2007-JAN-19 22:04:54.897
%      Start time, drdt = 2007-FEB-01 23:30:13.427
%      Stop time,  drdt = 2007-FEB-01 23:30:13.427
%      Start time, drdt = 2007-FEB-17 11:10:46.538
%      Stop time,  drdt = 2007-FEB-17 11:10:46.538
%      Start time, drdt = 2007-MAR-04 15:50:19.929
%      Stop time,  drdt = 2007-MAR-04 15:50:19.929
%      Start time, drdt = 2007-MAR-18 09:59:05.957
%      Stop time,  drdt = 2007-MAR-18 09:59:05.957
%
%      Relation condition: <
%      Start time, drdt = 2007-JAN-02 00:35:19.571
%      Stop time,  drdt = 2007-JAN-19 22:04:54.897
%      Start time, drdt = 2007-FEB-01 23:30:13.427
%      Stop time,  drdt = 2007-FEB-17 11:10:46.538
%      Start time, drdt = 2007-MAR-04 15:50:19.929
%      Stop time,  drdt = 2007-MAR-18 09:59:05.957
%
%      Relation condition: >
%      Start time, drdt = 2007-JAN-01 00:00:00.000
%      Stop time,  drdt = 2007-JAN-02 00:35:19.571
%      Start time, drdt = 2007-JAN-19 22:04:54.897
%      Stop time,  drdt = 2007-FEB-01 23:30:13.427
%      Start time, drdt = 2007-FEB-17 11:10:46.538
%      Stop time,  drdt = 2007-MAR-04 15:50:19.929
%      Start time, drdt = 2007-MAR-18 09:59:05.957
%      Stop time,  drdt = 2007-APR-01 00:00:00.000
%
%      Relation condition: LOCMIN
%      Start time, drdt = 2007-JAN-11 07:03:58.991
%      Stop time,  drdt = 2007-JAN-11 07:03:58.991
%      Start time, drdt = 2007-FEB-10 06:26:15.441
%      Stop time,  drdt = 2007-FEB-10 06:26:15.441
%      Start time, drdt = 2007-MAR-12 03:28:36.404
%      Stop time,  drdt = 2007-MAR-12 03:28:36.404
%
%      Relation condition: ABSMIN
%      Start time, drdt = 2007-JAN-11 07:03:58.991
%      Stop time,  drdt = 2007-JAN-11 07:03:58.991
%
%      Relation condition: LOCMAX
%      Start time, drdt = 2007-JAN-26 02:27:33.762
%      Stop time,  drdt = 2007-JAN-26 02:27:33.762
%      Start time, drdt = 2007-FEB-24 09:35:07.812
%      Stop time,  drdt = 2007-FEB-24 09:35:07.812
%      Start time, drdt = 2007-MAR-25 17:26:56.148
%      Stop time,  drdt = 2007-MAR-25 17:26:56.148
%
%      Relation condition: ABSMAX
%      Start time, drdt = 2007-MAR-25 17:26:56.148
%      Stop time,  drdt = 2007-MAR-25 17:26:56.148
%
%
%-Particulars
%
%   This routine determines if the caller-specified constraint
%   condition on the geometric event (range rate) is satisfied for
%   any time intervals within the confinement window `cnfine'. If one
%   or more such time intervals exist, those intervals are added
%   to the `result' window.
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
%   periods, within the confinement window, over which the
%   range rate function is monotone increasing and monotone
%   decreasing. Each of these time periods is represented by a SPICE
%   window. Having found these windows, all of the range rate
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
%   change of range rate will be sampled. Starting at
%   the left endpoint of an interval, samples will be taken at each
%   step. If a change of sign is found, a root has been bracketed; at
%   that point, the time at which the time derivative of the
%   range rate is zero can be found by a refinement process, for
%   example, using a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the range rate function is monotone:
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
%   attained and times when the range rate function is equal to a
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
%       if the step size is non-positive, an error is signaled by a
%       routine in the call tree of this routine.
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
%   3)  If the SPICE window `result' has insufficient capacity to
%       contain the number of intervals on which the specified
%       distance condition is met, an error is signaled by a routine
%       in the call tree of this routine.
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
%   8)  If `adjust' has a non-zero value when `relate' has any value other
%       than 'ABSMIN' or 'ABSMAX', an error is signaled by a routine
%       in the call tree of this routine.
%
%   9)  If either of the input body names do not map to NAIF ID
%       codes, an error is signaled by a routine in the call tree of
%       this routine.
%
%   10) If required ephemerides or other kernel data are not
%       available, an error is signaled by a routine in the call tree
%       of this routine.
%
%   11) If any of the input arguments, `target', `abcorr', `obsrvr',
%       `relate', `refval', `adjust', `step', `nintvls' or `cnfine',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   12) If any of the input arguments, `target', `abcorr', `obsrvr',
%       `relate', `refval', `adjust', `step', `nintvls' or `cnfine',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
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
%   -  In some cases the observer's state may be computed at times
%      outside of `cnfine' by as much as 2 seconds; data required to
%      compute this state must be provided by loaded kernels. See
%      -Particulars for details.
%
%   Kernel data are normally loaded once per program run, NOT every
%   time this routine is called.
%
%-Restrictions
%
%   1)  The kernel files to be used by this routine must be loaded
%       (normally using the Mice routine cspice_furnsh) before this
%       routine is called.
%
%   2)  This routine has the side effect of re-initializing the
%       range rate quantity utility package. Callers may themselves
%       need to re-initialize the range rate quantity utility
%       package after calling this routine.
%
%-Required_Reading
%
%   MICE.REQ
%   GF.REQ
%   SPK.REQ
%   NAIF_IDS.REQ
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
%       Edited the header to comply with NAIF standard. Added -Parameters,
%       -Exceptions, -Files, -Restrictions, -Literature_References and
%       -Author_and_Institution sections.
%
%       Updated header to describe use of expanded confinement window.
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
%   -Mice Version 1.0.1, 05-SEP-2012 (EDW)
%
%       Edit to comments to correct search description.
%
%       Corrected minor typo in header.
%
%       Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.0, 16-FEB-2010 (EDW)
%
%-Index_Entries
%
%   GF range rate search
%
%-&

function [result] = cspice_gfrr( target, abcorr, obsrvr,  relate, refval, ...
                                 adjust, step,   nintvls, cnfine )

   switch nargin

      case 9

         target  = zzmice_str(target);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfrr( `target`, `abcorr`, '    ...
                                     '`obsrvr`, `relate`, refval, adjust, ' ...
                                     'step, nintvls, cnfine )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfrr_c', target, abcorr, obsrvr, relate, ...
                                refval, adjust, step, nintvls,  ...
                                [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end
