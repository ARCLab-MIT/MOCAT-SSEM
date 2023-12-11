%-Abstract
%
%   CSPICE_GFOCLT determines time intervals when an observer sees one target
%   body occulted by, or in transit across, another.
%
%   The surfaces of the target bodies may be represented by triaxial
%   ellipsoids or by topographic data provided by DSK files.
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
%      occtyp   the string naming the type of occultation to find.
%
%               [1,c1] = size(occtyp); char = class(occtyp)
%
%               Note that transits are considered to be a type of
%               occultation.
%
%               Supported values and corresponding definitions are:
%
%                  'FULL'               denotes the full occultation
%                                       of the body designated by
%                                       `back' by the body designated
%                                       by `front', as seen from
%                                       the location of the observer.
%                                       In other words, the occulted
%                                       body is completely invisible
%                                       as seen from the observer's
%                                       location.
%
%                  'ANNULAR'            denotes an annular
%                                       occultation: the body
%                                       designated by `front' blocks
%                                       part of, but not the limb of,
%                                       the body designated by `back',
%                                       as seen from the location of
%                                       the observer.
%
%                  'PARTIAL'            denotes a partial,
%                                       non-annular occultation: the
%                                       body designated by `front'
%                                       blocks part, but not all, of
%                                       the limb of the body
%                                       designated by `back', as seen
%                                       from the location of the
%                                       observer.
%
%                  'ANY'                denotes any of the above three
%                                       types of occultations:
%                                       'PARTIAL', 'ANNULAR', or
%                                       'FULL'.
%
%                                       'ANY' should be used to search
%                                       for times when the body
%                                       designated by `front' blocks
%                                       any part of the body designated
%                                       by `back'.
%
%                                       The option 'ANY' must be used
%                                       if either the front or back
%                                       target body is modeled as
%                                       a point.
%
%               Case and leading or trailing blanks are not significant in
%               the string `occtyp'.
%
%      front    the string naming the target body that occults---that
%               is, passes in front of---the other.
%
%               [1,c2] = size(front); char = class(front)
%
%               Optionally, you may supply the integer NAIF ID code for the
%               body as a string. For example both 'MOON' and '301' are
%               legitimate strings that designate the Moon.
%
%               The `front' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      fshape   the string naming the geometric model used
%               to represent the shape of the front target body.
%
%               [1,c3] = size(fshape); char = class(fshape)
%
%               The supported options are:
%
%                 'ELLIPSOID'     Use a triaxial ellipsoid model
%                                 with radius values provided via the
%                                 kernel pool. A kernel variable
%                                 having a name of the form
%
%                                    'BODYnnn_RADII'
%
%                                 where nnn represents the NAIF
%                                 integer code associated with the
%                                 body, must be present in the kernel
%                                 pool. This variable must be
%                                 associated with three numeric
%                                 values giving the lengths of the
%                                 ellipsoid's X, Y, and Z semi-axes.
%
%                 'POINT'         Treat the body as a single point.
%                                 When a point target is specified,
%                                 the occultation type must be
%                                 set to 'ANY'.
%
%                 'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                     Use topographic data provided by DSK files to
%                     model the body's shape. These data must be
%                     provided by loaded DSK files.
%
%                     The surface list specification is optional. The
%                     syntax of the list is
%
%                        <surface 1> [, <surface 2>...]
%
%                     If present, it indicates that data only for the
%                     listed surfaces are to be used; however, data
%                     need not be available for all surfaces in the
%                     list. If absent, loaded DSK data for any surface
%                     associated with the target body are used.
%
%                     The surface list may contain surface names or
%                     surface ID codes. Names containing blanks must
%                     be delimited by double quotes, for example
%
%                        SURFACES = "Mars MEGDR 128 PIXEL/DEG"
%
%                     If multiple surfaces are specified, their names
%                     or IDs must be separated by commas.
%
%                     See the -Particulars section below for details
%                     concerning use of DSK data.
%
%               The combinations of the shapes of the target bodies
%               `front' and `back' must be one of:
%
%                  One ELLIPSOID, one POINT
%                  Two ELLIPSOIDs
%                  One DSK, one POINT
%
%               Case and leading or trailing blanks are not
%               significant in the string `fshape'.
%
%      fframe   the string naming the body-fixed, body-centered reference
%               frame associated with the front target body.
%
%               [1,c4] = size(fframe); char = class(fframe)
%
%               Examples of such names are 'IAU_SATURN' (for Saturn) and
%               'ITRF93' (for the Earth).
%
%               If the front target body is modeled as a point, `fframe'
%               should be left empty or blank.
%
%               The `fframe' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      back     the string naming the target body that is occulted
%               by---that is, passes in back of---the other.
%
%               [1,c5] = size(back); char = class(back)
%
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               The `back' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      bshape   the string naming the shape specification for the body
%               designated by `back'.
%
%               [1,c6] = size(bshape); char = class(bshape)
%
%               The supported options are those for `fshape'. See the
%               description of `fshape' above for details.
%
%      bframe   the string naming the body-fixed, body-centered
%               reference frame associated with the '`back'' target body.
%
%               [1,c7] = size(bframe); char = class(bframe)
%
%               Examples of such names are 'IAU_SATURN' (for Saturn)
%               and 'ITRF93' (for the Earth).
%
%               If the back target body is modeled as a point, `bframe'
%               should be left empty or blank.
%
%               The `bframe' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      abcorr   the string indicating the aberration corrections to to apply
%               to the state of the target body to account for one-way
%               light time.
%
%               [1,c8] = size(abcorr); char = class(abcorr)
%
%               Stellar aberration corrections are ignored if specified,
%               since these corrections don't improve the accuracy of the
%               occultation determination.
%
%               This routine accepts the same aberration corrections as does
%               the Mice routine cspice_spkezr. See the abcorr.req
%               for a detailed description of the aberration correction
%               options. For convenience, the options are listed below:
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
%      obsrvr   the name of the observing body.
%
%               [1,c9] = size(obsrvr); char = class(obsrvr)
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the observer is
%               Earth.
%
%               Case and leading or trailing blanks are not significant in
%               the string `obsrvr'.
%
%      step     the step size to use in the search.
%
%               [1,1] = size(step); double = class(step)
%
%               `step' must be shorter than any interval, within the
%               confinement window, over which the specified occultation
%               condition is met. In other words, `step' must be shorter
%               than the shortest occultation event the user wishes to
%               detect; `step' must also be shorter than the shortest time
%               interval between two occultation events that occur within
%               the confinement window (see below). However, `step' must not
%               be *too* short, or the search will take an unreasonable
%               amount of time.
%
%               The choice of `step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               `step' has units of TDB seconds.
%
%      cnfine   the SPICE window that confines the time period over which
%               the specified search is conducted.
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
%      nintvls  the maximum number of intervals to return in `result'.
%
%               [1,1] = size(nintvls); int32 = class(nintvls)
%
%               Note: this value should equal at least the number of expected
%               intervals. Recall two double precision values define
%               an interval.
%
%   the call:
%
%      result = cspice_gfoclt( occtyp, front,  fshape, fframe,             ...
%                              back,   bshape, bframe, abcorr,             ...
%                              obsrvr, step,   cnfine, nintvls)
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window `cnfine', on which the specified
%               constraint is satisfied.
%
%               [2n,1] = size(result); double = class(result)
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
%               is the convergence tolerance used for finding
%               endpoints of the intervals comprising the result
%               window. SPICE_GF_CNVTOL is used to determine when
%               binary searches for roots should terminate: when a
%               root is bracketed within an interval of length
%               SPICE_GF_CNVTOL, the root is considered to have been
%               found.
%
%               The accuracy, as opposed to precision, of roots found
%               by this routine depends on the accuracy of the input
%               data. In most cases, the accuracy of solutions will
%               be inferior to their precision.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Find occultations of the Sun by the Moon (that is, solar
%      eclipses)as seen from the center of the Earth over the month
%      December, 2001.
%
%      Use light time corrections to model apparent positions of Sun
%      and Moon. Stellar aberration corrections are not specified
%      because they don't affect occultation computations.
%
%      We select a step size of 3 minutes, which means we
%      ignore occultation events lasting less than 3 minutes,
%      if any exist.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfoclt_ex1.tm
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
%            pck00008.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00008.tpc',
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
%      function gfoclt_ex1()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfoclt_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2001 DEC 01 00:00:00 TDB',                ...
%                               '2002 JAN 01 00:00:00 TDB'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % Select a 3-minute step. We'll ignore any occultations
%         % lasting less than 3 minutes.
%         %
%         step    = 180.;
%
%         occtyp  = 'any';
%         front   = 'moon';
%         fshape  = 'ellipsoid';
%         fframe  = 'iau_moon';
%         back    = 'sun';
%         bshape  = 'ellipsoid';
%         bframe  = 'iau_sun';
%         obsrvr  = 'earth';
%         abcorr  = 'lt';
%
%         result = cspice_gfoclt( occtyp, front,  fshape, fframe,          ...
%                                 back,   bshape, bframe, abcorr,          ...
%                                 obsrvr, step,   cnfine, MAXWIN );
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
%               disp( ' ')
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
%      From : 2001-DEC-14 20:10:14.195952 (TDB)
%      To   : 2001-DEC-14 21:35:50.317994 (TDB)
%
%
%   2) Find occultations of Titan by Saturn or of Saturn by
%      Titan as seen from the center of the Earth over the
%      last three months of 2008. Search for every type
%      of occultation.
%
%      Use light time corrections to model apparent positions of
%      Saturn and Titan. Stellar aberration corrections are not
%      specified because they don't affect occultation computations.
%
%      We select a step size of 15 minutes, which means we
%      ignore occultation events lasting less than 15 minutes,
%      if any exist.
%
%      Use the SPK kernel below for providing the Titan ephemerides
%      and the meta-kernel from example 1 above.
%
%         sat427.bsp
%
%
%      Example code begins here.
%
%
%      function gfoclt_ex2()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%         OCCTYP  = {'FULL', 'ANNULAR', 'PARTIAL', 'ANY' };
%
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfoclt_ex1.tm' );
%         cspice_furnsh( 'sat427.bsp'   );
%
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2008 SEP 01 00:00:00 TDB',                ...
%                               '2009 JAN 01 00:00:00 TDB'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % Select a 15-minute step. We'll ignore any occultations
%         % lasting less than 15 minutes.
%         %
%         step    = 900.;
%
%         %
%         % The observation location is the Earth.
%         %
%         obsrvr  = 'earth';
%         shape   = 'ellipsoid';
%         abcorr = 'lt';
%
%         for i=1:numel(OCCTYP)
%
%            %
%            % For each type, do a search for both transits of
%            % Titan across Saturn and occultations of Titan by
%            % Saturn.
%            %
%            for j=1:2
%
%               if isequal(j,1)
%
%                  front  = 'TITAN';
%                  fframe = 'IAU_TITAN';
%                  back   = 'SATURN';
%                  bframe = 'IAU_SATURN';
%
%               else
%
%                  front  = 'SATURN';
%                  fframe = 'IAU_SATURN';
%                  back   = 'TITAN';
%                  bframe = 'IAU_TITAN';
%
%               end
%
%               result = cspice_gfoclt( OCCTYP(i), front, shape,  fframe,  ...
%                                       back,      shape, bframe, abcorr,  ...
%                                       obsrvr,    step,  cnfine, MAXWIN );
%
%
%               fprintf( 'Condition      : %s\n',   char(OCCTYP(i)) )
%               fprintf( 'Occultation of : %s\n',   back  )
%               fprintf( 'by             : %s\n\n', front )
%
%
%               %
%               % List the beginning and ending times in each interval
%               % if result contains data.
%               %
%
%               count = numel(result)/2;
%
%               if isequal(count,0)
%
%                  fprintf( 'Result window is empty.\n\n' )
%
%               else
%
%                  for k=1:count
%
%                     [left, right] = cspice_wnfetd( result, k );
%
%                     output = cspice_timout( [left,right], TIMFMT );
%
%                     if( isequal( left, right) )
%
%                        disp( ['Event time: ' output(1,:)] )
%
%                     else
%
%                        disp( ['From : ' output(1,:)] )
%                        disp( ['To   : ' output(2,:)] )
%                        disp( ' ')
%
%                     end
%
%                  end
%
%               end
%
%               %
%               % We've finished displaying the results of the
%               % current search.
%               %
%
%            end
%
%            %
%            % We've finished displaying the results of the
%            % searches using the current occultation type.
%            %
%
%         end
%
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
%      Condition      : FULL
%      Occultation of : SATURN
%      by             : TITAN
%
%      Result window is empty.
%
%      Condition      : FULL
%      Occultation of : TITAN
%      by             : SATURN
%
%      From : 2008-OCT-27 22:08:01.672540 (TDB)
%      To   : 2008-OCT-28 01:05:03.332576 (TDB)
%
%      From : 2008-NOV-12 21:21:59.270692 (TDB)
%      To   : 2008-NOV-13 02:06:05.034713 (TDB)
%
%      From : 2008-NOV-28 20:49:02.415745 (TDB)
%      To   : 2008-NOV-29 02:13:58.978005 (TDB)
%
%      From : 2008-DEC-14 20:05:09.258916 (TDB)
%      To   : 2008-DEC-15 01:44:53.517960 (TDB)
%
%      From : 2008-DEC-30 19:00:56.586894 (TDB)
%      To   : 2008-DEC-31 00:42:43.219312 (TDB)
%
%      Condition      : ANNULAR
%      Occultation of : SATURN
%      by             : TITAN
%
%      From : 2008-OCT-19 21:29:20.694709 (TDB)
%      To   : 2008-OCT-19 22:53:34.442728 (TDB)
%
%      From : 2008-NOV-04 20:15:38.652651 (TDB)
%      To   : 2008-NOV-05 00:18:59.130645 (TDB)
%
%      From : 2008-NOV-20 19:38:59.674044 (TDB)
%      To   : 2008-NOV-21 00:35:26.726756 (TDB)
%
%      From : 2008-DEC-06 18:58:34.093679 (TDB)
%      To   : 2008-DEC-07 00:16:17.653067 (TDB)
%
%      From : 2008-DEC-22 18:02:46.308376 (TDB)
%      To   : 2008-DEC-22 23:26:52.721881 (TDB)
%
%      Condition      : ANNULAR
%      Occultation of : TITAN
%      by             : SATURN
%
%      Result window is empty.
%
%      Condition      : PARTIAL
%      Occultation of : SATURN
%      by             : TITAN
%
%      From : 2008-OCT-19 20:44:30.377190 (TDB)
%      To   : 2008-OCT-19 21:29:20.694709 (TDB)
%
%      From : 2008-OCT-19 22:53:34.442728 (TDB)
%      To   : 2008-OCT-19 23:38:26.219866 (TDB)
%
%      From : 2008-NOV-04 19:54:40.368045 (TDB)
%      To   : 2008-NOV-04 20:15:38.652651 (TDB)
%
%      From : 2008-NOV-05 00:18:59.130645 (TDB)
%      To   : 2008-NOV-05 00:39:58.607160 (TDB)
%
%      From : 2008-NOV-20 19:21:46.714397 (TDB)
%      To   : 2008-NOV-20 19:38:59.674044 (TDB)
%
%      From : 2008-NOV-21 00:35:26.726756 (TDB)
%      To   : 2008-NOV-21 00:52:40.606954 (TDB)
%
%      From : 2008-DEC-06 18:42:36.120123 (TDB)
%      To   : 2008-DEC-06 18:58:34.093679 (TDB)
%
%      From : 2008-DEC-07 00:16:17.653067 (TDB)
%      To   : 2008-DEC-07 00:32:16.331200 (TDB)
%
%      From : 2008-DEC-22 17:47:10.796148 (TDB)
%      To   : 2008-DEC-22 18:02:46.308376 (TDB)
%
%      From : 2008-DEC-22 23:26:52.721881 (TDB)
%      To   : 2008-DEC-22 23:42:28.860689 (TDB)
%
%      Condition      : PARTIAL
%      Occultation of : TITAN
%      by             : SATURN
%
%      From : 2008-OCT-27 21:37:17.003994 (TDB)
%      To   : 2008-OCT-27 22:08:01.672540 (TDB)
%
%      From : 2008-OCT-28 01:05:03.332576 (TDB)
%      To   : 2008-OCT-28 01:35:49.235671 (TDB)
%
%      From : 2008-NOV-12 21:01:47.121213 (TDB)
%      To   : 2008-NOV-12 21:21:59.270692 (TDB)
%
%      From : 2008-NOV-13 02:06:05.034713 (TDB)
%      To   : 2008-NOV-13 02:26:18.211754 (TDB)
%
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 156 lines have been
%      provided.
%
%
%-Particulars
%
%   This routine provides a simple interface for conducting searches for
%   occultation events.
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when a specified type of
%   occultation occurs. The resulting set of intervals is returned as
%   a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%
%   The Search Process
%   ==================
%
%   The search for occultations is treated as a search for state
%   transitions: times are sought when the state of the `back' body
%   changes from "not occulted" to "occulted" or vice versa.
%
%   Step Size
%   =========
%
%   Each interval of the confinement window is searched as follows:
%   first, the input step size is used to determine the time separation
%   at which the occultation state will be sampled. Starting at the left
%   endpoint of the interval, samples of the occultation state will be
%   taken at each step. If a state change is detected, a root has been
%   bracketed; at that point, the "root"--the time at which the state
%   change occurs---is found by a refinement process, for example, via
%   binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the occultation state is constant:
%   the step size should be shorter than the shortest occultation
%   duration and the shortest period between occultations, within
%   the confinement window.
%
%   Having some knowledge of the relative geometry of the targets and
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
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie. This
%   refinement process terminates when the location of the root has been
%   determined to within an error margin called the "convergence
%   tolerance." The convergence tolerance used by this routine is set
%   via the parameter SPICE_GF_CNVTOL.
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't limit the accuracy of solutions found by this
%   routine. In general the accuracy of input data will be the limiting
%   factor.
%
%   To use a different tolerance value, a lower-level GF routine such as
%   gfocce_c must be called. Making the tolerance tighter than
%   SPICE_GF_CNVTOL is unlikely to be useful, since the results are
%   unlikely to be more accurate. Making the tolerance looser will speed
%   up searches somewhat, since a few convergence steps will be omitted.
%   However, in most cases, the step size is likely to have a much
%   greater effect on processing time than would the convergence
%   tolerance.
%
%
%   The Confinement Window
%   ======================
%
%   The simplest use of the confinement window is to specify a time
%   interval within which a solution is sought.
%
%   The confinement window also can be used to restrict a search to
%   a time window over which required data (typically ephemeris
%   data, in the case of occultation searches) are known to be
%   available.
%
%   In some cases, the confinement window be used to make searches
%   more efficient. Sometimes it's possible to do an efficient search
%   to reduce the size of the time period over which a relatively
%   slow search of interest must be performed. See the "CASCADE"
%   example program in gf.req for a demonstration.
%
%
%   Using DSK data
%   ==============
%
%      DSK loading and unloading
%      -------------------------
%
%      DSK files providing data used by this routine are loaded by
%      calling cspice_furnsh and can be unloaded by calling cspice_unload or
%      cspice_kclear. See the documentation of cspice_furnsh for limits on
%      numbers of loaded DSK files.
%
%      For run-time efficiency, it's desirable to avoid frequent
%      loading and unloading of DSK files. When there is a reason to
%      use multiple versions of data for a given target body---for
%      example, if topographic data at varying resolutions are to be
%      used---the surface list can be used to select DSK data to be
%      used for a given computation. It is not necessary to unload
%      the data that are not to be used. This recommendation presumes
%      that DSKs containing different versions of surface data for a
%      given body have different surface ID codes.
%
%
%      DSK data priority
%      -----------------
%
%      A DSK coverage overlap occurs when two segments in loaded DSK
%      files cover part or all of the same domain---for example, a
%      given longitude-latitude rectangle---and when the time
%      intervals of the segments overlap as well.
%
%      When DSK data selection is prioritized, in case of a coverage
%      overlap, if the two competing segments are in different DSK
%      files, the segment in the DSK file loaded last takes
%      precedence. If the two segments are in the same file, the
%      segment located closer to the end of the file takes
%      precedence.
%
%      When DSK data selection is unprioritized, data from competing
%      segments are combined. For example, if two competing segments
%      both represent a surface as a set of triangular plates, the
%      union of those sets of plates is considered to represent the
%      surface.
%
%      Currently only unprioritized data selection is supported.
%      Because prioritized data selection may be the default behavior
%      in a later version of the routine, the UNPRIORITIZED keyword is
%      required in the `fshape' and `bshape' arguments.
%
%
%      Syntax of the shape input arguments for the DSK case
%      ----------------------------------------------------
%
%      The keywords and surface list in the target shape arguments
%      `bshape' and `fshape' are called "clauses." The clauses may
%      appear in any order, for example
%
%         "DSK/<surface list>/UNPRIORITIZED"
%         "DSK/UNPRIORITIZED/<surface list>"
%         "UNPRIORITIZED/<surface list>/DSK"
%
%      The simplest form of the `method' argument specifying use of
%      DSK data is one that lacks a surface list, for example:
%
%         "DSK/UNPRIORITIZED"
%
%      For applications in which all loaded DSK data for the target
%      body are for a single surface, and there are no competing
%      segments, the above string suffices. This is expected to be
%      the usual case.
%
%      When, for the specified target body, there are loaded DSK
%      files providing data for multiple surfaces for that body, the
%      surfaces to be used by this routine for a given call must be
%      specified in a surface list, unless data from all of the
%      surfaces are to be used together.
%
%      The surface list consists of the string
%
%         "SURFACES = "
%
%      followed by a comma-separated list of one or more surface
%      identifiers. The identifiers may be names or integer codes in
%      string format. For example, suppose we have the surface
%      names and corresponding ID codes shown below:
%
%         Surface Name                              ID code
%         ------------                              -------
%         "Mars MEGDR 128 PIXEL/DEG"                1
%         "Mars MEGDR 64 PIXEL/DEG"                 2
%         "Mars_MRO_HIRISE"                         3
%
%      If data for all of the above surfaces are loaded, then
%      data for surface 1 can be specified by either
%
%         'SURFACES = 1'
%
%      or
%
%         'SURFACES = "Mars MEGDR 128 PIXEL/DEG"'
%
%      Double quotes are used to delimit the surface name because
%      it contains blank characters.
%
%      To use data for surfaces 2 and 3 together, any
%      of the following surface lists could be used:
%
%         'SURFACES = 2, 3'
%
%         'SURFACES = "Mars MEGDR  64 PIXEL/DEG", 3'
%
%         'SURFACES = 2, Mars_MRO_HIRISE'
%
%         'SURFACES = "Mars MEGDR 64 PIXEL/DEG", Mars_MRO_HIRISE'
%
%      An example of a shape argument that could be constructed
%      using one of the surface lists above is
%
%         'DSK/UNPRIORITIZED/SURFACES = "Mars MEGDR 64 PIXEL/DEG", 3'
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
%   3)  If name of either target or the observer cannot be translated
%       to a NAIF ID code, an error is signaled by a routine
%       in the call tree of this routine.
%
%   4)  If the radii of a target body modeled as an ellipsoid cannot
%       be determined by searching the kernel pool for a kernel
%       variable having a name of the form
%
%          'BODYnnn_RADII'
%
%       where nnn represents the NAIF integer code associated with
%       the body, an error is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If either of the target bodies `front' or `back' coincides with
%       the observer body `obsrvr', an error is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If the body designated by `front' coincides with that
%       designated by `back', an error is signaled by a routine
%       in the call tree of this routine.
%
%   7)  If either of the body model specifiers `fshape' or `bshape'
%       is not recognized, an error is signaled by a routine
%       in the call tree of this routine.
%
%   8)  If both of the body model specifiers `fshape' and `bshape'
%       specify point targets, an error is signaled by a
%       routine in the call tree of this routine.
%
%   9)  If one of the body model specifiers `fshape' and `bshape'
%       specifies a DSK model, and the other argument does not
%       specify a point target, an error is signaled by a routine in
%       the call tree of this routine.
%
%   10) If a target body-fixed reference frame associated with a
%       non-point target is not recognized, an error is signaled by a
%       routine in the call tree of this routine.
%
%   11) If a target body-fixed reference frame is not centered at the
%       corresponding target body, an error is signaled by a routine
%       in the call tree of this routine.
%
%   12) If the loaded kernels provide insufficient data to compute any
%       required state vector, an error is signaled by a routine in
%       the call tree of this routine.
%
%   13) If an error occurs while reading an SPK or other kernel file,
%       the error is signaled by a routine in the call tree
%       of this routine.
%
%   14) If a point target is specified and the occultation type is set
%       to a valid value other than 'ANY', an error is signaled by a
%       routine in the call tree of this routine.
%
%   15) If the output SPICE window `result' has insufficient capacity
%       to contain the number of intervals on which the specified
%       occultation condition is met, an error is signaled
%       by a routine in the call tree of this routine.
%
%   16) If the occultation type `occtyp' is invalid, an error is
%       signaled by a routine in the call tree of this routine.
%
%   17) If the aberration correction specification `abcorr' is invalid,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   18) If either `fshape' or `bshape' specifies that the target surface
%       is represented by DSK data, and no DSK files are loaded for
%       the specified target, an error is signaled by a routine in
%       the call tree of this routine.
%
%   19) If either `fshape' or `bshape' specifies that the target surface
%       is represented by DSK data, but the shape specification is
%       invalid, an error is signaled by a routine in the call tree
%       of this routine.
%
%   20) If any of the input arguments, `occtyp', `front', `fshape',
%       `fframe', `back', `bshape', `bframe', `abcorr', `obsrvr',
%       `step', `cnfine' or `nintvls', is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   21) If any of the input arguments, `occtyp', `front', `fshape',
%       `fframe', `back', `bshape', `bframe', `abcorr', `obsrvr',
%       `step', `cnfine' or `nintvls', is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: the calling application must load ephemeris data
%      for the targets, source and observer that cover the time
%      period specified by the window `cnfine'. If aberration
%      corrections are used, the states of the target bodies and of
%      the observer relative to the solar system barycenter must be
%      calculable from the available ephemeris data. Typically
%      ephemeris data are made available by loading one or more SPK
%      files via cspice_furnsh.
%
%   -  PCK data: bodies modeled as triaxial ellipsoids must have
%      semi-axis lengths provided by variables in the kernel pool.
%      Typically these data are made available by loading a text
%      PCK file via cspice_furnsh.
%
%   -  FK data: if either of the reference frames designated by
%      `bframe' or `fframe' are not built in to the SPICE system,
%      one or more FKs specifying these frames must be loaded.
%
%   The following data may be required:
%
%   -  DSK data: if either `fshape' or `bshape' indicates that DSK
%      data are to be used, DSK files containing topographic data
%      for the target body must be loaded. If a surface list is
%      specified, data for at least one of the listed surfaces must
%      be loaded.
%
%   -  Surface name-ID associations: if surface names are specified
%      in `fshape' or `bshape', the association of these names with
%      their corresponding surface ID codes must be established by
%      assignments of the kernel variables
%
%         NAIF_SURFACE_NAME
%         NAIF_SURFACE_CODE
%         NAIF_SURFACE_BODY
%
%      Normally these associations are made by loading a text
%      kernel containing the necessary assignments. An example
%      of such a set of assignments is
%
%         NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
%         NAIF_SURFACE_CODE += 1
%         NAIF_SURFACE_BODY += 499
%
%   -  CK data: either of the body-fixed frames to which `fframe' or
%      `bframe' refer might be a CK frame. If so, at least one CK
%      file will be needed to permit transformation of vectors
%      between that frame and the J2000 frame.
%
%   -  SCLK data: if a CK file is needed, an associated SCLK
%      kernel is required to enable conversion between encoded SCLK
%      (used to time-tag CK data) and barycentric dynamical time
%      (TDB).
%
%   Kernel data are normally loaded once per program run, NOT every
%   time this routine is called.
%
%-Restrictions
%
%   1)  The kernel files to be used by cspice_gfoclt must be loaded (normally
%       via the Mice routine cspice_furnsh) before cspice_gfoclt is called.
%
%-Required_Reading
%
%   MICE.REQ
%   DSK.REQ
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
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 2.1.0, 25-NOV-2021 (EDW) (JDR) (NJB)
%
%       Changed input argument name "room" to "nintvls" for consistency
%       with other routines.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Edited the header to comply with NAIF standard.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 2.0.0, 04-APR-2017 (EDW) (NJB)
%
%       Header update to reflect support for use of DSKs.
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.1.0, 12-MAY-2012 (EDW)
%
%       Renamed the argument "size" to "room". "size" is a Matlab function
%       name and it's seriously dumb to use a function name word as an
%       argument name.
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.0, 15-APR-2009 (EDW)
%
%-Index_Entries
%
%   GF occultation search
%
%-&

function [result] = cspice_gfoclt( occtyp, front,  fshape, fframe,         ...
                                   back,   bshape, bframe, abcorr,         ...
                                   obsrvr, step,   cnfine, nintvls )

   switch nargin

      case 12

         occtyp  = zzmice_str(occtyp);
         front   = zzmice_str(front);
         fshape  = zzmice_str(fshape);
         fframe  = zzmice_str(fframe);
         back    = zzmice_str(back);
         bshape  = zzmice_str(bshape);
         bframe  = zzmice_str(bframe);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         step    = zzmice_dp(step);
         cnfine  = zzmice_win(cnfine);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );

      otherwise

         error ( [ 'Usage: [result] = cspice_gfoclt( `occtyp`, '           ...
                                           '`front`, `fshape`, `fframe`, ' ...
                                           '`back`, `bshape`, `bframe`, '  ...
                                           '`abcorr`, `obsrvr`, step, '    ...
                                           'cnfine, nintvls )' ] )

   end

   %
   % Call the GF routine.  Add to `cnfine' the 6x1 space needed for
   % the control segment.
   %
   try

      [result] = mice('gfoclt_c',  occtyp, front,  fshape, fframe, ...
                                   back,   bshape, bframe,         ...
                                   abcorr, obsrvr,  step,          ...
                                   [zeros(6,1); cnfine], nintvls);
   catch spiceerr
      rethrow(spiceerr)
   end
