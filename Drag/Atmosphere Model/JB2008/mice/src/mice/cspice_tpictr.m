%-Abstract
%
%   CSPICE_TPICTR creates a time format picture suitable for use by the
%   routine cspice_timout from a given sample time string.
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
%      sample   a representative time string to use as a model to
%               format time strings (the string need not describe
%               an actual date - only format matters).
%
%               [1,c1] = size(sample); char = class(sample)
%
%                  or
%
%               [1,1] = size(sample); cell = class(sample)
%
%   the call:
%
%      [pictur, ok, errmsg] = cspice_tpictr(sample)
%
%   returns:
%
%      pictur   a format picture string suitable for use with the
%               SPICE routine cspice_timout.
%
%               [1,c2] = size(pictur); char = class(pictur)
%
%      ok       a boolean indicating whether `sample' parsed
%               without error, true, or some parse error occurred, false.
%
%               [1,1] = size(ok); logical = class(ok)
%
%      errmsg   a string containing the explanation of
%               the parse error.
%
%               [1,c3] = size(errmsg); char = class(errmsg)
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Given a sample with the format of the UNIX date string,
%      create a SPICE time picture for use in cspice_timout.
%
%
%      Example code begins here.
%
%
%      function tpictr_ex1()
%
%         sample = 'Thu Oct 1 11:11:11 PDT 1111';
%
%         %
%         % Make the call. 'ok' returns false is an error occurred.
%         % The error description returns in the err variable.
%         %
%         [pictur, ok, errmsg] = cspice_tpictr( sample );
%
%         %
%         % If a false error flag, print the picture; if
%         % a true error flag, print the error message.
%         %
%         if ( ok )
%            disp( pictur )
%         else
%            disp( errmsg )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7
%
%
%-Particulars
%
%   Although the routine cspice_timout provides SPICE users with a great
%   deal of flexibility in formatting time strings, users must
%   master the means by which a time picture is constructed
%   suitable for use by cspice_timout.
%
%   This routine allows SPICE users to supply a sample time string
%   from which a corresponding time format picture can be created,
%   freeing users from the task of mastering the intricacies of
%   the routine cspice_timout.
%
%   Note that cspice_timout can produce many time strings whose patterns
%   can not be discerned by this routine. When such outputs are
%   called for, the user must consult cspice_timout and construct the
%   appropriate format picture "by hand." However, these exceptional
%   formats are not widely used and are not generally recognizable
%   to an uninitiated reader.
%
%-Exceptions
%
%   1)  All problems with the inputs are reported via `ok' and `errmsg'.
%
%   2)  If a format picture can not be created from the sample
%       time string, `pictur' is returned as a blank string.
%
%   3)  If the input argument `sample' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `sample' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   TIME.REQ
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
%   -Mice Version 1.2.0, 25-AUG-2021 (EDW) (JDR)
%
%       Changed the output argument name "error" to "errmsg" for
%       consistency with other routines.
%
%       Edited the header to comply with NAIF standard.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.1, 09-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.1.0, 10-MAY-2011 (EDW)
%
%       "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.1, 31-MAR-2010 (EDW)
%
%       Renamed error message argument 'error' to 'errmsg'.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   Use a sample time string to produce a time format picture
%
%-&

function [pictur, ok, errmsg] = cspice_tpictr(sample)

   switch nargin
      case 1

         sample = zzmice_str(sample);

      otherwise

         error ( [ 'Usage: [`pictur`, ok, `errmsg`] = ', ...
                                        'cspice_tpictr( `sample` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [pictur, ok, errmsg] =  mice('tpictr_c', sample );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      ok = zzmice_logical(ok);
   catch spiceerr
      rethrow(spiceerr)
   end


