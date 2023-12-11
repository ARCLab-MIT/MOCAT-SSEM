%-Abstract
%
%   CSPICE_BODS2C translates a string containing a body name or ID code to the
%   corresponding integer code.
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
%      name     containing the name or ID code of a body or  object, such as
%               a planet, satellite, comet, asteroid, barycenter, DSN station,
%               spacecraft, or instrument.
%
%               [n,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,n] = size(name); cell = class(name)
%
%               If `name' contains the name of a body or object, that name
%               must be "known" to the SPICE system, whether through
%               hard-coded registration or run-time registration in the SPICE
%               kernel pool.
%
%               Case and leading and trailing blanks in `name' are not
%               significant. However when a name is made up of more than one
%               word, they must be separated by at least one blank. That is,
%               all of the following strings are equivalent names:
%
%                  'JUPITER BARYCENTER'
%                  'Jupiter Barycenter'
%                  'JUPITER BARYCENTER   '
%                  'JUPITER    BARYCENTER'
%                  '   JUPITER BARYCENTER'
%
%               However, 'JUPITERBARYCENTER' is not equivalent to the names
%               above.
%
%               If `name' is a string representation of an integer, for
%               example
%
%                  '399'
%
%               the string will be translated to the equivalent integer datum.
%               The input integer need not be one recognized by the SPICE
%               system: the integer need not be a built-in NAIF ID code, nor
%               need it be associated with a name via run-time registration.
%
%   the call:
%
%      [code, found] = cspice_bods2c( name )
%
%   returns:
%
%      code     containing the SPICE code(s) for `name' if `name' contains
%               the name of a body or object as determined by the SPICE
%               name-code mapping subsystem.
%
%               [1,n] = size(code); int32 = class(code)
%
%               If the input argument `name' represents an integer, the same
%               integer is returned in `code'. If neither mapping exists
%               `code' returns as 0.
%
%      found    true if `name' has a translation or represents an integer
%               within the bounds of representable integers as defined by the
%               Mice routines cspice_intmax and cspice_intmin.
%
%               [1,n] = size(found); logical = class(found)
%
%               `found' and `code' return with the same vectorization
%               measure, N, as `name'.
%
%-Parameters
%
%   MAXL        is the maximum allowable length of a body name. The
%               current value of this parameter is 36.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Apply the cspice_bods2c call to several body names to retrieve
%      their associated NAIF IDs included in the default SPICE ID-name
%      lists, a name not included in that list, a string representing a
%      positive integer and another representing a negative integer.
%
%      Example code begins here.
%
%
%      function bods2c_ex1()
%         %
%         % Retrieve the NAIF ID associated to a body name.
%         %
%         disp( 'Scalar:' )
%         name         = 'Hyperion';
%         [code,found] = cspice_bods2c( name );
%
%         %
%         % Output the mapping if it exists.
%         %
%         if ( found )
%            txt = sprintf( 'String %s maps to ID %i', ...
%                            name, code );
%            disp(txt)
%         end
%
%         disp(' ')
%
%         %
%         % Create an array of strings. Include one string not an integer
%         % and unknown to the SPICE system.
%         %
%         disp( 'Vector:' )
%         name          = strvcat( 'Cassini'   , '399',  ...
%                                  'Planet Bob', 'MARS', ...
%                                  '-123456'   , '987654' );
%         [code, found] = cspice_bods2c( name );
%
%         n_elements = size(code,2);
%
%         %
%         % Loop over the output array.
%         %
%         for n=1:n_elements
%
%            %
%            % Check for a valid name/ID mapping.
%            %
%            if( found(n))
%               txt = sprintf( 'String %s maps to ID %i', ...
%                               name(n,:), code(n) );
%               disp(txt)
%            else
%               txt = sprintf( 'Unknown string ID %s', name(n,:) );
%               disp(txt)
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Scalar:
%      String Hyperion maps to ID 607
%
%      Vector:
%      String Cassini    maps to ID -82
%      String 399        maps to ID 399
%      Unknown string ID Planet Bob
%      String MARS       maps to ID 499
%      String -123456    maps to ID -123456
%      String 987654     maps to ID 987654
%
%
%-Particulars
%
%   A sister version of this routine exists named mice_bods2c that returns
%   the output arguments as fields in a single structure.
%
%   cspice_bods2c is one of five related subroutines,
%
%      cspice_bods2c      Body string to code
%      cspice_bodc2s      Body code to string
%      cspice_bodn2c      Body name to code
%      cspice_bodc2n      Body code to name
%      cspice_boddef      Body name/code definition
%
%   cspice_bods2c, cspice_bodc2s, cspice_bodn2c, and cspice_bodc2n
%   perform translations between body names and their corresponding
%   integer ID codes which are used in SPICE files and routines.
%
%   cspice_bods2c is a slightly more general version of cspice_bodn2c:
%   support for strings containing ID codes in string format enables a caller
%   to identify a body using a string, even when no name is associated with
%   that body.
%
%   cspice_bodc2s is a general version of cspice_bodc2n; the routine returns
%   either the name assigned in the body ID to name mapping or a string
%   representation of the `code' value if no mapping exists.
%
%   cspice_boddef assigns a body name to ID mapping. The mapping has
%   priority in name-to-ID and ID-to-name translations.
%
%   Refer to naif_ids.req for the list of name/code associations built
%   into SPICE, and for details concerning adding new name/code
%   associations at run time by loading text kernels.
%
%-Exceptions
%
%   1)  If there is any problem with the body name-ID mapping kernel
%       variables present in the kernel pool, an error is signaled by
%       a routine in the call tree of this routine.
%
%   2)  Body name strings are upper-cased, their leading and trailing
%       blanks removed, and embedded blanks are compressed out, after
%       which they get truncated to the maximum body name length MAXL.
%       Therefore, two body names that differ only after that maximum
%       length are considered equal.
%
%   3)  If the input argument `name' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `name' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Body-name mappings may be defined at run time by loading text
%   kernels containing kernel variable assignments of the form
%
%      NAIF_BODY_NAME += ( <name 1>, ... )
%      NAIF_BODY_CODE += ( <code 1>, ... )
%
%   See naif_ids.req for details.
%
%-Restrictions
%
%   1)  See exception <2>.
%
%-Required_Reading
%
%   MICE.REQ
%   NAIF_IDS.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Updated `found' output
%       argument description in -I/O.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 12-MAR-2012 (EDW) (SCK)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.1, 16-MAY-2009 (EDW)
%
%       Edit to -Particulars section to document the cspice_bodc2s routine.
%       Extended argument descriptions in the -I/O section.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   body name to code
%
%-&

function [code,found] = cspice_bods2c(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_code_, _found_] = cspice_bods2c(_`name`_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [bods2c] = mice('bods2c_s',name);
      code     = reshape( [bods2c.code], 1, [] );
      found    = reshape( [bods2c.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end




