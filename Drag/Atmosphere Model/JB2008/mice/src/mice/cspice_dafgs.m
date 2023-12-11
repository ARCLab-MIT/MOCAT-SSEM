%-Abstract
%
%   CSPICE_DAFGS returns the summary of the current DAF array,
%   i.e. the array found by the previous call to cspice_daffna
%   or cspice_daffpa.
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
%      nd       a scalar integer defining the size of the return
%               double precision array.
%
%               [1,1] = size(nd); int32 = class(nd)
%
%      ni       a scalar integer defining the size of the return
%               integer array.
%
%               [1,1] = size(ni); int32 = class(ni)
%
%               For an SPK file, `nd' always equals 2, `ni' always equals 6.
%               The precise contents of the vectors depend on the type of DAF
%               but the final two elements of the `ic' (integer) vector always
%               contains the initial and final addresses respectively of the
%               array.
%
%   the call:
%
%      [dc, ic] = cspice_dafgs( nd, ni )
%
%   returns:
%
%      dc       the array of double precision components of the summary.
%
%               [1,nd] = size(dc); double = class(dc)
%
%      ic       the array of integer components of the summary.
%
%               [1,ni] = size(ic); int32 = class(ic)
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
%   1) Create a simple program to output the double precision and integer
%      values stored in an SPK's segments' descriptors. This program opens
%      a DAF for read, performs a forward search for the DAF arrays,
%      prints the segment descriptor for each array found, then closes
%      the DAF.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         de421.bsp
%
%
%      Example code begins here.
%
%
%      function dafgs_ex1()
%
%         %
%         % Local constants.
%         %
%         kernel = 'de421.bsp';
%
%         %
%         % Open a DAF for read. Return a 'handle' referring
%         % to the file.
%         %
%         handle = cspice_dafopr( kernel );
%
%         %
%         % Define the summary parameters appropriate
%         % for an SPK file.
%         %
%         ND = 2;
%         NI = 6;
%
%         %
%         % Begin a forward search on the file.
%         %
%         cspice_dafbfs( handle );
%
%         %
%         % Search until a DAF array is found.
%         %
%         found = cspice_daffna;
%
%         %
%         % Loop while the search finds subsequent DAF arrays.
%         %
%         while found
%
%            [dc, ic ] = cspice_dafgs( ND, NI );
%
%            fprintf( 'Doubles:  ' )
%            fprintf( '%f   ', dc )
%            fprintf( '\n' )
%
%            fprintf( 'Integers: ' )
%            fprintf( '%d   ', ic )
%            fprintf( '\n\n' )
%
%
%            %
%            % Check for another segment.
%            %
%            found = cspice_daffna;
%
%         end
%
%         %
%         % Safely close the DAF.
%         %
%         cspice_dafcls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 1   0   1   2   641   310404
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 2   0   1   2   310405   423048
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 3   0   1   2   423049   567372
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 4   0   1   2   567373   628976
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 5   0   1   2   628977   674740
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 6   0   1   2   674741   715224
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 7   0   1   2   715225   750428
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 8   0   1   2   750429   785632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 9   0   1   2   785633   820836
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 10   0   1   2   820837   944040
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 301   3   1   2   944041   1521324
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 399   3   1   2   1521325   2098608
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 199   1   1   2   2098609   2098620
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 299   2   1   2   2098621   2098632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 499   4   1   2   2098633   2098644
%
%
%      Note, the specific contents of `ic' and `dc' depend on the
%      type of DAF.
%
%      Note, the final entries in the integer array contain the segment
%      start/end indexes. The output indicates the search proceeded
%      from the start of the file (low value index) towards the end
%      (high value index).
%
%-Particulars
%
%   A single call to cspice_dafgs equates to the CSPICE calls:
%
%      dafgs_c( sum );
%      dafus_c( sum, nd, ni, dc, ic );
%
%   without use of the `sum' variable.
%
%-Exceptions
%
%   1)  If this routine is called when no search is in progress in the
%       the current DAF, the error SPICE(DAFNOSEARCH) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the DAF for which the "current" array's summary is to be
%       returned has actually been closed, an error is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If no array is current in the current DAF, the error
%       SPICE(NOCURRENTARRAY) is signaled by a routine in the call tree
%       of this routine. There is no current array when a search is
%       started by cspice_dafbfs or cspice_dafbbs, but no calls to
%       cspice_daffna or cspice_daffpa have been made yet, or whenever
%       cspice_daffna or cspice_daffpa return the value False.
%
%   4)  If `nd' is zero or negative, no double precision components
%       are returned.
%
%   5)  If `ni' is zero or negative, no integer components are returned.
%
%   6)  If any of the input arguments, `nd' or `ni', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   7)  If any of the input arguments, `nd' or `ni', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
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
%   DAF.REQ
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
%   -Mice Version 1.1.0, 09-JUL-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Modified code example to hardcode the input DAF file.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 21-MAR-2014 (EDW)
%
%-Index_Entries
%
%   get DAF summary
%
%-&

function [dc, ic] = cspice_dafgs( nd, ni)

   switch nargin
      case 2

         nd  = zzmice_int(nd);
         ni  = zzmice_int(ni);

      otherwise

         error ( 'Usage: [dc, ic] = cspice_dafgs( nd, ni)' )

   end

   %
   % Call the MEX library.
   %
   try
      [dc, ic] = mice( 'dafgs_c', nd, ni );
   catch spiceerr
      rethrow(spiceerr)
   end

