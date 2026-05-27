      SUBROUTINE STRLOWER( STRING )

C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Convert all characters in the range [A-Z] to lowercase.
C
C !INPUT PARAMETERS:
C     STRING        On input, string variable.
C
C !OUTPUT PARAMETERS:
C     STRING        On output, string with all characters in the
C                   range [A-Z] converted to lowercase.
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C     Original version by Liam.Gumley@ssec.wisc.edu
C
C !END
C--------------------------------------------------------------------

      IMPLICIT NONE

c ... Input arguments

      CHARACTER*(*) string

c ... Local variables

      INTEGER i, code

c ... Convert characters in the range [A-Z] to lower case

      do i = 1, len( string )
        code = ichar( string( i : i ) )
        if ( code .ge. 65 .and. code .le. 90 ) then
          string( i : i ) = char( code + 32 )
        endif
      end do

      END
