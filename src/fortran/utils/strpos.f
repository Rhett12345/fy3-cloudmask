      INTEGER FUNCTION STRPOS( STRING, SUBSTRING )
      
C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Find a substring within a string.
C
C !INPUT PARAMETERS:
C     STRING        String variable
C     SUBSTRING     The substring to be searched for within STRING
C 
C !OUTPUT PARAMETERS:
C     STRPOS        If SUBSTRING occurs in STRING, the character
C                   position of the match. Otherwise STRPOS=-1.
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

!      IMPLICIT NONE

c ... Arguments

      CHARACTER*(*) string, substring

c ... Local variables
      
      INTEGER string_len, substring_len
      INTEGER i

c ... Set return value

      strpos = -1

c ... Get string lengths

      string_len = len( string )
      substring_len = len( substring )
      
c ... If string or substring is empty, return

      if ( string_len .eq. 0 .or. substring_len .eq. 0 ) return

c ... If substring is longer than string, return
      
      if ( substring_len .gt. string_len ) return
      
c ... Search string for substring

      do i = 1, string_len - substring_len + 1
        if( string( i : i + substring_len - 1 ) .eq. substring ) then
          strpos = i
          return
        endif
      end do
            
      END
