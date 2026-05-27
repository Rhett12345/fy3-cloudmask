      integer function strlen(str)

      implicit none

C-----------------------------------------------------------------------
C !F77
C
C !Description: This function determines the byte position of the last
C               non-blank character in a string buffer.  This position 
C               is referred to as the string length and it is returned 
C               as the function value.  If the string buffer contains 
C               all blank characters, a string length of zero is 
C               returned. 
C
C !Input Parameters:
C     character*(*) str           The string buffer
C
C !Output Parameters: None
C
C !Revision History:
C
C
c Revision 1.1  1997/11/05  23:34:47  rhucek
c Initial revision
c
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C !Internals:
C
C    integer        buflen            The length of the string buffer
C    integer        ipos              Byte location of string character
C     
C !END
C----------------------------------------------------------------------------

C-----Declaration of PARAMETERs 
      character*1 BLANK
      PARAMETER (BLANK = ' ')

C-----Declaration of function arguments
      character*(*) str

C-----Declaration of local variables
      integer buflen,ipos

      
      If (str .EQ. BLANK) Then
         strlen = 0
      Else

c----initialize position of last non-blank character to end of buffer
         buflen = len(str)
         ipos = buflen

         Do While( str(ipos:ipos).EQ.BLANK .AND. ipos.GT.1 )
            ipos = ipos - 1   
         End Do 
              
         strlen = ipos
      EndIf

      Return
      End 
