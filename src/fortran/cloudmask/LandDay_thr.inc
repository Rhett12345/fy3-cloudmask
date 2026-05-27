!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime land processing path.
!
!!Revision History:
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!     The four parameters represent:
!        1 - the low clear confidence threshold (0 confidence)
!        2 - the actual straight threshold value (50% confidence)
!            (pass/no pass)
!        3 - the high clear confidence threshold (100% confidence)
!        4 - the s curve power used to interpolate confidences
!
!!END******************************************************

! ... Daytime land thresholds

      real dlvrat(4)   
      real dlref1(4)   
      real dlref3(4)   
      real dl11_4lo(4) 
      real dl11_12hi(1)
      real dlco2(4)    
      real dlh20(4)    
      real dltci(2)

      common /dayland_thresholds/ dlvrat,dlref1,dlref3,dl11_4lo, &
         dl11_12hi,dlco2,dlh20,dltci
