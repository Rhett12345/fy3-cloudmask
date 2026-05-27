!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime polar land processing path.
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

! ... Daytime polar land thresholds

      real pdlvrat(4)   
      real pdlref1(4)   
      real pdlref3(4)   
      real pdl11_4lo(4) 
      real pdl11_12hi(1)
      real pdlh20(4)    
      real pdltci(2)

      common /dpland_thresholds/ pdlvrat,pdlref1,pdlref3,pdl11_4lo,  &
        pdl11_12hi,pdlh20,pdltci
