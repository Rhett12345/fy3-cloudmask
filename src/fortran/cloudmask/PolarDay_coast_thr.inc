!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime polar coastal land processing path.
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

!     Daytime polar coastal land thresholds.

      real pdlref1_t2(4)   
      real pdlref3_t2(4)   
      real pdl11_4lo_t2(4) 
      real pdl11_12hi_t2(1)
      real pdlh20_t2(4)    
      real pdltci_t2(2)

      common /dpcoast_thresholds/ pdlref1_t2,pdlref3_t2,   &
         pdl11_4lo_t2,pdl11_12hi_t2,pdlh20_t2,pdltci_t2
