!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime coastal land processing path.
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

! ... Daytime coastal land thresholds

      real dlref1_t2(4)   
      real dlref3_t2(4)   
      real dl11_4lo_t2(4) 
      real dl11_12hi_t2(1)
      real dlco2_t2(4)    
      real dlh20_t2(4)    
      real dltci_t2(2)

      common /daycoast_thresholds/ dlref1_t2,dlref3_t2,dl11_4lo_t2,  &
             dl11_12hi_t2,dlco2_t2,dlh20_t2,dltci_t2
