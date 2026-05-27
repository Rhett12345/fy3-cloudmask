!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!               the daytime coastal desert processing path
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

! ... Day time coastal desert thresholds.

      real lds11_12hi_c(4)
      real ldsref2_c(4)   
      real ldsref3_c(4)   
      real lds11_4lo_c(4) 
      real lds11_4hi_c(4) 
      real ldsco2_c(4)    
      real ldsh20_c(4)
      real ldstci_c(2)

      common /desertcoast_thresholds/ lds11_12hi_c,ldsref2_c,  &
         ldsref3_c,lds11_4lo_c,lds11_4hi_c,ldsco2_c,ldsh20_c,  &
         ldstci_c
