!!F77-IN!*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime desert processing path.
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

!     Daytime desert thresholds.

      real lds11_12hi(4)
      real ldsref2(4)   
      real ldsref3(4)   
      real lds11_4lo(4) 
      real lds11_4hi(4) 
      real ldsco2(4)    
      real ldsh20(4)
      real ldstci(2)

      common /desertday_thresholds/ lds11_12hi,ldsref2,ldsref3,   &
          lds11_4lo,lds11_4hi,ldsco2,ldsh20,ldstci
