!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime ocean processing path.
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

!     Daytime ocean thresholds

      real dobt11(4)   
      real do11_12hi(1)
      real dovrathi(4) 
      real dovratlo(4) 
      real doref2(4)   
      real doref3(4)   
      real do11_4lo(4) 
      real doh20(4)    
      real doco2(4)    
      real dotci(2)

      common /dayocean_thresholds/ dobt11,do11_12hi,dovrathi,   &
        dovratlo,doref2,doref3,do11_4lo,doh20,doco2,dotci
