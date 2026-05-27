!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime polar ocean processing path.
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

!     Daytime polar ocean thresholds

      real pdobt11(4)   
      real pdo11_12hi(1)
      real pdovrathi(4) 
      real pdovratlo(4) 
      real pdoref2(4)   
      real pdoref3(4)   
      real pdo11_4lo(4) 
      real pdoh20(4)    
      real pdotci(2)

      common /dpocean_thresholds/ pdobt11,pdo11_12hi,pdovrathi,  &
        pdovratlo,pdoref2,pdoref3,pdo11_4lo,pdoh20,pdotci
