!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime polar ocean processing path.
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Added 'pno86_73' and 'pno_11var', removed 'pnosig'.
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
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

!     Nighttime polar ocean thresholds

      real pnobt11(4)   
      real pno11_12hi(4)
      real pno11_4lo(4) 
      real pnoh20(4)    
      real pno86_73(4)
      real pno_11var(4)

      common /pnocean_thresholds/ pnobt11,pno11_12hi,pno11_4lo,  &
         pnoh20,pno86_73,pno_11var
