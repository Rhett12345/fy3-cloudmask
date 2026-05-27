!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime ocean processing path.
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Added 'no86_73' and 'no_11var', removed 'nosig'.
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

!     Nighttime ocean thresholds

      real nobt11(4)   
      real no11_12hi(4)
      real no11_4lo(4) 
      real noco2(4)    
      real noh20(4)    
      real no86_73(4)
      real no_11var(4)

      common /niteocean_thresholds/ nobt11,no11_12hi,no11_4lo,  &
        noco2,noh20,no86_73,no_11var
