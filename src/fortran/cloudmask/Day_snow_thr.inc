!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime snow processing path.
!
!!Revision History:
! Added 'ds11_12hi','ds11_12adj'
! 06/04 Collection 5   R. Frey
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

!     Daytime snow thresholds

      real dsco2(4) 
      real dsref3(4)
      real ds4_11(4)
      real ds4_11hel(4)
      real dsh20(4) 
      real dstci(2)
      real ds11_12hi(1)
      real ds11_12adj(1)

      common /daysnow_thresholds/ dsco2,dsref3,ds4_11,dsh20,   &
              ds4_11hel,dstci,ds11_12hi,ds11_12adj
