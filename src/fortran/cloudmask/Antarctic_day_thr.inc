!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime Antarctica processing path.
!
!!Revision History:
!
! Added 'ant11_12hi', 'ant11_12adj'
! 06/04 Collection 5   R. Frey
! Added 'ant4_11' variables
! 10/04 Collection 5b  R. Frey
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

! ... Daytime polar snow thresholds

      real anth20(4) 
      real ant4_11l(4)
      real ant4_11h(4)
      real ant4_11m1(4)
      real ant4_11m2(4)
      real ant4_11m3(4)
      real bt_11_bnds4(4)

      common /antarctic_thresholds/ anth20,             &
              ant4_11l,ant4_11h,ant4_11m1,ant4_11m2,    &
              ant4_11m3,bt_11_bnds4
