!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime polar snow processing path.
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Added threshold variables for 11-12 um thin cirrus test,
! 11-3.9 um, 7.3-11 cloud tests
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

!     Nighttime polar snow thresholds

      real pnsh20(4)   
      real pn65_11(1)
      real pn13_11(1)
      real pns11_12hi(1)
      real pn7_11(1)
      real pn_11_4l(4)
      real pn_11_4h(4)
      real bt_11_bounds(4)
      real bt_11_bnds2(4)
      real pn_11_4m1(4)
      real pn_11_4m2(4)
      real pn_11_4m3(4)
      real pn_4_12l(4)
      real pn_4_12h(4)
      real pn_4_12m1(4)
      real pn_4_12m2(4)
      real pn_4_12m3(4)
      real pn_7_11l(4)
      real pn_7_11h(4)
      real pn_7_11m1(4)
      real pn_7_11m2(4)
      real pn_7_11m3(4)
      real pn_7_11lw(4)
      real pn_7_11hw(4)
      real pn_7_11m1w(4)
      real pn_7_11m2w(4)
      real pn_7_11m3w(4)
      real pn11_12adj(1)

      common /pnsnow_thresholds/ pnsh20,                                &
         pn65_11,pn13_11,pn_11_4l,pn_11_4h,pn_11_4m1,bt_11_bounds,      &
         pn_11_4m2,pn_4_12l,pn_4_12h,pn_4_12m1,pn_4_12m2,pn_4_12m3,     &
         pn_11_4m3,pn_7_11l,pn_7_11h,pn_7_11m1,pn_7_11m2,               &
         pn_7_11m3,bt_11_bnds2,pn7_11,pns11_12hi,                       &
         pn_7_11lw,pn_7_11hw,pn_7_11m1w,pn_7_11m2w,pn_7_11m3w,          &
         pn11_12adj
