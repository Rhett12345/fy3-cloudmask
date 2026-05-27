!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime polar snow processing path.
!
!!Revision History:
! 06/04 Collection 5  R. Frey:
! Added 'dps11_12hi', 'dfs11_12adj'
! 10/04 Collection 5b R. Frey:
! Added 'dps4_11' variables.
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

      real dpsref1(4)
      real dpsref3(4)
      real dpsh20(4) 
      real dpstci(2)
      real dps11_12hi(1)
      real dps11_12adj(1)
      real dps4_11l(4)
      real dps4_11h(4)
      real dps4_11m1(4)
      real dps4_11m2(4)
      real dps4_11m3(4)
      real bt_11_bnds3(4)

      common /dpsnow_thresholds/ dpsref1,dpsref3,           &
        dpsh20,dpstci,dps11_12hi,dps11_12adj,dps4_11l,      &
        dps4_11h,dps4_11m1,dps4_11m2,dps4_11m3,bt_11_bnds3
