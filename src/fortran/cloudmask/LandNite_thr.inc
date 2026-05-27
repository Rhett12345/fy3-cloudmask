!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime land processing path.
!
!!Revision History:
! Added 'nl4_12lo'
! 06/04 Collection 5   R. Frey
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

! ... Nighttime land thresholds

      real nl4_12hi(4)
      real nl4_12lo(4)
      real nlco2(4)   
      real nlh20(4)   
      real nl11_12hi(4)
      real nl7_11s(4)
      real nl_11_4l(4)
      real nl_11_4h(4)
      real nl_11_4m(4)
      real bt_diff_bounds(2)

      common /niteland_thresholds/ nl4_12hi,nlco2,nlh20,      &
         nl11_12hi,nl_11_4l,nl_11_4h,bt_diff_bounds,nl7_11s,  &
         nl_11_4m,nl4_12lo
