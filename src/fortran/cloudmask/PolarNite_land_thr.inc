!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime polar land processing path.
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Removed 'pnl4_12hi'.
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

!     Nighttime polar land thresholds

      real pnlh20(4)   
      real pnl11_12hi(4)

      common /npland_thresholds/ pnlh20, pnl11_12hi
