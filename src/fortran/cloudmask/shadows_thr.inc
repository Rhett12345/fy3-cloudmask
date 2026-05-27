!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!               the shadow detection subroutine
!
!!Revision History:
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!
!!END******************************************************

! ... Shadow threshold

      real shavrat(1)
      real shadnir(2)
      real shad124(1)

      common /shadow_thresholds/ shavrat,shadnir,shad124
