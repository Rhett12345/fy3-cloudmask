!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the non-cloud obstruction detection subroutine
!
!!Revision History:
!
! Added 'ncrat', 'ncvrat', 'ncsig'.
! 08/20/03 raf
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!
!!END******************************************************

! ... Thresholds

      real nc_bt37(1)
      real nc37_11(1)
      real nc21(1)
      real nc11_12(1)
      real ncrat(1)
      real ncvrat(1)
      real ncsig(1)

      common /ncobs_thresholds/ nc_bt37,nc37_11,nc21,nc11_12,  &
                               ncrat,ncvrat,ncsig
