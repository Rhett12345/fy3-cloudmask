!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the snow detection subroutine
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

! ... Thresholds

      real sm_bt11(1)
      real sm_ref2(1)
      real sm_ref3(1)
      real sm_co2(1)
      real sm85_11(1)
      real sm37_11(1)
      real sm37_11hel(1)
      real sm_ndsi(1)
      real sm_mnir(1)

      common /snow_thresholds/ sm_bt11,sm_ref2,sm_ref3,     &
        sm_co2,sm85_11,sm37_11,sm37_11hel,sm_ndsi,sm_mnir
