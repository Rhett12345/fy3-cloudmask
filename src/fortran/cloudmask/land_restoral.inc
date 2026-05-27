!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the land_restoral processing (chk_land.f).
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

!     Land restoral thresholds.

      real ldsr5_4_thr(1)
      real ldr5_4_thr(1)
      real ld20m22(1)
      real ld22m31(1)
      real ldsbt11(3)
      real ldsbt11bd(3)
      real lnbt11(3)

      common /land_restoral_thresholds/ ldsr5_4_thr,ldr5_4_thr, &
             ld20m22,ld22m31,ldsbt11,ldsbt11bd,lnbt11
