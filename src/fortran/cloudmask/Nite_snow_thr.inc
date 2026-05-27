!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime snow processing path.
!
!!Revision History:
! Added 'ns11_12hi','ns11_12adj'
! 06/04 Collection 5    R. Frey
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

!     Nighttime snow thresholds

      real ns11_12hi(1)
      real ns4_12hi(4)
      real ns11_4lo(4)
      real nsco2(4)   
      real nsh20(4)   
      real n65_11(1)  
      real ns11_12adj(1)

      common /nitesnow_thresholds/ ns4_12hi,ns11_4lo,nsco2,nsh20,  &
         n65_11,ns11_12hi,ns11_12adj
