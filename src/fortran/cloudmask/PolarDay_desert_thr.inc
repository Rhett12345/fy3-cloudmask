!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the daytime polar desert processing path
!
!!Revision History:
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

! ... Daytime polar desert thresholds.

      real pds11_12hi(4)
      real pdsref2(4)   
      real pdsref3(4)   
      real pds11_4lo(4) 
      real pds11_4hi(4) 
      real pdsh20(4)
      real pdstci(2)

      common /dpdesert_thresholds/ pds11_12hi,pdsref2,pdsref3,   &
         pds11_4lo,pds11_4hi,pdsh20,pdstci
