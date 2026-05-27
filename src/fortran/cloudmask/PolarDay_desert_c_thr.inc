!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!               the daytime polar coastal desert processing path
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

! ... Daytime polar coastal desert thresholds.

      real pds11_12hi_c(4)
      real pdsref2_c(4)   
      real pdsref3_c(4)   
      real pds11_4lo_c(4) 
      real pds11_4hi_c(4) 
      real pdsh20_c(4)
      real pdstci_c(2)

      common /dpdesertc_thresholds/ pds11_12hi_c,pdsref2_c,pdsref3_c,   &
          pds11_4lo_c,pds11_4hi_c,pdsh20_c,pdstci_c
