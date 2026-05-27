!!F77-INC*******************************************************
!
!!Description: Include file for initializing thresholds used in
!              the nighttime snow processing path.
!
!!Revision History:
! Added 'ns11_12hi','ns11_12adj'
! 06/04 Collection 5    Min Min
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

!     pfmft and nfmft thresholds
      
      real pfmft_11maxthre(1)
      real pfmft_btd_min(1)
      real pfmft_land(4)
      real pfmft_ocean(4)
      real pfmft_snow(4)
      real pfmft_cold(4)   
      real nfmft_land(4)   
      real nfmft_ocean(4)  
      real nfmft_snow(4)
      real nfmft_desert(4)
      real nfmft_maxthre(1)
      
      common /p_n_fmft_thresholds/ pfmft_land,pfmft_ocean,pfmft_snow,pfmft_cold,  &
         nfmft_land,nfmft_ocean,nfmft_snow,nfmft_desert,nfmft_maxthre, pfmft_11maxthre, &
         pfmft_btd_min