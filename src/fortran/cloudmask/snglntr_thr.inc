!!F77-INC*******************************************************
!
!!Description: Included file containing thresholds for sun glint
!               areas used in the cloud mask MOD35 processing
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
!!END**************************************************

! ... Sun Glint Thresholds

      real snglntv(2)  
      real snglntvcl(2)
      real snglntvch(2)
   
      real sg_tbdfl(1)
      real sg_tbdfh(1)
      real snglrat(1)

      real snglnt0(4)
      real snglnt10(4)
      real snglnt20(4)
      real snglnt_bounds(4)

      common /snglnt_thresholds/ snglntv,snglntvcl,snglntvch,sg_tbdfl,sg_tbdfh, &
            snglnt0,snglnt10,snglnt20,snglnt_bounds,snglrat
