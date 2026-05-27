      subroutine shadows(pxldat,shadow,visusd,qa_bits)

      implicit none
      save

!!F77 ************************************************************
!!Description:
! determines the presence of shadows and sets the appropriate bit flag.
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! visusd        Logical variable indicating whether visible data
!               was used or not
! 
!!Output Parameters:
! shadow        logical variable indicating shadow is present if .true.
! qa_bits       Byte array holding current pixel qa results
!
!!Revision History:
!
!!Team-unique Header:
!
!    This software is developed by the MODIS Science Data Support Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-32373.
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!
!!END****************************************************************
!
      include 'shadows_thr.inc'
      include 'global.inc'

! ... scalar arguments ..
      logical shadow,visusd
! ...
! ... array arguments ..
      real pxldat(inband)
      byte qa_bits(10)
! ...
! ... local scalars ..
      real masv88,masv66,masv945,vrat,masv124
      integer debug,h_output

! ... Common statement for debug purposes
!      common / bug / debug, h_output

! ...
      masv66 = pxldat(3)
      masv88 = pxldat(4)
      masv945 = pxldat(18)
!      masv124 = pxldat(5) ! 1.03 replace 1.24

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Within shadows testing routine '',/)')
!      endif
! ................................................................

 
!     Reflectance at 0.945 um must be ge 12% and "visible ratio" gt the
!     ocean threshold or else we're seeing a shadow.  The test on the 
!     ratio is to assure that the scene is not clear-sky conditions over
!     a sub-grid scale water body. 
      if (visusd) then
        if (nint(masv66) .ne. nint(bad_data) .and. &
            nint(masv88) .ne. nint(bad_data) .and. &
           nint(masv945) .ne. nint(bad_data) ) then !.and. &
!           nint(masv124) .ne. nint(bad_data)) then

          vrat = (masv88 - masv66) / (masv66 + masv88)

!         Set qa bits which show that we actually did test for shadows
          call set_qa_bit(qa_bits,10)

! ...     debug statement ............................................
!          if (debug .gt. 0) then
!             write(h_output,'(10x,'' Shadows test '',5f10.2,/)') 
!     *           masv66,masv88,vrat,masv945,masv124
!          endif
! ....................................................................

          if(masv945 .lt. shadnir(1) .and. vrat .gt. shavrat(1) .and.   &
             masv945 .gt. shadnir(2) ) then !.and. masv124 .lt. shad124(1)) then
            shadow = .true.
          else
            shadow = .false.
          end if

! ....    debug statement ............................................
!          if (debug .gt. 0) then
!            if (shadow) then
!              write(h_output,'(10x,'' Shadow found ''/)') 
!            else 
!              write(h_output,'(10x,'' Shadow not found ''/)') 
!            endif
!          endif
! ....................................................................
        endif
      endif

      return
      end
