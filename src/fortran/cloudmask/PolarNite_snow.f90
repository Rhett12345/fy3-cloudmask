      subroutine PolarNite_snow(pxldat,vza,hi_elev,land,antarctic,   &
                                testbits,qa_bits,nmtests,confdnc,btclr)

      implicit none
      save
 
!--------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over polar snow 
!      surfaces during nighttime hours.
!
!      For polar night over snow/ice the cloud test groups are:
!          Group 1: High thick cloud
!                   6.75 micron bt test (not is use with mas)
!
!          Group 2: Low cloud - thick
!                   11-12 micron bt test
!                   7.3-11 micron bt test
!                   11-3.9 micron bt tests
!
!          Group 5: High cloud - thin
!                   3.9-12 micron bt test
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Viewing zenith angle for current pixel
! hi_elev       Indicates elevation > 2000 meters
! land          Indicates land surface beneath snow
! antarctic     Indicates south of -60 latitude
!
!!Output Parameters:
! testbits      six byte integer containing bit results
! qa_bits       ten byte array containing QA bit results
! nmtests       Acutal number of tests applied in this subroutine
! confdnc       product of all applied individual confidences
!
!!Revision History:
! 06/04 Collection 5   R. Frey:
! Added 11-12 um thin cirrus test, 7.3-11 um cloud test
! 10/04 Collection 5b  R. Frey:
! Removed 11-12 test from Antarctica
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!    Externals:
!       Subroutines conf_test,set_bit,set_qa_bit,tview,get_pn_thresholds
!
!!END
!--------------------------------------------------------------------

      include 'global.inc'
      include 'PolarNite_snow_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real confdnc, vza
      integer nmtests
      logical hi_elev,land,antarctic
! ...
! ... array arguments ..
      real pxldat(inband),btclr(7),tv11_12
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,mas4_12,masir11,masir12,masir13,masir4,masir7,mas7_11,        &
           c3,masir65,mas11_4,c4,c5,cmin1,cmin2,cmin5,fac,dtr,cosvza,pi,    &
           pre_confdnc,groups,locut,hicut,midpt,power,schi,diftemp,         &
           dfthrsh,masdf1,c6
      integer nptests,debug,h_output,kk

!     local arrays
      integer ngtests(3)

      real,parameter :: Rel_equality_EPS = 0.000001
       
! ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview,get_pn_thresholds

! ... Common statement for debug purposes
!      common / bug / debug, h_output

! ... initialize variables
      pi = acos(-1.0)
      dtr = pi/180.0

! ... nmtests counts the number of tests applied to this pixel
      nmtests = 0

! ... nptests counts the number of tests passed
      nptests = 0

! ... set confidence to 1.0 to begin with
      confdnc = 1.0

! ... place band values into individual variables for easy
! ... identification
!      masir4 = pxldat(22)
!      masir65 = pxldat(27)
!      masir7 = pxldat(28)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)
!      masir13 = pxldat(33)

      masir4 = pxldat(21) ! 4.05 replace 3.959
!      masir65 = pxldat(27)
      masir7 = pxldat(22)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
!      masir13 = pxldat(33)
! ...
! ...
      mas4_12 = 0.0
      mas11_4 = 0.0
      groups = 0.0
      pre_confdnc = 0.0
      fac = 0

! ... the ! suffix variables represent individual test confidences
	  c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin5 = 1.0

! ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine PolarNite_snow '',
!     +                   /)')
!      endif
! ................................................................


!     **** GROUP 1 TESTS *************************************

! === PFMFT test disabled (btclr requires NWP RTM) ===
!       if (nint(masir11) .ne. nint(bad_data) .and.   &
!           nint(masir12) .ne. nint(bad_data) .and.   &
!           (masir11 < pfmft_11maxthre(1)) .and.   &
! !          (masir11-masir12) < pfmft_btd_min(1) ) then
!           (btclr(5)-btclr(6)) > pfmft_btd_min(1) ) then          !jincheng
! 		  nmtests = nmtests + 1
!         if ((masir11 > 270.0) .and. (btclr(5) > 270.0)) then
!             tv11_12 = (masir11 - masir12) -  &
!                       (btclr(5) - btclr(6)) *(masir11 - 260.0) / &
!                       (btclr(5) - 260.0)
!         else
!             tv11_12 = (masir11 - masir12)
!         endif
!         call set_qa_bit(qa_bits,14)
!         !if (masir11.gt.dlco2(2)) then
!           call set_bit(testbits,14)
!           nptests = nptests + 1
!         !end if
! !        call conf_test(tv11_12,pfmft_ocean(1),pfmft_ocean(3),pfmft_ocean(4),   &
! !                       pfmft_ocean(2),1,c1)
!         call conf_test(tv11_12,pfmft_snow(1),pfmft_snow(3),pfmft_snow(4),   &              !jincheng
!                        pfmft_snow(2),1,c1)
!         cmin1 = min(cmin1,c1)
!         ngtests(1) = ngtests(1) + 1
! === PFMFT test disabled end ===

! === NFMFT test disabled (btclr requires NWP RTM) ===
!       if (nint(masir11) .ne. nint(bad_data) .and.   &
!           nint(masir12) .ne. nint(bad_data) .and.   &
!           (masir11-masir12) <= nfmft_maxthre(1) ) then
!         nmtests = nmtests + 1
!         !tv11_12 = (btclr(5) - btclr(6)) - (masir11 - masir12)
!         tv11_12 =  (masir11 - masir12) - (btclr(5) - btclr(6))
!         call set_qa_bit(qa_bits,15)
!         !if (masir11.gt.dlco2(2)) then
!           call set_bit(testbits,15)
!           nptests = nptests + 1
!         !end if
! !        call conf_test(tv11_12,nfmft_ocean(1),nfmft_ocean(3),nfmft_ocean(4),   &
! !                       nfmft_ocean(2),1,c2)
!         call conf_test(tv11_12,nfmft_snow(1),nfmft_snow(3),nfmft_snow(4),   &            !jincheng
!                        nfmft_snow(2),1,c2)
!         cmin1 = min(cmin1,c2)
!         ngtests(1) = ngtests(1) + 1
! === NFMFT test disabled end ===
      
!     H20 vapor channel (6.7 micron) high cloud test
!      if ( .not. (antarctic .and. land) ) then
!        if (nint(masir65) .ne. nint(bad_data)) then
!          nmtests = nmtests + 1
!          call set_qa_bit(qa_bits,15)
!          if (masir65.gt.pnsh20(2)) then
!            call set_bit(testbits,15)
!            nptests = nptests + 1
!          end if
!          call conf_test(masir65,pnsh20(1),pnsh20(3),pnsh20(4),
!     *                   pnsh20(2),1,c2)
!          cmin1 = min(cmin1,c2)
!          ngtests(1) = ngtests(1) + 1
!        endif
!      endif
!
! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2,l5)') masir65,pnsh20(1),
!     +          pnsh20(2),pnsh20(3),pnsh20(4),antarctic
!      endif
! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
!
!
!     ****  GROUP 2 TESTS  ***********************************

! ... 11-12um brightness temperature difference test
! ... for thin cirrus).

      if ( .not. (antarctic .and. land) ) then

        if (nint(masir11) .ne. nint(bad_data) .and.  &
            nint(masir12) .ne. nint(bad_data) .and.  &
            vza .gt. 0.0) then

          masdf1 = masir11 - masir12
! ...     calculate secant of viewing zenith angle.
          cosvza = cos(vza*dtr)
          if (abs(cosvza).gt.Rel_equality_EPS) then
            schi = 1.0/cosvza
          else
            schi = 99.0
          end if

! ...     interpolate look-up table values of 11 - 12 micron bt
! ...     difference thresholds (function of viewing zenith
! ...     and 11 micron brightness temperature).
          call tview(1,schi,masir11,diftemp)

! ...     if a threshold was determined by apollo, then use this
! ...     as the thin cirrus test, otherwise use a standard threshold
          if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
            dfthrsh = pns11_12hi(1)
          else
!           Add adjustment for snow cover.
            dfthrsh = diftemp + pn11_12adj(1)
          end if

!...      Set flags if test passed
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,18)
          if (masdf1.le.dfthrsh) then
            call set_bit(testbits,18)
            nptests = nptests + 1
          end if

          locut = dfthrsh
          midpt = dfthrsh - (0.3 * dfthrsh)
          if(masir11 .lt. 270.0) then
            hicut = midpt - (0.2 * dfthrsh)
          else
            hicut = midpt - 1.25
          end if

          call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c6)
          cmin2 = min(cmin2,c6)
          ngtests(2) = ngtests(2) + 1

        endif

! ..... debug statement ..........................................
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''APOLLO masdf1: '',8f10.2)') masdf1,
!     +          pns11_12hi(1),pn11_12adj(1),
!     +          masir11,schi,dfthrsh,locut,hicut
!        endif
! ................................................................

      endif

! ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir4) .ne.  nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
        mas11_4 = masir11 - masir4
        !mas11_4 = masir11 - (masir4 - 2.0)  ! revised by minmin 20190109

        call get_pn_thresholds(masir11,bt_11_bounds,pn_11_4l,pn_11_4m1,     &
                               pn_11_4m2,pn_11_4m3,pn_11_4h,locut,hicut,    &
                               midpt,power)

        if (mas11_4 .le. midpt) then
          call set_bit(testbits,19)
          nptests = nptests + 1
        end if
        call conf_test(mas11_4,locut,hicut,power,midpt,1,c3)
        cmin2 = min(cmin2,c3)
        ngtests(2) = ngtests(2) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''mas11_4: '',6f10.2)')mas11_4,locut,
!     +            hicut,midpt,power,masir11
!      endif
! ................................................................

! ... 7.3 minus 11 micron cloud test
      if (nint(masir11) .ne. nint(bad_data) .and.    &
          nint(masir7) .ne.  nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,23)
        mas7_11 = masir7 - masir11

        if(land) then
!         Use land (snow) thresholds.
          call get_pn_thresholds(masir11,bt_11_bnds2,pn_7_11l,pn_7_11m1,      &
                                 pn_7_11m2,pn_7_11m3,pn_7_11h,locut,hicut,    &
                                 midpt,power)
        else
!         Use water (ice) thresholds.
          call get_pn_thresholds(masir11,bt_11_bnds2,pn_7_11lw,pn_7_11m1w,    &
                                 pn_7_11m2w,pn_7_11m3w,pn_7_11hw,locut,hicut, &
                                 midpt,power)
        end if

        if (mas7_11 .gt. midpt) then
          call set_bit(testbits,23)
          nptests = nptests + 1
        end if
        call conf_test(mas7_11,locut,hicut,power,midpt,1,c5)
    !    cmin2 = min(cmin2,c5)
    !    ngtests(2) = ngtests(2) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''mas7_11: '',7f10.2)')mas7_11,locut,
!     +            hicut,midpt,power,masir11,c5
!      end if
! ................................................................

! *******     END OF GROUP 2 TESTS ****************************




! *******    START OF GROUP 5 TESTS  **************************
! ... 4-12um brightness temperature difference test
! ... for thin cirrus).
      if (nint(masir12) .ne. nint(bad_data) .and.     &
          nint(masir4) .ne.  nint(bad_data) .and.     &
          nint(masir11) .ne.  nint(bad_data)) then

        if( (.not. hi_elev) ) then

          mas4_12 = masir4 - masir12
          !mas4_12 = (masir4 - 2.0) - masir12 ! revised by minmin 20190109
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,17)

          call get_pn_thresholds(masir11,bt_11_bounds,pn_4_12l,       &
                                 pn_4_12m1,pn_4_12m2,pn_4_12m3,       &
                                 pn_4_12h,locut,hicut,midpt,power)

          if (mas4_12 .le. midpt) then
            nptests = nptests + 1
            call set_bit(testbits,17)
          end if
          call conf_test(mas4_12,locut,hicut,power,midpt,1,c4)
          cmin5 = min(cmin5,c4)
          ngtests(3) = ngtests(3) + 1

        end if
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!      write(h_output,'(1x,''mas4_12: '',6f10.2)')mas4_12,locut,
!     +            hicut,midpt,power,masir11
!      endif
! ................................................................
! ********    END OF GROUP 5 TESTS  *****************************


!     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin5

!     Next, make sure you have all groups covered
      groups = 0
      do kk = 1,3
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if(groups .gt. 0) then
        fac = 1.0 / groups
      else
        fac = 0.0
      end if
!     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac

! ................................................................

!     Three last tests.  If the 6.5 micron brightness temperature is
!     greater than the 11 micron value, then there is an inversion
!     and is very, very likely to be clear.  Same deal for the 13
!     minus 11 and the 7.2 minus 11 difference.
!      if (nint(masir11) .ne. nint(bad_data) .and.
!          nint(masir65) .ne.  nint(bad_data)) then
!        call set_qa_bit(qa_bits,26)
!        if (masir65 - masir11 .gt. pn65_11(1)) then 
!           confdnc = 1.0
!           call set_bit(testbits,26)
!        endif
!      endif

!      if(confdnc .lt. 1.0) then
!
!        if (nint(masir11) .ne. nint(bad_data)
!     +                   .and.
!     +      nint(masir13) .ne.  nint(bad_data)) then
!          call set_qa_bit(qa_bits,26)
!          if ((masir13 - masir11) .gt. pn13_11(1)) then 
!            confdnc = 1.0
!           call set_bit(testbits,26)
!          endif
!        endif
!
!      end if

!      if(confdnc .lt. 1.0) then
!
!        if (nint(masir11) .ne. nint(bad_data) .and.   &
!            nint(masir7) .ne.  nint(bad_data)) then
!          call set_qa_bit(qa_bits,26)
!          if ((masir7 - masir11) .gt. pn7_11(1)) then 
!            confdnc = 1.0
!            call set_bit(testbits,26)
!          endif
!        endif
!
!      end if

! ... debug statement ............................................
!      if (debug .gt. 0) then
!          write(h_output,'(10x,''Final inversion tests:'',/,4f10.2,/)')
!     +                masir65,masir13,masir11,confdnc
!      endif
! ...................................................................


! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',9f8.5/,f8.5)') c2,c3,c4,c5,c6,
!     +         cmin1,cmin2,cmin5,fac,confdnc
!      endif
! ................................................................

      return
      end
