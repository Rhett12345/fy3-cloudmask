      subroutine PolarNite_land(pxldat,vza,desert,hi_elev,sfctmp,    &
                               eco_type,testbits,qa_bits,nmtests,    &
                               confdnc,btclr,is_cold_sfc)  

      implicit none
      save

!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over polar land
!      surfaces during nightime hours.
!
!      For nighttime polar land the groups are:
!          Group 1: High thick cloud
!                   6.75 micron bt test 
!                   surface temperature test
!
!          Group 2: Low cloud - thick
!                   11-12 micron bt test
!                   11-4 micron bt test
!                   7.3-11 micron bt test
!
!          Group 5: High cloud - thin
!                   3.7-12 micron bt test
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           viewing zenith angle
! desert        flag indicating desert processing
! hi_elev       flag indicating high elevation processing (> 2000 m)
! sfctmp        Surface air temperature from model data
! eco_type      Ecosystem index

!
!!Output Parameters:
! testbits      four-byte integer containing bit results
! qa_bits       10 byte array containing QA bit results
! nmtests       Number of tests actually applied in this routine
! confdnc       product of all applied individual confidences
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Added surface temperature cloud test (GDAS sfc air vs 11 micron bt)
! Added 7.3-11 micron cloud test
! Implemented new 11-12 micron thin cirrus test (Key version)
! Modified 3.9-12 micron test (dynamic thresholds based on 11 micron bt)
! 10/04 Collection 5  R. Frey
! Added 11-12 and 3.9-11 um BTD conditions on choice of LST threshold
! Changed basic LST threshold from 10K to 12K
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!Design Notes:
!
! Externals:
!       Subroutines: conf_test,set_bit,set_qa_bit,tview,get_pn_thresholds
!
!!END
!-----------------------------------------------------------------------

      include 'global.inc'
      include 'PolarNite_land_thr.inc'
      include 'PolarNite_snow_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real confdnc,vza,sfctmp
      integer nmtests
      logical desert,hi_elev
      byte eco_type
! ...
! ... array arguments ..
      real pxldat(inband),btclr(7),tv11_12
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,mas4_12,masir11,masir12,masir4,masir7,mas7_11,             &
           c3,masir65,mas11_4,c4,cmin1,cmin2,cmin5,groups,               &
           fac,pre_confdnc,c5,c6,masdf1,schi,cosvza,dtr,pi,diftemp,      &
           dfthrsh,locut,hicut,midpt,power,a,                            &
           lst_thrsh,corr,sfcdif,c7,masdf2
      integer nptests,debug,h_output,kk

! ... local arrays
      integer ngtests(3)
! ...
      real, parameter :: Rel_equality_EPS = 0.000001
      real, parameter :: max_vza = 65.49
      
! ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview,get_pn_thresholds

! ... Common statement for debug purposes
!      common / bug / debug, h_output

! ...
! ... initialize variables
      pi = acos(-1.0)
      dtr = pi/180.0
! ... nmtests counts the number of tests applied to this pixel
      nmtests = 0
! ... nptests counts the number of tests passed
      nptests = 0
! ... confidence to 1.0 to begin with
      confdnc = 1.0

! ... place band values into individual variables for easy
! ... identification
!      masir4 = pxldat(22)
!      masir65 = pxldat(27)
!      masir7 = pxldat(28)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)

      masir4 = pxldat(20)  ! 3.8 replace 3.959
!      masir65 = pxldat(27)
      masir7 = pxldat(22)
      masir11 = pxldat(24)
      masir12 = pxldat(25)

! ...
      mas4_12 = 0.0
      mas11_4 = 0.0
      masdf1 = 0.0
      schi = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      c7 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin5 = 1.0

! ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine LandNite '',/)')
!      endif
! ................................................................

!   pfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11 < pfmft_11maxthre(1)) .and.   &
!          (masir11-masir12) < pfmft_btd_min(1) ) then
          (btclr(5)-btclr(6)) > pfmft_btd_min(1) ) then          !jincheng
		  nmtests = nmtests + 1
        if ((masir11 > 270.0) .and. (btclr(5) > 270.0)) then
            tv11_12 = (masir11 - masir12) -  &
                      (btclr(5) - btclr(6)) *(masir11 - 260.0) / &
                      (btclr(5) - 260.0)
        else
            tv11_12 = (masir11 - masir12)
        endif
        call set_qa_bit(qa_bits,14)
        !if (masir11.gt.dlco2(2)) then
          call set_bit(testbits,14)
          nptests = nptests + 1
        !end if
        if (is_cold_sfc == 1) then
           call conf_test(tv11_12,pfmft_cold(1),pfmft_cold(3),pfmft_cold(4),   &
                          pfmft_cold(2),1,c1)        
        else
           call conf_test(tv11_12,pfmft_land(1),pfmft_land(3),pfmft_land(4),   &
                          pfmft_land(2),1,c1)
        endif
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

!   nfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11-masir12) <= nfmft_maxthre(1) ) then
        nmtests = nmtests + 1
        !tv11_12 = (btclr(5) - btclr(6)) - (masir11 - masir12)
        tv11_12 =  (masir11 - masir12) - (btclr(5) - btclr(6))
        call set_qa_bit(qa_bits,15)
        !if (masir11.gt.dlco2(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        !end if
        call conf_test(tv11_12,nfmft_land(1),nfmft_land(3),nfmft_land(4),   &
                       nfmft_land(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif
      
!     **** GROUP 1 TESTS *************************************

!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65 .gt. pnlh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,pnlh20(1),pnlh20(3),pnlh20(4),
!     *                pnlh20(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,pnlh20(1),
!     +          pnlh20(2),pnlh20(3),pnlh20(4)
!      endif
! ................................................................


! ... Surface Temperature Test

      if ( nint(masir11) .ne. nint(bad_data) .and. (.not. hi_elev) .and.  &
           nint(masir12) .ne. nint(bad_data) .and. eco_type .ne. 8) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        masdf1 = masir11 - masir12
        !masdf2 = masir11 - masir4
        masdf2 = masir11 - (masir4-2.0)  ! revised by minmin 20190108 to correct masir4
        
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

        if(desert) then
          lst_thrsh = 20.0
        else if(masdf1 .ge. 0.0 .or. (masdf1 .lt. 0.0 .and. (masdf2 .le. -0.5 .or. masdf2 .ge. 1.0))) then
          lst_thrsh = 12.0
        else
          lst_thrsh = 20.0
        end if

        if(masdf1 .ge. 1.0) then
          midpt = lst_thrsh + (2.0 * nint(masdf1))
        else
          midpt = lst_thrsh
        end if
        a = vza / max_vza
        corr = (a**4) * 3.0
        midpt = midpt + corr
        locut = midpt + 2.0
        hicut = midpt - 2.0

        sfcdif = sfctmp - masir11

        if( sfcdif .lt. midpt ) then
          call set_bit(testbits,27)
          nptests = nptests + 1
        end if

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c7)
        cmin1 = min(cmin1,c7)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''sfctmp: '',9f9.3)') masdf1,masdf2,
!     +   sfctmp,masir11,sfcdif,locut,midpt,hicut,c7
!      endif

! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
!
! 
!     ****  GROUP 2 TESTS  ***********************************
! ... 11-12um brightness temperature difference test
! ... for thin cirrus).
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          vza .gt. 0.0) then
      
        masdf1 = masir11 - masir12
! ...   added apollo viewing angle/av4t regressed threshold.
! ...   calculate secant of viewing zenith angle.
        cosvza = cos(vza*dtr)
        if (abs(cosvza).gt.Rel_equality_EPS) then
          schi = 1.0/cosvza
        else
          schi = 99.0
        end if

! ...   interpolate look-up table values of 11 - 12 micron bt
! ...   difference thresholds (function of viewing zenith
! ...   and 11 micron brightness temperature).
        call tview(1,schi,masir11,diftemp)

! ...   if a threshold was determined by apollo, then use this
! ...   as the thin cirrus test, otherwise use a standard threshold
        if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
          dfthrsh = pnl11_12hi(1)
        else
!         Add 0.2 for likely snow cover.
          dfthrsh = diftemp + 0.2
        end if

!...    Set flags if test passed
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

        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c5)
        cmin2 = min(cmin2,c5)
        ngtests(2) = ngtests(2) + 1

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',5f10.2)') masdf1,
!     +          pnl11_12hi(1),dfthrsh,locut,hicut
!      endif
! ................................................................

! ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir4) .ne.  nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
       ! mas11_4 = masir11 - masir4
         mas11_4 = masir11 - (masir4-2.0)   ! revised by minmin 20180109
    

!       Use polar night snow thresholds here.  Logic is that land
!       surfaces poleward of 60 latitude are snow covered most of
!       the year and many times ancillary snow map is incomplete.

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
!      end if
! ................................................................

! ... 7.3 minus 11 micron cloud test 
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir7) .ne.  nint(bad_data)) then

!       Check 11 um brightness temperature.  This to guard against
!       false cloud retrievals during polar summer.
        if(masir11 .lt. 270.0) then

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,23)
          mas7_11 = masir7 - masir11

!         Use polar night snow thresholds here.  Logic is that land
!         surfaces poleward of 60 latitude are snow covered most of
!         the year and many times ancillary snow map is incomplete.

          call get_pn_thresholds(masir11,bt_11_bnds2,pn_7_11l,pn_7_11m1,       &
                                 pn_7_11m2,pn_7_11m3,pn_7_11h,locut,hicut,     &
                                 midpt,power)

          if (mas7_11 .gt. midpt) then
            call set_bit(testbits,23)
            nptests = nptests + 1
          end if
          call conf_test(mas7_11,locut,hicut,power,midpt,1,c6)
 !         cmin2 = min(cmin2,c6)
 !         ngtests(2) = ngtests(2) + 1

        end if
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''mas7_11: '',6f10.2)')mas7_11,locut,
!     +            hicut,midpt,power,masir11
!      end if
! ................................................................

! *******     END OF GROUP 2 TESTS ****************************



! *******    START OF GROUP 5 TESTS  **************************
! ... 4-12um brightness temperature difference test
! ... for thin cirrus)
      if (nint(masir12) .ne. nint(bad_data) .and.  &
          nint(masir4) .ne.  nint(bad_data)) then
        !mas4_12 = masir4 - masir12  !mas11_4 = masir11 - (masir4-2.0)   ! revised by minmin 20180109
        mas4_12 = (masir4-2.0) - masir12 
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,17)

!       Use polar night snow thresholds here.  Logic is that land
!       surfaces poleward of 60 latitude are snow covered most of
!       the year and many times ancillary snow map is incomplete.

        call get_pn_thresholds(masir11,bt_11_bounds,pn_4_12l,pn_4_12m1,      &
                               pn_4_12m2,pn_4_12m3,pn_4_12h,locut,hicut,     &
                               midpt,power)

        if (mas4_12 .le. midpt) then
          nptests = nptests + 1
          call set_bit(testbits,17)
        end if
        call conf_test(mas4_12,locut,hicut,power,midpt,1,c4)
        cmin5 = min(cmin5,c4)
        ngtests(3) = ngtests(3) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!      write(h_output,'(1x,''mas4_12: '',6f10.2)')mas4_12,locut,
!     +            hicut,midpt,power,masir11
!      endif
! ................................................................
! ********    END OF GROUP 5 TESTS  *****************************

!     Determine final confidence based on group values
      pre_confdnc = max(cmin1, 0.1) * max(cmin2, 0.1) * max(cmin5, 0.1)

!     Next, make sure you have all groups covered
      groups = 0
      do kk = 1,3
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if (groups .gt. 0) then
        fac = 1.0 / groups
!       Find final pixel confidence as nth root of group tests
        confdnc = pre_confdnc**fac
      else
        confdnc = 1.0
      end if


! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',5i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',8f8.5/,2f8.5)') c2,c3,c4,c5,
!     +         c7,cmin1,cmin2,cmin5,fac,confdnc
!      endif
! ................................................................

      return
      end
