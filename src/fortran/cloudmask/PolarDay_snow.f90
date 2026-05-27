      subroutine PolarDay_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,   &
                               testbits,qa_bits,nmtests,confdnc,btclr)

      implicit none
      save
 
!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over polar snow 
!      surfaces during daylight hours.
!
!      For daytime land type 1 the groups are:
!          Group 1: High thick cloud
!                   6.75 micron bt test (not is use with mas)
!
!          Group 2: Low cloud - thick
!                   11-4 micron bt tests
!                   11-12 micron bt tests
!        
!          Group 4: Thin cirrus test
!                   1.38 micron reflectance test 
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Viewing zenith angle for current pixel
! visusd        Logical variable indicating whether vis data used or not
! cirrus_vis 	Logical variable flagging thin cirrus contaminated
!               scenes in the visible
! hi_elev       Logical variable indicating elevations > 2000 meters.

!
!!Output Parameters:
! testbits      six byte integer containing bit results
! qa_bits       ten byte integer containing qa bit results
! confdnc       product of all applied individual confidences
!
!!Revision History:
! 06/04 Collection 5  R. Frey:
! Added 11-12 um thin cirrus test (J. Key version)
! 10/04 Collection 5b R. Frey:
! Added 11 um BT-dependent 3.9-11 BTD test thresholds; removed static
! thresholds.
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!Design Notes:
!    Externals:
!       Subroutines conf_test,set_bit,clear_bit,set_qa_bit,tview
!
!!END
!-----------------------------------------------------------------------

      include 'global.inc'
      include 'PolarDay_snow_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical visusd,cirrus_vis,hi_elev
! ...
! ... array arguments ..
      real pxldat(inband),btclr(7),tv11_12
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,c3,c4,mas11_4,cmin1,cmin2,cmin4,pi,dtr,masdf1,cosvza,schi,   &
           masir11,masir12,masir4,masv188,masir65,groups,                  &
           fac,pre_confdnc,diftemp,dfthrsh,locut,hicut,c6,midpt,power
      integer nptests,kk,debug,h_output

! ... local arrays
      integer ngtests(3)
      real sn4_11(4)

      real, parameter :: Rel_equality_EPS = 0.000001

! ... external subroutines ..
      external conf_test,set_bit,clear_bit,set_qa_bit,tview

!     Common statement for debug purposes
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
!      masv188 = pxldat(26)
!     masir4 = pxldat(22)
!      masir65 = pxldat(27)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)

      masv188 = pxldat(19)
      masir4 = pxldat(20) ! 3.8 replace 3.959
!      masir65 = pxldat(27)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
      
! ...
      mas11_4 = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c6 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin4 = 1.0

! ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine PolarDay_snow '',
!     +                   /)')
!      endif
! ................................................................

!   pfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11 < pfmft_11maxthre(1)) .and.   &
!          (masir11-masir12) < pfmft_btd_min(1) ) then
          (btclr(5)-btclr(6)) > pfmft_btd_min(1) ) then          !jincheng
!		  nmtests = nmtests + 1
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
!        print*,'tv11_12=',tv11_12
        call conf_test(tv11_12,pfmft_snow(1),pfmft_snow(3),pfmft_snow(4),   &
                       pfmft_snow(2),1,c1)
!        cmin1 = min(cmin1,c1)
!        cmin1 = 1.0  ! added by minmin 20180322
!        ngtests(1) = ngtests(1) + 1
      endif
!      print*,'c1 pfmft = ',cmin1
!   nfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11-masir12) <= nfmft_maxthre(1) ) then
!        nmtests = nmtests + 1
        !tv11_12 = (btclr(5) - btclr(6)) - (masir11 - masir12)
        tv11_12 =  (masir11 - masir12) - (btclr(5) - btclr(6))
        call set_qa_bit(qa_bits,15)
        !if (masir11.gt.dlco2(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        !end if
        !print*,'tv11_12_2=',tv11_12
        call conf_test(tv11_12,nfmft_snow(1),nfmft_snow(3),nfmft_snow(4),   &
                       nfmft_snow(2),1,c2)
!                print*,'c2 nfmft = ',c2
!        cmin1 = min(cmin1,c2)
!        cmin1 = 1.0  ! added by minmin 20180322
 !       ngtests(1) = ngtests(1) + 1
      endif
      

 
! ... perform tests.  
 
!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65 .gt. dpsh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,dpsh20(1),dpsh20(3),dpsh20(4),
!     *               dpsh20(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,dpsh20(1),
!     +          dpsh20(2),dpsh20(3),dpsh20(4)
!      endif
! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
 
 
 
!     ****  GROUP 2 TESTS  ***********************************

! ... 11-12um brightness temperature difference test
! ... for thin cirrus).
      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir12) .ne. nint(bad_data) .and.  &
                   vza .gt. 0.0) then

        masdf1 = masir11 - masir12
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
          dfthrsh = dps11_12hi(1)
        else
!         Add adjustment for snow cover.
          dfthrsh = diftemp + dps11_12adj(1)
        end if

!...    Set flags if test passed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf1 .le. dfthrsh) then
          call set_bit(testbits,18)
          nptests = nptests + 1
        end if
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - (0.3 * dfthrsh)
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',8f10.2)') masdf1,
!     +          dps11_12hi(1),dps11_12adj(1),
!     +          masir11,schi,dfthrsh,locut,hicut
!      endif
! ................................................................

! ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) .and.   &
            nint(masir4) .ne.  nint(bad_data)) then

          if(masir11 .gt. 230.0) then

            nmtests = nmtests + 1
            call set_qa_bit(qa_bits,19)
            mas11_4 = masir4 - masir11

            call get_pn_thresholds(masir11,bt_11_bnds3,dps4_11l,dps4_11m1,       &
                                   dps4_11m2,dps4_11m3,dps4_11h,locut,hicut,     &
                                   midpt,power)

            if (mas11_4 .le. midpt) then
              call set_bit(testbits,19)
              nptests = nptests + 1
            end if
            !print*,'thres=',locut,hicut,power,midpt
            call conf_test(mas11_4,locut,hicut,power,midpt,1,c3)
            cmin2 = min(cmin2,c3)
            ngtests(2) = ngtests(2) + 1
			
          end if
        end if

! ..... debug statement ............................................
!        if (debug .gt. 0) then
!         write(h_output,'(1x,''mas11_4: '',6f9.3)')masir11,mas11_4,
!     +         locut,midpt,hicut,c3
!        endif
! ..................................................................
      end if

! *******     END OF GROUP 2 TESTS ****************************
!
!
! ***********   START OF GROUP 4 TESTS  *************************
! ... near infrared high cloud test
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (masv188 .le. dpsref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,dpsref3(1),dpsref3(3),dpsref3(4),    &
                         dpsref3(2),1,c4)
          cmin4 = min(cmin4,c4)
         ! print*,'c1,c2,c3,c6,c4=',c1,c2,c3,c6,c4
          ngtests(3) = ngtests(3) + 1
        end if

! ...   debug statement ............................................
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''masv188: '',6f10.4)')masv188,dpsref3(1),
!     +                dpsref3(2),dpsref3(3),dpsref3(4)
!        endif
! ................................................................
      endif

! ************   END OF GROUP 4 TESTS   ****************************
!
!     Check to see if thin cirrus bit should be set
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if(masv188 .lt. dpstci(1) .and. masv188 .ge. dpstci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
! ...     debug statement ............................................
!          if (debug .gt. 0) then
!            write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)')masv188,   
!     +                    dpstci(1),dpstci(2)
!          endif
! ................................................................
        endif
      endif
!
!     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin4
!     print*,'c1 c2 c4',cmin1 , cmin2 , cmin4,nmtests

!     Next, make sure you have all groups covered
      groups = 0
      do kk = 1,3
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if(groups .gt. 0) fac = 1.0 / groups
!     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',8f8.5/,2f8.5)') c2,c3,c4,c6,
!     +         cmin1,cmin2,cmin4,fac,confdnc
!      endif
! ................................................................

      return
      end
