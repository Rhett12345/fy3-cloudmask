subroutine Day_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,   &
                    testbits,qa_bits,nmtests,confdnc,btclr)

      implicit none
      save
 
!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over snow 
!      surfaces during daylight hours.
!
!      For daytime snow the groups are:
!          Group 1: High thick cloud
!                   13.9 micron bt test
!                   6.75 micron bt test 
!
!          Group 2: Low cloud - thick
!                   11-4 micron bt test
!                   11-12 micron bt test
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
! hi_elev       Logical variable indicating high elevation (> 2000 meters)

!
!!Output Parameters:
! testbits      six byte integer containing bit results
! qa_bits       ten byte integer containing qa bit results
! confdnc       product of all applied individual confidences
!
!!Revision History:
!
! Added 11-12 thin cirrus test 
! 06/04 Collection 5   R. Frey
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!Design Notes:
!    Externals:
!       Subroutines conf_test,set_bit,clear_bit,set_qa_bit
!
!!END
!-----------------------------------------------------------------------

      include 'global.inc'
      include 'Day_snow_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical visusd,cirrus_vis,hi_elev
! ...
! ... array arguments ..
      real pxldat(inband), btclr(7),tv11_12
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,c3,c4,c6,mas11_4,cmin1,cmin2,cmin4,locut,hicut,       &
           masir11,masir13,masir4,masv188,masir65,masir12,             &
           groups,fac,pre_confdnc,schi,cosvza,masdf1,pi,dtr,diftemp,   &
           dfthrsh
      integer nptests,kk,debug,h_output

! ... local arrays
      integer ngtests(3)
      real sn4_11(4)
      
      real,parameter :: Rel_equality_EPS = 0.000001

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
      !masv188 = pxldat(26)
      !masir4 = pxldat(22)
      !masir65 = pxldat(27)
      !masir11 = pxldat(31)
      !masir12 = pxldat(32)
      !masir13 = pxldat(35)
      masv188 = pxldat(19)
!      masir4 = pxldat(23)
      masir4 = pxldat(20)               !3.8 replace 3.9   jincheng
      !masir65 = pxldat(27)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
      !masir13 = pxldat(35)
! ...
      mas11_4 = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin4 = 1.0

! ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine Day_snow '',
!     +                   /)')
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
        call conf_test(tv11_12,pfmft_snow(1),pfmft_snow(3),pfmft_snow(4),   &
                       pfmft_snow(2),1,c1)
        cmin1 = min(cmin1,c1)  
       ! cmin1 = 1  ! added by minmin
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
!        call conf_test(tv11_12,nfmft_land(1),nfmft_land(3),nfmft_land(4),   &
!                       nfmft_land(2),1,c2)
         call conf_test(tv11_12,nfmft_snow(1),nfmft_snow(3),nfmft_snow(4),   &        ! jincheng
                       nfmft_snow(2),1,c2)
        cmin1 = min(cmin1,c2)
       ! cmin1 = 1  ! added by minmin
        ngtests(1) = ngtests(1) + 1
      endif
       
! ... perform tests.  note that some tests are not used
! ... in sunglint conditions.
!
! ... co2 high cloud test
!      if (nint(masir13) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,14)
!        if (masir13 .gt. dsco2(2)) then
!          call set_bit(testbits,14)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir13,dsco2(1),dsco2(3),dsco2(4),  &
!                       dsco2(2),1,c1)
!        cmin1 = min(cmin1,c1)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,dsco2(1),
!     +          dsco2(2),dsco2(3),dsco2(4)
!      endif
! ................................................................


!     H20 vapor channel (6.7 micron) high cloud test    ![no 6.5um channel]
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65 .gt. dsh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,dsh20(1),dsh20(3),dsh20(4),
!     *               dsh20(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,dsh20(1),
!     +          dsh20(2),dsh20(3),dsh20(4)
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
          dfthrsh = ds11_12hi(1)
        else
!         Add adjustment for snow cover.
          dfthrsh = diftemp + ds11_12adj(1)
        end if

!...    Set flags if test passed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf1.le.dfthrsh) then
          call set_bit(testbits,18)
          nptests = nptests + 1
        end if
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - (0.3 * dfthrsh)
!       hicut = dfthrsh - 1.25
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif
 
! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',8f10.2)') masdf1,
!     +          ds11_12hi(1),ds11_12adj(1),
!     +          masir11,schi,dfthrsh,locut,hicut
!      endif
! ................................................................

! ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) .and.  &
            nint(masir4) .ne.  nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19)
          mas11_4 = masir4 - masir11

          if(hi_elev) then
            sn4_11(1) = ds4_11hel(1)
            sn4_11(2) = ds4_11hel(2)
            sn4_11(3) = ds4_11hel(3)
            sn4_11(4) = ds4_11hel(4)
          else
            sn4_11(1) = ds4_11(1)
            sn4_11(2) = ds4_11(2)
            sn4_11(3) = ds4_11(3)
            sn4_11(4) = ds4_11(4)
          end if

          if (mas11_4 .le. sn4_11(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if
          call conf_test(mas11_4,sn4_11(1),sn4_11(3),sn4_11(4), &
                         sn4_11(2),1,c3)
          cmin2 = min(cmin2,c3)
          ngtests(2) = ngtests(2) + 1
        end if

! ..... debug statement ............................................
!        if (debug .gt. 0) then
!         write(h_output,'(1x,''mas11_4: '',5f10.2)')mas11_4,sn4_11(1),
!     +            sn4_11(2),sn4_11(3),sn4_11(4)
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
          if (masv188 .le. dsref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,dsref3(1),dsref3(3),dsref3(4),  &
                         dsref3(2),1,c4)
          cmin4 = min(cmin4,c4)
          ngtests(3) = ngtests(3) + 1
        end if

! ...   debug statement ............................................
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''masv188: '',6f10.4)')masv188,dsref3(1),
!     +                dsref3(2),dsref3(3),dsref3(4)
!        endif
! ................................................................
      endif

! ************   END OF GROUP 4 TESTS   ****************************
!
!     Check to see if thin cirrus bit should be set
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if(masv188 .lt. dstci(1) .and. masv188 .ge. dstci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
! ...     debug statement ............................................
!          if (debug .gt. 0) then
!            write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)')masv188,
!     +                    dstci(1),dstci(2)
!          endif
! ................................................................
        endif
      endif
!
!     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin4

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

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',9f8.5/,2f8.5)') c1,c2,c3,c4,c6,
!     +         cmin1,cmin2,cmin4,fac,confdnc
!      endif
! ................................................................

      return
      end
