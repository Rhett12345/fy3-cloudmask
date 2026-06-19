subroutine LandDay(pxldat,vza,visusd,vrused,cirrus_vis,   &
                   hi_elev,testbits,qa_bits,nmtests,confdnc,btclr,is_cold_sfc)
 
!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over land 
!      surfaces during daylight hours.
!
!      For daytime land the groups are:
!          Group 1: High thick cloud
!                   13.9 micron bt test (masir13) 
!                   6.75 micron bt test 
!
!          Group 2: Low cloud - thick
!                   11-12 micron bt tests
!                   11-4 micron bt tests
!        
!          Group 3: Thick cloud
!                   .66 micron reflectance test (masv66)
!                   .87/.66 micron reflectance ratio test
! 
!          Group 4: Thin cirrus test
!                   1.38 micron reflectance test 
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Current pixel viewing angle
! visusd        Logical variable indicating whether vis data used or not
! vrused        Logical variable defining when reflectance ratio
!               test can be used.
! cirrus_vis 	Logical variable flagging thin cirrus contaminated
!               scenes in the visible
! hi_elev       Logical variable indicating high elevation (> 2000 meters)
!
!!Output Parameters:
! testbits      6 byte array containing cloud mask bit results
! qa_bits       10 byte array containing QA bit results
! confdnc       product of all applied individual confidences
!
!!Revision History:
! Changed implementation of 11-12 um thin cirrus test (Key version)
! 06/04 Collection 5    R. Frey
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!Design Notes:
!
!   Externals:
!      Subroutines: 
!      conftest,tview,set_qa_bit,set_bit,clear_bit
!
!!End
!---------------------------------------------------------------------

      implicit none
      save

      include 'global.inc'
      include 'LandDay_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      integer(kind=1) :: is_cold_sfc
      logical visusd,vrused,cirrus_vis,hi_elev
! ...
! ... array arguments ..
      real pxldat(inband)
      real btclr(7)
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,c3,c4,c5,cosvza,dfthrsh,diftemp,dtr,mas11_4,masdf1,    &
           masir11,masir12,masir13,masir4,masv188,masv66,               &
           masv88,pi,schi,vrat,c6,masir65,c7,cmin1,                     &
           cmin2,cmin3,cmin4,groups,fac,pre_confdnc,                    &
           eta,etad,etan,locut,hicut,s1,s2, tv11_12
      real, parameter :: Rel_equality_EPS = 0.000001
 
      integer nptests,debug,h_output,kk
! ...
! ... local arrays ..
      integer ngtests(4)
! ...
! ... external subroutines ..
      external conf_test,tview,set_bit,clear_bit,set_qa_bit

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
!      masv66 = pxldat(1)
!      masv88 = pxldat(2)
!      masv188 = pxldat(26)
!      masir4 = pxldat(22)
!      masir65 = pxldat(27)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)
!      masir13 = pxldat(35)

      masv66 = pxldat(3)
      masv88 = pxldat(4)
      masv188 = pxldat(19)
      masir4 = pxldat(20)  ! 3.8 replace 3.959
!       masir4 = 0.30*pxldat(21) + 0.70*pxldat(20) 
!      masir65 = pxldat(27)
!      masir65 = pxldat(22) + 0.5 ! added by minmin
      masir11 = pxldat(24)
      masir12 = pxldat(25)
!      masir13 = pxldat(35)

! ... Initialization
      masdf1 = 0.0
      cosvza = 0.0
      schi = 0.0
      diftemp = 0.0
      dfthrsh = 0.0
      vrat = 0.0
      mas11_4 = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      c7 = 0.0
! ... cmin are the group confidences
      cmin1 = 1.0
      cmin2 = 1.0
      cmin3 = 1.0
      cmin4 = 1.0

! ... initialize group number holder
      do 10 kk = 1 , 4
         ngtests(kk) = 0 
  10  continue

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine LandDay '',/)') 
!      endif
! ................................................................

 
!     **** GROUP 1 TESTS *************************************
! ... co2 high cloud test
!      if (nint(masir13) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,14)
!        if (masir13.gt.dlco2(2)) then
!          call set_bit(testbits,14)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir13,dlco2(1),dlco2(3),dlco2(4),   &
!                       dlco2(2),1,c1)
!        cmin1 = min(cmin1,c1)
!        ngtests(1) = ngtests(1) + 1
!      endif

!   pfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11 < pfmft_11maxthre(1)) .and.   &
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
        
         !  print*,'tv11_12=',tv11_12,c1
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif
      !print*,'vv1=',tv11_12,pfmft_land(1),pfmft_land(3),pfmft_land(4),pfmft_land(2)
      
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
     ! print*,'vv2=',tv11_12, nfmft_land(1),nfmft_land(3),nfmft_land(4), nfmft_land(2)
! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,dlco2(1),
!     +          dlco2(2),dlco2(3),dlco2(4)
!      endif
! ................................................................


!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65.gt.dlh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,dlh20(1),dlh20(3),dlh20(4),dlh20(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif
      
      
! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,dlh20(1),
!     +          dlh20(2),dlh20(3),dlh20(4)
!      endif
! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
 
 !print*,masir11,masir12,masir4,masv66,masv88,masv188
 
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
          dfthrsh = dl11_12hi(1)
        else
          dfthrsh = diftemp
        end if

!...    Set flags if test passed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf1.le.dfthrsh) then
          call set_bit(testbits,18)
          nptests = nptests + 1
        end if
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - 1.25
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c3)
        cmin2 = min(cmin2,c3)
        ngtests(2) = ngtests(2) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',4f10.2)') masdf1,
!     +          dfthrsh,locut,hicut
!      endif
! ................................................................


! ... 11 minus 4 micron BTDIF fog and low cloud test.
! ... for now placing in the SWIR bit place holder
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) .and.  &
            nint(masir4) .ne.  nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19)
          mas11_4 = masir11 - masir4
          if (mas11_4.ge.dl11_4lo(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if
          call conf_test(mas11_4,dl11_4lo(1),dl11_4lo(3),dl11_4lo(4),  &
                         dl11_4lo(2),1,c4)
         cmin2 = min(cmin2,c4)
         ngtests(2) = ngtests(2) + 1
        endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!         write(h_output,'(1x,''mas11_4: '',5f10.2)')mas11_4,dl11_4lo(1),
!     +            dl11_4lo(2),dl11_4lo(3),dl11_4lo(4)
!      endif
! ................................................................
      end if
! *******     END OF GROUP 2 TESTS ****************************
 
 
 
! ********  START OF GROUP 3 TESTS ****************************
! ... visible (channel 1) reflectance threshold test.
      if (visusd) then
        if (nint(masv66) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20)
          if (masv66.le.dlref1(2)) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(masv66,dlref1(1),dlref1(3),dlref1(4),dlref1(2),1,c5)
          cmin3 = min(cmin3,c5)
          ngtests(3) = ngtests(3) + 1
        end if
!       print*,'vvs',masv66,dlref1(1),dlref1(2),dlref1(3),cmin3 
! ...   debug statement ............................................
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''masv66: '',5f10.2)') masv66,dlref1(1),
!     +            dlref1(2),dlref1(3),dlref1(4)
!        endif
! ................................................................
      end if

! ... visible channel ratio test (channel 2 / channel 1)
! ... Changed to implement GEMI test instead of straight ratio
      if (visusd .and. vrused) then
        if (nint(masv66) .ne. nint(bad_data) .and.  &
            nint(masv88) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,21)
! ...     Scale values by 100 to make consistent with MAS version
          s1 = masv66 * 100.
          s2 = masv88 * 100.
          etan = 2.0 * (s2-s1) + 1.5*s2 + 0.5*s1
          etad = s2 + s1 + 0.5
          eta = etan / etad
          vrat=eta * (1.0-0.25*eta) - ((s1-0.125) / (1.0-s1))
          if(vrat .gt. dlvrat(2)) then
            nptests = nptests + 1
            call set_bit(testbits,21)
          end if
          call conf_test(vrat,dlvrat(1),dlvrat(3),dlvrat(4),dlvrat(2),1,c6)
          cmin3 = min(cmin3,c6)
          ngtests(3) = ngtests(3) + 1
        end if

! ...   debug statement ............................................
!        if (debug .gt. 0) then
!           write(h_output,'(1x,''GEMI: '',7f10.2)')vrat,masv88,
!     +                masv66,dlvrat(1),dlvrat(2),dlvrat(3),dlvrat(4)
!        endif
! ................................................................
      end if

! ******       END OF GROUP 3 TESTS   ****************************
!
!
!
!
! ***********   START OF GROUP 4 TESTS  *************************
! ... near infrared high cloud test
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (masv188.le.dlref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,dlref3(1),dlref3(3),dlref3(4),   &
                         dlref3(2),1,c7)
          cmin4 = min(cmin4,c7)
          ngtests(4) = ngtests(4) + 1
        endif

! ...   debug statement ............................................
!        if (debug .gt. 0) then
!           write(h_output,'(1x,''masv188: '',6f10.4)')masv188,dlref3(1),
!     +                dlref3(2),dlref3(3),dlref3(4)
!        endif
! ................................................................
      end if
! ************   END OF GROUP 4 TESTS   ****************************
 
 
!     Check to see if thin cirrus bit should be set
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then 
          call set_qa_bit(qa_bits,9)
          if (masv188 .lt. dltci(1) .and. masv188 .ge. dltci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
            !cirrus_vis = .false. ! added by minmin
          endif
! ...     debug statement ............................................
!          if (debug .gt. 0) then
!             write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)')masv188,
!     +                         dltci(1),dltci(2)
!          endif
! ................................................................
        endif
      endif

!     Determine intermediate confidence based on group values

!     Next, make sure you have all groups covered
      groups = 0
      pre_confdnc = 1.0
      do kk = 1,4
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
          if (kk .eq. 1) pre_confdnc = pre_confdnc * cmin1
          if (kk .eq. 2) pre_confdnc = pre_confdnc * cmin2
          if (kk .eq. 3) pre_confdnc = pre_confdnc * cmin3
          if (kk .eq. 4) pre_confdnc = pre_confdnc * cmin4
        end if
      enddo
      if (groups .gt. 0) then
        fac = 1.0 / groups
!       Find final pixel confidence as nth root of group tests
        confdnc = pre_confdnc**fac
      else
        confdnc = 1.0
      end if

!       print*,'LandDay============================'      
!       print*,'cmin1=',cmin1                           !jincheng
!       print*,'cmin2=',cmin2
!       print*,'cmin3=',cmin3
!       print*,'cmin4=',cmin4
       !print*,'t11-12',tv11_12,masv188,visusd,'confdnc = ',cmin1,cmin2,cmin3,cmin4,'c-',confdnc,fac,ngtests(1:4)


! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',8f8.5/,5f8.5)') c1,c2,c3,c4,c5,
!     +         c6,c7,cmin1,cmin2,cmin3,cmin4,fac,confdnc
!      endif
! ................................................................

      return
      end
