subroutine LandDay_desert_c(pxldat,vza,visusd,cirrus_vis,  &
                            hi_elev,testbits,qa_bits,      &
                            nmtests,confdnc,btclr,is_cold_sfc)
 
      implicit none
      save

!--------------------------------------------------------------------
!!F77 
!
!!Description:
!
!      Routine for performing clear sky tests over coastal desert 
!      surfaces during daylight hours.
!
!      For daytime coastal desert the groups are:
!          Group 1: High thick cloud
!                   13.9 micron bt test (masir13) 
!                   6.75 micron bt test 
!
!          Group 2: Mainly Low cloud - thick
!                   11-4 micron bt tests
!                   11-12 Thin cirrus tests
!        
!          Group 3: Thick cloud
!                   .87 micron reflectance test (masv87)
! 
!          Group 4: Thin cirrus test
!                   1.38 micron reflectance test 
!
!!Input Parameters:
! real    pxldat Array containing reflectance or brightness temperatures
!                for all bands for a single pixel
! real    vza    Current pixel viewing angle
! logical visusd Logical variable indicating whether vis data used 
! logical cirrus_vis Logical variable flagging thin cirrus contaminated
!               scenes in the visible
! logical hi_elev Logical variable indicating elevation > 2000 meters
!
!!Output Parameters:
! byte    testbits six byte integer containing bit results
! byte    qa_bits  ten byte array containing QA bit results
! integer nmtests  Counts number of inidividual tests applied
! real    confdnc  product of all applied individual confidences
!
!----------------------------------------------------------------------

      include 'global.inc'
      include 'LandDay_desert_c_thr.inc'
      include 'pfmft_nfmft_thr.inc'

! ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical visusd,cirrus_vis,hi_elev
! ...
! ... array arguments ..
      integer(kind=1) :: is_cold_sfc
      real pxldat(inband),btclr(7),tv11_12
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,c3,c4,c5,cosvza,dfthrsh,diftemp,dtr,mas11_4,masdf1,   &
           masir11,masir12,masir13,masir4,masv188,masv88,pi,schi,      &
           masir65,c6,cmin1,cmin2,cmin3,cmin4,pre_confdnc,             &
           groups,fac,locut,hicut
      real, parameter :: Rel_equality_EPS = 0.000001

      integer nptests,kk,debug,h_output
! ...
! ... local arrays ..
      real hicuta(2),locuta(2),midpta(2)
      integer ngtests(4)
! ...
! ... external subroutines ..
      external conf_test,conf_test_2val,tview,set_bit,clear_bit,set_qa_bit

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
!      masv88 = pxldat(2)
!      masv188 = pxldat(26)
!      masir4 = pxldat(22)
!      masir65 = pxldat(27)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)
!      masir13 = pxldat(35)

      masv88 = pxldat(4)
      masv188 = pxldat(19)
      masir4 = pxldat(20) ! 3.8 replace 3.959
!      masir4 = 0.30*pxldat(21) + 0.70*pxldat(20) 
!      masir65 = pxldat(27)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
!      masir13 = pxldat(35)

! ...
      masdf1 = 0.0
      cosvza = 0.0
      schi = 0.0
      diftemp = 0.0
      dfthrsh = 0.0
      mas11_4 = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
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
!        write(h_output,'(10x/,''Subroutine LandDay_desert_c '',
!     +                      /)')
!      endif
! ................................................................

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
        call conf_test(tv11_12,nfmft_desert(1),nfmft_desert(3),nfmft_desert(4),   &
                       nfmft_desert(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif
      
! ... perform tests.  note that some tests are not used
! ... in sun glint conditions.
!
!     **** GROUP 1 TESTS *************************************
! ... co2 high cloud test
!      if (nint(masir13) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,14)
!        if (masir13 .gt. ldsco2_c(2)) then
!          call set_bit(testbits,14)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir13,ldsco2_c(1),ldsco2_c(3),ldsco2_c(4),  &
!                       ldsco2_c(2),1,c1)
!        cmin1 = min(cmin1,c1)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,ldsco2_c(1),
!     +          ldsco2_c(2),ldsco2_c(3),ldsco2_c(4)
!      endif
! ................................................................


!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65 .gt. ldsh20_c(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,ldsh20_c(1),ldsh20_c(3),ldsh20_c(4),   &
!                       ldsh20_c(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,ldsh20_c(1),
!     +          ldsh20_c(2),ldsh20_c(3),ldsh20_c(4)
!      endif
! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
!
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
          dfthrsh = lds11_12hi_c(2)
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
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) .and.  &
            nint(masir4) .ne.  nint(bad_data)) then
          if(masir11 .le. 320.0) then
            nmtests = nmtests + 1
            mas11_4 = masir11 - masir4
            call set_qa_bit(qa_bits,19)
            if (mas11_4 .ge. lds11_4lo_c(2) .and.    &
                mas11_4 .le. lds11_4hi_c(2)) then
              call set_bit(testbits,19)
              nptests = nptests + 1
            end if
            locuta(1) = lds11_4lo_c(1)
            locuta(2) = lds11_4hi_c(1)
            hicuta(1) = lds11_4lo_c(3)
            hicuta(2) = lds11_4hi_c(3)
            midpta(1) = lds11_4lo_c(2)
            midpta(2) = lds11_4hi_c(2)
  
            call conf_test_2val(mas11_4,locuta,hicuta,1.0,midpta,2,c4)
            cmin2 = min(cmin2,c4)
            ngtests(2) = ngtests(2) + 1
          endif
        end if

! ..... debug statement ............................................
!        if (debug .gt. 0) then
!        write(h_output,'(1x,''mas11_4a: '',5f10.2)')mas11_4,
!     +        lds11_4lo_c(1),lds11_4lo_c(2),lds11_4lo_c(3),lds11_4lo_c(4)
!        write(h_output,'(1x,''mas11_4b: '',5f10.2)')mas11_4,
!     +        lds11_4hi_c(1),lds11_4hi_c(2),lds11_4hi_c(3),lds11_4hi_c(4)
!        endif
! ..................................................................
      endif


! *******     END OF GROUP 2 TESTS ****************************
!
!
!
! ********  START OF GROUP 3 TESTS ****************************
! ... visible (.88 micron) reflectance threshold test.
      if (visusd) then
        if (nint(masv88) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20)
          if (masv88 .le. ldsref2_c(2)) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(masv88,ldsref2_c(1),ldsref2_c(3),ldsref2_c(4),  &
                         ldsref2_c(2),1,c5)
          cmin3 = min(cmin3,c5)
          ngtests(3) = ngtests(3) + 1
        end if

!  ...  debug statement ............................................
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''masv88: '',5f10.2)') masv88,ldsref2_c(1),
!     +            ldsref2_c(2),ldsref2_c(3),ldsref2_c(4)
!        endif
! ..................................................................
      endif
! ******       END OF GROUP 3 TESTS   ****************************
 
 
 
 
! ***********   START OF GROUP 4 TESTS  *************************
! ... near infrared high cloud test
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (masv188 .le. ldsref3_c(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,ldsref3_c(1),ldsref3_c(3),ldsref3_c(4),   &
                         ldsref3_c(2),1,c6)
          cmin4 = min(cmin4,c6)
          ngtests(4) = ngtests(4) + 1
        end if
 
! ...   debug statement ............................................
!        if (debug .gt. 0) then
!           write(h_output,'(1x,''masv188: '',6f10.4)')masv188,
!     +           ldsref3_c(1),ldsref3_c(2),ldsref3_c(3),ldsref3_c(4)
!        endif
! ................................................................
      endif
! ************   END OF GROUP 4 TESTS   ****************************
 
 

!     Check to see if thin cirrus bit should be set
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (masv188.lt.ldstci_c(1) .and. masv188.ge.ldstci_c(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
        endif
      endif
!
!     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4

!     Next, make sure you have all groups covered
      groups = 0
      do kk = 1,4
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if (groups .gt. 0) fac = 1.0 / groups
!     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',8f8.5/,4f8.5)') c1,c2,c3,c4,c5,
!     +         c6,cmin1,cmin2,cmin3,cmin4,fac,confdnc
!      endif
! ................................................................

      return
      end
