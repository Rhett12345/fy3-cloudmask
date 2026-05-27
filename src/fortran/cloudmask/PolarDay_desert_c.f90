      subroutine PolarDay_desert_c(pxldat,vza,visusd,cirrus_vis,hi_elev,  &
                                   testbits,qa_bits,nmtests,confdnc,      &
                                   btclr,is_cold_sfc)
 
      implicit none
      save

!--------------------------------------------------------------------
!!F77 
!
!!Description:
! Performs clear-sky spectral tests for polar coastal desert surfaces 
! during daylight conditions.  
!
! Each spectral test is placed in one of five test groups. The groups
! represented in this routine are:
!
!          Group 1: High, thick clouds using IR brightness temperatures
!                   6.7 micron test 
!
!          Group 2: Low, thick clouds and thin cirrus clouds using
!                   brightness temperature differences
!                   11-4 micron test
!                   11-12 thin cirrus test
!        
!          Group 3: Thick clouds using reflectance information
!                   .86 micron reflectance test
! 
!          Group 4: Thin cirrus clouds using NIR reflectances
!                   1.38 micron reflectance test 
!
!          Group 5: None.
!
! A "confidence of clear sky" is computed for each spectral test. 
! Confidences from single-threshold tests are calculated in subroutine
! 'conf_test'.  Those from double-threshold tests (where clear-sky 
! radiance data falls in a range between two thresholds or lies on
! either side of them) are generated in 'conf_test_2val'.
!
! The minimum confidence value in each group is defined as the group
! confidence.  Final confidence is defined as the nth root of the
! product of the group confidences, where n is the number of groups.
!
! A "qa bit" (in array 'qa_bits') is set for each test performed and
! a "test bit" is set (in array 'testbits') for each clear-sky test
! passed.
!
!!Input Parameters:
! real       pxldat       Reflectances/brightness temperatures in all
!                         bands for current 1km pixel
! real       vza          Viewing zenith angle for current pixel
! logical    visusd       Indicates visible data use
! logical    hi_elev      Flag indicating elevation > 2000 meters
!
!!Output Parameters:
! logical    cirrus_vis   Thin cirrus indicated in current pixel
! byte       testbits     Byte array containing bit flags
! byte       qa_bits      Byte array containing QA bit flags
! integer    nmtests      Counts number of inidividual tests applied
! real       confdnc      Final clear-sky confidence
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Implemented new 11-12 um thin cirrus test (J. Key version)
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
! 1996/4/3 13:55:25
! K. Strabala (kathys@ssec.wisc.edu)
! Original version.
!
!!Design Notes:
!
!   Externals:
!      Subroutines:  
!      conf_test_2val,conf_test,tview,set_bit,clear_bit
!      set_qa_bit
!
!!End
!----------------------------------------------------------------------

!     Declarations.

!     Include files.
      include 'global.inc'
      include 'PolarDay_desert_c_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
!     Parameter statements.
      real, parameter :: Rel_equality_EPS = 0.000001

!     Scalar arguments.
      real confdnc,vza
      integer nmtests
      logical visusd,cirrus_vis,hi_elev
        
!     Array arguments. 
      real pxldat(inband),btclr(7),tv11_12
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)
         
!     Local scalars.
      real c1,c2,c3,c4,c5,c6,cosvza,dfthrsh,diftemp,dtr,m31_22,m31_32,    &
           m31,m32,m22,m26,m02,pi,schi,m27,cmin1,cmin2,cmin3,cmin4,    &
           pre_confdnc,groups,fac,locut,hicut,masir11,masir12
      integer nptests,kk,debug,h_output

!     Local arrays.
      real hicuta(2),locuta(2),midpta(2)
      integer ngtests(5)
        
!     External routines.
      external conf_test,conf_test_2val,tview,set_bit,clear_bit,set_qa_bit

!     Common statement for debug variables.
!      common / bug / debug, h_output
        
!----------------------------------------------------------------------

!     Initialize variables.
      pi = acos(-1.0)
      dtr = pi/180.0

!     'nmtests' counts the number of tests applied to this pixel.
      nmtests = 0

!     'nptests' counts the number of tests which found no evidence
!     of cloud.
      nptests = 0

!     Place reflectance and brightness temperature values into easy-to-
!     identify variables.
!      m02 = pxldat(2)
!      m26 = pxldat(26)
!      m22 = pxldat(22)
!      m27 = pxldat(27)
!      m31 = pxldat(31)
!      m32 = pxldat(32)

      m02 = pxldat(4)
      m26 = pxldat(19)
      m22 = pxldat(21) ! 4.05 replace 3.959
!      m27 = pxldat(27)
      m31 = pxldat(24)
      m32 = pxldat(25)
      masir11 = m31
      masir12 = m32
      
!     Initialize test group confidences.
      cmin1 = 1.0
      cmin2 = 1.0
      cmin3 = 1.0
      cmin4 = 1.0

!     Initialize array containing number of tests in each test group.
      do 10 kk = 1 , 5
         ngtests(kk) = 0
  10  continue

!----------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Subroutine PolarDay_desert_c '',
!     +                      /)')
!      endif

!----------------------------------------------------------------------

!     Begin clear-sky tests.

!----------------------------------------------------------------------

!     GROUP 1 TESTS 

!   pfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11 < pfmft_11maxthre(1)) .and.   &
 !         (masir11-masir12) < pfmft_btd_min(1) ) then
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
      
      
!     Water vapor channel (6.7 micron) high cloud test.

!      if (nint(m27) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (m27 .gt. pdsh20_c(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(m27,pdsh20_c(1),pdsh20_c(3),pdsh20_c(4),   &
!                       pdsh20_c(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

!----------------------------------------------------------------------

!     debug statement
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''6.7um test: '',6f10.2)') m27,pdsh20_c(1),
!     +          pdsh20_c(2),pdsh20_c(3),c2,cmin1
!      endif

!----------------------------------------------------------------------
 
!     GROUP 2 TESTS

!     11-12um brightness temperature difference test for thin cirrus.

      if (nint(m31) .ne. nint(bad_data) .and.   &
          nint(m32) .ne. nint(bad_data)) then

        m31_32 = m31 - m32

!       Calculate secant of viewing zenith angle.
        cosvza = cos(vza*dtr)
        if (abs(cosvza) .gt. Rel_equality_EPS) then
          schi = 1.0/cosvza
        else
          schi = 99.0
        end if

!       Interpolate cloud threshold table of 11-12 um differences
!       which are defined as a function of viewing zenith angle and
!       11 um brightness temperature.
        call tview(1,schi,m31,diftemp)

!       If a valid threshold was found, use it as the test threshold,
!       otherwise use a static threshold.
        if ((diftemp .lt. 0.1) .or. (abs(schi-99.0) .lt. 0.0001)) then
          dfthrsh = pds11_12hi_c(2)
        else
          dfthrsh = diftemp
        end if

!       Perform thin cirrus test.
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (m31_32 .le. dfthrsh) then
          call set_bit(testbits,18)
        nptests = nptests + 1
        end if
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - (0.3 * dfthrsh)
        call conf_test(m31_32,locut,hicut,1.0,dfthrsh,1,c3)
        cmin2 = min(cmin2,c3)
        ngtests(2) = ngtests(2) + 1

      endif

!----------------------------------------------------------------------

!     debug statement
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO m31_32: '',6f10.2)') m31_32,
!     +          dfthrsh,locut,hicut,c3,cmin2
!      endif

!----------------------------------------------------------------------
  
!     11-4um brightness temperature difference test for fog and low
!     clouds.

      if (visusd) then
        if (nint(m31) .ne. nint(bad_data) .and.  &
            nint(m22) .ne. nint(bad_data)) then
          if(m31 .le. 320.0) then

            nmtests = nmtests + 1
            m31_22 = m31 - m22
            call set_qa_bit(qa_bits,19)
            if (m31_22 .ge. pds11_4lo_c(2) .and.  &
                m31_22 .le. pds11_4hi_c(2)) then
              call set_bit(testbits,19)
              nptests = nptests + 1
            end if

            locuta(1) = pds11_4lo_c(1)
            locuta(2) = pds11_4hi_c(1)
            hicuta(1) = pds11_4lo_c(3)
            hicuta(2) = pds11_4hi_c(3)
            midpta(1) = pds11_4lo_c(2)
            midpta(2) = pds11_4hi_c(2)
  
            call conf_test_2val(m31_22,locuta,hicuta,1.0,midpta,2,c4)
            cmin2 = min(cmin2,c4)
            ngtests(2) = ngtests(2) + 1

          endif
        end if

!----------------------------------------------------------------------

!       debug statement
!        if (debug .gt. 0) then
!        write(h_output,'(1x,''11-4um a: '',4f10.2)')m31_22,
!     +        pds11_4lo_c(1),pds11_4lo_c(2),pds11_4lo_c(3)
!        write(h_output,'(1x,''11-4um b: '',6f10.2)')m31_22,
!     +        pds11_4hi_c(1),pds11_4hi_c(2),pds11_4hi_c(3),c4,cmin2
!        endif

!----------------------------------------------------------------------

      endif

!----------------------------------------------------------------------

!     GROUP 3 TESTS

!     0.86 micron reflectance threshold test.

      if (visusd) then
        if (nint(m02) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20)
          if (m02 .le. pdsref2_c(2)) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(m02,pdsref2_c(1),pdsref2_c(3),pdsref2_c(4),  &
                         pdsref2_c(2),1,c5)
          cmin3 = min(cmin3,c5)
          ngtests(3) = ngtests(3) + 1
        end if

!----------------------------------------------------------------------

!       debug statement
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''.86um test : '',6f10.2)') m02,
!     +            pdsref2_c(1),pdsref2_c(2),pdsref2_c(3),c5,cmin3
!        endif

!----------------------------------------------------------------------

      endif

!     GROUP 4 TESTS 

!     Near-infrared high cloud test.

      if (visusd .and. (.not. hi_elev) ) then
        if (nint(m26) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (m26 .le. pdsref3_c(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(m26,pdsref3_c(1),pdsref3_c(3),pdsref3_c(4),   &
                         pdsref3_c(2),1,c6)
          cmin4 = min(cmin4,c6)
          ngtests(4) = ngtests(4) + 1
        end if
 
!----------------------------------------------------------------------

!       debug statement
!        if (debug .gt. 0) then
!           write(h_output,'(1x,''1.38um test: '',6f10.4)')m26,
!     +           pdsref3_c(1),pdsref3_c(2),pdsref3_c(3),c6,cmin4
!        endif

!----------------------------------------------------------------------

      endif

!----------------------------------------------------------------------

!     Thin cirrus test.

      if (visusd .and. (.not. hi_elev) ) then
        if (nint(m26) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (m26 .lt. pdstci_c(1) .and. m26 .ge. pdstci_c(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
        endif
      endif
 
!----------------------------------------------------------------------

!     Determine final confidence based on group values.
      pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4

!     Next, find the number of test groups used for the pixel.
      groups = 0
      do kk = 1,5
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      fac = 1.0
      if (groups .gt. 0) fac = 1.0 / groups

!     Final pixel confidence is nth root of group confidence product.
      confdnc = pre_confdnc**fac

!----------------------------------------------------------------------

!     debug statement
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',7f8.5/,4f8.5)') c2,c3,c4,c5,
!     +         c6,cmin1,cmin2,cmin3,cmin4,fac,confdnc
!      endif

!----------------------------------------------------------------------

      return
      end
