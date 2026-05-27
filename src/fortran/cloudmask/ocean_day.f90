      subroutine ocean_day(pxldat,vza,snglnt,visusd,cirrus_vis,sfctmp,   &
                           refang,sh_ocean,testbits,qa_bits,nmtests,     &
                           confdnc,btclr)

       implicit none
       save 

!---------------------------------------------------------------------
!!F77 
!
!!Description:
! Performs clear-sky spectral tests for water surfaces during daylight
! conditions.
!
! Each spectral test is placed in one of five test groups. The groups
! represented in this routine are:
!
!      
!          Group 1: High thick cloud
!                   11 micron bt test 
!                   13.9 micron bt test  
!                   6.75 micron bt test 
!
!          Group 2: Low cloud - thick
!                   8-11 micron and 11-12 micron bt tests
!                   11-4 micron bt tests
!                  
!          Group 3: Thick cloud
!                   .87 micron reflectance test
!                   .87/.66 micron reflectance ratio test
! 
!          Group 4: Thin cirrus test
!                   1.38 micron reflectance test 
!
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
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Current pixel viewing angle
! snglnt        Logical variable indicating sun glint contamination
! visusd        Logical variable indicating whether vis data used or not
! cirrus_vis 	Logical variable flagging thin cirrus contaminated
! 		scenes in the visible
! sfctmp        SST from ancillary data
! refang        Reflectance angle
! sh_ocean      Logical flag indicating ocean depths < 50 m or within 5 km
!               of shoreline
!
!!Output Parameters:
! testbits      6 byte array containing bit results
! qa_bits       10 byte array contining qa bit results
! nmtests 	 Counts number of inidividual tests applied
! confdnc       product of all applied individual confidences
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Implemented new 11-12 um thin cirrus test (Key version)
! 10/04 Collection 5b R. Frey
! Added SST test.
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!    Externals:
!        Subroutines conf_test_2val,conf_test,tview,set_bit,
!                    clear_bit,check_bits,set_qa_bit
!        Functions rega,regb
!
!!END
!-------------------------------------------------------------------

!     Declarations.

      include 'global.inc'
      include 'ocean_day_thr.inc'
      include 'snglntr_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
!     Scalar arguments. 
      real confdnc,vza,refang,sfctmp
      logical visusd,snglnt,cirrus_vis,sh_ocean
      integer nmtests

!     Array arguments. 
      real pxldat(inband),btclr(7),tv11_12,masir11,masir12
      byte testbits(6),qa_bits(10)
        
!     Local scalars.
      real c1,c2,c3,c4,c6,c7,c8,c9,c11,cosvza,dfthrsh,diftemp,    &
           dtr,r24_21,r24_25,r23_24,r24,r25,m35,r21,r23,r19,      &
           pi,schi,vrat,m27,Rel_equality_EPS,r03,r04,             &
           cmin1,cmin2,cmin3,cmin4,tri_thres,fac,groups,          &
           pre_confdnc,locut,hicut,midpt,power,max_vza,a,corr,    &
           sfcdif,c10,sst_thrsh
      integer nptests,rtn,debug,h_output,kk

!     Local arrays.
      real hicuta(2),locuta(2),midpta(2)
      integer ngtests(5)
        
!     Parameter statements.
      parameter (Rel_equality_EPS = 0.000001)
      parameter(max_vza = 65.49)

!     External functions.
      real trispc 
      external trispc
         
!     External subroutines.
      external conf_test,tview,set_bit,clear_bit,conf_test_2val,set_qa_bit,get_sg_thresholds
!    
!     Intrinsic functions.
      intrinsic acos,cos

!     Common statement for debug purposes
!      common / bug / debug, h_output

!-------------------------------------------------------------------

!     Initialize variables.

      pi = acos(-1.0)
      dtr = pi/180.0

!     'nmtests' counts the number of tests applied to this pixel.
      nmtests = 0

!     'nptests' counts the number of tests which found no evidience 
!     of cloud.
      nptests = 0

!     Place reflectance and brightness temperature values into easy-to-
!     identify variables.
!      m01g = pxldat(1)
!      m02g = pxldat(2)
!      m26g = pxldat(26)
!      m22g = pxldat(22)
!      m27g = pxldat(27)
!      m29g = pxldat(29)
!      m31g = pxldat(31)
!      m32g = pxldat(32)
!      m35g = pxldat(35)
      r03 = pxldat(3)  ! 0.65
      r04 = pxldat(4)  ! 0.86
      r19 = pxldat(19) ! 1.38
      r21 = pxldat(21) ! 3.959 (c22 of modis) replaced by MERSI_II 4.05
!      m27 = pxldat(27) ! 6.5 
      r23 = pxldat(23) ! 8.5 
      r24 = pxldat(24) ! 11
      r25 = pxldat(25) ! 12
      masir11 = r24
      masir12 = r25
!      m35 = pxldat(35)

!     Initialize test group confidences.
      cmin1 = 1.0
      cmin2 = 1.0
      cmin3 = 1.0
      cmin4 = 1.0
 
!     Initialize array containing number of tests in each test group.
      do 10 kk = 1 , 5
         ngtests(kk) = 0 
  10  continue

!-------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine ocean_day '',/)')
!      endif

!-------------------------------------------------------------------

!     Begin clear sky tests.

!-------------------------------------------------------------------

!     GROUP 1 TESTS

!     11 micron brightness temperature threshold test.
      !print*, 'ocean day',dobt11(1),dobt11(3),dobt11(4),dobt11(2)
      if (nint(r24) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,13)
        if (r24 .ge. dobt11(2)) then
          call set_bit(testbits,13)
          nptests = nptests + 1
        end if
        call conf_test(r24,dobt11(1),dobt11(3),dobt11(4),dobt11(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

!-------------------------------------------------------------------

!     debug statement
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''r24: '',5f10.2)') r24,dobt11(1),
!     +          dobt11(2),dobt11(3),dobt11(4)
!      endif

!-------------------------------------------------------------------

!   pfmft test
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          (masir11 < pfmft_11maxthre(1)) .and.   &
          (masir11-masir12) < pfmft_btd_min(1) ) then
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
     !     call conf_test(tv11_12,pfmft_ocean(1),pfmft_ocean(3),pfmft_ocean(4),   &
     !                    pfmft_ocean(2),1,c2)
     !   cmin1 = min(cmin1,c2)          ! annotation by minmin (to close this threshold)
     !   ngtests(1) = ngtests(1) + 1
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
      !  call conf_test(tv11_12,nfmft_ocean(1),nfmft_ocean(3),nfmft_ocean(4),   &
      !                 nfmft_ocean(2),1,c3)
      !  cmin1 = min(cmin1,c3)           ! annotation by minmin (to close this threshold)
      !  ngtests(1) = ngtests(1) + 1
      endif
       
!     co2 high cloud test

!      if (nint(m35) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,14)
!        if (m35 .gt. doco2(2)) then
!          call set_bit(testbits,14)
!          nptests = nptests + 1
!        end if
!        call conf_test(m35,doco2(1),doco2(3),doco2(4),
!     *                doco2(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

!-------------------------------------------------------------------
 
!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''m35: '',5f10.2)') m35,doco2(1),
!     +          doco2(2),doco2(3),doco2(4)
!      endif

!-------------------------------------------------------------------

!     H20 vapor channel (6.7 micron) high cloud test 

!      if (nint(m27) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (m27 .gt. doh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(m27,doh20(1),doh20(3),doh20(4),
!     *                doh20(2),1,c3)
!        cmin1 = min(cmin1,c3)
!        ngtests(1) = ngtests(1) + 1
!      endif

!-------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''m27: '',5f10.2)') m27,doh20(1),
!     +          doh20(2),doh20(3),doh20(4)
!      endif

!-------------------------------------------------------------------

! ... SST test

      if ( (nint(r24) .ne. nint(bad_data)) .and.    &
           (nint(r25) .ne. nint(bad_data)) ) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

        if(sh_ocean) then
          sst_thrsh = 10.0
        else
          sst_thrsh = 6.0
        end if

        r24_25 = r24 - r25
        if(r24_25 .ge. 1.0) then
          midpt = sst_thrsh + (2.0 * nint(r24_25))
        else
          midpt = sst_thrsh
        end if

        a = vza / max_vza
        corr = (a**4) * 3.0
        midpt = midpt + corr
        locut = midpt + 1.0
        hicut = midpt - 2.0

        sfcdif = sfctmp - r24

        if( sfcdif .lt. midpt ) then
          call set_bit(testbits,27)
          nptests = nptests + 1
        end if

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c10)
        cmin1 = min(cmin1,c10)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

!-----------------------------------------------------------------

!     debug statement
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''sfctmp: '',9f10.3)') sfctmp,r24,
!     +           sfcdif,a,corr,locut,midpt,hicut,c10
!      endif

!-----------------------------------------------------------------

!     GROUP 2 TESTS

!     tri-spectral tests - 8, 11 and 12 micron BTDIF's

      if (nint(r24) .ne. nint(bad_data) .and. &
          nint(r25) .ne. nint(bad_data) .and. &
          nint(r23)  .ne. nint(bad_data)) then
 
         r23_24 = r23 - r24
         r24_25 = r24 - r25

!       Get clear sky 8-11 micron clear sky thresholds based
!       upon 11-12 difference and compare to global regressions
!       determined from global HIRS data
        tri_thres = trispc(r24_25) 
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (r23_24 .lt. tri_thres) then
          nptests = nptests + 1
          call set_bit(testbits,18)
        end if
        locut = tri_thres + .5
        hicut = tri_thres - .5
        call conf_test(r23_24,locut,hicut,1.0,tri_thres,1,c4)
!        cmin2 = min(cmin2,c4)            ! revised by wuxiao
!        ngtests(2) = ngtests(2) + 1
      endif
 
!-------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''trispec: '',5f10.2)') r24_25,r23_24,
!     +          locut,tri_thres,hicut
!      endif

!-------------------------------------------------------------------

!       11-12um brightness temperature difference test
!       for thin cirrus.

      if (nint(r24) .ne. nint(bad_data) .and.   &
          nint(r25) .ne. nint(bad_data) .and.   &
          vza .gt. 0.0) then

        r24_25 = r24 - r25

!       calculate secant of viewing zenith angle.
        cosvza = cos(vza*dtr)
        if (abs(cosvza).gt.Rel_equality_EPS) then
          schi = 1.0/cosvza
        else
          schi = 99.0
        end if
 
!       Interpolate look-up table values of 11 - 12 micron bt
!       difference thresholds (function of viewing zenith
!       and 11 micron brightness temperature).
        call tview(1,schi,r24,diftemp)
 
        if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
          dfthrsh = do11_12hi(1)
        else
          dfthrsh = diftemp
        end if
 
!       Since the IR BTDIF bit has possibly been set already,
!       change the bit only if the current test failed.
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (r24_25 .le. dfthrsh) then
          nptests = nptests + 1
        else
          call check_bits(testbits,18,rtn)
          if (rtn .eq. 1) then
            call clear_bit(testbits,18)
          end if
        endif
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - 1.25
        call conf_test(r24_25,locut,hicut,1.0,dfthrsh,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif

!-------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO r24_25: '',4f10.2)') r24_25,
!     +          dfthrsh,locut,hicut
!      endif

!-------------------------------------------------------------------

!     11 minus 4 micron BTDIF fog and low cloud test.

      if (visusd .and. .not. snglnt) then
        if (nint(r24) .ne. nint(bad_data) .and.  &
            nint(r21) .ne. nint(bad_data)) then

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19) 
          r24_21 = r24 - r21

          if (r24_21 .ge. do11_4lo(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if

          call conf_test(r24_21,do11_4lo(1),do11_4lo(3),do11_4lo(4),  &
                         do11_4lo(2),1,c7)

          cmin2 = min(cmin2,c7)
          ngtests(2) = ngtests(2) + 1

        endif

!-------------------------------------------------------------------

!       debug statement 
!        if (debug .gt. 0) then
!         write(h_output,'(1x,''r24_21: '',5f10.2)')r24_21,do11_4lo(1),
!     +            do11_4lo(2),do11_4lo(3),do11_4lo(4)
!        endif

!-------------------------------------------------------------------

      end if

!     GROUP 3 TESTS 

!     NIR reflectance threshold test.

      if (visusd) then
        if (nint(r04) .ne. nint(bad_data)) then

!         Take into account sunglint problems
          if(snglnt) then
            call get_sg_thresholds(refang,locut,hicut,midpt,power)
          else
            locut = doref2(1)
            hicut = doref2(3)
            midpt = doref2(2)
            power = doref2(4)
          end if

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20) 
          if (r04.le.midpt) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(r04,locut,hicut,power, &
                         midpt,1,c8)
          cmin3 = min(cmin3,c8)
          ngtests(3) = ngtests(3) + 1
        endif

!-------------------------------------------------------------------

!       debug statement 
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''r04: '',6f10.4)') r04,locut,
!     +            hicut,midpt,power,refang
!        endif

!-------------------------------------------------------------------

      end if

!     Visible channel ratio test 

      if (visusd) then
        if (nint(r03) .ne. nint(bad_data) .and. &
            nint(r04) .ne. nint(bad_data)) then

!         Account for sun glint contamination
          if (snglnt) then
            locuta(1) = snglntvcl(1)
            locuta(2) = snglntvcl(2)
            hicuta(1) = snglntvch(1)
            hicuta(2) = snglntvch(2)
            midpta(1) = snglntv(1)
            midpta(2) = snglntv(2)
          else
            locuta(1) = dovratlo(1)
            locuta(2) = dovrathi(1)
            hicuta(1) = dovratlo(3)
            hicuta(2) = dovrathi(3)
            midpta(1) = dovratlo(2)
            midpta(2) = dovrathi(2)
          end if

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,21) 
          vrat = r04 / r03
          if (vrat .lt. midpta(1) .or. vrat .gt. midpta(2)) then
            call set_bit(testbits,21)
            nptests = nptests + 1
          end if
          call conf_test_2val(vrat,locuta,hicuta,1.0,midpta,2,c9)
          cmin3 = min(cmin3,c9)
          ngtests(3) = ngtests(3) + 1
        endif

!-------------------------------------------------------------------

!       debug statement 
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''vrat: '',7f10.2)') vrat,locuta(1),
!     +            locuta(2),hicuta(1),hicuta(2),midpta(1),midpta(2)
!        endif

!-------------------------------------------------------------------

      end if

!     GROUP 4 TESTS 

! ... Near-infrared high cloud test.

      if (visusd) then
        if (nint(r19) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16) 
          if (r19 .le. doref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(r19,doref3(1),doref3(3),doref3(4),  &
                         doref3(2),1,c11)
          cmin4 = min(cmin4,c11)
          ngtests(4) = ngtests(4) + 1
        endif

!-------------------------------------------------------------------

!       debug statement 
!        if (debug .gt. 0) then
!           write(h_output,'(1x,''r19: '',6f10.4)')r19,doref3(1),
!     +                doref3(2),doref3(3),doref3(4)
!        endif

!-------------------------------------------------------------------

      end if

!     Thin cirrus test.

      if (visusd) then
        if (nint(r19) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (r19 .lt. dotci(1) .and. r19 .ge. dotci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif

!-------------------------------------------------------------------

!         debug statement 
!          if (debug .gt. 0) then
!             write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)') r19,
!     +             dotci(1),dotci(2)
!          endif

!-------------------------------------------------------------------

        endif
      endif
 
!-------------------------------------------------------------------

!     Determine initial confidence based on group values.
      pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4
      
!     Find the number of test groups used for current pixel.
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
      !confdnc = cmin4
!-------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',7f8.5/,8f8.5)') c1,c2,c3,c4,
!     +         c6,c7,c8,c9,c11,cmin1,cmin2,cmin3,cmin4,fac,confdnc
!      endif

!-------------------------------------------------------------------

      return
      end
