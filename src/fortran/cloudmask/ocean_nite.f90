      subroutine ocean_nite(indat,kele,pxldat,vza,sfctmp,sh_ocean,   &
                            uniform,testbits,qa_bits,nmtests,confdnc,btclr)

!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over water
!      surfaces during nightime hours.
!
!      For nighttime ocean the groups are:
!          Group 1: High thick cloud
!                   11 micron bt test 
!                   13.9 micron bt test  
!                   6.75 micron bt test 
!                   SST test
!
!          Group 2: Btdif tests
!                   8-11 micron and 11-12 micron bt tests
!                   11-4 micron bt test
!                   8.6-7.3 micron bt test
!                   11 micron variability test
!
!!Input Parameters:
! indat         Array containing reflectance, brightness temperatures for
!                 three complete scan lines
! kele          Position of first pixel in current context
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           viewing zenith angle in degrees
! sfctmp        SST for current pixel
! sh_ocean      Logical flag indicating ocean depths < 50 m or within 5 km
!               of shoreline
! uniform       Logical variable indicating uniform conditions over context
!
!!Output Parameters:
! testbits      6-byte array containing bit results
! qa_bits       10-byte array containing qa bit results
! nmtests       Number of tests applied for this pixel
! confdnc       product of all applied individual confidences
!
!!Revision History:
! 06/04 Collection 5  R. Frey:
! Added SST test
! Implemented new 11 um variability test
! Added 8.6-7.3 test
! Implemented new 11-12 um thin cirrus test (J. Key version)
! 10/04 Collection 5b R. Frey:
! Added shallow ocean condition on choice of SST threshold.
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!Design Notes:
!    Externals:
!       Subroutines conf_test,set_bit,set_qa_bit,tview,check_bits,
!                   clear_bit,chk_spatial2
!       Function trispc
!
!!END
!---------------------------------------------------------------------

      implicit none
      save 

      include 'global.inc'
      include 'ocean_nite_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real vza,sfctmp,confdnc
      integer nmtests,kele
      logical sh_ocean,uniform
! ...
! ... array arguments ..
      real pxldat(inband),indat(necntx,nlcntx,inband),btclr(7),tv11_12
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..

      real c1,c2,c3,c4,mas11_4,masdf1,                            &
           masdf2,masir11,masir12,masir13,masir4,masir8,          &
           masir65,c6,cmin1,cmin2,groups,masir73,                 &
           fac,pre_confdnc,tri_thres,c7,diftemp,dfthrsh,schi,     &
           cosvza,dtr,pi,dwvs,midpt,locut,hicut,c9,               &
           c10,sfcdif,a,corr,c11,np,sst_thrsh
      integer nptests,rtn,kk,debug,h_output,npix

      real, parameter :: Rel_equality_EPS = 0.000001
      real, parameter :: max_vza = 65.49
! ...
!     local arrays
      integer ngtests(2)

! ... external functions ..
      real trispc
      external trispc
! ...
! ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview,check_bits, &
               clear_bit,chk_spatial2

! ... Common statement for debug purposes
!      common / bug / debug, h_output

! ... compute degrees-to-radians conversion
      pi = acos(-1.0)
      dtr = pi/180.0
! ... nmtests counts the number of tests applied to this pixel
      nmtests = 0
! ... nptests counts the number of tests passed
      nptests = 0
!     ngtests counts the number of tests applied per group
      ngtests(1) = 0
      ngtests(2) = 0
! ... confidence to 1.0 to begin with
      confdnc = 1.0

! ... place band values into individual variables for easy
! ... identification
!      masir4 = pxldat(22)
!      masir65 = pxldat(27)
!      masir73 = pxldat(28)
!      masir8 = pxldat(29)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)
!      masir13 = pxldat(35)

      masir4 = pxldat(21) ! 4.05 replace 3.959 
!      masir65 = pxldat(27)
      masir73 = pxldat(22)
      masir8 = pxldat(23)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
!      masir13 = pxldat(35)

! ...
      masdf2 = 0.0
      masdf1 = 0.0
      mas11_4 = 0.0
      dfthrsh = 0.0
      schi = 0.0
      diftemp = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c6 = 0.0
      c7 = 0.0
      c9 = 0.0
      c10 = 0.0
      c11 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
 

! ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine ocean_nite '',/)')
      endif
! ................................................................

!     **** GROUP 1 TESTS *************************************
!     11 micron brightness temperature threshold test
      if (nint(masir11) .ne. nint(bad_data)) then
    !    nmtests = nmtests + 1
        call set_qa_bit(qa_bits,13)
! ...   compare to daytime ocean threshold, set clear bit if passed
        if (masir11 .ge. nobt11(2)) then
          call set_bit(testbits,13)
          nptests = nptests + 1
        end if
! ...   calculate confidence compared to low and high confidence cutoffs
        call conf_test(masir11,nobt11(1),nobt11(3),nobt11(4),nobt11(2),1,c1)
        cmin1 = min(cmin1,c1)          ! annotation by minmin (to close this threshold)
        ngtests(1) = ngtests(1) + 1    ! annotation by minmin (to close this threshold)
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir11: '',5f10.2)') masir11,
!     +          nobt11(1),nobt11(2),nobt11(3),nobt11(4)
!      endif
! ................................................................

! === PFMFT test disabled (btclr requires NWP RTM) ===
!       if (nint(masir11) .ne. nint(bad_data) .and.   &
!           nint(masir12) .ne. nint(bad_data) .and.   &
!           (masir11 < pfmft_11maxthre(1)) .and.   &
! !          (masir11-masir12) < pfmft_btd_min(1) ) then
!           (btclr(5)-btclr(6)) > pfmft_btd_min(1) ) then          !jincheng
! !		  nmtests = nmtests + 1
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
!           call conf_test(tv11_12,pfmft_ocean(1),pfmft_ocean(3),pfmft_ocean(4),   &
!                          pfmft_ocean(2),1,c2)
! !        cmin1 = min(cmin1,c2)
! !        ngtests(1) = ngtests(1) + 1
! === PFMFT test disabled end ===
      
! === NFMFT test disabled (btclr requires NWP RTM) ===
!       if (nint(masir11) .ne. nint(bad_data) .and.   &
!           nint(masir12) .ne. nint(bad_data) .and.   &
!           (masir11-masir12) <= nfmft_maxthre(1) ) then
! !        nmtests = nmtests + 1
!         !tv11_12 = (btclr(5) - btclr(6)) - (masir11 - masir12)
!         tv11_12 =  (masir11 - masir12) - (btclr(5) - btclr(6))
!         call set_qa_bit(qa_bits,15)
!         !if (masir11.gt.dlco2(2)) then
!           call set_bit(testbits,15)
!           nptests = nptests + 1
!         !end if
!         call conf_test(tv11_12,nfmft_ocean(1),nfmft_ocean(3),nfmft_ocean(4),   &
!                        nfmft_ocean(2),1,c3)
!      !   cmin1 = min(cmin1,c3)      ! annotation by minmin (to close this threshold)
!      !   ngtests(1) = ngtests(1) + 1
! === NFMFT test disabled end ===
      
! ... co2 high cloud test
!      if (nint(masir13) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,14)
!        if (masir13.gt.noco2(2)) then
!          call set_bit(testbits,14)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir13,noco2(1),noco2(3),noco2(4),noco2(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,noco2(1),
!     +          noco2(2),noco2(3),noco2(4)
!      endif
! ................................................................


!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65.gt.noh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,noh20(1),noh20(3),noh20(4),noh20(2),1,c3)
!        cmin1 = min(cmin1,c3)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,noh20(1),
!     +          noh20(2),noh20(3),noh20(4)
!      endif
! ................................................................

! ... SST test

      if ( (nint(masir11) .ne. nint(bad_data)) .and.   &
           (nint(masir12) .ne. nint(bad_data)) ) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

        if(sh_ocean) then
          sst_thrsh = 10.0
        else
          sst_thrsh = 6.0
        end if

        masdf1 = masir11 - masir12
        if(masdf1 .ge. 1.0) then
          midpt = sst_thrsh + (2.0 * nint(masdf1))
        else
          midpt = sst_thrsh
        end if

        a = vza / max_vza
        corr = (a**4) * 3.0
        midpt = midpt + corr
        locut = midpt + 1.0
        hicut = midpt - 2.0

        sfcdif = sfctmp - masir11
        
        if( sfcdif .lt. midpt ) then
          call set_bit(testbits,27)
          nptests = nptests + 1
        end if

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c10)
        cmin1 = min(cmin1,c10)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''sfctmp: '',9f10.3)') sfctmp,masir11,
!     +           sfcdif,a,corr,locut,midpt,hicut,c10
!      endif
! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
 
 
!     ****  GROUP 2 TESTS  ***********************************
! ... tri-spectral tests - 8, 11 and 12 micron BTDIF's
! ... calculate 8 minus 11 and 11 minus 12 micron BTDIFs
      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir12) .ne. nint(bad_data) .and.  &
          nint(masir8)  .ne. nint(bad_data)) then

        masdf2 = masir8 - masir11
        masdf1 = masir11 - masir12

!       Get clear sky 8-11 micron clear sky thresholds based
!       upon 11-12 difference and compare to global regressions
!       determined from global HIRS data
        tri_thres = trispc(masdf1)
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf2.lt.tri_thres) then
          nptests = nptests + 1
          call set_bit(testbits,18)
        end if
        locut = tri_thres + .5
        hicut = tri_thres - .5
        call conf_test(masdf2,locut,hicut,1.0,tri_thres,1,c4)
!        cmin2 = min(cmin2,c4)           ! revised by wuxiao
!        ngtests(2) = ngtests(2) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''trispec: '',6f10.2)') masdf1,masdf2,
!     +          tri_thres,locut,hicut,c4
!      endif
! ................................................................


      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir12) .ne. nint(bad_data) .and.  &
          vza .gt. 0.0) then

        masdf1 = masir11 - masir12

! ...   11-12um brightness temperature difference test
! ...   for thin cirrus).
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
           dfthrsh = no11_12hi(1)
        else
           dfthrsh = diftemp
        end if
 
! ...   Since the IR BTDIF testbit was already potentially set,
! ...   then only change the bit if the current test failed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (masdf1.le.dfthrsh) then
          nptests = nptests + 1
        else
          call check_bits(testbits,18,rtn)
          if (rtn .eq. 1) then
            call clear_bit(testbits,18)
          end if
        endif
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - 1.25
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c7)
        cmin2 = min(cmin2,c7)
        ngtests(2) = ngtests(2) + 1
      endif
!print*,'c2-1',cmin2
! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',6f10.2)') masdf1,
!     +          no11_12hi(1),dfthrsh,locut,hicut,c7
!      endif
! ................................................................


! ... 11 minus 4 micron BTDIF fog and low cloud test.
! ... for now placing in the SWIR bit place holder
      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir4) .ne. nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
        mas11_4 = masir11 - masir4
        if (mas11_4.le.no11_4lo(2)) then
          call set_bit(testbits,19)
          nptests = nptests + 1
        end if
        call conf_test(mas11_4,no11_4lo(1),no11_4lo(3),no11_4lo(4),no11_4lo(2),1,c6)
!        cmin2 = min(cmin2,c6)        ! revised by wuxiao
!        ngtests(2) = ngtests(2) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''mas11_4: '',6f10.2)')mas11_4,no11_4lo(1),
!     +          no11_4lo(2),no11_4lo(3),no11_4lo(4),c6
!      endif
! ................................................................

! ... Water vapor cloud test.

      if( nint(masir73) .ne. nint(bad_data) .and. &
          nint(masir8) .ne. nint(bad_data)) then

        dwvs = masir8 - masir73

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,29)

        if ( dwvs .gt. no86_73(2) ) then
          nptests = nptests + 1
          call set_bit(testbits,29)
        end if

        call conf_test(dwvs,no86_73(1),no86_73(3),no86_73(4),no86_73(2),1,c9)
        cmin2 = min(cmin2,c9)
        ngtests(2) = ngtests(2) + 1

      endif
!print*,'c2-2',c9,dwvs
! ... debug statement ............................................
!      if (debug .gt. 0) then
!         write(h_output,'(1x,''dwvs: '',5f6.2)') dwvs,no86_73(1),
!     +      no86_73(2),no86_73(3),c9
!      endif
! ................................................................

! ... Variability test

      if (uniform) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,30)

        call chk_spatial2(indat,kele,npix)
        np = npix * 1.0

        if ( np .gt. no_11var(2) ) then
           nptests = nptests + 1
           call set_bit(testbits,30)
        end if

        call conf_test(np,no_11var(1),no_11var(3),no_11var(4),no_11var(2),1,c11)
        cmin2 = min(cmin2,c11)
        ngtests(2) = ngtests(2) + 1
!print*,'c2-3',c11
! ..... debug statement ..........................................
!        if (debug .gt. 0) then
!          write(h_output,'(1x,''var: '',5f10.3)') np,no_11var(1),
!     +            no_11var(2),no_11var(3),c11
!        endif
! ................................................................

      end if

! *******     END OF GROUP 2 TESTS ****************************
!
!
!     Determine final confidence based on group values
!cmin2 = 1.
!cmin1 = 1.
      pre_confdnc = cmin1 * cmin2 
      groups = 0.0
      do kk = 1,2
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if (groups .gt. 0) then
        fac = 1.0 / groups
      else
        fac = 0.0
      end if
      !print*,'test',cmin1,cmin2,fac,uniform,pre_confdnc**fac
      confdnc = pre_confdnc**fac

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',9f8.5/,4f8.5)') c1,c2,c3,c4,
!     +         c6,c7,c9,c10,c11,cmin1,cmin2,fac,confdnc
!      endif
! ................................................................
!
      return
      end
