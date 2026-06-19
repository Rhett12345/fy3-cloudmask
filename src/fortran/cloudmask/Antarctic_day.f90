      subroutine Antarctic_day(pxldat,visusd,testbits,qa_bits,   &
                               nmtests,confdnc,btclr,is_cold_sfc)

      implicit none
      save
 
!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over Antarctic snow 
!      surfaces during daylight hours.
!
!      For daytime land type 1 the groups are:
!          Group 1: High thick cloud
!                   6.7 micron bt test
!
!          Group 2: Low cloud - thick
!                   11-4 micron bt tests
!        
!          Group 4: Thin cirrus test
!                   1.38 micron reflectance test 
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! visusd        Logical variable indicating whether vis data used or not
!
!
!!Output Parameters:
! testbits      six byte integer containing bit results
! qa_bits       ten byte integer containing qa bit results
! nmtests       number spectral tests performed
! confdnc       product of all applied individual confidences
!
!!Revision History:
!
! Added 11-12 um thin cirrus test
! 06/04 Collection 5    R. Frey
! Removed 11-12 um thin cirrus test.
! Added 11 um BT-dependent 3.0-11 um BTD test.  Replaces test with static
! thresholds.
! 10/04 Collection 5b   R. Frey
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
      include 'Antarctic_day_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real confdnc
      integer nmtests
      logical visusd
! ...
! ... array arguments ..
      real pxldat(inband),btclr(7),tv11_12
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,c3,mas11_4,cmin1,cmin2,        &
           masir11,masir4,locut,hicut,       &
           masir65,groups,fac,pre_confdnc,midpt,power,masir12
      integer nptests,kk,debug,h_output

! ... local arrays
      integer ngtests(3)
     
      real, parameter :: Rel_equality_EPS = 0.000001

! ... external subroutines ..
      external conf_test,set_bit,clear_bit,set_qa_bit

!     Common statement for debug purposes
!      common / bug / debug, h_output

! ...
! ... nmtests counts the number of tests applied to this pixel
      nmtests = 0

! ... initialize variables

! ... nptests counts the number of tests passed
      nptests = 0
! ... set confidence to 1.0 to begin with
      confdnc = 1.0
! ... place band values into individual variables for easy
! ... identification
!      masir4 = pxldat(22)  !
!      masir65 = pxldat(27)
!      masir11 = pxldat(31)
      masir4 = pxldat(20)  ! 4.05 replace 3.959
!      masir65 = pxldat(27)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
! ...
      mas11_4 = 0.0

! ... the ! suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0

! ... initialize group number holder
      do 10 kk = 1 , 2
         ngtests(kk) = 0
  10  continue

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing subroutine Antarctic_day '',
!     +                   /)')
!      endif
! ................................................................

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',8f10.2)') masdf1,
!     +          dps11_12hi(1),dps11_12adj(1),
!     +          masir11,schi,dfthrsh,locut,hicut
!      endif
! ................................................................

! === PFMFT test disabled (btclr requires NWP RTM) ===
!       if (nint(masir11) .ne. nint(bad_data) .and.   &
!           nint(masir12) .ne. nint(bad_data) .and.   &
!           (masir11 < pfmft_11maxthre(1)) .and.   &
! !          (masir11-masir12) < pfmft_btd_min(1) ) then
!           (btclr(5)-btclr(6)) > pfmft_btd_min(1) ) then          !jincheng
! !        nmtests = nmtests + 1
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
!         if (is_cold_sfc == 1) then
!            call conf_test(tv11_12,pfmft_cold(1),pfmft_cold(3),pfmft_cold(4),   &
!                           pfmft_cold(2),1,c1)        
!         else
!            call conf_test(tv11_12,pfmft_land(1),pfmft_land(3),pfmft_land(4),   &
!                           pfmft_land(2),1,c1)
!         endif
!         cmin1 = min(cmin1,c1)
!         ngtests(1) = ngtests(1) + 1
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
!         call conf_test(tv11_12,nfmft_land(1),nfmft_land(3),nfmft_land(4),   &
!                        nfmft_land(2),1,c2)
!         cmin1 = min(cmin1,c2)
!         ngtests(1) = ngtests(1) + 1
! === NFMFT test disabled end ===
       
! ... perform tests.  
 
!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65 .gt. anth20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,anth20(1),anth20(3),anth20(4),
!     *               anth20(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,anth20(1),
!     +          anth20(2),anth20(3),anth20(4)
!      endif
! ................................................................
!     *****  END OF GROUP 1 TESTS  ***************************
 
 
 
!     ****  GROUP 2 TESTS  ***********************************

! ... 11 minus 3.9 micron BTDIF fog and low cloud test.
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) .and.  &
            nint(masir4) .ne.  nint(bad_data)) then

          if(masir11 .gt. 230.0) then

            nmtests = nmtests + 1
            call set_qa_bit(qa_bits,19)

            mas11_4 = masir4 - masir11

            call get_pn_thresholds(masir11,bt_11_bnds4,ant4_11l,ant4_11m1,      &
                                   ant4_11m2,ant4_11m3,ant4_11h,locut,hicut,    &
                                   midpt,power)

            if (mas11_4 .le. midpt) then
              call set_bit(testbits,19)
              nptests = nptests + 1
            end if
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
 
!     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2

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
!        write(h_output,'(1x,''tests '',3i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',6f8.5)') c2,c3,
!     +         cmin1,cmin2,fac,confdnc
!      endif
! ................................................................

      return
      end
