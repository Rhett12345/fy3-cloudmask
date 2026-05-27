      subroutine LandNite(pxldat,plat,vza,coast,desert,hi_elev,sh_lake,   &
                          sfctmp,eco_type,testbits,qa_bits,nmtests,       &
                          confdnc,ptwp,btclr,is_cold_sfc)

      use names_module, only: fylat_sensor_id, fylat_nwp_opt
      implicit none
      save

!---------------------------------------------------------------------
!!F77 
!
!!Description:
!      Routine for performing clear sky tests over land
!      surfaces during nightime hours.
!
!      For nighttime land the groups are:
!          Group 1: High thick cloud
!                   13.9 micron bt test (masir13) 
!                   6.75 micron bt test 
!                   Surface temperature test 
!
!          Group 2: Low cloud - thick
!                   11-12 micron bt tests
!                   11-4 micron bt tests
!                   7.3-11 micron bt test (thick mid-level)
!
!          Group 5: High cloud - thin
!                   3.7-12 micron bt test
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           viewing zenith angle
! coast         flag indicating coast processing
! desert        flag indicating desert processing
! hi_elev       flag indicating high elevation processing (> 2000 m)
! sh_lake       Logical flag indicating shallow inland lakes
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
!
! 07-15-02 raf:
! Modified 3.9-11 um test so that the test threshold is a function of
! the 11-12 um BTD.  Thresholds are calculated in get_nl_thresholds.f
! Added the 7.3-11 um mid-level cloud test.
! Added 2K to 3.9-11 um test threshold for coastal areas.
! 06/04 Collection 5   raf:
! Added surface temperature test    
! Implemented new 11-12 um thin cirrus test (J. Key version)
! 10/04 Collection 5   raf:
! Added 11-12 and 3.9-11 um BTD conditions to choice of LST thresholds.
! Changed basic LST threshold from 10K to 12K.
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!
! Externals:
!       Subroutines: conf_test,set_bit,set_qa_bit,tview,get_nl_thresholds
!
!!END
!-----------------------------------------------------------------------

      include 'global.inc'
      include 'LandNite_thr.inc'
      include 'pfmft_nfmft_thr.inc'
      
! ...
! ... scalar arguments ..
      real confdnc,vza,plat,sfctmp,ptwp
      integer nmtests
      logical coast,desert,hi_elev,sh_lake
      byte eco_type
! ...
! ... array arguments ..
      real pxldat(inband), btclr(7),tv11_12
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,mas4_12,masir11,masir12,masir13,masir4,           &
           c3,masir65,mas11_4,c4,cmin1,cmin2,cmin5,groups,         &
           fac,pre_confdnc,c5,masdf1,schi,cosvza,dtr,pi,diftemp,   &
           dfthrsh,locut,hicut,masir73,mas7_11,c6,                 &
           power,midpt,a,c9,sfcdif,corr,lst_thrsh,masdf2
      integer nptests,debug,h_output,kk
      ! added by minmin
      real masir37, mas37_12, delta_t

! ... local arrays
      integer ngtests(3)
! ...
      real, parameter :: Rel_equality_EPS = 0.000001
      real, parameter :: max_vza = 65.49
      integer i4
! ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview,get_nl_thresholds

! ... Common statement for debug purposes
!      common / bug / debug, h_output

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
!      masir73 = pxldat(28)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)
!      masir13 = pxldat(35)
      masir37 = pxldat(20)
!      masir4 = 0.30*pxldat(21) + 0.70*pxldat(20)  ! 4.05 replace 3.959
      masir4 = pxldat(20)           !3.8 replace 3.959   jincheng
!      masir65 = pxldat(27)
      masir73 = pxldat(22)
      masir11 = pxldat(24)
      masir12 = pxldat(25)
!      masir13 = pxldat(35)

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
      c9 = 0.0
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
  !      cmin1 = min(cmin1,c1)
  !      ngtests(1) = ngtests(1) + 1
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
    !    cmin1 = min(cmin1,c2)
    !    ngtests(1) = ngtests(1) + 1
      endif
      
!     **** GROUP 1 TESTS *************************************
! ... co2 high cloud test
!      if (nint(masir13) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,14)
!        if (masir13.gt.nlco2(2)) then
!          call set_bit(testbits,14)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir13,nlco2(1),nlco2(3),nlco2(4),
!     *                nlco2(2),1,c1)
!        cmin1 = min(cmin1,c1)
!        ngtests(1) = ngtests(1) + 1
!      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,nlco2(1),
!     +          nlco2(2),nlco2(3),nlco2(4)
!      endif
! ................................................................

!     H20 vapor channel (6.7 micron) high cloud test
!      if (nint(masir65) .ne. nint(bad_data)) then
!        nmtests = nmtests + 1
!        call set_qa_bit(qa_bits,15)
!        if (masir65.gt.nlh20(2)) then
!          call set_bit(testbits,15)
!          nptests = nptests + 1
!        end if
!        call conf_test(masir65,nlh20(1),nlh20(3),nlh20(4),
!     *                nlh20(2),1,c2)
!        cmin1 = min(cmin1,c2)
!        ngtests(1) = ngtests(1) + 1
!      endif
     
    
!     channel (12.0 -3.7 micron) high cloud test from NPP/VIIRS ATBD   [added by minmin]
     i4=0
     if (i4==1) then 
      if (nint(masir37) .ne. nint(bad_data) .and. &
          masir37       .gt. 230.0          .and. &   
          nint(masir12) .ne. nint(bad_data) .and. &
          .not. sh_lake .and. ptwp .lt. 6.0) then
          
        mas37_12 = masir37 - masir12  
        nmtests = nmtests + 1       
        ! 4.5 4.0 3.5  cloudy un clear
        call set_qa_bit(qa_bits,15)
!        if (masir65.gt.nlh20(2)) then
        if (mas37_12.gt.4.0) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
!        call conf_test(masir65,nlh20(1),nlh20(3),nlh20(4),nlh20(2),1,c2)
        call conf_test(mas37_12,4.5,3.5,1.0,4.0,1,c2)

        ! adjust confidence for modis method
        !  clear Pro Clear  Pro Cloudy   Cloudy
        !   90%   90%-50%     50%-0%       0%
        !   99%   99%-95%     95%-66%      66%
        !if (c2 > 0.90 .and. c2 < 1.0) then 
        !    c2 = 0.99
        !endif
        !if (c2 > 0.50 .and. c2 <= 0.90 ) then 
        !    c2 = 0.95 + 0.04*((c2-0.50)/0.40)
        !endif       
        !if (c2 > 0.0 .and. c2 <= 0.50 ) then 
        !    c2 = 0.66 + 0.31*(c2/0.50)
        !endif   
        !if (c2 == 0.0 .or. c2 == 1.0) then 
        !    c2 = c2
        !endif   
        cmin1 = min(cmin1,c3)

        ngtests(1) = ngtests(1) + 1
      endif
    endif
! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,nlh20(1),
!     +          nlh20(2),nlh20(3),nlh20(4)
!      endif

! ................................................................

! ... Surface Temperature Test
      i4=0
      if ( nint(masir11) .ne. nint(bad_data) .and. (.not. hi_elev) .and.   &
           nint(masir12) .ne. nint(bad_data) .and. eco_type .ne. 8) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        masdf1 = masir11 - masir12
       ! masdf2 = masir11 - masir4
        masdf2 = masir11 - (masir4-1.5)  ! revised by minmin 20190108 to correct masir4
!print*,desert,masir11,masir12,masir4,masdf1,masdf2
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

       ! if(desert) then !! revised by minmin 20190108
       !   lst_thrsh = 20.0
       ! else if(masdf1 .ge. -0.2 .or. (masdf1 .lt. -0.2 .and. (masdf2 .le. -0.5 .or. masdf2 .ge. 1.0))) then ! revised by minmin 20190108
       !   lst_thrsh = 12.0  
       ! else
       !   lst_thrsh = 20.0
       ! end if
       
       
       ! revised by minmin 20190108  !due to grapes higher lst
        delta_t = 0.0
        if (fylat_nwp_opt == 6) then ! grapes 0.25*0.25
            delta_t = 0.0
        endif
        
        if(desert) then ! revised by minmin 20190108  !due to grapes higher lst
          lst_thrsh = 20.0 + delta_t
        else if(masdf1 .ge. -0.2 .or. (masdf1 .lt. -0.2 .and. (masdf2 .le. -0.5 .or. masdf2 .ge. 1.0))) then ! revised by minmin 20190108
        !else if(masdf1 .ge. -0.2 ) then
          lst_thrsh = 12.0 + delta_t 
        else
          lst_thrsh = 20.0 + delta_t
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

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c9)
        cmin1 = min(cmin1,c9)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''sfctmp: '',9f9.3)') masdf1,masdf2,
!     +   sfctmp,masir11,sfcdif,locut,midpt,hicut,c9
!      end if

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
          dfthrsh = nl11_12hi(1)
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

        locut = dfthrsh
        midpt = dfthrsh - (0.3 * dfthrsh)
        if(masir11 .lt. 270.0) then
          if(abs(plat) .le. 30.0) then
            hicut = midpt - 1.25
          else
            a = (90.0 - abs(plat)) / 60.0
            hicut = -0.1 - ((a**4) * 1.15)
          end if
        else
          hicut = midpt - 1.25
        end if

        call conf_test(masdf1,locut,hicut,1.0,midpt,1,c5)
        cmin2 = min(cmin2,c5)
        ngtests(2) = ngtests(2) + 1

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''APOLLO masdf1: '',5f10.2)') masdf1,
!     +          nl11_12hi(1),dfthrsh,locut,hicut
!      endif
! ................................................................

! ... 11 minus 4 micron BTDIF fog and low cloud test.

   !if (i4==1) then
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir4) .ne.  nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
        mas11_4 = masir11 - masir4
      !  mas11_4 = masir11 - (masir4-1.5)   ! corrected by minmin  20190109

        masdf1 = masir11 - masir12
        call get_nl_thresholds(masdf1,locut,hicut,midpt,power)

        if(sh_lake) then
          locut = locut + 2.0 
          midpt = midpt + 2.0 
          hicut = hicut + 2.0 
        end if

        if (mas11_4 .le. midpt) then
          call set_bit(testbits,19)
          nptests = nptests + 1
        end if

        call conf_test(mas11_4,locut,hicut,power,midpt,1,c3)

        cmin2 = min(cmin2,c3)
        ngtests(2) = ngtests(2) + 1

      endif
   !endif
! ... debug statement ............................................
!      if (debug .gt. 0) then
!         write(h_output,'(1x,''mas11_4: '',l4,6f10.2)')coast,mas11_4,
!     +       masdf1,locut,midpt,hicut,c3
!      endif
! ................................................................


! ... 7.3-11um brightness temperature difference test
! ... for thick, mid-level clouds
    !if (i4==1) then
      if (nint(masir11) .ne. nint(bad_data)  .and.  &
          nint(masir73) .ne.  nint(bad_data) .and.  &
          nint(masir4) .ne. nint(bad_data)) then

        mas11_4 = masir11 - masir4

        if(mas11_4 .le. -2.0) then
       !  if(mas11_4 .le. -1.5) then  ! corrected by minmin 
         
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,23)
          mas7_11 = masir73 - masir11

          if ( mas7_11 .le. nl7_11s(2) ) then
            nptests = nptests + 1
            call set_bit(testbits,23)
          end if

          call conf_test(mas7_11,nl7_11s(1),nl7_11s(3),nl7_11s(4),   &
                         nl7_11s(2),1,c6)

         ! cmin2 = min(cmin2,c6)
         ! ngtests(2) = ngtests(2) + 1

        end if

      endif
   !endif
! ... debug statement ............................................
!      if (debug .gt. 0) then
!         write(h_output,'(1x,''mas7_11: '',6f10.2)')mas7_11,masir11,
!     +            nl7_11s(1),nl7_11s(2),nl7_11s(3),c6
!      endif
! ................................................................


! *******     END OF GROUP 2 TESTS ****************************



! *******    START OF GROUP 5 TESTS  **************************
! ... 4-12um brightness temperature difference test
! ... for thin cirrus)
      if (nint(masir12) .ne. nint(bad_data) .and.  &
          nint(masir4) .ne.  nint(bad_data)) then
        mas4_12 = masir4 - masir12
        !mas4_12 = (masir4 - 1.5) - masir12  ! revised by minmin 20190109
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,17)
        if (mas4_12.le.nl4_12hi(2)) then
          nptests = nptests + 1
          call set_bit(testbits,17)
        end if
        call conf_test(mas4_12,nl4_12hi(1),nl4_12hi(3),nl4_12hi(4),   &
                       nl4_12hi(2),1,c4)
        cmin5 = min(cmin5,c4)
        ngtests(3) = ngtests(3) + 1
      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!         write(h_output,'(1x,''mas4_12: '',5f10.2)')mas4_12,nl4_12hi(1),
!     +            nl4_12hi(2),nl4_12hi(3),nl4_12hi(4)
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
      if (groups .gt. 0) fac = 1.0 / groups
!     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac


! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(1x,''tests '',5i10)') nmtests,nptests,ngtests
!        write(h_output,'(1x,''confdnc '',10f8.5/,f8.5)') c1,c2,c3,c4,c5,
!     +      c6,cmin1,cmin2,cmin5,fac,confdnc
!      endif
! ................................................................

      return
      end
