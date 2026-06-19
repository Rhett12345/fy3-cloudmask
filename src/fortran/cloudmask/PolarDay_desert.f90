      subroutine PolarDay_desert(pxldat,vza,visusd,cirrus_vis,hi_elev,  &
                                 testbits,qa_bits,nmtests,confdnc,      &
                                 btclr,is_cold_sfc)

      implicit none
      save

!--------------------------------------------------------------------
!!F77
!
!!Description:
!      Routine for performing clear sky tests over polar desert
!      surfaces during daylight hours.
!
!      For daytime polar desert the groups are:
!          Group 1: 11-12um BTD thin cirrus test
!          Group 2: 11-4um BTD fog/low cloud test
!          Group 3: Visible reflectance tests (0.86um + 1.38um)
!
!!Input Parameters:
! real    pxldat   Array containing reflectance or brightness temperatures
!                  for all bands for a single pixel
! real    vza      Current pixel viewing angle
! logical visusd   Logical variable indicating whether vis data used
! logical hi_elev  Flag indicating elevation > 2000 meters
!
!!Output Parameters:
! byte    testbits six byte integer containing bit results
! byte    qa_bits  ten byte array containing QA bit results
! integer nmtests  Counts number of inidividual tests applied
! real    confdnc  product of all applied individual confidences
! logical cirrus_vis Logical variable flagging thin cirrus contaminated
!                scenes in the visible
!
!----------------------------------------------------------------------

      include 'global.inc'
      include 'PolarDay_desert_thr.inc'

! ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical visusd,cirrus_vis,hi_elev
! ...
! ... array arguments ..
      real pxldat(inband),btclr(7)
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)
! ...
! ... local scalars ..
      real c1,c2,c3,mas11_4,masdf1,   &
           masir11,masir12,masir4,masv188,masv88,   &
           cmin1,cmin2,cmin3,pre_confdnc,          &
           groups,fac
      integer nptests,kk
! ...
! ... local arrays ..
      integer ngtests(3)
! ...
! ... external subroutines ..
      external conf_test,set_bit,clear_bit,set_qa_bit

! ... initialize variables
      nmtests = 0
      nptests = 0
      confdnc = 1.0

! ... place band values into individual variables
      masv88 = pxldat(4)
      masv188 = pxldat(19)
      masir4 = pxldat(20) ! 3.8um
      masir11 = pxldat(24)
      masir12 = pxldat(25)

! ... initialize group confidences
      cmin1 = 1.0
      cmin2 = 1.0
      cmin3 = 1.0

! ... initialize group test counts
      do kk = 1, 3
         ngtests(kk) = 0
      enddo

!----------------------------------------------------------------------
!     GROUP 1: 11-12um BTD thin cirrus test
!----------------------------------------------------------------------
      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir12) .ne. nint(bad_data)) then
        masdf1 = masir11 - masir12

        call set_qa_bit(qa_bits,18)
        if (masdf1.le.pds11_12hi(2)) then
          nmtests = nmtests + 1
          call set_bit(testbits,18)
          nptests = nptests + 1
          ngtests(1) = ngtests(1) + 1
        end if
        call conf_test(masdf1,pds11_12hi(1),pds11_12hi(3),   &
                       pds11_12hi(4),pds11_12hi(2),1,c1)
        cmin1 = min(cmin1,c1)
      endif

!----------------------------------------------------------------------
!     GROUP 2: 11-4um BTD fog/low cloud test
!----------------------------------------------------------------------
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) .and.   &
            nint(masir4) .ne.  nint(bad_data)) then
          mas11_4 = masir11 - masir4
          call set_qa_bit(qa_bits,19)
          if (mas11_4.ge.pds11_4lo(2)) then
            nmtests = nmtests + 1
            call set_bit(testbits,19)
            nptests = nptests + 1
            ngtests(2) = ngtests(2) + 1
          end if

          call conf_test(mas11_4,pds11_4lo(1),pds11_4lo(3),   &
                         pds11_4lo(4),pds11_4lo(2),1,c2)
          call conf_test(mas11_4,pds11_4hi(1),pds11_4hi(3),   &
                         pds11_4hi(4),pds11_4hi(2),1,c3)
          cmin2 = min(cmin2,c2,c3)
        end if
      endif

!----------------------------------------------------------------------
!     GROUP 3: Visible reflectance tests (0.86um + 1.38um)
!----------------------------------------------------------------------
      if (visusd) then
!       0.86um test
        if (nint(masv88) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,20)
          if (masv88 .le. pdsref2(2)) then
            nmtests = nmtests + 1
            call set_bit(testbits,20)
            nptests = nptests + 1
            ngtests(3) = ngtests(3) + 1
          end if
          call conf_test(masv88,pdsref2(1),pdsref2(3),pdsref2(4),   &
                         pdsref2(2),1,c3)
          cmin3 = min(cmin3,c3)
        end if

!       1.38um test
        if (.not. hi_elev) then
          if (nint(masv188) .ne. nint(bad_data)) then
            call set_qa_bit(qa_bits,16)
            if (masv188 .le. pdsref3(2)) then
              nmtests = nmtests + 1
              call set_bit(testbits,16)
              nptests = nptests + 1
              ngtests(3) = ngtests(3) + 1
            end if
            call conf_test(masv188,pdsref3(1),pdsref3(3),pdsref3(4),   &
                           pdsref3(2),1,c3)
            cmin3 = min(cmin3,c3)
          end if
        end if
      endif

!     Thin cirrus bit check
      if (visusd .and. (.not. hi_elev) ) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (masv188 .lt. pdstci(1) .and. masv188 .ge. pdstci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
        endif
      endif

!----------------------------------------------------------------------
!     Determine confidence based on group values
!----------------------------------------------------------------------
      groups = 0
      pre_confdnc = 1.0
      do kk = 1,3
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
          if (kk .eq. 1) pre_confdnc = pre_confdnc * cmin1
          if (kk .eq. 2) pre_confdnc = pre_confdnc * cmin2
          if (kk .eq. 3) pre_confdnc = pre_confdnc * cmin3
        end if
      enddo
      if (groups .gt. 0) then
        fac = 1.0 / groups
        confdnc = pre_confdnc**fac
      else
        confdnc = 1.0
      end if

      return
      end
