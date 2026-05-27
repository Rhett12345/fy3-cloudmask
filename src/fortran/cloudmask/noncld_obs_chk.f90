      subroutine noncld_obs_chk(indat,pxldat,confdnc,kele,          &
                                line_edge,klin,qa_bits,testbits,    &
                                smoke)

      implicit none
      save

!!F77 ************************************************************
!!Description:
! ... Routine which checks for the presence of a non-cloud obstruction.
! ... Currently checks for smoke.
!
!!Input Parameters:
! indat         Array containing nlcltx number of lines of
!               reflectance or BT values for each channel
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! confdnc       clear-sky confidence for current pixel
! kele          First element pixel in the region (context)
! maxele        Maximum number of pixels in a scan
! line_edge     Logical variable indicating first or last line of
!               data in a granule
! klin          Counter indicating number of lines processed
!
!!Output Parameters:
! qa_bits       Byte array containing qa results
! testbits      Byte array containing test results
! smoke         logical variable indicating whether smoke is present
!
!!Revision History:
! 10/04  Collection 5b    R. Frey
! Changed call to get_regstd.
!
!!Team-unique Header:
!
!    This software is developed by the MODIS Science Data Support Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-32373.
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!Design Notes:
!
!!END****************************************************************
!
      include 'global.inc'
      include 'noncld_obs_chk.inc'

! ... scalar arguments ..
      logical smoke,line_edge
      real confdnc
      integer kele,maxele,klin

! ... array arguments
      real pxldat(inband),indat(necntx,nlcntx,inband)
      byte qa_bits(10),testbits(6)

! ... local scalars ..
      logical thk_smoke,fire,bit_test
      real masv21,masv66,masir3,masir11,tdif,coef,masir12,tdiff,    &
           masv55,masv47,masv86,smkrat,vrat,sigma,masir3f,mean
      integer band,debug,h_output,j

! ... local arrays
      integer bitno(4),rtn(4),rtnqa(4)

! ... Common statement for debug purposes
!      common / bug / debug, h_output

!     Routine which checks for the possible presence of smoke.

!     Set bit numbers to test.
      data bitno /15,16,18,19/

!     Initializations.
      fire = .false.
      thk_smoke = .false.

!      masv47 = pxldat(3)
!      masv86 = pxldat(2)
!      masv66 = pxldat(1)
!      masv55 = pxldat(4)
!      masv21 = pxldat(7)
!      masir3 = pxldat(20)
!      masir3f = pxldat(21)
!      masir11 = pxldat(31)
!      masir12 = pxldat(32)

      masv47 = pxldat(1)
      masv86 = pxldat(4)
      masv66 = pxldat(3)
      masv55 = pxldat(2)
      masv21 = pxldat(7)
      masir3 = pxldat(20)
!      masir3f = pxldat(21)
      masir11 = pxldat(24)
      masir12 = pxldat(25)

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Within nco bit testing routine '',/)')
!      endif
! ................................................................

!     Check clear sky tests.
      bit_test = .true.
      do j = 1,4
        call check_qa_bits(qa_bits,bitno(j),rtnqa(j))
        call check_bits(testbits,bitno(j),rtn(j))
        if(rtnqa(j) .eq. 1) then
          if(rtn(j) .eq. 0) then
            bit_test = .false.
          end if
        else
          bit_test = .false.
        end if
      enddo

!     Do not perform fire or smoke tests if tests reported in bits
!     15, 16, 18, or 19 reported cloud.

      if(bit_test) then

!       Set bit to say that we did attempt nco smoke test
        call set_qa_bit(qa_bits,8)

!       Check for fires (hot spots). ???

!        if(masir3f .ne. nint(bad_data) .and. masir11 .ne. nint(bad_data)) then
!          tdif = masir3f - masir11
!          if (masir3f .gt. nc_bt37(1) .and. tdif .gt. nc37_11(1)) then
!            fire = .true.
!          end if 
!        end if 
        
        ! added by min min  [we should change thresholds]
!        if(masir3 .ne. nint(bad_data) .and. masir11 .ne. nint(bad_data)) then
!          tdif = masir3 - masir11
!          if (masir3 .gt. nc_bt37(1) .and. tdif .gt. nc37_11(1)) then
!            fire = .true.
!          end if 
!        end if 

!       Test for thick smoke. 

        if (nint(masv66) .ne. nint(bad_data) .and.   &
            nint(masv21) .ne. nint(bad_data) .and.   &
            nint(masv47) .ne. nint(bad_data) .and.   &
            nint(masv47) .ne. nint(bad_data) .and.   &
            nint(masir3) .ne. nint(bad_data) .and.   &
            nint(masir11) .ne. nint(bad_data) ) then

          if((masv21*100.0) .lt. nc21(1)) then

            coef = 6.0 + masv21*100.0
            if((masv66*100.0) .gt. coef) then

              smkrat = masv47 / masv66
              if(smkrat .ge. ncrat(1)) then


                vrat = masv86 / masv66
                if(vrat .ge. ncvrat(1)) then

                  band = 1
                  !call get_regstd(indat,kele,maxele,line_edge,klin,   &
                  !                band,sigma,mean)
                   call get_regstd(indat,kele,line_edge,klin,band,sigma,mean)
                  if(sigma .le. ncsig(1)) then

                    thk_smoke = .true.

                  end if
                end if
              end if
            end if
          end if
        end if

      end if

      if(thk_smoke .or. fire) smoke = .true.

! ...   debug statement ............................................
!        if (debug .gt. 0) then
!          write(h_output,'(10x,'' NCO test '',/,11f8.2,3l4/)')
!     *    masir3,masir11,tdif,masv21,coef,masv66,masv47,smkrat,masv86,
!     *    vrat,sigma,bit_test,fire,thk_smoke
!        endif
! ..................................................................

!     Perform dust test.
      if (nint(masir11) .ne. nint(bad_data) .and.  &
          nint(masir12) .ne. nint(bad_data)) then

        if(confdnc .gt. 0.67) then

!         Set bit to say that we attempted the nco dust test
          call set_qa_bit(qa_bits,28)

          tdiff = masir11 - masir12

          if(tdiff .lt. nc11_12(1)) then
            call clear_bit(testbits,28)
          end if

        end if

      end if

! ..................................................................

! ....debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x,'' Dust test '',/,4f8.2/)')
!     *    masir11,masir12,tdiff,nc11_12(1)
!      endif
! ................................................................

      return
      end
