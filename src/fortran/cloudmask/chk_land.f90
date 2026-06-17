      subroutine chk_land(pxldat,eco_type,desert,tbadj,confdnc,qa_bits,testbits)

      implicit none
      save

!-----------------------------------------------------------------------
!!F77
!
!!Description:
!     Perform final clear-sky confidence check on desert pixels.
!     If confidence of clear sky is low but temp is warm and if the
!     IR clear sky tests all passed, then use the 11 um brightness
!     temperature to assign a final confidence.
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! tbadj         11 um brightness temperature threshold adjustment for
!               deserts
! eco_type      Byte variable containing ecosystem index for current pixel
! desert        Logical flag indicating desert processing path
!
!!Output Parameters:
! testbits      Byte array containing cloud mask results
! qa_bits       Byte array containing qa bit results
! confdnc       Current pixel unobstructed confidence
!
!!Revision History:
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!END
!-----------------------------------------------------------------------

      include 'global.inc'
      include 'land_restoral.inc'

!     scalar arguments
      real confdnc,tbadj
      byte eco_type
      logical desert

!     array arguments
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

!     local scalars
      integer j,h_output,debug
      real m31, m22, m20, m5, m4, m5_4, m5_4_thr, md1, md2
      logical bit_test

!     local arrays
      integer bitno(5),rtn(5),rtnqa(5)
      real hds11(3)

!     external subroutines
      external check_bits,check_qa_bits,set_qa_bit

!     Common statement for debug purposes
!      common / bug / debug, h_output

!     Set bit numbers to test in "final test".
      data bitno /14,15,16,18,19/

!-----------------------------------------------------------------------

!      m5 = pxldat(5)
!      m4 = pxldat(4)
!      m31 = pxldat(31)
!      m22 = pxldat(22)
!      m20 = pxldat(20)

      m5 = pxldat(5)   ! 1.03 replace 1.24
      m4 = pxldat(2)
      m31 = pxldat(24)
      m22 = pxldat(20) ! 3.8 replace 3.959
!      m22 = pxldat(20)*0.7+pxldat(21)*0.3 
      m20 = pxldat(20)
      
!     Check IR clear sky tests.
      bit_test = .true.
      do j = 1,5
        call check_qa_bits(qa_bits,bitno(j),rtnqa(j))
        call check_bits(testbits,bitno(j),rtn(j))
        if(rtnqa(j) .eq. 1) then
          if(rtn(j) .eq. 0) then
            bit_test = .false.
          end if
        end if
      enddo

      if(bit_test) then

        if (nint(m31) .ne. nint(bad_data)) then

!         Get elevation-adjusted 11 micron brightness temperature threshold.

          call set_qa_bit(qa_bits,26)

          if(eco_type .eq. 8) then
            hds11(1) = ldsbt11bd(1) - tbadj
            hds11(2) = ldsbt11bd(2) - tbadj
            hds11(3) = ldsbt11bd(3) - tbadj
          else
            hds11(1) = ldsbt11(1) - tbadj
            hds11(2) = ldsbt11(2) - tbadj
            hds11(3) = ldsbt11(3) - tbadj
          end if

!         Check for hot scene.
          if(m31 .gt. hds11(1)) then

!           Assign confidence level based on 11 micron Tbb.
            if(m31 .gt. hds11(3)) then
!             Assign pixel to probably clear, set bit #26.
              confdnc = 1.0
              call set_bit(testbits,26)
            else if(m31 .gt. hds11(2)) then
!             Assign pixel to probably clear, set bit #26.
              confdnc = 0.96
              call set_bit(testbits,26)
            else
!             Assign pixel to uncertain, do not set bit #26.
              confdnc = 0.95
            end if

          end if

!-----------------------------------------------------------------------

!         debug statement
!          if (debug .gt. 0) then
!            write(h_output,'(''Land restoral test 1: '',/,6f10.2,l5,f10.5,
!     *            /)') ldsbt11(1),tbadj,hds11(1),hds11(2),hds11(3),
!     *            m31,bit_test,confdnc
!          endif

!-----------------------------------------------------------------------

        end if

        if (confdnc .le. 0.95) then

          if (nint(m20) .ne. nint(bad_data) .and.  &
              nint(m22) .ne. nint(bad_data) .and.  &
              nint(m31) .ne. nint(bad_data) .and.  &
              nint(m5) .ne. nint(bad_data)  .and.  &
              nint(m4) .ne. nint(bad_data)) then

            if(desert) then
              m5_4_thr = ldsr5_4_thr(1)
            else
              m5_4_thr = ldr5_4_thr(1)
            end if

            m5_4 = m5 / m4                   ! jincheng
            md1 = m20 - m22
            md2 = m22 - m31

!            if (md1 .lt. ld20m22(1) .and. md2 .lt. ld22m31(1) .and.   &            !jincheng
!                m5_4 .gt. m5_4_thr ) then
            if (  md2 .lt. ld22m31(1) .and.   &
                m5_4 .gt. m5_4_thr ) then
              confdnc = 0.96
              call set_bit(testbits,26)
            end if

!-----------------------------------------------------------------------

!          debug statement
!           if (debug .gt. 0) then
!             write(h_output,'(''Land restoral test 2: '',/,7f12.3,
!     *            /)') m5_4,md1,md2,ld20m22(1),ld22m31(1),m5_4_thr,confdnc
!           endif

!-----------------------------------------------------------------------


          end if

        end if

      end if

!-----------------------------------------------------------------------

      return
      end
