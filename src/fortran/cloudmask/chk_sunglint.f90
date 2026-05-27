      subroutine chk_sunglint(indat,pxldat,kele,confdnc,klin,  &
                              line_edge,qa_bits,testbits)


!---------------------------------------------------------------------
!!F77 
!
!!Description:
!
!     Performs clear sky restoral tests in sun-glint conditions.
!
!!Input parameters:
! indat         Array containing nlcntx lines of data
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for current pixel
! kele          Current granule element number being processed
! confdnc       Current pixel unobstructed confidence
!
!!Output Parameters:
! qa_bits       Byte array containing qa bits
! testbits      Byte array containing test results
!
!!Revision History:
! $Id: chk_sunglint.f,v 1.1.2.6 2004/10/25 17:13:35 raf Exp $
! 06/04 Collection 5  R. Frey
! Updated argument list in call to spatial_var
! 10/04 Collection 5  R. Frey
! Added band 2 "sigma*mean" clear-sky restoral test
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD35.
!
!!END
!---------------------------------------------------------------------

      save
   
      include 'global.inc'
      include 'snglntr_thr.inc'

!     Scalar arguments
      integer kele,klin,maxele
      real confdnc
      logical line_edge

!     Array arguments
      real indat(necntx,nlcntx,inband),pxldat(inband)
      byte testbits(6),qa_bits(10)

!     Local scalars
      integer h_output,debug,varslt,rtn,ipt
      real modir37,modir11,d37_11,modv895,modv935,modv443,rat,sigma,mean
      logical irclr

      integer band

!     Local Arrays
      real diff(var_band,8),rgdata(nlcntx,necntx,var_band)

!     External subroutines
      external spatial_var,get_regdif,check_bits,set_qa_bit,set_bit,get_regstd
      


! ... Common statement for debug purposes
!      common / bug / debug, h_output

!---------------------------------------------------------------------

!     Get brightness temperature differences between pixel of interest
!     and the 8 surrounding it.
      call get_regdif(indat,kele,rgdata,diff)

!     Check variation in the region.
      call spatial_var(diff,ipt,varslt)
        
      if(varslt .eq. 1) then

!       Region is uniform in the 11 micron IR window.
!       Determine the logical flag 'irclr' - true if ir cloud tests below
!       have all been passed. APOLLO test makes final decision for bit 18.
        irclr = .true.
        call check_bits(testbits,13,rtn)
        if(rtn .eq. 0) irclr = .false.
!        call check_bits(testbits,14,rtn)  ! revised by wuxiao
!        if(rtn .eq. 0) irclr = .false.
!        call check_bits(testbits,15,rtn)
!        if(rtn .eq. 0) irclr = .false.
!        call check_bits(testbits,18,rtn)
!        if(rtn .eq. 0) irclr = .false.
        call check_bits(testbits,27,rtn)
        if(rtn .eq. 0) irclr = .false.

        if(irclr) then

!          modir37 = pxldat(20)
!          modir11 = pxldat(31)
!          modv895 = pxldat(17)
!          modv935 = pxldat(18)
!          modv443 = pxldat(9)

          modir37 = pxldat(20)
          modir11 = pxldat(24)
          modv895 = pxldat(16)
          modv935 = pxldat(17)
          modv443 = pxldat(9)

          if(modir37 .ne. nint(bad_data) .and.   &
             modir11 .ne. nint(bad_data) .and.   &
             modv895 .ne. nint(bad_data) .and.   &
             modv935 .ne. nint(bad_data)) then

!           Set bit which indicates this test was applied.
            call set_qa_bit(qa_bits,26)

            d37_11 = modir37 - modir11

!           Set bit if tests passed.
            if(d37_11 .ge. sg_tbdfl(1) ) then

              confdnc = 0.67

              rat = modv895 / modv935

              if(rat.gt.snglrat(1) .and. modv443.ne.nint(bad_data)) then

                call set_bit(testbits,26)
                confdnc = 0.96

              else

                band = 2
              call get_regstd(indat,kele,line_edge,klin,band,sigma,mean)

                if(mean .ne. bad_data) then
                  if( (sigma * mean) .lt. 0.001) then
                      call set_bit(testbits,26)
                      confdnc = 0.96
                  end if
                end if

              end if

            end if

          end if

! ......  debug statement ............................................
!          if(debug .gt. 0) then
!            write(h_output,'(10x,'' Sun-glint clear sky restoral: '')')
!            write(h_output,'(10x,'' 3.7-11 um difference '',i5,f8.2)') varslt,d37_11
!            write(h_output,'(10x,'' Channel 17-18 ratio  '',f8.2)') rat
!            write(h_output,'(10x,'' Sigma test  '',3f10.5)') mean,sigma,sigma*mean
!            write(h_output,'(10x,'' Confidence after sun-glint CSR tests
!     +         is  '',f10.2,/)') confdnc
!          end if
! ....................................................................

        end if

      end if

      return
      end
