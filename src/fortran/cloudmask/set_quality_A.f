      subroutine set_quality_A(nmtests,nbands,lsf,qa_bits)

      implicit none
      save

c--------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for setting output qa bit flags according
c     to final number of spectral tests and bands used
c     to determine the cloud mask
c
c!Input Parameters:
c nmtests       Number of tests applied to this pixel
c nbands	Number of bands successfully read for this pixel
c lsf           Current pixel land/sea flag
c h_eco1        Ecosystem file handle number
c
c!Output Parameters:
c qa_bits       Byte array containing qa bit results
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c--------------------------------------------------------------------

c ... scalar arguments
      integer nmtests,nbands,h_eco1,lsf

c ... array arguments
      byte qa_bits(10)

c ... local scalars
      integer debug,h_output

c ... external subroutines
      external set_qa_bit

c ... Common statement for debug purposes
!      common / bug / debug, h_output

c ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Within set_quality_A routine '',/)')
!      endif
c ...............................................................


c ... set number of spectral tests applied bits
      !if (nmtests .gt. 6) then
      if (nmtests .gt. 4) then
         call set_qa_bit(qa_bits,50)
         call set_qa_bit(qa_bits,51)
      !else if (nmtests .gt. 3) then
      else if (nmtests .gt. 2) then
         call set_qa_bit(qa_bits,51)
      else if (nmtests .gt. 0) then
         call set_qa_bit(qa_bits,50)
      end if

c ... set number of bands with good data read
      if (nbands .gt. 14) then
         call set_qa_bit(qa_bits,48)
         call set_qa_bit(qa_bits,49)
      else if (nbands .gt. 7) then
         call set_qa_bit(qa_bits,49)
      else if (nbands .gt. 0) then
         call set_qa_bit(qa_bits,48)
      end if

c ... Set ecosystem file bit
!      if (h_eco1 .eq. -5555) call set_qa_bit(qa_bits,64)
      call set_qa_bit(qa_bits,64)

c ... Set the Land/Sea Mask file
      if (lsf .eq. -1) then
        call set_qa_bit(qa_bits,70)
        call set_qa_bit(qa_bits,71)
      endif

c ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x,'' qa bits:    nmtest   nbands'',2i10,/)')
!     *        nmtests,nbands
!      endif
c ................................................................

      return
      end
