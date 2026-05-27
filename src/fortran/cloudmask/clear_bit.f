      subroutine clear_bit(testbits,bit_num)

      implicit none
      save

c-------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine for clearing a single bit within the 48-bit
c     cloud mask output array.
c
c!Input Parameters:
c bit_num       Bit number in 48 bit cloud mask
c
c!Output Parameters:
c testbits      four-byte integer containing bit results
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c-------------------------------------------------------------------

c     scalar arguments
      integer bit_num
c     argument arrays
      byte testbits(6)

c     local scalars 
      integer itest,iword,ipos

c ... intrinsic functions ..
      intrinsic ibclr

c     Determine which word (1-6) of the 1-byte array contains the 
c     bit of interest.

      iword = (bit_num / 8) + 1

c     Determine the position of the bit within the current
c     8-bit segment (1-byte word).

      ipos = bit_num - ((iword-1) * 8)

      itest = testbits(iword)
      itest = ibclr(itest,ipos)

      testbits(iword) = itest

      return
      end
