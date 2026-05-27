      subroutine check_qa_bits(qa_bits,bit_num,rtn)

      implicit none
      save

c----------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine for checking if a bit has been set or not
c     out of a 10 byte array variable
c
c!Input parameters:
c qa_bits       10-byte integer containing qa bit results
c bit_num       Bit number in 48 bit cloud mask
c
c!Output Parameters:
c rtn           Result: 1 = bit set, 0 = bit not set
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------
c
c     scalar arguments
      integer bit_num,rtn
c     scalar arrays
      byte qa_bits(10)

c     local scalars 
      integer iword,itest,ipos

c ... intrinsic functions ..
      intrinsic btest

c     initialize result
      rtn = 0

c     Determine which word (1-10) of the 1-byte array contains the 
c     bit of interest.

      iword = (bit_num / 8) + 1

c     Determine the position of the bit within the current
c     8-bit segment (1-byte word).

      ipos = bit_num - ((iword-1) * 8)

      itest = qa_bits(iword)
  
      if (btest(itest,ipos)) then
         rtn = 1
      else
         rtn = 0
      endif

      return
      end
