      SUBROUTINE SET_QA_BIT(QA_BITS,BIT_NUM)

c--------------------------------------------------------------------
c!F77
c
c!Description:
c     Routine for setting a single bit within the 80-bit
c     cloud mask qa output array.
c
c!Input parameters:
c bit_num       Bit number in 48 bit cloud mask
c
c!Output Parameters:
c qa_bits       10 byte array containing QA bit results
c
c!Revision History:
c
c!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c--------------------------------------------------------------------

      implicit none

      save

c     scalar arguments
      integer bit_num
c     argument arrays
      byte qa_bits(10)
c
c     local scalars
      integer iword,itest,ipos

c ... intrinsic functions ..
      intrinsic ibset

c     Determine which word (1-10) of the 1-byte array contains the
c     bit of interest.

      iword = (bit_num / 8) + 1

c     Determine the position of the bit within the current
c     8-bit segment (1-byte word).

      ipos = bit_num - ((iword-1) * 8)

      itest = qa_bits(iword)
      itest = ibset(itest,ipos)

      qa_bits(iword) = itest

      END

