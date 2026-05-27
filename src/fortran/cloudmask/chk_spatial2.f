      subroutine chk_spatial2(indat,kele,npix)


c---------------------------------------------------------------------
C!F77 
c
c!Description:
c
c     Routine for checking spatial variability of ocean scenes over
c     a 3X3 pixel region.
c
c!Input parameters:
c indat         Array containing nlcntx lines of data
c kele          Current granule element number being processed
c
c!Output Parameters:
c npix          Number pixels (max of 8) which satisfy variability test
c
c!Revision History:
c $Id: chk_spatial2.f,v 1.1.2.1 2004/05/27 15:55:24 raf Exp $
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'

c     Scalar arguments
      integer kele,npix

c     Array arguments
      real indat(necntx,nlcntx,inband)

c     Local scalars
      integer varslt,h_output,debug

c     Local Arrays
      real diff(var_band,8),rgdata(nlcntx,necntx,var_band)

c     External subroutines
      external get_regdif,spatial_var

c     Common statement for debug purposes
!      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Get brightness temperature differences between pixel of interest
c     and the 8 surrounding it.
      call get_regdif(indat,kele,rgdata,diff)

c     Check variation in the region.
      call spatial_var(diff,npix,varslt)
        
      return
      end
