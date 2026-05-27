      subroutine chk_spatial_var(indat,kele,confdnc,qa_bits,
     *                           testbits)


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
c confdnc       Current pixel unobstructed confidence
c
c!Output Parameters:
c qa_bits       Byte array containing qa bits
c testbits      Byte array containing test results
c
c!Revision History:
c $Id: chk_spatial_var.f,v 1.1.2.2 2004/05/26 15:17:27 raf Exp $
c 06/04 Collection 5  R. Frey:
c Updated calling arguments in call to routine spatial_var
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'

c     Scalar arguments
      integer kele
      real confdnc

c     Array arguments
c      real indat(npixel,nlcntx,inband)
      real indat(necntx,nlcntx,inband)
      byte testbits(6),qa_bits(10)

c     Local scalars
      integer varslt,h_output,debug,ipt

c     Local Arrays
      real diff(var_band,8),rgdata(nlcntx,necntx,var_band)

c     External subroutines
      external get_regdif,spatial_var,set_qa_bit,set_bit

c     Common statement for debug purposes
!      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x,'' Checking spatial variability'',/)')
!      endif

c---------------------------------------------------------------------

c     Get brightness temperature differences between pixel of interest
c     and the 8 surrounding it.
      call get_regdif(indat,kele,rgdata,diff)

c     Check variation in the region.
      call spatial_var(diff,ipt,varslt)
        
c     Set bit which indicates this test was applied.
      call set_qa_bit(qa_bits,25)

c     Set bit if test passed (found no evidence of cloud).
      if (varslt .eq. 1) then
        call set_bit(testbits,25)
      end if

c     Increase the confidence if spatial variability test showed
c     uniform conditions.
      if ((varslt .eq. 1) .and. (confdnc .gt. .66)) then
        confdnc = 0.96

c ......debug statement ................................................
!        if (debug .gt. 0) then
!          write(h_output,'(10x,'' Confidence Bumped up through spatial 
!     +          uniformity testing  '',f10.2,/)') confdnc
!        endif
c ......................................................................

      else if ((varslt .eq. 1) .and. (confdnc .le. 0.66)) then
        confdnc = 0.67

c ......debug statement ................................................
!        if (debug .gt. 0) then
!          write(h_output,'(10x,'' Confidence Bumped up through spatial
!     +        uniformity testing  '',f10.2,/)') confdnc
!        endif
c ......................................................................
      else if (varslt .eq. 0) then

c ......debug statement ................................................
!          if (debug .gt. 0) then
!            write(h_output,'(10x,'' Confidence Not Bumped up through   
!     +           spatial uniformity testing  '',f10.2,/)') confdnc
!          endif
c ......................................................................

      endif

      return
      end
