      subroutine spatial_var(diff,ipt,result)

      implicit none
      save

c---------------------------------------------------------------------
C!F77
c
c!DESCRIPTION:
c     routine for performing spatial variability tests and
c     setting the appropriate bit flag.
c
c!Input Parameters:
c diff          Array of surrounding pixel reflectance or brightness
c               temperature differences from center pixel value
c
c!Output Parameters:
c ipt           number pixels satisfying variability test (8 max)
c results	results of spatial variability test (1=uniform)
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Updated argument list.
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.  
c 
c!END
c---------------------------------------------------------------------

      include 'global.inc'

      integer masir11
      parameter (masir11=2)
c ...
c ... scalar arguments ..
      integer ipt,result
c ...
c ... array arguments ..
      real diff(var_band,8)
c ...
c ... local scalars ..
      integer i
      INCLUDE 'spatial_var_thr.inc'
c ...
c ... Have chose to only test the 11 micron uniformity for now
      ipt = 0
c ... compare surrounding bt differences to threshold
      do 200 i = 1,8
        if (abs(diff(masir11,i)) - bad_data .gt. 0.1)  then
          if (diff(masir11,i).le.dovar11(1)) then
            ipt = ipt + 1
          end if
        endif
  200 continue

c ... if all surrounding pixel differences were less than the
c ... threshold value, scene is declared to be uniform.
      if (ipt.eq.8) then
        result = 1
      else
        result = 0
      end if

      return
      end
