      subroutine get_regstd(indat,kele,line_edge,klin,band,
     *                      sigma,mean)

      implicit none
      save

c---------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c     Computes standard deviation of (center) pixel of interest and the     
c     surrounding eight pixel values for indicated channel.   
c
c
c!Input parameters:
c indat         Array containing nlcltx number of lines of
c               reflectance or BT values for each channel
c kele          First element pixel in the region (context)
c maxele        Maximum number of pixels in a scan
c line_edge     Logical variable indicating first or last line of
c               data in a granule
c klin          Counter indicating number of lines processed
c band          MODIS band number (1-36)
c
c!Output Parameters:
c sigma         Standard deviation of indicated MODIS channel
c mean          Mean of pixel reflectances or brightness temperatures
c
c!Revision History:
c 10/04 Collection 5b   R. Frey
c Added calculation and output of mean value.
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD35.
c
c!END
c---------------------------------------------------------------------------

      INCLUDE 'global.inc'

c     scalar arguments
      integer band,kele,maxele,klin
      real sigma,mean
      logical line_edge

c     array arguments
      real indat(necntx,nlcntx,inband)

c     local scalars
      integer imv,ide,i,j,k,nl,debug,h_output
      double precision n,num,den,sqsum,sumsq,sum,tbb,sig,mn

c     local arrays
      real rgdata(nlcntx,necntx)

c ... Common statement for debug purposes
!      common / bug / debug, h_output

c ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Within routine get_regstd.f'',/)')
!      endif
c ................................................................

c     Initialize variables.

      imv = 0
      ide = 0
      do 10 i = 1 , nlcntx
        do 20 j = 1 , necntx
          rgdata(i,j) = 0.0
   20   continue
   10 continue

c ... put pixel data into context data array.

      imv = ((necntx-1)/2) + 1

      if(line_edge .and. klin .eq. 1) then
        nl = 1
      else if(line_edge .and. klin .gt. 1) then
        nl = 2
      else
        nl = nlcntx
      end if

      do 100 i = 1,nl
        do 200 j = 1,necntx
          !ide = kele + (j-imv)
          ide = j
          !if(ide .gt. 0 .and. ide .le. maxele) then
          if(ide .gt. 0 ) then
            if (indat(ide,i,band).gt.0.0 .and.
     +          indat(ide,i,band).lt.1000000.0) then
              rgdata(i,j) = indat(ide,i,band)
            else
              rgdata(i,j) = bad_data
            end if
          else
            rgdata(i,j) = bad_data
          end if
  200   continue
  100 continue

c     Get standard deviation central pixel.

      n = 0.d0
      sum = 0.d0
      sumsq = 0.d0
      do 400 i = 1,nl
        do 500 j = 1,necntx
	  tbb = dble(rgdata(i,j))
	  if(abs(tbb-bad_data) .gt. 0.1) then
	    n = n + 1.d0
	    sum = sum + tbb          
	    sumsq = sumsq + (tbb * tbb)
          end if
 500    continue
 400  continue 

      if(n .gt. 1.d0) then
        sqsum = sum * sum
        num = (n * sumsq) - sqsum
        den = n * (n-1)
        sig = dsqrt(num / den)
        mn = sum / n
      else if(n .eq. 1.d0) then
        sig = 0.d0
        mn = sum
      else
        sig = dble(bad_data)
        mn = dble(bad_data)
      end if 
      sigma = sngl(sig)
      mean = sngl(mn)

c ....debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(3f12.4)') rgdata
!        write(h_output,'(10x,'' Standard dev. calc. '',/,l2,5i6,8f12.4/)')
!     *    line_edge,klin,nl,kele,maxele,band,sqsum,sumsq,n,num,den,
!     *    sig,sigma,mean
!      endif
c ................................................................

      return
      end
