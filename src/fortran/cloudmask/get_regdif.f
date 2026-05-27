      subroutine get_regdif(indat,kele,rgdata,diff)
      implicit none
      save

c---------------------------------------------------------------------
c!F77
c
c!DESCRIPTION: 
c     computes reflectance or brightness temperature differences
c     between the (center) pixel of interest and the surrounding
c     eight pixel values.  also outputs the input data for all 9
c     pixel positions. This represents the spatial variability 
c     tests for now.
c
c!Input parameters:
c indat         Array containing nlcltx number of lines of
c               reflectance or BT values for each channel
c kele          First element pixel in the region (context)
c
c!Output Parameters:
c rgdata        Array containing current pixel and surrounding pixel
c               reflectance and brightness temperature information
c diff          Array of surrounding pixel reflectance or brightness
c               temperature differences from center pixel value
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------------

      INCLUDE 'global.inc'
c ... scalar arguments ..
      integer kele
c ...
c ... array arguments ..
      real diff(var_band,8),rgdata(nlcntx,necntx,var_band)
      real indat(necntx,nlcntx,inband)
c     +     indat(npixel,nlcntx,inband)
c ...
c ... local scalars ..
      real a,b
      integer i,i1,i2,ide,imv,j,j1,j2,k,kk
c ...
c ... local arrays ..
      integer i1loc(8),i2loc(8),j1loc(8),j2loc(8),band(var_band)
c ...
c ... data statements ..
c ... order for differences from central value to be taken 
      data i1loc/1,1,1,2,2,3,3,3/
      data i2loc/2,2,2,2,2,2,2,2/
      data j1loc/1,2,3,1,3,1,2,3/
      data j2loc/2,2,2,2,2,2,2,2/

c ... bands to do spatial variability tests on
!      data band /1,31/
      data band /3,24/
      
c ...
c ... initialize variables
      i1 = 0
      i2 = 0
      imv = 0
      ide = 0
      j1 = 0
      j2 = 0

c ... initialize arrays     
      do 10 i = 1 , nlcntx
        do 20 j = 1 , necntx
          do 30 k = 1 , var_band
            rgdata(i,j,k) = 0.0
   30     continue
   20   continue
   10 continue
c
      do 40 i = 1 , var_band
        do 50 j = 1 , 8
          diff(i,j) = 0.0
   50   continue
   40 continue
  
 
c ... put pixel data into context data array.
      imv = ((necntx-1)/2) + 1
      do 100 i = 1,nlcntx
        do 200 j = 1,necntx
c          ide = kele + (j-imv)
           ide = j
          do 300 k = 1,var_band
            if (indat(ide,i,band(k)).gt.0.0 .and. 
     +        indat(ide,i,band(k)).lt.1000000.0) then
              rgdata(i,j,k) = indat(ide,i,band(k)) 
            else
              rgdata(i,j,k) = bad_data
            end if
  300     continue
  200   continue
  100 continue
 
c ... get differences about pixel of interest for all channels.
      do 400 k = 1,var_band
        do 500 kk = 1,8
          i1 = i1loc(kk)
          j1 = j1loc(kk)
          i2 = i2loc(kk)
          j2 = j2loc(kk)
          a = rgdata(i1,j1,k)
          b = rgdata(i2,j2,k)
          if (abs(a-bad_data).gt.0.1 .AND. abs(b-bad_data).gt.0.1) then
            diff(k,kk) = a - b
          else
            diff(k,kk) = bad_data
          end if
  500   continue
  400 continue
 
      return
      end
