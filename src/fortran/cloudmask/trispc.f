      real function trispc(tdf1)

      implicit none
      save

c!F77 ************************************************************
c
c!Description:
c ... regression function relating the pw values to the 8-11um btdif.
c ... The difference is related to the amount of pw in the atmosphere
c ... due to the weak water vapor lines present in this spectral
c ... region.  the regressions were determined from actual hirs
c ... data, and provides a threshold for clear which can be compared
c ... to actual 8-11 differences
c
c!Input Parameters:
c tdf1      11-12 micron BTDIF observed
c
c!Output Parameters: None
c 
c!Revision History:
c
c!Team-unique Header:
c
c    This software is developed by the MODIS Science Data Support Team
c    for the National Aeronautics and Space Administration,
c    Goddard Space Flight Center, under contract NAS5-32373.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!Design Notes:
c   Returns the function value trispc as the 8-11 micron calculated
c    clear value.
c
c!END****************************************************************
c
c ...
c ... scalar arguments ..
      real tdf1
c ...
c ... local scalars ..
c ... variables: a's regression coefficients
      real a1,a2,a3,a4,x
c ...
c ... data statements ..
c ... coefficients for hirs data
c     data a1 /-1.7681/, a2 /-3.729/ ,a3 /1.054/, a4 /-0.102/
      data a1 /2.7681/, a2 /-3.729/ ,a3 /1.054/, a4 /-0.102/

      x = tdf1
c     trispc = a1 + a2*tdf1 + a3*(tdf1*tdf1) + a4*(tdf1**3) 
      trispc = a1 + x*(a2 + x*(a3 + x*a4))
c ...
      return
      end
