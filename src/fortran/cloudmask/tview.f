      subroutine tview(key,xmu,bt11,corr)

      implicit none
      save

C
c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION: subroutine which computes the bi-dimensional linear or
c              quadratic interpolation (lagrange form) of a given value
c              to a table of scan angle and 11um brightness temperature
C              dependent 11-12 micron brightness temperature 
c              differences.  These values were taken from APOLLO.
c 
c!Input Parameters:
c key           linear (1) or quadratic (2) interpolation               
c xmu           secant of the viewing zenith angle
c bt11          pixel 11 micron brightness temperature
c
c!Output Parameters:
c corr          computed 11-12 micron BTDIF threshold for 
c               thin cirrus detection
c!Revision History:
c Extended temperature range (J. Key version)
c 06/04 Collection 5     R. Frey
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c Revision 01.00  1995/7/17 02:35:00
c K. Strabala (kathys@ssec.wisc.edu)
c Initial delivery of software.  Modified to comply with ESDIS and
c MODIS standards.  It was a pain.
c Revised to include extension of thresholds to cold temperatures (J. Key).
c R. Frey  2003/10/31 10:55
c
c!END
c-----------------------------------------------------------------------

c ... scalar arguments ..
      real bt11,corr,xmu
      integer key
c ...
c ... local scalars ..
      real lt0,lt1,lt2,lu0,lu1,lu2,p,p0,p1,p2,t,t0,t1,t2,u,u0,u1,u2
      integer i,i0,i1,i2,ii,j,j0,j1,j2,jj,Iflag,debug,h_output
c ...
c ... local arrays ..
      real tab(5,13),ttab(13),utab(5)

!      common / bug / debug, h_output
c ...
c ... data statements ..
 
      data utab/2.00,1.75,1.5,1.25,1.00/
      data ttab /190.,200.,210.,220.,230.,240.,250.,
     *           260.,270.,280.,290.,300.,310./
      data (tab(i,1),i=1,5)/0.559,0.542,0.520,0.491,0.450/
      data (tab(i,2),i=1,5)/0.424,0.416,0.405,0.391,0.370/
      data (tab(i,3),i=1,5)/0.286,0.294,0.305,0.319,0.340/
      data (tab(i,4),i=1,5)/0.137,0.162,0.194,0.238,0.300/
      data (tab(i,5),i=1,5)/0.123,0.156,0.199,0.257,0.340/
      data (tab(i,6),i=1,5)/0.198,0.240,0.294,0.367,0.470/
      data (tab(i,7),i=1,5)/0.333,0.366,0.409,0.467,0.550/
      data (tab(i,8),i=1,5)/0.696,0.704,0.715,0.729,0.750/
      data (tab(i,9),i=1,5)/1.217,1.184,1.141,1.083,1.000/
      data (tab(i,10),i=1,5)/3.184,2.926,2.591,2.140,1.500/
      data (tab(i,11),i=1,5)/5.178,4.854,4.433,3.866,3.060/
      data (tab(i,12),i=1,5)/8.269,7.885,7.389,6.720,5.770/
      data (tab(i,13),i=1,5)/12.452,11.985,11.381,10.567,9.410/
c ...
c ... initialize variables
      lt0 = 0.0
      lt1 = 0.0
      lt2 = 0.0
      lu0 = 0.0
      lu1 = 0.0
      lu2 = 0.0
      p = 0.0
      p0 = 0.0
      p1 = 0.0
      p2 = 0.0
      t = 0.0
      t0 = 0.0
      t1 = 0.0
      t2 = 0.0
      u = 0.0
      u1 = 0.0
      u2 = 0.0
      i0 = 0
      i1 = 0
      i2 = 0
      j0 = 0
      j1 = 0
      j2 = 0
      jj = 0
 
c ...
c ... bounds check
      u = xmu
      t = bt11
      if (u.gt.utab(1)) u = utab(1)
      if (u.lt.utab(5)) u = utab(5)
      if (t.lt.ttab(1)) t = ttab(1)
      if (t.gt.ttab(13)) t = ttab(13)

      Iflag = 0
      do 1 i = 2, 5
          ii = i
          if (u.ge.utab(i) .AND. Iflag.eq.0) then
             Iflag = 1
             if (key.eq.1) then
                i0 = ii - 1
                i1 = ii
             else
                if (ii.eq.5) then
                   i0 = ii - 2
                   i1 = ii - 1
                   i2 = ii
                else
                   i0 = ii - 1
                   i1 = ii
                   i2 = ii + 1
                end if
             end if
          end if
    1 continue

      Iflag = 0
      do 6 j = 2, 13
          jj = j
          if (t.le.ttab(j) .AND. Iflag.eq.0) then
             Iflag = 1
             if (key.eq.1) then
                 j0 = jj - 1
                 j1 = jj
             else
                if (jj.eq.13) then
                   j0 = jj - 2
                   j1 = jj - 1
                   j2 = jj
                else
                   j0 = jj - 1
                   j1 = jj
                   j2 = jj + 1
                end if
             end if
          end if
    6 continue

c ... branch on scheme type
      IF (key.eq.1) THEN

c ...   linear scheme
c ...   designate index values

        u0 = utab(i0)
        u1 = utab(i1)
        t0 = ttab(j0)
        t1 = ttab(j1)

c ...   lagrange polynomials
        lu0 = (u-u1)/ (u0-u1)
        lu1 = (u-u0)/ (u1-u0)
        lt0 = (t-t1)/ (t0-t1)
        lt1 = (t-t0)/ (t1-t0)

c ...   interpolating polynomials for the first dimension
        p0 = tab(i0,j0)*lu0 + tab(i1,j0)*lu1
        p1 = tab(i0,j1)*lu0 + tab(i1,j1)*lu1

c ...   interpolating polynomial for second dimension
        p = p0*lt0 + p1*lt1
c ...
        corr = p

c ..... Debug statement.  ........................................
!        if(debug .eq. 4) then
!          write(h_output,'(''tview '',7f8.3)') u0,u1,t0,t1,p0,p1,p
!        end if
c ................................................................

        return

      END IF
c ...
c ... quadratic scheme
c ... designate index values
      u0 = utab(i0)
      u1 = utab(i1)
      u2 = utab(i2)
      t0 = ttab(j0)
      t1 = ttab(j1)
      t2 = ttab(j2)
c ...
c ... lagrange polynomials
      lu0 = (u-u1)* (u-u2)/ (u0-u1)/ (u0-u2)
      lu1 = (u-u0)* (u-u2)/ (u1-u0)/ (u1-u2)
      lu2 = (u-u0)* (u-u1)/ (u2-u0)/ (u2-u1)
      lt0 = (t-t1)* (t-t2)/ (t0-t1)/ (t0-t2)
      lt1 = (t-t0)* (t-t2)/ (t1-t0)/ (t1-t2)
      lt2 = (t-t0)* (t-t1)/ (t2-t0)/ (t2-t1)
c ...
c ...
c ... interpolating polynomials for the first dimension
      p0 = tab(i0,j0)*lu0 + tab(i1,j0)*lu1 + tab(i2,j0)*lu2
      p1 = tab(i0,j1)*lu0 + tab(i1,j1)*lu1 + tab(i2,j1)*lu2
      p2 = tab(i0,j2)*lu0 + tab(i1,j2)*lu1 + tab(i2,j2)*lu2
c ...
c ... interpolating polynomial for second dimension
      p = p0*lt0 + p1*lt1 + p2*lt2
c ...
      corr = p
      return
      end
