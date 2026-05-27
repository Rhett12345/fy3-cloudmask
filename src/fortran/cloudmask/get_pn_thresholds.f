      subroutine get_pn_thresholds(bt_11,bt_bnds,th_low,th_mid1,
     *                             th_mid2,th_mid3,th_hi,locut,hicut,
     *                             midpt,power)

      implicit none
      save

c---------------------------------------------------------------------
C!F77
c
c!Description:
c     Routine for setting cloud test thresholds for 11 minus 3.9 um
c     test in polar night conditions.
c
c!Input parameters:
c bt_11         11 um brightness temperature for current pixel
c bt_bnds       11 um brightness temperature bounds for current test
c th_low        Test thresholds for 'bt_11' lt 'bt_bnds(1)'.
c th_mid1       Test thresholds for 'bt_11' between 'bt_bnds(1)' and 
c                    'bt_bnds(2)'.
c th_mid2       Test thresholds for 'bt_11' between 'bt_bnds(2)' and
c                    'bt_bnds(3)'.
c th_mid3       Test thresholds for 'bt_11' between 'bt_bnds(3)' and
c                    'bt_bnds(4)'.
c th_hi         Test thresholds for 'bt_11' gt 'bt_bnds(4)'.
c
c!Output Parameters:
c locut         Zero confidence value of clear sky confidence interval
c hicut         1.0 confidence value of clear sky confidence interval
c midpt         0.5 confidence value of clear sky confidence interval
c power         Power of clear sky confidence curve
c
c!Revision History:
c Original version - 03/07/2002     
c R. Frey            MODIS Science Group at UW-Madison
c Revised          - 10/30/2003
c R. Frey
c 06/04 Collection 5  R. Frey
c Added 'th_mid1','th_mid2','th_mid3'.
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

c     Array arguments.
      real bt_bnds(4),th_low(4),th_hi(4),th_mid1(4)
      real th_mid2(4),th_mid3(4)

c     Scalar arguments.
      real bt_11,locut,hicut,midpt,power

c     Local scalars.
      real lo_tmp,hi_tmp,lo_tmp_thr,hi_tmp_thr,a,conf_range

c     Compute clear sky confidence interval boundaries and mid-point.

      if(bt_11 .lt. bt_bnds(1)) then

        hicut = th_low(3)
        locut = th_low(1)
        midpt = th_low(2)
        power = th_low(4)

      else if(bt_11 .gt. bt_bnds(4)) then

        hicut = th_hi(3)
        locut = th_hi(1)
        midpt = th_hi(2)
        power = th_hi(4)

      else

        if(bt_bnds(2) .eq. 0.0 .and. bt_bnds(3) .eq. 0.0) then

          lo_tmp = bt_bnds(1)
          hi_tmp = bt_bnds(4)
          lo_tmp_thr = th_mid1(1)
          hi_tmp_thr = th_mid1(2)
          power = th_mid1(4)
          conf_range = th_mid1(3)

        else if(bt_11 .ge. bt_bnds(1) .and. bt_11 .lt. bt_bnds(2)) then

          lo_tmp = bt_bnds(1)
          hi_tmp = bt_bnds(2)
          lo_tmp_thr = th_mid1(1)
          hi_tmp_thr = th_mid1(2)
          power = th_mid1(4)
          conf_range = th_mid1(3)

        else if(bt_11 .ge. bt_bnds(2) .and. bt_11 .lt. bt_bnds(3)) then

          lo_tmp = bt_bnds(2)
          hi_tmp = bt_bnds(3)
          lo_tmp_thr = th_mid2(1)
          hi_tmp_thr = th_mid2(2)
          power = th_mid2(4)
          conf_range = th_mid2(3)

        else if(bt_11 .ge. bt_bnds(3) .and. bt_11 .le. bt_bnds(4)) then

          lo_tmp = bt_bnds(3)
          hi_tmp = bt_bnds(4)
          lo_tmp_thr = th_mid3(1)
          hi_tmp_thr = th_mid3(2)
          power = th_mid3(4)
          conf_range = th_mid3(3)

        end if

        a = (bt_11 - lo_tmp) / (hi_tmp - lo_tmp)
        midpt = lo_tmp_thr + a*(hi_tmp_thr - lo_tmp_thr)
        hicut = midpt - conf_range
        locut = midpt + conf_range
       
      end if

      return
      end
