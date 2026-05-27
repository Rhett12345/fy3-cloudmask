module water_module

!!-----------------------------------------------------------------------
!! !F90                                                                  
!!
!! !Description: 
!!    cloud mask algorithm for water 
!!
!! !Input parameters
!!    none
!! 
!! !Output parameters
!!    none
!!
!! !Author's information
!!    Author: Min Min/Wu Xiao/Zheng Zhaojun/Liu Ruixia/Zhang Miao/Yang Changjun/Qiu Hong
!!    E-mail: minmin@cma.gov.cn
!!    Tel   : 86-010-68406763
!!    National Satellite Meteorological Center, CMA 
!!  
!! !END
!!----------------------------------------------------------------------

use names_module
use data_arrays_module
use constant
use thresholds_read_module
!use cloudmask_data_arrays

implicit none
!+++++++++++++++++++ step 1: Variables +++++++++++++++++++++++++++++++++

             
contains
!+++++++++++++++++++ step 2: Subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!subroutine water_day(pxldat,vza,snglnt,visusd,refang,cirrus_vis,
!                     sfctmp,hi_elev,uniform,ice,maxele,klin,
!                     line_edge,sh_ocean,indat,kele,nmtests,
!                     testbits,qa_bits,confdnc)
subroutine water_day(pxldat,vza,snglnt,visusd,refang,cirrus_vis,     &
                     sfctmp,hi_elev,uniform,ice,snow,kele,klin,      &
                     line_edge,sh_ocean,indat,nmtests,               &
                     testbits,qa_bits,confdnc,btclr)
!      implicit none
!      save

!---------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for setting appropriate flags and processing path
!     for daytime observations over water surfaces.
!
!     If the confidence determined from spectral tests is
!     uncertain, then other tests may be applied. 
!
!!Input parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Current pixel viewing angle
! snglnt        Logical variable flagging sunglint pixels
! visusd        Logical variable indicating whether vis data used or not
! refang        Reflectance angle for current pixel
! cirrus_vis    Logical variable flagging thin cirrus contaminated
!               scenes in the visible
! sfctmp        SST from ancillary data
! hi_elev       Logical variable indicating elevations > 2000 meters
! uniform       Logical variable indicating uniformity of context
! ice           Logical variable flagging ice backgrounds
! sh_ocean      Logical variable indicating shallow ocean
! indat         Array containing nlcntx lines of data
! kele          Current granule element number being processed
!
!!Output Parameters:
! nmtests       Number of tests applied to this pixel
! testbits      Byte array containing cloud mask results
! qa_bits       Byte array containing qa bit results
! confdnc       Current pixel unobstructed confidence
!
!---------------------------------------------------------------------

      include 'global.inc'

!     scalar arguments
      real vza,confdnc,refang,sfctmp
      integer kele,nmtests,maxele,klin
      logical visusd,snglnt,uniform,ice,cirrus_vis,hi_elev,sh_ocean,line_edge
      logical snow

!     array arguments
      real pxldat(inband),indat(necntx,nlcntx,inband),btclr(7)
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)

!     local scalars
!      integer debug,h_output

!     external subroutines
      !external ocean_day,Day_snow,chk_spatial_var,chk_sunglint,chk_shallow_water

!     Common statement for debug purposes
!      common / bug / debug, h_output

!---------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Using water_day processing path '',/)')
!      endif

!---------------------------------------------------------------------

!     Further define processing path.

      if (ice .or. snow) then
!        Processing for ice-covered scenes.

         call Day_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,testbits, &
                       qa_bits,nmtests,confdnc,btclr)

      else 
  
!        Normal ocean processing
         call ocean_day(pxldat,vza,snglnt,visusd,cirrus_vis,sfctmp,   &
                        refang,sh_ocean,testbits,qa_bits,nmtests,     &
                        confdnc,btclr)

      endif

!---------------------------------------------------------------------

!     debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,'' Water day premlinary confidence '',
!     +        f10.2)') confdnc
!        write(h_output,'(10x/,'' uniform,snglnt,sh_ocean,ice '',
!     +        4l5/)') uniform,snglnt,sh_ocean,ice
!      endif

!---------------------------------------------------------------------

!     Perform clear sky confidence confirmation tests.

!     Under certain conditions, apply the spatial variability test.
      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform) then
         call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

!     Perform clear sky restoral tests in sun-glint regions.
      if(snglnt .and. uniform .and. confdnc.lt.0.95) then
!         call chk_sunglint(indat,pxldat,kele,confdnc,maxele,klin,    &
!                           line_edge,qa_bits,testbits)
         call chk_sunglint(indat,pxldat,kele,confdnc,klin,    &
                           line_edge,qa_bits,testbits)
      end if

!     Perform clear sky restoral tests in shallow ocean conditions.
      if(sh_ocean .and. (.not. ice)) then
         call chk_shallow_water(pxldat,confdnc,qa_bits,testbits)
      end if

!---------------------------------------------------------------------

!      return
end subroutine water_day
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
subroutine water_nite(pxldat,vza,uniform,ice,snow,indat,sfctmp,sh_ocean, &
                      kele,nmtests,testbits,qa_bits,confdnc,        &
                      btclr)

!      implicit none
!      save

!---------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for setting appropriate flags and processing path
!     for nighttime observations over water surfaces.
!
!     If the confidence determined from the water background is
!     uncertain, than a spatial variability test is applied if
!     the pixels within the processing context are all from the
!     same ecosystem (not including a coastline).
!
!!Input parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! uniform       Logical variable flagging context of similar ecosystem
! ice           Locical variable flagging ice backgrounds
! indat         Array containing nlcntx lines of data
! sfctmp        Surface temperature for current pixel
! sh_ocean      Logical flag indicating ocean depths < 50 m or within
!               5 km of shorlines
! kele          Current granule element number being processed
! vza           Viewing zenith angle
!
!!Output Parameters:
! nmtests       Number of tests applied for this pixel
! testbits      Byte array containing cloud mask results
! qa_bits       10-byte array containing qa bit results
! confdnc       Current pixel unobstructed confidence
!
!!Revision History:
! Removed 11 micron standard deviation calculation
! Added 'lnd' and 'sfctmp' for Collection 5 processing.
! 10/04  Collection 5b  R. Frey
! Updated calling argument list for ocean_nite.
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!END
!---------------------------------------------------------------------
!
      include 'global.inc'
    
!     scalar arguments
      real vza,confdnc,sfctmp
      integer kele,nmtests
      logical uniform,ice,sh_ocean
      logical snow  ! added by minmin 20190110
!
!     array arguments
      real pxldat(inband),indat(necntx,nlcntx,inband),btclr(7)
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)

!     local scalars
      integer debug,h_output
      logical lnd

!     local arrays

!     external subroutines
!      external ocean_nite,Nite_snow,chk_spatial_var

!     Common statement for debug purposes
!      common / bug / debug, h_output

!---------------------------------------------------------------------
!
!     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Using water_nite processing path '',/)')
!      endif

!---------------------------------------------------------------------

      if (ice .or. snow) then

!        Ice processing

         lnd = .false.
         call Nite_snow(pxldat,vza,lnd,testbits,qa_bits,nmtests,confdnc,btclr)

      else
!        Normal nighttime ocean processing

         call ocean_nite(indat,kele,pxldat,vza,sfctmp,sh_ocean,   &
                         uniform,testbits,qa_bits,nmtests,confdnc,&
                         btclr)

      endif

!---------------------------------------------------------------------

!     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,'' Water nite preliminary confidence '',
!     +        f10.2,/)') confdnc
!      endif

!---------------------------------------------------------------------

!     Perform clear sky confidence confirmation tests.

      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform) then
        call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

!---------------------------------------------------------------------

!      return
end subroutine water_nite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!-------------------------- end MODULE ---------------------------------

end module water_module
