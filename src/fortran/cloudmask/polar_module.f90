module polar_module

!!-----------------------------------------------------------------------
!! !F90                                                                  
!!
!! !Description: 
!!    cloud mask algorithm for land
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
subroutine polar_day(pxldat,vza,snglnt,visusd,refang,vrused,               &
                     cirrus_vis,land,ice,snow,desert,coast,eco_type,       &
                     uniform,hi_elev,kele,indat,nmtests,testbits,tbadj,    &
                     antarctic,sh_ocean,sfctmp,qa_bits,confdnc,            &
                     btclr,is_cold_sfc)

!      implicit none
!      save

!----------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for providing conditional input parameter pertaining to 
!     polar daytime processing.
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Current pixel viewing angle
! snglnt        Logical variable flagging sunglint pixels
! visusd        Logical variable indicating whether vis data used or not
! visusd        Logical variable indicating when reflect ratio test
!               can be implemented
! refang        reflectance angle for current pixel
! cirrus_vis    Logical variable flagging thin cirrus contaminated
!               scenes in the visible
! land          Logical variable flagging land backgrounds
! ice           Logical variable flagging ice backgrounds
! snow          Logical variable flagging snow backgrounds
! desert        Logical variable flagging desert backgrounds
! coast         Logical variable flagging coast backgrounds
! uniform       Logical variable flagging contexts with similar surface properties
! hi_elev       Logical flag indicating elevations > 2000 meters.
! kele          Current pixel being processed
! indat         array containing 'nlcntx' lines of data
! tbadj         11 um brightness temperature threshold adjustment for 
!               deserts (based on elevation)
! eco_type      Byte variable containing ecosystem index for current pixel
! antarctic     Logical flag indicating Antarctic region (< -60 latitude)
! sh_ocean      Logical variable indicating shallow ocean
! sfctmp        SST for current pixel
!
!!Output Parameters:
! nmtests       Number of tests applied to this pixel
! testbits      Byte array containing cloud mask results
! qa_bits       Byte array containing qa bit results
! confdnc       Current pixel unobstructed confidence
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Updated calling arguments to PolarDay_snow and Antarctic_day.
! 10/04 Collection 5b R. Frey
! Updated calling arguments to PolarDay_ocean and Antarctic_day.
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!END
!----------------------------------------------------------------------

      include 'global.inc'

!     scalar arguments
      real vza,confdnc,tbadj,refang,sfctmp
      logical visusd,land,snow,snglnt,ice,coast,desert,vrused,cirrus_vis, &
              uniform,hi_elev,antarctic,sh_ocean
      integer nmtests,kele
      byte eco_type

!     array arguments
      real pxldat(inband),indat(necntx,nlcntx,inband),btclr(7)
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)

!     local scalars
      integer debug,h_output

!     external subroutines
      external PolarDay_snow,PolarDay_land,PolarDay_ocean,PolarDay_coast,    &
               PolarDay_desert,PolarDay_desert_c,chk_land,chk_coast,         &
               chk_spatial_var,chk_shallow_water

!     Common statement for debug purposes
!      common / bug / debug, h_output

!----------------------------------------------------------------------

!     Debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Using polar_day processing path'',/)')
!      endif

!----------------------------------------------------------------------

!     Polar processing further defined based on ecosystem map, land/sea
!     tag file, and snow/ice data.

      if(antarctic .and. land) then
!        Antarctica
         call Antarctic_day(pxldat,visusd,testbits,qa_bits,nmtests,confdnc,btclr,is_cold_sfc)
         
         if (nmtests == 0) then   ! added by minmin to save antarctica cloudmask test
            call PolarDay_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,  &
                               testbits,qa_bits,nmtests,confdnc,      &
                               btclr)
         endif 

      else if(snow .or. ice) then
!        snow or ice covered surfaces

         call PolarDay_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,  &
                            testbits,qa_bits,nmtests,confdnc,      &
                            btclr)

      else if (land) then
!        Land surfaces

         if (desert .and. coast) then

            call PolarDay_desert_c(pxldat,vza,visusd,cirrus_vis,hi_elev,   &
                                   testbits,qa_bits,nmtests,confdnc,       &
                                   btclr,is_cold_sfc)

         else if (coast) then

            call PolarDay_coast(pxldat,vza,visusd,cirrus_vis,hi_elev,      &
                                testbits,qa_bits,nmtests,confdnc,          &
                                btclr,is_cold_sfc)

         else if (desert) then

            call PolarDay_desert(pxldat,vza,visusd,cirrus_vis,hi_elev,     &
                                 testbits,qa_bits,nmtests,confdnc,         &
                                 btclr,is_cold_sfc)

         else 

            call PolarDay_land(pxldat,vza,visusd,vrused,cirrus_vis,        &
                               hi_elev,testbits,qa_bits,nmtests,           &
                               confdnc,btclr,is_cold_sfc)

         endif

      else

!        water surface

         call PolarDay_ocean(pxldat,vza,snglnt,visusd,cirrus_vis,         &
                             refang,sfctmp,sh_ocean,testbits,qa_bits,     &
                             nmtests,confdnc,btclr)

      end if

!----------------------------------------------------------------------

!     Debug statement 
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,'' Polar day premlinary confidence '',
!     +        f10.2,3l5/)') confdnc,land,coast,sh_ocean
!      endif

!----------------------------------------------------------------------

!     Perform clear sky confidence confirmation tests.

!     Under certain conditions, apply spatial variability test.
      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform .and. (.not. land)) then
        call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

!     Check land.
      if(land .and. (.not. (snow .or. ice))) then
        if(confdnc .le. 0.95) then
          call chk_land(pxldat,eco_type,desert,tbadj,confdnc,qa_bits,testbits)
        end if
      end if

!     Check coastal regions.
      if(coast .and. (.not. (snow .or. ice))) then
        call chk_coast(pxldat,confdnc,qa_bits,testbits)
      end if

!     Perform "final test" in shallow ocean conditions.
      if(sh_ocean .and. (.not. ice)) then
         call chk_shallow_water(pxldat,confdnc,qa_bits,testbits)
      end if

!      return
end subroutine polar_day
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine polar_nite(pxldat,vza,land,ice,snow,desert,hi_elev,         &
                      sfctmp,eco_type,nmtests,testbits,uniform,        & 
                      indat,kele,antarctic,sh_ocean,qa_bits,confdnc,   &
                      btclr,is_cold_sfc)

!----------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for providing conditional input parameters pertaining
!     to polar nighttime processing.
!
!!Input parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Viewing zenith angle
! land          Logical variable flagging land backgrounds
! ice           Logical variable flagging ice backgrounds
! snow          Logical variable flagging snow backgrounds
! desert        Logical variable flagging desert backgrounds
! hi_elev       Logical variable flagging high elevations (> 2000 meters)
! sfctmp        SST from Reynolds Blended data set
! eco_type      Ecosystem index from Olson ecosystem data set
! uniform       Logical variable flagging uniform background
! kele          Current element number being processed
! indat         Array containing 'nlcntx' lines of radiance data
! antarctic     Logical variable indicating pixel is south of -60 lat
! sh_ocean      Logical variable indicating ocean depths < 50 m or within
!               5 km of shoreline
!
!!Output Parameters:
! nmtests       Number of tests applied to this pixel
! testbits      Byte array containing cloud mask results
! qa_bits       Byte array containing qa bit results
! confdnc       Current pixel unobstructed confidence
!
!!END
!----------------------------------------------------------------------

      include 'global.inc'

!     scalar arguments
      real confdnc,sfctmp
      logical land,snow,ice,uniform,desert,hi_elev,antarctic,sh_ocean
      integer nmtests,kele
      byte eco_type

!     array arguments
      real pxldat(inband),vza,indat(necntx,nlcntx,inband),btclr(7)
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)

!     local scalars
      integer debug,h_output

!     local arrays

!     external subroutines
      external PolarNite_snow,PolarNite_land,PolarNite_ocean,chk_spatial_var

!     Common statement for debug purposes
!      common / bug / debug, h_output

!----------------------------------------------------------------------

!     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Using polar_nite processing path '',/)')
!      endif

!----------------------------------------------------------------------

      if (snow .or. ice) then

!        snow or ice covered surfaces
         call PolarNite_snow(pxldat,vza,hi_elev,land,antarctic,testbits,   &
                             qa_bits,nmtests,confdnc,btclr)

      else if (land) then

!        Land surfaces
         call PolarNite_land(pxldat,vza,desert,hi_elev,sfctmp,             &
                             eco_type,testbits,qa_bits,nmtests,confdnc,    &
                             btclr,is_cold_sfc)

      else

!        Water surfaces.
         call PolarNite_ocean(indat,kele,pxldat,vza,sfctmp,sh_ocean,       &
                              uniform,testbits,qa_bits,nmtests,confdnc,    &
                              btclr)

      end if

!----------------------------------------------------------------------

!     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,'' Polar Nite preliminary confidence '',
!     +        f10.2,/)') confdnc
!      endif

!----------------------------------------------------------------------

!     Perform clear sky confidence confirmation tests. 

      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform .and. (.not. land)) then
        call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

!----------------------------------------------------------------------

end subroutine polar_nite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------- end MODULE ---------------------------------

end module polar_module
