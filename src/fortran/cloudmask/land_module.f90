module land_module

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
subroutine land_day(pxldat,vza,visusd,vrused,cirrus_vis,             &
                    desert,coast,snow,ice,hi_elev,tbadj,eco_type,    &
                    testbits,qa_bits,nmtests,confdnc,btclr,is_cold_sfc, &
                    indat_11um, ielem, iline)

!      implicit none
!      save

!-----------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for setting appropriate flags and processing path
!     for daytime observations over land surfaces.
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Current pixel viewing angle
! visusd        Logical variable indicating whether vis data used or not
! vrused        Logical variable indicating when vrat test can be used
! cirrus_vis    Logical variable flagging thin cirrus contaminated
!               scenes in the visible
! hi_elev       Logical variable indicating high elevation (> 2000 meters)
! tbadj         11 um brightness temperature threshold adjustment for 
!               deserts
! desert        Logical variable indicating desert ecosystems
! coast         Logical variable indicating coast ecosystems
! snow          Logical variable flagging snow backgrounds
! ice           Logical variable flagging ice backgrounds
! eco_type      Byte variable containing ecosystem index for current pixel
!
!!Output Parameters:
! testbits      Byte array containing cloud mask results
! qa_bits       Byte array containing qa bit results
! nmtests       Number of tests actually applied to the given pixel
! confdnc       Current pixel unobstructed confidence
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Updated calling arguments to Day_snow.
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!END
!-----------------------------------------------------------------------

      include 'global.inc'

!     scalar arguments
      real vza,confdnc,tbadj
      integer nmtests
      logical visusd,vrused,snow,ice,desert,coast,cirrus_vis,hi_elev
      byte eco_type
      integer ielem, iline

!     array arguments
      real pxldat(inband), btclr(7)
      real indat_11um(3,3)
      byte testbits(6),qa_bits(10)

!     local scalars
      integer h_output,debug
      integer(kind=1) :: is_cold_sfc

!     external subroutines
      external LandDay_desert,LandDay,Day_snow,LandDay_desert_c, &
               LandDay_coast,chk_land,chk_coast

!     Common statement for debug purposes
!      common / bug / debug, h_output

!-----------------------------------------------------------------------

!     debug statement 
!     if (debug .gt. 0) then
!        write(h_output,'(10x/,''Using land_day processing path '',/)')
!      endif

!-----------------------------------------------------------------------

!     Land processing further defined based on ecosystem map, land/sea
!     tag and snow/ice data.

      if (snow .or. ice) then

         call Day_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,   &
                       testbits,qa_bits,nmtests,confdnc,btclr)

      else if (desert .and. coast) then

         call LandDay_desert_c(pxldat,vza,visusd,cirrus_vis,   &
                               hi_elev,testbits,qa_bits,       &
                               nmtests,confdnc,btclr,is_cold_sfc)

      else if (coast) then

         call LandDay_coast(pxldat,vza,visusd,cirrus_vis,      &
                            hi_elev,testbits,qa_bits,nmtests,  &
                            confdnc,btclr,is_cold_sfc)
 
      else if (desert) then


         call LandDay_desert(pxldat,vza,visusd,cirrus_vis,     &
                             hi_elev,testbits,qa_bits,nmtests, &
                             confdnc,btclr,is_cold_sfc)
 
      else
 
         call LandDay(pxldat,vza,visusd,vrused,cirrus_vis,     &
                      hi_elev,testbits,qa_bits,nmtests,confdnc,& 
                      btclr,is_cold_sfc)
         !print*,'confdnc2=',confdnc
      endif

!-----------------------------------------------------------------------

!     Perform clear sky confidence confirmation tests.

!     Check non-snow covered land regions.

      if(.not. (snow .or. ice)) then
        if(confdnc .le. 0.95) then
          call chk_land(pxldat,eco_type,desert,tbadj,confdnc,qa_bits,testbits)
        end if
      end if

!     Check coastal regions.
      if(coast .and. (.not. (snow .or. ice))) then
        call chk_coast(pxldat,confdnc,qa_bits,testbits)
      end if

!-----------------------------------------------------------------------

!      return
end subroutine land_day
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine land_nite(pxldat,plat,vza,ice,snow,coast,tbadj,desert,   &
                     hi_elev,sh_lake,sfctmp,eco_type,nmtests,       &
                     testbits,qa_bits,confdnc,ptwp,                 &
                     btclr,is_cold_sfc,                             &
                     indat_11um, ielem, iline)

!      implicit none
!      save

!----------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for setting appropriate flags and processing path
!     for nighttime observations over land surfaces.
!
!!Input Parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           viewing zenith angle
! ice           Logical variable flagging ice backgrounds
! snow          Logical variable flagging snow backgrounds
! coast         Logical variable indicating coast processing
! tbadj         11 um brightness temperature elevation adjustment
! desert        Logical variable flagging arid ecosystem
! hi_elev       Logical variable flagging high elevation (> 2000 m)
! sh_lake       Logical flag indicating shallow inland lakes
! sfctmp        Surface air temperature from model output 
! eco_type      Ecosystem index
!
!!Output Parameters:
! nmtests       Number of tests applied to this pix
! testbits      Byte array containing cloud mask results
! qa_bits       Byte array containing qa bit results
! confdnc       Current pixel unobstructed confidence
!
!!Revision History:
! 06/04 Collection 5  R. Frey
! Updated calling argument lists for Nite_snow and LandNite.
! Added 'desert','hi_elev','sfctmp','eco_type','lnd', removed 'plat'.
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-35.
!
!!END
!----------------------------------------------------------------------

      include 'global.inc'

!     scalar arguments
      real confdnc,vza,tbadj,plat,sfctmp,ptwp
      logical snow,ice,coast,desert,hi_elev,sh_lake
      integer nmtests
      byte eco_type
      integer ielem, iline

!     array arguments
      real pxldat(inband) , btclr(7)
      real indat_11um(3,3)
      integer(kind=1) :: is_cold_sfc
      byte testbits(6),qa_bits(10)

!     local scalars
      integer debug,h_output
      logical lnd

!     external subroutines
      external LandNite,Nite_snow,chk_land_nite

!     Common statement for debug purposes
!      common / bug / debug, h_output

!----------------------------------------------------------------------

!     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Using land_nite processing path '',/)')
!      endif

!----------------------------------------------------------------------

      if (snow .or. ice) then
!        snow or ice covered surfaces

         lnd = .true.
         call Nite_snow(pxldat,vza,lnd,testbits,qa_bits,nmtests,confdnc,btclr)
      
      else
!        Standard processing
        
         call LandNite(pxldat,plat,vza,coast,desert,hi_elev,sh_lake,    &
                       sfctmp,eco_type,testbits,qa_bits,nmtests,        &
                       confdnc,ptwp,btclr,is_cold_sfc)
 
      endif

!----------------------------------------------------------------------

!     Debug statement.
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,'' Land nite confidence '',
!     +        f10.2,/)') confdnc
!      endif

!----------------------------------------------------------------------

!     Perform clear-sky restoral tests.

      if( .not. (snow .or. ice)) then
        if(confdnc .le. 0.95) then
          call chk_land_nite(pxldat,tbadj,confdnc,qa_bits,testbits)
        end if
      end if

!      return
end subroutine land_nite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------- end MODULE ---------------------------------

end module land_module
