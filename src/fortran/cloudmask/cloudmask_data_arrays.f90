module cloudmask_data_arrays

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    fylat fy3/MERSI data arrays for cloudmask algorithm.
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C    none
!C
!C !Author's information
!C    Author: Min Min
!C    E-mail: minmin@cma.gov.cn
!C    Tel   : 86-010-68406763
!C    National Satellite Meteorological Center 
!C  
!C !end
!C----------------------------------------------------------------------

!use data_arrays_module

implicit none

logical :: line_edge,ele_edge,polar,land,day,night,ice,snglnt,visusd,water,bad_value,       &
           coast,desert,vrused,snow,bad_geo,map_ice,map_snow,ndsi_snow,                     &
           hi_elev,antarctic,sh_ocean,sg_bad_data,sh_lake,                                  &
           New_Zealand,Greenland,process,cirrus_ir,cirrus_vis,no_250,                       &
           uniform,shadow,smoke

integer         :: lsf, nmtests, nbands, nbad_1km, nbad_250
real(kind=4)    :: confdnc,precip_water,vza,plat,plon,sfctmp,pmsl,u_wind,v_wind,refang

real(kind=4), dimension(:,:), pointer :: out_pwater             !lyj
real(kind=4), dimension(:,:), pointer :: out_sfctmp             !lyj
integer, dimension(:,:), pointer :: out_polar                   !lyj
integer, dimension(:,:), pointer :: out_day                     !lyj
integer, dimension(:,:), pointer :: out_night                   !lyj
integer, dimension(:,:), pointer :: out_land                    !lyj
integer, dimension(:,:), pointer :: out_water                   !lyj
integer, dimension(:,:), pointer :: out_coast                   !lyj
integer, dimension(:,:), pointer :: out_snglnt                  !lyj
integer, dimension(:,:), pointer :: out_snow                    !lyj
integer, dimension(:,:), pointer :: out_ice                     !lyj
integer, dimension(:,:), pointer :: out_desert                  !lyj
integer, dimension(:,:), pointer :: out_uniform                 !lyj
integer, dimension(:,:), pointer :: out_shadow                  !lyj

byte testbits(6), qa_bits(10)

!,cube_eco(npixel,scans_cube),
     !+     contx_eco(npixel,nlcntx),v250_band(nx,ny,vis_band),
     !+     v1km_band(nx,ny,inband),qa_bits(10),qa_bitarray(npixel,10)
             

end module cloudmask_data_arrays
