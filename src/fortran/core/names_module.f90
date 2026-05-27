module names_module


!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is to define names for FY3/MERSI-II product code/
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
!C    National Satellite Meteorological Center, CMA 
!C  
!C !END
!C----------------------------------------------------------------------

implicit none

!+++++++++ 1. senser option+++++++++++++++++++++ +++++++++++++++++++++++
integer(kind=1)     :: fylat_sensor_id            ! number of sensor [1=fy3d/MERSI_II(modis to mersi_II); 2=fy3d/MERSI_II(real mersi_II); 3= fy3d/MERSI_II(viirs to mersi_II); 21=fy3d/MERSI_II]
integer(kind=1)     :: fylat_nwp_opt              ! option of nwp 
integer(kind=1)     :: fylat_rtm_opt              ! option of rtm [0=none; 1=PFAAST]
character(len=1000) :: code_root_path 
character(len=1000) :: L1b_data_path 
character(len=1000) :: nwp_data_path 
character(len=1000) :: oisst_data_path
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++ 2. file names for i/o +++++++++++++++++++++++++++++++++++++++
character(len=1000) :: fy3_mersi_GEO_data  
character(len=1000) :: fy3_mersi_L1b_data       
character(len=1000) :: fy3_mersi_CLM_data
character(len=1000) :: fy3_mersi_CLA_data
character(len=1000) :: fy3_mersi_CLP_data
character(len=1000) :: fy3_mersi_CTP_data
character(len=1000) :: fy3_mersi_COT_data
character(len=1000) :: fy3_mersi_CON_data
character(len=1000) :: fy3_mersi_SST_data
character(len=1000) :: fy3_intermediate         !lyj
character(len=1000) :: nwp_grib_data1
character(len=1000) :: nwp_grib_data2
character(len=1000) :: oisst_data
integer(kind=1)     :: cloudmask_id,         & ! cloud mask id
                       cloudamount_id,       & ! cloud amount id
                       cloudphase_id,        & ! cloud phase and type id
                       cloudtopz_id,         & ! cloud top height
                       cloudtau_day_id,      & ! cloud optical properties at daytime
                       cloudtau_night_id,    & ! cloud optical properties at nighttime
                       cloudtypeII_id,       & ! cloud type II
                       surface_sst_id          ! sst
integer(kind=1)     :: write_inter_id          ! write intermediate result
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++ 3. algorithm option +++++++++++++++++++++++++++++++++++++++++
type, public :: algorithm_option
   ! 1. cloud
   integer(kind=1) :: cloudmask_index             ! 1.cloud mask
   integer(kind=1) :: cloudamount_index             ! 1.cloud mask
   integer(kind=1) :: cloudphase_index            ! 2.cloud phase
   integer(kind=1) :: cloudtopz_index             ! 3.cloud top height
   integer(kind=1) :: cloudtau_day_index          ! 4.cloud optical and microphysical [daytime]
   integer(kind=1) :: cloudtau_night_index        ! 5.cloud optical and microphysical [nighttime]  
   integer(kind=1) :: cloudtypeII_index           ! 6.cloud type II
   integer(kind=1) :: surface_sst_index           ! 7.sea surface temperature 
end type 
type(algorithm_option), public :: fylat_alg_opt   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++ 4. ancillary data +++++++++++++++++++++++++++++++++++++++++++
character(len=14),parameter :: gogel_name='goge1_2_img.v1'                      ! 1. 1KM ECOSYSTEM 
character(len=34),parameter :: IGBP_name ='IGBP.EcoMap.NtoS.2004.149.v004.hdf'  ! 2. IGBP ECOSYSTEM   
character(len=19),parameter :: ecosystem_name  ='fylat_ecosystem.hdf'           ! 3. Ecosystem 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module names_module
