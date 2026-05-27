module data_arrays_module

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    fylat fy3/MERSI data arrays.
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

implicit none

integer(kind=4), public :: i, j, k

!---------- 1. fylat mersi data public variables -------------------------
!$$$$$ 1.1. two nearest nwp DATA times
type, public :: nwpt
   character(LEN=4),dimension(2)  :: year
   character(LEN=2),dimension(2)  :: month
   character(LEN=2),dimension(2)  :: day
   character(LEN=2),dimension(2)  :: hour
end type nwpt
type(nwpt), public :: nwptime

!$$$$$ 1.2. julian times of satellite and nwp DATA
type, public :: julian_time
   real(kind=8)  :: sate
   real(kind=8)  :: nwp1
   real(kind=8)  :: nwp2
end type julian_time
type(julian_time), public :: jutime

integer(kind=1), dimension(:,:), pointer      :: cube_eco

!$$$$$ 1.3 surface type 
!logical polar,day,night,land,water,snglnt,snow,ice,uniform,shadow,coast,desert
!-------------------------------------------------------------------------


!---------- 2. fylat fy3 mersi data variables ----------------------------
type, public :: fylat_fy3_mersi_geo
    real(kind=4), dimension(:,:), pointer      :: lon 	
    real(kind=4), dimension(:,:), pointer      :: lat 
    real(kind=4), dimension(:,:), pointer      :: SolarZenith 	
    real(kind=4), dimension(:,:), pointer      :: SolarAzimuth 
    real(kind=4), dimension(:,:), pointer      :: SensorZenith
    real(kind=4), dimension(:,:), pointer      :: SensorAzimuth 
    real(kind=4), dimension(:,:), pointer      :: RelAzimuth
    real(kind=4), dimension(:,:), pointer      :: GlintAngle  	
    real(kind=4), dimension(:,:), pointer      :: dem  	
    integer(kind=1), dimension(:,:), pointer   :: lsm   !landseamask
    integer(kind=1), dimension(:,:,:), pointer :: flag  !flag
    real(kind=4)                               :: sun_earth_distance ! distance between sun to earth
    real(kind=4), dimension(:,:), pointer      :: Cos_Satzen         ! cos[satellite zenith angle]   
    real(kind=4), dimension(:,:), pointer      :: Cos_Solzen         ! cos[solar zenith angle]   
    real(kind=4), dimension(:,:), pointer      :: Scatzen            ! scattering zenith angle
end type fylat_fy3_mersi_geo
type(fylat_fy3_mersi_geo), public  :: geo


type, public :: fylat_fy3_mersi_L1b
    integer(kind=4)                           :: year			! # year
    integer(kind=4)                           :: month			! # month
    integer(kind=4)                           :: day			! # day
    integer(kind=4)                           :: hour			! # hour
    integer(kind=4)                           :: mint			! # minute
    integer(kind=4)                           :: nday			! # number of day in one year
    integer(kind=4)                           :: nLine			! # lines 
    integer(kind=4)                           :: nElem			! # pixels 
    integer(kind=4)                           :: nChan			! # fylat/MERSI-II channels
    integer(kind=4)                           :: nvis			! # visible channel number
    integer(kind=4)                           :: nir   		    ! # ir channel number
    real(kind=4), dimension(4,19)             :: vis_cal_coef
    real(kind=4), dimension(4,6)              :: ir_cal_coef
    real(kind=4), dimension(:,:,:), pointer   :: ref_vis            ! %
    real(kind=4), dimension(:,:,:), pointer   :: rad_ir             ! mW/(m2 sr cm)
    real(kind=4), dimension(:,:,:), pointer   :: tbb_ir             ! K
    integer(kind=1), dimension(30)            :: chan_flag		    ! # chan_flag
    real(kind=4), dimension(30)               :: midwave		    ! # chan_mid wavelength
    real(kind=4), dimension(30)               :: midwnum		    ! # chan_mid wavenumber
    real(kind=4), dimension(30)               :: planck_const1		! # planck_constant1
    real(kind=4), dimension(30)               :: planck_const2		! # planck_constant2
    real(kind=4), dimension(30)               :: planck_const1_wlen ! # planck_constant1_wlen
    real(kind=4), dimension(30)               :: planck_const2_wlen	! # planck_constant2_wlen
    real(kind=4), dimension(30)               :: solar_const		! # solar constant
    real(kind=4), dimension(30)               :: a                  ! # planck_calibration_coefficient1
    real(kind=4), dimension(30)               :: b	                ! # planck_calibration_coefficient2
    integer(kind=1), dimension(:,:), pointer  :: snow_mask   	! # snow mask 
    integer(kind=1), dimension(:,:), pointer  :: eco          	! # ecosystem map
    real(kind=4),    dimension(:,:), pointer  :: sst            ! sate y pos in nwp data 
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss38    ! surface emissivity at channel 3.8 um
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss40    ! surface emissivity at channel 4.0 um
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss73    ! surface emissivity at channel 7.3 um
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss86    ! surface emissivity at channel 8.6 um
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss11    ! surface emissivity at channel 11.0 um
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss12    ! surface emissivity at channel 12.0 um
    real(kind=4),    dimension(:,:), pointer  :: sfc_emiss13    ! surface emissivity at channel 13.3 um
    real(kind=4),    dimension(:,:), pointer  :: ws_albedo66    ! surface albedo at channel 0.66
    real(kind=4),    dimension(:,:), pointer  :: ws_albedo87    ! surface albedo at channel 0.86
    real(kind=4),    dimension(:,:), pointer  :: ws_albedo124   ! surface albedo at channel 1.24
    real(kind=4),    dimension(:,:), pointer  :: ws_albedo164   ! surface albedo at channel 1.64
    real(kind=4),    dimension(:,:), pointer  :: ws_albedo213   ! surface albedo at channel 2.13
    real(kind=4),    dimension(:,:), pointer  :: rad_clr38      ! toa clear radiance at channel 7 
    real(kind=4),    dimension(:,:), pointer  :: rad_clr40      ! toa clear radiance at channel 8
    real(kind=4),    dimension(:,:), pointer  :: rad_clr73      ! toa clear radiance at channel 9
    real(kind=4),    dimension(:,:), pointer  :: rad_clr86      ! toa clear radiance at channel 10
    real(kind=4),    dimension(:,:), pointer  :: rad_clr11      ! toa clear radiance at channel 11
    real(kind=4),    dimension(:,:), pointer  :: rad_clr12      ! toa clear radiance at channel 12
    real(kind=4),    dimension(:,:), pointer  :: rad_clr13      ! toa clear radiance at channel 13
    real(kind=4),    dimension(:,:), pointer  :: bt_clr38       ! toa clear bt at channel 7   [unit:K]
    real(kind=4),    dimension(:,:), pointer  :: bt_clr40       ! toa clear bt at channel 8   [unit:K]
    real(kind=4),    dimension(:,:), pointer  :: bt_clr73       ! toa clear bt at channel 9   [unit:K]
    real(kind=4),    dimension(:,:), pointer  :: bt_clr86       ! toa clear bt at channel 10  [unit:K]
    real(kind=4),    dimension(:,:), pointer  :: bt_clr11       ! toa clear bt at channel 11  [unit:K]
    real(kind=4),    dimension(:,:), pointer  :: bt_clr12       ! toa clear bt at channel 12  [unit:K]
    real(kind=4),    dimension(:,:), pointer  :: bt_clr13       ! toa clear bt at channel 13  [unit:K]
    integer(kind=4), dimension(:,:), pointer  :: x_nwp          ! sate x pos in nwp data
    integer(kind=4), dimension(:,:), pointer  :: y_nwp          ! sate y pos in nwp data  
    integer(kind=4), dimension(:,:), pointer  :: ivza           ! view zenith angle associated with clear sky calculations
    integer(kind=4), dimension(:,:), pointer  :: isfc           ! index of layer nearest to sfc 
    real(kind=4),    dimension(:,:), pointer  :: zsfc
end type fylat_fy3_mersi_L1b
type(fylat_fy3_mersi_L1b), public  :: sat
!--------------------------------------------------------------------------


!---------- 3. fylat nwp data variables -----------------------------------
type, public :: nwpdata ! original nwp data [ncep reanalysis data]
  real(kind=4), dimension(:,:), pointer    :: lon
  real(kind=4), dimension(:,:), pointer    :: lat
  real(kind=4), dimension(:,:), pointer    :: plev_nointerp
  real(kind=4), dimension(:,:,:), pointer  :: psfc
  real(kind=4), dimension(:,:,:), pointer  :: pmsl
  real(kind=4), dimension(:,:,:), pointer  :: tsfc
  real(kind=4), dimension(:,:,:), pointer  :: zsfc
  real(kind=4), dimension(:,:,:), pointer  :: albedo
  real(kind=4), dimension(:,:,:), pointer  :: t_sigma
  real(kind=4), dimension(:,:,:), pointer  :: rh_sigma
  real(kind=4), dimension(:,:,:), pointer  :: u_sigma
  real(kind=4), dimension(:,:,:), pointer  :: v_sigma
  real(kind=4), dimension(:,:,:), pointer  :: tpw
  real(kind=4), dimension(:,:,:), pointer  :: weasd
  real(kind=4), dimension(:,:,:), pointer  :: o3col
  real(kind=4), dimension(:,:,:), pointer  :: ttropo
  real(kind=4), dimension(:,:,:,:), pointer:: tlev
  real(kind=4), dimension(:,:,:,:), pointer:: zlev
  real(kind=4), dimension(:,:,:,:), pointer:: o3lev
  real(kind=4), dimension(:,:,:,:), pointer:: rhlev
  real(kind=4), dimension(:,:,:,:), pointer:: clwlev
  real(kind=4), dimension(:,:,:,:), pointer:: ulev
  real(kind=4), dimension(:,:,:,:), pointer:: vlev
end type nwpdata
type(nwpdata), save, public :: nwpo ! the neighbouring two nwp data


type, public :: nwpinterp   ! the interpolated nwp data
  integer(kind=4)                        :: nlon = 360
  integer(kind=4)                        :: nlat = 181
  integer(kind=4)                        :: nlon05 = 720
  integer(kind=4)                        :: nlat05 = 361
  integer(kind=4)                        :: nlevels = 26
  integer(kind=4)                        :: nlevels_rh = 21
  real(kind=4), dimension(:,:), pointer  :: lon
  real(kind=4), dimension(:,:), pointer  :: lat
  real(kind=4), dimension(:,:), pointer  :: psfc
  real(kind=4), dimension(:,:), pointer  :: pmsl
  real(kind=4), dimension(:,:), pointer  :: tsfc
  real(kind=4), dimension(:,:), pointer  :: zsfc
  real(kind=4), dimension(:,:), pointer  :: albedo
  real(kind=4), dimension(:,:), pointer  :: t_sigma
  real(kind=4), dimension(:,:), pointer  :: rh_sigma
  real(kind=4), dimension(:,:), pointer  :: u_sigma
  real(kind=4), dimension(:,:), pointer  :: v_sigma
  real(kind=4), dimension(:,:), pointer  :: tpw
  real(kind=4), dimension(:,:), pointer  :: weasd
  real(kind=4), dimension(:,:), pointer  :: o3col
  real(kind=4), dimension(:,:), pointer  :: ttropo
  real(kind=4), dimension(:,:,:), pointer:: plev
  real(kind=4), dimension(:,:,:), pointer:: tlev
  real(kind=4), dimension(:,:,:), pointer:: zlev
  real(kind=4), dimension(:,:,:), pointer:: o3lev
  real(kind=4), dimension(:,:,:), pointer:: rhlev
  real(kind=4), dimension(:,:,:), pointer:: wlev
  real(kind=4), dimension(:,:,:), pointer:: clwlev
  real(kind=4), dimension(:,:,:), pointer:: ulev
  real(kind=4), dimension(:,:,:), pointer:: vlev
end type nwpinterp
type(nwpinterp), save, public :: nwp26   ! the interpolated nwp data at satellite time

TYPE, PUBLIC :: nwpinterp2   ! the interpolated nwp DATA
  INTEGER(KIND=4)                           :: nlon = 360
  INTEGER(KIND=4)                           :: nlat = 181
  INTEGER(KIND=4)                           :: nlon05 = 720
  INTEGER(KIND=4)                           :: nlat05 = 361
  INTEGER(KIND=4)                           :: nlon25 = 1440
  INTEGER(KIND=4)                           :: nlat25 = 721
  INTEGER(KIND=4)                           :: nlevels = 31
  INTEGER(KIND=4)                           :: nlevels_rh = 31
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: lon
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: lat
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: psfc
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: pmsl
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: tsfc
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: zsfc
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: albedo
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: t_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: rh_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: u_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: v_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: tpw
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: weasd
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: o3col
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: ttropo
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: plev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: tlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: zlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: o3lev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: rhlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: wlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: clwlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: ulev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: vlev
END TYPE nwpinterp2
TYPE(nwpinterp2), SAVE, PUBLIC :: nwp31   ! the interpolated nwp DATA at satellite time

TYPE, PUBLIC :: nwpinterp3   ! the interpolated nwp DATA
  INTEGER(KIND=4)                           :: nlon = 360
  INTEGER(KIND=4)                           :: nlat = 181
  INTEGER(KIND=4)                           :: nlon05 = 720
  INTEGER(KIND=4)                           :: nlat05 = 361
  INTEGER(KIND=4)                           :: nlon25 = 1440
  INTEGER(KIND=4)                           :: nlat25 = 721
  INTEGER(KIND=4)                           :: nlevels = 41
  INTEGER(KIND=4)                           :: nlevels_rh = 41
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: lon
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: lat
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: psfc
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: pmsl
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: tsfc
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: zsfc
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: albedo
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: t_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: rh_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: u_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: v_sigma
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: tpw
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: weasd
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: o3col
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)   :: ttropo
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: plev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: tlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: zlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: o3lev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: rhlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: wlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: clwlev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: ulev
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: vlev
END TYPE nwpinterp3
TYPE(nwpinterp3), SAVE, PUBLIC :: nwp41   ! the interpolated nwp DATA at satellite time


type, public :: nwpinterp_T639   ! the interpolated nwp data
  integer(kind=4)                        :: nlon = 1280
  integer(kind=4)                        :: nlat = 641
  integer(kind=4)                        :: nlevels = 36
  integer(kind=4)                        :: nlevels_rh = 36
  real(kind=4), dimension(:,:), pointer  :: lon
  real(kind=4), dimension(:,:), pointer  :: lat
  real(kind=4), dimension(:,:), pointer  :: psfc
  real(kind=4), dimension(:,:), pointer  :: pmsl
  real(kind=4), dimension(:,:), pointer  :: tsfc
  real(kind=4), dimension(:,:), pointer  :: zsfc
  real(kind=4), dimension(:,:), pointer  :: albedo
  real(kind=4), dimension(:,:), pointer  :: t_sigma
  real(kind=4), dimension(:,:), pointer  :: rh_sigma
  real(kind=4), dimension(:,:), pointer  :: u_sigma
  real(kind=4), dimension(:,:), pointer  :: v_sigma
  real(kind=4), dimension(:,:), pointer  :: tpw
  real(kind=4), dimension(:,:), pointer  :: weasd
  real(kind=4), dimension(:,:), pointer  :: o3col
  real(kind=4), dimension(:,:), pointer  :: ttropo
  real(kind=4), dimension(:,:,:), pointer:: plev
  real(kind=4), dimension(:,:,:), pointer:: tlev
  real(kind=4), dimension(:,:,:), pointer:: zlev
  real(kind=4), dimension(:,:,:), pointer:: o3lev
  real(kind=4), dimension(:,:,:), pointer:: rhlev
  real(kind=4), dimension(:,:,:), pointer:: wlev
  real(kind=4), dimension(:,:,:), pointer:: clwlev
  real(kind=4), dimension(:,:,:), pointer:: ulev
  real(kind=4), dimension(:,:,:), pointer:: vlev
end type nwpinterp_T639
type(nwpinterp_T639), save, public :: nwp36   ! the interpolated nwp data at satellite time


TYPE, PUBLIC :: nwpinterp_grapes_gfs   ! the interpolated nwp DATA
  INTEGER(kind=4)                           :: nlon = 1440
  INTEGER(kind=4)                           :: nlat = 720
  INTEGER(kind=4)                           :: nlevels = 40
  INTEGER(kind=4)                           :: nlevels_rh = 30
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: lon
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: lat
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: psfc
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: pmsl
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: tsfc
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: zsfc
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: albedo
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: t_sigma
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: rh_sigma
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: u_sigma
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: v_sigma
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: tpw
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: weasd
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: o3col
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: ttropo
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: plev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: tlev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: zlev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: o3lev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: rhlev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: wlev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: clwlev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: ulev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:,:):: vlev
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: u10m
  REAL(kind=4),ALLOCATABLE,DIMENSION(:,:)  :: v10m
END TYPE nwpinterp_grapes_gfs
TYPE(nwpinterp_grapes_gfs), SAVE, PUBLIC :: nwp40   ! the interpolated nwp DATA at satellite time

type, public :: nwpdat 
  real(kind=4)                            :: a
  real(kind=4)                            :: lon
  real(kind=4)                            :: lat
  real(kind=4)                            :: psfc
  real(kind=4)                            :: pmsl
  real(kind=4)                            :: tsfc
  real(kind=4)                            :: tsfc_uni
  real(kind=4)                            :: zsfc
  real(kind=4)                            :: albedo
  real(kind=4)                            :: t_sigma
  real(kind=4)                            :: rh_sigma
  real(kind=4)                            :: u_sigma
  real(kind=4)                            :: v_sigma
  real(kind=4)                            :: tpw
  real(kind=4)                            :: weasd
  real(kind=4)                            :: o3col
  real(kind=4)                            :: ttropo
  real(kind=4), dimension(:), pointer     :: plev 
  real(kind=4), dimension(:), pointer     :: tlev 
  real(kind=4), dimension(:), pointer     :: zlev 
  real(kind=4), dimension(:), pointer     :: o3lev 
  real(kind=4), dimension(:), pointer     :: wlev 
! real(kind=4), dimension(:), pointer     :: clwlev 
  real(kind=4), dimension(:), pointer     :: rhlev 
  real(kind=4), dimension(:), pointer     :: tpwlev
  real(kind=4), dimension(:), pointer     :: ulev 
  real(kind=4), dimension(:), pointer     :: vlev 
  integer(kind=4), dimension(:), pointer  :: inversion_lev
  integer(kind=4)                         :: sfc_level
  integer(kind=4)                         :: tropo_level
  integer(kind=4)                         :: strato_level
  integer(kind=4)                         :: ninversion
end type nwpdat

type, public :: nwpprofile101      ! the interpolated nwp data to 101 layers
  integer(kind=4)                         :: rtm_nvzen
  integer(kind=4)                         :: nlon = 360
  integer(kind=4)                         :: nlat = 181
  integer(kind=4)                         :: nlon05 = 720
  integer(kind=4)                         :: nlat05 = 361
  integer(kind=4)                         :: nlon_T639 = 1280
  integer(kind=4)                         :: nlat_T639 = 641
  integer(kind=4)                         :: n2lon25 = 1440
  integer(kind=4)                         :: n2lat25 = 720
  integer(kind=4)                         :: nlevels
  real(kind=4)                            :: first_lat2       != 90.0 for grib1 -90.0 for grib2
  real(kind=4)                            :: first_lat = 90.0 != 90.0 for grib1 -90.0 for grib2
  real(kind=4)                            :: first_lon = -180.0
  real(kind=4)                            :: dlat = 1.0
  real(kind=4)                            :: dlon = 1.0
  real(kind=4)                            :: dlat05 = 0.5
  real(kind=4)                            :: dlon05 = 0.5
  real(kind=4)                            :: dlat_T639 = 0.2815
  real(kind=4)                            :: dlon_T639 = 0.2823
  real(kind=4)                            :: dlat0p25 = 0.25
  real(kind=4)                            :: dlon0p25 = 0.25
  REAL(kind=4)                            :: dlat25 != 0.25
  REAL(kind=4)                            :: dlon25 != 0.25
  
  INTEGER(kind=4)                         :: num_lon != 360
  INTEGER(kind=4)                         :: num_lat != 181
  REAL(kind=4)                            :: first_lon25 != 360
  REAL(kind=4)                            :: first_lat25 != 181
  INTEGER(kind=4)                         :: nlon25 = 1440
  INTEGER(kind=4)                         :: nlat25 = 720
  INTEGER(kind=4)                         :: nlon25_gfs = 1440
  INTEGER(kind=4)                         :: nlat25_gfs = 721
  INTEGER(kind=4)                         :: nlon25_ne != 721
  INTEGER(kind=4)                         :: nlat25_ne != 360


  type(nwpdat), dimension(:,:), pointer   :: dat
end type nwpprofile101
type(nwpprofile101), save, public :: nwp   ! the interpolated nwp data to 101 layers
!--------------------------------------------------------------------------

!---------- 4. fylat product output arrays --------------------------------
!+++ cloud mask 
integer(kind=1), parameter :: cm_byte_dim = 6
integer(kind=1), parameter :: cm_qa_dim   = 10
byte, dimension(:,:,:), pointer ::  cm_bitarray, cm_qa_bitarray  ! cloud mask
integer(kind=1), dimension(:,:,:), pointer :: cm_tmp  ! 0-3 

!+++ cloud amount
integer(kind=4)                          :: ix_5km, iy_5km
integer(kind=1), dimension(:,:), pointer :: cloud_amount     ! cloud amount
integer(kind=1), dimension(:,:), pointer :: cloud_amount_qa  ! cloud amount quality 
integer(kind=4), dimension(:,:), pointer :: lat_5km     ! 
integer(kind=4), dimension(:,:), pointer :: lon_5km

!+++ cloud phase and type
type, public :: CLP_out      ! the interpolated nwp data to 101 layers
  !-A local variable used to point at the cldtype output variable
  integer(kind=1), dimension(:,:), pointer :: Cldtype
  integer(kind=1), dimension(:,:), pointer :: Cldtype_Tmpy 
  
  !-A local variable used to point at the cldphase output variable
  integer(kind=1), dimension(:,:), pointer :: Cldphase 
  
  !-A local variable used to point at the cldphase QF variable
  integer(kind=1), dimension(:,:,:), pointer :: Cldphase_Qf 
  
  !-A local variable used to point at the cldphase Qpi variable
  integer(kind=1), dimension(:,:,:), pointer :: Cldphase_Qpi 
  
  !-Needed for gradient filter output
  integer(kind=4), dimension(:,:), pointer :: Xgrad_Emiss14 
  
  !-Needed for gradient filter output
  integer(kind=4), dimension(:,:), pointer :: Ygrad_Emiss14
  
  !-Needed for gradient filter output
  integer(kind=4), dimension(:,:), pointer :: Num_Steps_Gradient
  
  !-Needed for gradient filter input
  integer(kind=1), dimension(:,:), pointer :: Cldphase_Lrc_Mask
  
  !-Ch 10 (7.3 um) cloud emissivity array for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Emiss_Chn10_Tot
  
  !-Ch 10 (7.3 um) cloud emissivity array for single layered clouds at LRC
  real(kind=4), dimension(:,:), pointer :: Emiss_Chn10_Tot_Lrc
  
  !-Ch 14 (11 um) cloud emissivity array for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Emiss_Chn14_Tot 
  
  !-Ch 14 (11 um) cloud emissivity array for multi-layered clouds
  real(kind=4), dimension(:,:), pointer :: Emiss_Chn14_Tot_Multi 
  
  !-Beta (7.3/11 um) for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta7311_Tot 
  
  !-Beta (7.3/11 um) for multi-layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta7311_Tot_Multi 
  
  !-Beta (12/11 um) for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta1112_Tot 
  
  !-Beta (12/11 um) for single layered clouds at LRC
  real(kind=4), dimension(:,:), pointer :: Cldbeta1112_Tot_Lrc 
  
  !-Beta (12/11 um) for multi-layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta1112_Tot_Multi 
  
  !-Beta (12/11 um) opaque assumption for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta1112_Opaque 
  
  !-Beta (12/11 um) opaque assumption for multi-layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta1112_Opaque_Multi 
  
  !-Beta (8.5/11 um) for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta8511_Tot 
  
  !-Beta (8.5/11 um) for multi-layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta8511_Tot_Multi 
  
  !-Beta (8.5/11 um) opaque assumption for single layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta8511_Opaque 
  
  !-Beta (8.5/11 um) opaque assumption for single layered clouds at LRC
  real(kind=4), dimension(:,:), pointer :: Cldbeta8511_Opaque_Lrc 
  
  !-Beta (8.5/11 um) opaque assumption for multi-layered clouds
  real(kind=4), dimension(:,:), pointer :: Cldbeta8511_Opaque_Multi 
  
  !-The temperature an optically thick cloud would have assuming a near opaque 7.3 um emissivity
  real(kind=4), dimension(:,:), pointer :: Opaque_Cld_Temp_Chn10 
  
  !-The temperature an optically thick cloud would have assuming a near opaque 11 um emissivity
  real(kind=4), dimension(:,:), pointer :: Opaque_Cld_Temp_Chn14 
end type CLP_out
type(CLP_out), save, public :: clp ! the interpolated nwp data to 101 layers

!+++ cloud top height properties
type, public :: CTP_out
  real(kind=4), dimension(:,:), pointer :: cldt            !cloud temperature (K)
  real(kind=4), dimension(:,:), pointer :: cldemiss        !cloud emissivity at 11 microns
  real(kind=4), dimension(:,:), pointer :: cldp            !cloud pressure (hPa)
  real(kind=4), dimension(:,:), pointer :: cldz            !cloud height (km)
!  real(kind=4), dimension(:,:), pointer :: cod_ir          !cloud IR optical depth
  real(kind=4), dimension(:,:), pointer :: cod_vis         !cloud visible optical depth
  real(kind=4), dimension(:,:), pointer :: cldbeta1112     !beta value for 11 and 12 microns
!  real(kind=4), dimension(:,:), pointer :: reff            !effective particle radius
!  real(kind=4), dimension(:,:), pointer :: ash_loading     !ash loading in ton/m^2
!  integer(kind=1), dimension(:,:,:), pointer :: qf         !retrieval quality flags
end type CTP_out
type(CTP_out), save, public :: ctp

!+++ cloud microphysical and optical properties at daytime
type, public :: COT_out
  real(kind=4), dimension(:,:), pointer :: cod_vis               !cloud temperature (K)
  real(kind=4), dimension(:,:), pointer :: cldreff               !cloud emissivity at 11 microns
  real(kind=4), dimension(:,:), pointer :: cldlwp                !cloud pressure (hPa)
  real(kind=4), dimension(:,:), pointer :: cldiwp                !cloud height (km)
  integer(kind=1), dimension(:,:,:), pointer :: qcflg_cotd         !retrieval quality flags
end type COT_out
type(COT_out), save, public :: cot

!+++ sea surface temperature
type, public :: sfc_out
  real(kind=4), dimension(:,:), pointer :: sst                    !sst(K)
  integer(kind=1), dimension(:,:), pointer :: qcflg_sst         !retrieval quality flags
end type sfc_out
type(sfc_out), save, public :: sfc
!--------------------------------------------------------------------------

!---------- 5. fylat ir rtm arrays ----------------------------------------
real(kind=4), public, parameter :: RTM_VZA_BINSIZE = 0.01

! 5.1 RTM structure definition
type, public :: rtm_prof
  integer(kind=1) :: flag = 0
  real(kind=4)                        :: satzen
  real(kind=4)                        :: bt_clr38, bt_clr40, bt_clr73, bt_clr86, bt_clr11, bt_clr12
  real(kind=4)                        :: rad_clr38, rad_clr40, rad_clr73, rad_clr86, rad_clr11, rad_clr12
  real(kind=4), dimension(:), pointer :: rtm_util
  real(kind=4), dimension(:), pointer :: trans_atm_clr38, trans_atm_clr40,   &
                                         trans_atm_clr73, trans_atm_clr86,   &
                                         trans_atm_clr11, trans_atm_clr12
  real(kind=4), dimension(:), pointer :: rad_atm_clr38, rad_atm_clr40,   &
                                         rad_atm_clr73, rad_atm_clr86,   &
                                         rad_atm_clr11, rad_atm_clr12
  real(kind=4), dimension(:), pointer :: cloud_prof38,  cloud_prof40,  &
                                         cloud_prof73,  cloud_prof86,  &
                                         cloud_prof11, cloud_prof12
end type rtm_prof

! 5.2 rtm parameters defination  
type, public :: rtm_params
  integer(kind=1) :: flag = 0
  type (rtm_prof), dimension(:), allocatable :: d
end type rtm_params

! 5.3 Declare a 2D pointer of "rtm_params" structures
type (rtm_params), public, dimension(:, :), allocatable :: rtm
!----------------------------------------------------------------------------------

!------------------- 6. planck fast table -----------------------------------------
!integer(kind=4), parameter, public :: ir_nchan_max = 6 ! for fy3c/mersi-ii
real(kind=4), parameter, public    :: planck_max_T = 360.0
real(kind=4), parameter, public    :: planck_min_T = 159.0
real(kind=4), parameter, public    :: planck_delta_T = 1.0
integer(kind=4), parameter, public :: nplanck = (planck_max_T - planck_min_T)/planck_delta_T
  
type, public :: planck_table
  real(kind=4), dimension(nplanck)    :: T_planck
  real(kind=4), dimension(nplanck,30) :: B_table
end type planck_table
  
type (planck_table), save, public :: rutil
!----------------------------------------------------------------------------------

!---------- 7. 101 layers US standard profile arrays ----------------------
real(kind=4), dimension(101) :: zstd, pstd, tstd, wstd, ostd 
! zstd=height(km), pstd=pressure(hPa), tstd=temperature(K), 
! wstd=water vapor quality mix ratio(g/kg),dimension ostd=o3 mix ratio(ppmv)

data zstd     / 84.3537,  77.4029,  71.9723,  67.3967,  63.4045,  &
      59.8409,  56.5871,  53.6288,  50.8840,  48.3846,  46.0977,  & 
      44.0340,  42.1639,  40.4546,  38.8843,  37.4368,  36.0908,  & 
      34.8416,  33.6707,  32.5763,  31.5408,  30.5617,  29.6312,  & 
      28.7396,  27.8896,  27.0751,  26.2924,  25.5438,  24.8256,  & 
      24.1320,  23.4634,  22.8192,  22.1970,  21.5955,  21.0137,  & 
      20.4505,  19.9045,  19.3741,  18.8584,  18.3564,  17.8676,  & 
      17.3911,  16.9267,  16.4733,  16.0310,  15.5992,  15.1777,  & 
      14.7658,  14.3631,  13.9693,  13.5839,  13.2067,  12.8373,  & 
      12.4756,  12.1215,  11.7746,  11.4346,  11.1006,  10.7715,  & 
      10.4464,  10.1250,   9.8064,   9.4904,   9.1771,   8.8667,  & 
       8.5591,   8.2544,   7.9524,   7.6530,   7.3563,   7.0622,  & 
       6.7706,   6.4815,   6.1949,   5.9108,   5.6292,   5.3499,  & 
       5.0731,   4.7985,   4.5263,   4.2563,   3.9887,   3.7232,  & 
       3.4600,   3.1989,   2.9401,   2.6833,   2.4286,   2.1761,  & 
       1.9257,   1.6772,   1.4309,   1.1865,   0.9441,   0.7036,  & 
       0.4649,   0.2277,   0.0000,  -0.2201,  -0.4373,  -0.6517/ 

data pstd     / 0.0050,    0.0161,    0.0384,    0.0769,    0.1370,   &
     0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,   &
     1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,   &
     5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,   &
    14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,   &
    29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,   &
    51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,   &
    83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,   &
   125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,   &
   180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,   &
   247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,   &
   328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,   &
   424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,   &
   535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,   &
   661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,   &
   802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,   &
   958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/

data tstd                     / 190.19, 203.65, 215.30, 226.87, 237.83, &     
        247.50, 256.03, 263.48, 267.09, 270.37, 266.42, 261.56, 256.40, &
        251.69, 247.32, 243.27, 239.56, 236.07, 232.76, 230.67, 228.71, &     
        227.35, 226.29, 225.28, 224.41, 223.61, 222.85, 222.12, 221.42, &
        220.73, 220.07, 219.44, 218.82, 218.23, 217.65, 217.18, 216.91, &
        216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, & 
        216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.71, &
        216.71, 216.72, 216.81, 217.80, 218.77, 219.72, 220.66, 222.51, &
        224.57, 226.59, 228.58, 230.61, 232.61, 234.57, 236.53, 238.48, &
        240.40, 242.31, 244.21, 246.09, 247.94, 249.78, 251.62, 253.45, &
        255.26, 257.04, 258.80, 260.55, 262.28, 264.02, 265.73, 267.42, &
        269.09, 270.77, 272.43, 274.06, 275.70, 277.32, 278.92, 280.51, &
        282.08, 283.64, 285.20, 286.74, 288.25, 289.75, 291.22, 292.68/ 

data wstd                    /  0.001,  0.001,  0.002,  0.003,  0.003,  &   
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  &   
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  &   
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  &   
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  &   
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  &   
        0.003,  0.003,  0.004,  0.004,  0.005,  0.005,  0.007,  0.009,  &   
        0.011,  0.012,  0.014,  0.020,  0.025,  0.030,  0.035,  0.047,  &   
        0.061,  0.075,  0.089,  0.126,  0.162,  0.197,  0.235,  0.273,  &   
        0.310,  0.356,  0.410,  0.471,  0.535,  0.601,  0.684,  0.784,  &   
        0.886,  0.987,  1.094,  1.225,  1.353,  1.519,  1.686,  1.852,  &   
        2.036,  2.267,  2.496,  2.721,  2.947,  3.170,  3.391,  3.621,  &   
        3.848,  4.084,  4.333,  4.579,  4.822,  5.061,  5.296,  5.528/ 

data ostd                  /  0.47330,0.27695,0.28678,0.51816,0.83229,  & 
      1.18466,1.69647,2.16633,3.00338,3.76287,4.75054,5.61330,6.33914,  & 
      7.03675,7.50525,7.75612,7.81607,7.69626,7.56605,7.28440,7.01002,  & 
      6.72722,6.44629,6.17714,5.92914,5.69481,5.47387,5.26813,5.01252,  & 
      4.68941,4.35141,4.01425,3.68771,3.37116,3.06407,2.77294,2.50321,  & 
      2.24098,1.98592,1.74840,1.54451,1.34582,1.17824,1.02513,0.89358,  & 
      0.78844,0.69683,0.62654,0.55781,0.50380,0.45515,0.42037,0.38632,  & 
      0.35297,0.32029,0.28832,0.25756,0.22739,0.19780,0.16877,0.14901,  & 
      0.13190,0.11511,0.09861,0.08818,0.07793,0.06786,0.06146,0.05768,  & 
      0.05396,0.05071,0.04803,0.04548,0.04301,0.04081,0.03983,0.03883,  & 
      0.03783,0.03685,0.03588,0.03491,0.03395,0.03368,0.03349,0.03331,  & 
      0.03313,0.03292,0.03271,0.03251,0.03190,0.03126,0.03062,0.02990,  & 
      0.02918,0.02850,0.02785,0.02721,0.02658,0.02596,0.02579,0.02579/ 
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
end module data_arrays_module

