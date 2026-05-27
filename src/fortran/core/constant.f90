module constant

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is to define contants in the fylat.
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C
!C !Author's information
!C    Author: Min Min
!C    E-mail: minmin@cma.gov.cn
!C    Tel   : 86-010-68406763
!C    National Satellite Meteorological Center 
!C  
!C !End
!C----------------------------------------------------------------------

implicit none

!+++++++++++++++++++ 1. public variables ++++++++++++++++++++++++++++++++
integer(kind=1), parameter, public  :: success = 0, warning = 1, error = 2, failure =3
integer(kind=1) :: errorlevel

integer(kind=1), parameter, public :: SUCCEED = 0, FAIL = -1
!type, public :: symbol_struct
!
!    integer(kind=1) :: NO                    ! 0
!    integer(kind=1) :: YES                   ! 1
!
!    integer(kind=1) :: NO_SNOW               ! 1
!    integer(kind=1) :: SEA_ICE               ! 2
!    integer(kind=1) :: SNOW                  ! 3   
!    
!end type symbol_struct
!type(symbol_struct), public, parameter :: sym = &
!    symbol_struct(0,1, &                             !no/yes
!                  1,2,3)  	  	                     !snow and ice mask
                  
!+++++++++++++++++++ 2. missing value +++++++++++++++++++++++++++++++++++
real(kind=8),    parameter, public      :: missing_value_real8 = -999.0
real(kind=4),    parameter, public      :: missing_value_real4 = -999.0
integer(kind=1), parameter, public      :: missing_value_int1 = -128
integer(kind=2), parameter, public      :: missing_value_int2 = -32768
integer(kind=4), parameter, public      :: missing_value_int4 = -999
!real(kind=4),    parameter, public      :: bad_data = -999.0

!+++++++++++++++++++ 3. symbol struct +++++++++++++++++++++++++++++++++++
type, public :: symbol_struct
    
    integer(kind=1) :: NO_SCALE              ! 0
    integer(kind=1) :: LINEAR_SCALE          ! 1
    integer(kind=1) :: LOG_SCALE             ! 2
    integer(kind=1) :: SQRT_SCALE            ! 3
    
    integer(kind=1) :: CLOUDY                ! 0
    integer(kind=1) :: PROB_CLOUDY           ! 1
    integer(kind=1) :: PROB_CLEAR            ! 2
    integer(kind=1) :: CLEAR                 ! 3
    
    integer(kind=1) :: CLEAR_TYPE            ! 0
    integer(kind=1) :: FOG_TYPE              ! 1
    integer(kind=1) :: WATER_TYPE            ! 2
    integer(kind=1) :: SUPERCOOLED_TYPE      ! 3
    integer(kind=1) :: MIXED_TYPE            ! 4
    integer(kind=1) :: TICE_TYPE             ! 5
    integer(kind=1) :: CIRRUS_TYPE           ! 6
    integer(kind=1) :: OVERLAP_TYPE          ! 7
    integer(kind=1) :: OVERSHOOTING_TYPE     ! 8
    integer(kind=1) :: UNKNOWN_TYPE          ! 9
    
    integer(kind=1) :: CLEAR_PHASE           ! 0
    integer(kind=1) :: WATER_PHASE           ! 1
    integer(kind=1) :: SUPERCOOLED_PHASE     ! 2
    integer(kind=1) :: MIXED_PHASE           ! 3
    integer(kind=1) :: ICE_PHASE             ! 4
    integer(kind=1) :: UNKNOWN_PHASE         ! 5
    
    integer(kind=1) :: NO_AEROSOL_AERO       ! 0
    integer(kind=1) :: MOSTLY_ASH_AERO       ! 1
    integer(kind=1) :: ASH_ICE_AERO          ! 2
    integer(kind=1) :: ASH_WATER_AERO        ! 3
    integer(kind=1) :: DUST_AERO             ! 4
    integer(kind=1) :: SMOKE_AERO            ! 5
    integer(kind=1) :: SULFATE_AERO          ! 6
    integer(kind=1) :: UNKNOWN_AERO          ! 7
    
    integer(kind=1) :: NO_SPACE              ! 0
    integer(kind=1) :: SPACE                 ! 1
    
    integer(kind=1) :: NO                    ! 0
    integer(kind=1) :: YES                   ! 1
    
    integer(kind=1) :: WATER_SFC             ! 0
    integer(kind=1) :: EVERGREEN_NEEDLE_SFC  ! 1
    integer(kind=1) :: EVERGREEN_BROAD_SFC   ! 2
    integer(kind=1) :: DECIDUOUS_NEEDLE_SFC  ! 3
    integer(kind=1) :: DECIDUOUS_BROAD_SFC   ! 4
    integer(kind=1) :: MIXED_FORESTS_SFC     ! 5
    integer(kind=1) :: WOODLANDS_SFC         ! 6
    integer(kind=1) :: WOODED_GRASS_SFC      ! 7
    integer(kind=1) :: CLOSED_SHRUBS_SFC     ! 8
    integer(kind=1) :: OPEN_SHRUBS_SFC       ! 9
    integer(kind=1) :: GRASSES_SFC           ! 10
    integer(kind=1) :: CROPLANDS_SFC         ! 11
    integer(kind=1) :: BARE_SFC              ! 12
    integer(kind=1) :: URBAN_SFC             ! 13
    
    integer(kind=1) :: NO_DESERT             ! 0
    integer(kind=1) :: NIR_DESERT            ! 1
    integer(kind=1) :: BRIGHT_DESERT         ! 2
    
    integer(kind=1) :: SHALLOW_OCEAN         ! 0
    integer(kind=1) :: LAND                  ! 1
    integer(kind=1) :: COASTLINE             ! 2
    integer(kind=1) :: SHALLOW_INLAND_WATER  ! 3
    integer(kind=1) :: EPHEMERAL_WATER       ! 4
    integer(kind=1) :: DEEP_INLAND_WATER     ! 5
    integer(kind=1) :: MODERATE_OCEAN        ! 6
    integer(kind=1) :: DEEP_OCEAN            ! 7
    
    integer(kind=1) :: NO_VOLCANO            ! 0
    integer(kind=1) :: VERY_CLOSE_VOLCANO    ! 1
    integer(kind=1) :: CLOSE_VOLCANO         ! 2
    
    integer(kind=1) :: NO_COAST              ! 0
    integer(kind=1) :: COAST_1KM             ! 1
    integer(kind=1) :: COAST_2KM             ! 2
    integer(kind=1) :: COAST_3KM             ! 3
    integer(kind=1) :: COAST_4KM             ! 4
    integer(kind=1) :: COAST_5KM             ! 5
    integer(kind=1) :: COAST_6KM             ! 6
    integer(kind=1) :: COAST_7KM             ! 7
    integer(kind=1) :: COAST_8KM             ! 8
    integer(kind=1) :: COAST_9KM             ! 9
    integer(kind=1) :: COAST_10KM            ! 10
    
    integer(kind=1) :: NO_SNOW               ! 1
    integer(kind=1) :: SEA_ICE               ! 2
    integer(kind=1) :: SNOW                  ! 3
    
    integer(kind=4) :: SUCCESS               ! 0
    integer(kind=4) :: FAILURE               ! 1
    integer(kind=4) :: INFORMATION           ! 2
    integer(kind=4) :: WARNING               ! 3
    integer(kind=4) :: EOF                   ! 4
    integer(kind=4) :: UNDEFINED             ! 5
    integer(kind=4) :: EXISTS                ! 6
    integer(kind=4) :: EXIT                  ! 7
    
    integer(kind=1) :: SNOW_NOT_AVAILABLE    ! 1
    integer(kind=1) :: NWP_SNOW              ! 2
    integer(kind=1) :: IMS_SNOW              ! 3
    
    integer(kind=1) :: CONSTANT_SFC_EMISS    ! 1
    integer(kind=1) :: TABLE_SFC_EMISS       ! 2
    integer(kind=1) :: SEEBOR_SFC_EMISS      ! 3
    
    integer(kind=1) :: CONSTANT_SFC_ALB      ! 1
    integer(kind=1) :: TABLE_SFC_ALB         ! 2
    integer(kind=1) :: MODIS_SFC_ALB         ! 3
    
    integer(kind=4) :: LITTLE_ENDIAN         ! 0
    integer(kind=4) :: BIG_ENDIAN            ! 1
    
    integer(kind=1) :: WATER_GEN             ! 0
    integer(kind=1) :: COAST_GEN             ! 1
    integer(kind=1) :: LAND_GEN              ! 2
    integer(kind=1) :: DESERT_GEN            ! 3
    integer(kind=1) :: SNOW_GEN              ! 4
    
    INTEGER(kind=1) :: QC_GOOD			    ! 0
    INTEGER(kind=1) :: QC_CYCLE_VZA          ! 1
    INTEGER(kind=1) :: QC_CYCLE_SZA          ! 2
    INTEGER(kind=1) :: QC_CYCLE_NOCLOUD      ! 3
    INTEGER(kind=1) :: QC_CYCLE_CLOUDTYPE    ! 4
    INTEGER(kind=1) :: QC_CYCLE_TCLOUD       ! 5
    INTEGER(kind=1) :: QC_MINERR_WATER_0     ! 6  no retrieval for water phase
    INTEGER(kind=1) :: QC_MINERR_ICE_0       ! 7  no retrieval for ice   phase
        
END TYPE symbol_struct
  
TYPE(symbol_struct), public, parameter :: sym = &
     symbol_struct(0,1,2,3, &                         !scaling method
                   0,1,2,3, &                         !cloud mask
                   0,1,2,3,4,5,6,7,8,9, &             !cloud type
                   0,1,2,3,4,5, &                     !cloud phase
                   0,1,2,3,4,5,6,7, &                 !aerosol mask
                   0,1, &                             !space    
                   0,1, &                             !no/yes
                   0,1,2,3,4,5,6,7,8,9,10,11,12,13, & !surface type
                   0,1,2, &                           !desert mask
                   0,1,2,3,4,5,6,7, &                 !land/water mask
                   0,1,2, &                           !volcano mask
                   0,1,2,3,4,5,6,7,8,9,10, &          !coast mask
                   1,2,3, &                           !snow mask
                   0,1,2,3,4,5,6,7, &                 !system flags
                   1,2,3, &                           !snow mask type
                   1,2,3, &                           !surface emissivity type
                   1,2,3, &                           !surface albedo type
                   0,1, &                             !endian config
                   0,1,2,3,4,&                        !generic surface type for output
                   0,1,2,3,4,5,6,7)  	  	          !quality flags for nighttime cldopt


!+++++++++++++++++++ 4. fundmental value ++++++++++++++++++++++++++++++++
  real(kind=8), parameter, private :: ONE = 1.0
  real(kind=8), parameter, private :: TWO = 2.0
  
  !#----------------------------------------------------------------------------#
  !#                       -- mpav - Earth Radius in meters  --                 #
  !#----------------------------------------------------------------------------#
  real(kind=8), parameter, public :: R_EARTH             = 6378206.4

  !#----------------------------------------------------------------------------#
  !#                -- IRRATIONAL NUMBERS AND ASSOCIATED BITS --                #
  !#----------------------------------------------------------------------------#

  ! PI
  real(kind=8), parameter, public :: PI             = 3.141592653589793238462643
  !mpav - added DTOR
  real(kind=8), parameter, public :: DTOR           = PI/180.0
  real(kind=8), parameter, public :: RTOD           = 180.0/PI

  ! E
  real(kind=8), parameter, public :: E              = 2.71828182845904523560287
  real(kind=8), parameter, public :: E_RECIPROCAL   = 0.367879441171442321595523
  real(kind=8), parameter, public :: E_SQUARED      = 7.389056098930650227230427
  real(kind=8), parameter, public :: E_LOG10        = 0.434294481903251827651129
  
!+++++++++++++++++++ 4. univeral constants ++++++++++++++++++++++++++++++
  ! ----------------------------------------------
  ! Speed of light
  ! Symbol:c,  Units:m/s,  Rel.Uncert.(ppm): exact
  ! ----------------------------------------------
  real(kind=8), parameter, public :: SPEED_OF_LIGHT = 2.99792458e+08

  ! --------------------------------------------------
  ! Permeability of vacuum
  ! Symbol:mu0,  Units:N/A^2,  Rel.Uncert.(ppm): exact
  ! --------------------------------------------------
  real(kind=8), parameter, public :: PERMEABILITY = PI * 4.0e-07

  ! -----------------------------------------------------
  ! Permittivity of vacuum
  ! Symbol:epsilon0,  Units:F/m,  Rel.Uncert.(ppm): exact
  ! -----------------------------------------------------
  real(kind=8), parameter, public :: PERMITTIVITY =                ONE                  / &
  !                                             ------------------------------------
                                                ( PERMEABILITY * SPEED_OF_LIGHT**2 )

  ! ---------------------------------------------
  ! Planck constant
  ! Symbol:h,  Units:Js,  Rel.Uncert.(ppm): 0.078
  ! ---------------------------------------------
  real(kind=8), parameter, public :: PLANCK_CONSTANT = 6.62606876e-34

  ! ----------------------------------------------------
  ! Gravitational constant
  ! Symbol:G,  Units:m^3/kg/s^2,  Rel.Uncert.(ppm): 1500
  ! ----------------------------------------------------
  real(kind=8), parameter, public :: GRAVITATIONAL_CONSTANT = 6.673e-11



  !#----------------------------------------------------------------------------#
  !#                          -- CONVERSION FACTORS --                          #
  !#----------------------------------------------------------------------------#

  ! ---------------------------------------------
  ! Electron volt
  ! Symbol:eV,  Units:J,  Rel.Uncert.(ppm): 0.039
  ! ---------------------------------------------
  real(kind=8), parameter, public :: ELECTRON_VOLT = 1.602176462e-19

  ! ---------------------------------------------
  ! Unified atomic mass unit
  ! Symbol:u,  Units:kg,  Rel.Uncert.(ppm): 0.079
  ! ---------------------------------------------
  real(kind=8), parameter, public :: UNIFIED_ATOMIC_MASS_UNIT = 1.66053873e-27

  ! ----------------------------------------------
  ! Standard atmosphere
  ! Symbol:P0,  Units:Pa,  Rel.Uncert.(ppm): exact
  ! ----------------------------------------------
  real(kind=8), parameter, public :: STANDARD_ATMOSPHERE = 101325.0

  ! ----------------------------------------------------------------------
  ! Standard temperature
  ! Symbol:T0,  Units:Kelvin,  Rel.Uncert.(ppm): exact
  !
  ! Note that the unit of thermodynamic temperature, the Kelvin, is the
  ! fraction 1/273.16 of the thermodynamic temperature of the triple point
  ! of water. The standard temperature is the ice point of water, NOT the
  ! triple point, hence the 0.01K difference.
  ! ----------------------------------------------------------------------
  real(kind=8), parameter, public :: STANDARD_TEMPERATURE = 273.15

  ! ------------------------------------------------
  ! Standard gravity
  ! Symbol:g,  Units:m/s^2,  Rel.Uncert.(ppm): exact
  ! ------------------------------------------------
  real(kind=8), parameter, public :: STANDARD_GRAVITY = 9.80665



  !#----------------------------------------------------------------------------#
  !#                        -- PHYSICOCHEMICAL CONSTANTS --                     #
  !#----------------------------------------------------------------------------#

  ! -----------------------------------------------------
  ! Avogadro constant
  ! Symbol:N(A),  Units:mole^-1,  Rel.Uncert.(ppm): 0.079
  ! -----------------------------------------------------
  real(kind=8), parameter, public :: AVOGADRO_CONSTANT = 6.02214199e+23


  ! -------------------------------------------------
  ! Molar gas constant
  ! Symbol:R,  Units:J/mole/K,  Rel.Uncert.(ppm): 1.7
  ! -------------------------------------------------
  real(kind=8), parameter, public :: MOLAR_GAS_CONSTANT = 8.314472

  ! --------------------------------------------
  ! Boltzmann constant
  ! Symbol:k,  Units:J/K,  Rel.Uncert.(ppm): 1.7
  !
  !         R
  !   k = ------
  !        N(A)
  !
  !     = 1.3806503(24)e-23
  !
  ! --------------------------------------------
  real(kind=8), parameter, public :: BOLTZMANN_CONSTANT = MOLAR_GAS_CONSTANT / &
  !                                                   ------------------
                                                       AVOGADRO_CONSTANT

  ! ------------------------------------------------------
  ! Stefan-Boltzmann constant
  ! Symbol:sigma,  Units:W/m^2/K^4,  Rel.Uncert.(ppm): 7.0
  !
  !             PI^2
  !             ----.k^4
  !              60                     h
  !   sigma = ------------   ( hbar = ----- )
  !            hbar^3.c^2              2PI
  !
  !         = 5.670400(40)e-08
  !
  ! I just placed the value here due to the mathematical
  ! gymnastics required to calculate it directly.
  ! ------------------------------------------------------
  real(kind=8), parameter, public :: STEFAN_BOLTZMANN_CONSTANT = 5.670400e-08

  ! -------------------------------------------------------
  ! First Planck function constant
  ! Symbol:c1,  Units:W.m^2.sr^-1,  Rel.Uncert.(ppm): 0.078
  !
  !   c1 = 2.h.c^2
  !
  !      = 1.191042722(93)e-16
  !
  ! -------------------------------------------------------
  real(kind=8), parameter, public :: C_1 = TWO * PLANCK_CONSTANT * SPEED_OF_LIGHT**2

  ! ---------------------------------------------
  ! Second Planck function constant
  ! Symbol:c2,  Units:K.m,  Rel.Uncert.(ppm): 1.7
  !
  !         h.c
  !   c2 = -----
  !          k
  !
  !      = 1.4387752(25)e-02
  !
  ! ---------------------------------------------
  real(kind=8), parameter, public :: C_2 = PLANCK_CONSTANT * SPEED_OF_LIGHT / &
  !                                    ----------------------------------
                                               BOLTZMANN_CONSTANT

  ! -----------------------------------------------------------------
  ! Molar volume of an ideal gas at standard temperature and pressure
  ! Symbol:Vm,  Units:m^3/mol,  Rel.Uncert.(ppm): 1.7
  !
  !         R.T0
  !   Vm = ------
  !          P0
  !
  !      = 2.2413996(39)e-02
  !
  ! -----------------------------------------------------------------
  real(kind=8), parameter, public :: STP_MOLAR_VOLUME = ( MOLAR_GAS_CONSTANT * STANDARD_TEMPERATURE ) / &
  !                                                 ---------------------------------------------
                                                                 STANDARD_ATMOSPHERE

  ! ------------------------------------------------------------------
  ! Loschmidt constant: The number density of one mole of an ideal gas
  ! at standard temperature and pressure
  ! Symbol:n0,  Units:m^-3,  Rel.Uncert.(ppm): 1.7
  !
  !         N(A).P0
  !   n0 = ---------
  !          R.T0
  !
  !         N(A)
  !      = ------     .....(1)
  !          Vm
  !
  !      = 2.6867775(47)e+25
  !
  ! Alternatively, using the ideal gas law directly, we know,
  !
  !   P.V = n.k.T     .....(2)
  !
  ! For V = 1m^3 (unit volume), and P = P0, T = T0, then eqn.(2)
  ! becomes,
  !
  !   P0 = n0.k.T0
  !
  ! which rearranges to
  !
  !          P0  
  !   n0 = ------     .....(3)
  !         k.T0 
  !
  ! Equation (1) rather than eqn(3) is used here.
  ! ------------------------------------------------------------------
  real(kind=8), parameter, public :: LOSCHMIDT_CONSTANT = AVOGADRO_CONSTANT / &
  !                                                   -----------------
                                                      STP_MOLAR_VOLUME
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module constant
