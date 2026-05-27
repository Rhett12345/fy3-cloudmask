module planck_module

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is a planck code for FY3/MERSI-II product code/
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
!C !end
!C----------------------------------------------------------------------

use names_module
use data_arrays_module
use constant
use numerical
!use platform_module

implicit none


contains
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~ function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function fylat_planck_rad2tbb(rad, band, units) result(bright)

!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute brightness temperature for a fylat infrared band
!
!
!!INPUT parameterS:
!    RAD (real)      Planck radiance (units are determined by UNITS)
!    BAND (LONG)     MODIS IR band number (20-25, 27-36)
!    UNITS (LONG)    Flag defining radiance units
!                    0 => milliWatts per square meter per
!                         steradian per inverse centimeter
!                    1 => Watts per square meter per
!                         steradian per micron
!
!!OUTPUT parameterS:
!    BRIGHT  Brightness temperature (Kelvin)
!                  Note that a value of -1.0 is returned if
!                  RAD .LE. 0.0, or BAND is not in range 20-25, 27-36.
!
!!end
!-----------------------------------------------------------------------

      implicit none

! ... Include files
!      include 'platform_name.inc'
      
! ... Arguments
      real(kind=4)    :: rad
      integer(kind=4) :: band, units

!c ... Local variables
! 1 = modis to fy3 mersi_ii
      real(kind=4) :: cwn_terra(6), tcs_terra(6), tci_terra(6)
      real(kind=4) :: cwn_aqua(6), tcs_aqua(6), tci_aqua(6)
! 2 = viirs to fy3 mersi_ii
      real(kind=4) :: cwn_viirs(6), tcs_viirs(6), tci_viirs(6)
! 21 = real fy3d mersi_ii
      real(kind=4) :: cwn_fy3d(6), tcs_fy3d(6), tci_fy3d(6)      
! 22 = real fy3e mersi_ii
      real(kind=4) :: cwn_fy3e(6), tcs_fy3e(6), tci_fy3e(6)      
      
      real(kind=4) :: cwn, tcs, tci, bright
      integer index

! ... External functions
!      real bright_m, brite_m
!      external bright_m, brite_m
            
! ... Data statements

!-----------------------------------------------------------------------

!     TERRA MODIS DETECTOR-AVERAGED SPECTRAL RESPONSE
!-----------------------------------------------------------------------
!-----
! fylat_sensor_id = 1 / convert modis to mersi II
! ... Effective central wavenumbers (inverse centimeters)
! modis
!c      20v,  21,  22,   23v,
!c      24,   25,  27,   28v,
!c      29v,  30,  31v,  32v,
!c      33,   34,  35,   36,
      data cwn_terra/               &
       2.641767E+03, 2.465422E+03,  &
       1.362741E+03, 1.173198E+03,  &
       9.081998E+02, 8.315149E+02/

! ... Temperature correction slopes (no units)
      data tcs_terra/               &  
       9.993487E-01, 9.998701E-01,  &
       9.994937E-01, 9.995643E-01,  &
       9.995880E-01, 9.997388E-01/

! ... Temperature correction intercepts (Kelvin)
      data tci_terra/               & 
       4.744530E-01, 8.856134E-02,  & 
       2.037728E-01, 1.559624E-01,  &
       1.176660E-01, 6.856633E-02/

!-----------------------------------------------------------------------

!     AQUA MODIS DETECTOR-AVERAGED SPECTRAL RESPONSE
!     (LIAM GUMLEY 2003/06/05)

!     BAND 20 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
!     BAND 21 TEMPERATURE RANGE WAS  180.00 K TO  400.00 K
!     BAND 22 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
!     BAND 23 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
!     BAND 24 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 25 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 27 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 28 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 29 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
!     BAND 30 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
!     BAND 31 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
!     BAND 32 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
!     BAND 33 TEMPERATURE RANGE WAS  180.00 K TO  330.00 K
!     BAND 34 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 35 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K
!     BAND 36 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K
!     BANDS
! modis
!c      20v,  21,  22,   23v,
!c      24,   25,  27,   28v,
!c      29v,  30,  31v,  32v,
!c      33,   34,  35,   36,
! fylat
!       20,  21,  
!       22,  23,
!       24,  25


!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_aqua/                &
       2.647418E+03, 2.462446E+03,  &
       1.361638E+03, 1.169637E+03,  &
       9.076808E+02, 8.308397E+02/

!c ... Temperature correction slopes (no units)
      data tcs_aqua/                &
       9.993438E-01, 9.998729E-01,  &
       9.994894E-01, 9.995439E-01,  &
       9.995483E-01, 9.997404E-01/

!c ... Temperature correction intercepts (Kelvin)
      data tci_aqua/                &
       4.792821E-01, 8.659482E-02,  &
       2.053504E-01, 1.628724E-01,  &
       1.290129E-01, 6.810679E-02/
       
!-----
! fylat_sensor_id = 2 / convert npp/viirs to mersi II
!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_viirs/               &
       2.647418E+03, 2.462446E+03,  &
       1.361638E+03, 1.169637E+03,  &
       9.076808E+02, 8.308397E+02/

!c ... Temperature correction slopes (no units)
      data tcs_viirs/               &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01/

!c ... Temperature correction intercepts (Kelvin)
      data tci_viirs/               &
       0.000000E-00, 0.000000E-00,  &
       0.000000E-00, 0.000000E-00,  &
       0.000000E-00, 0.000000E-00/
       
!-----
! fylat_sensor_id = 21 / real fy3d/mersi II
!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_fy3d/                &
!       2.647418E+03, 2.462446E+03,  &
!       1.361638E+03, 1.169637E+03,  &
!       9.076808E+02, 8.308397E+02/
       2.6434359E+03, 2.471654E+03, & 
        1.382621E+03, 1.168182E+03, &
        9.333640E+02, 8.369410E+02/

!c ... Temperature correction slopes (no units)
      data tcs_fy3d/                &
       0.9992917440, 0.9994814177,  &
       0.9989956900, 0.9997135336,  &
       0.9980397975, 0.9983777125/


!c ... Temperature correction intercepts (Kelvin)
      data tci_fy3d/                  &
       0.50718071650,  0.3493280160,  &
       0.40925130837,  0.1014073981,  &
       0.57633464244,  0.4317181810/      
!-----
! fylat_sensor_id = 22 / real fy3e/mersi II
!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_fy3e/                &
       2.647418E+03, 2.462446E+03,  &
       1.361638E+03, 1.169637E+03,  &
       9.076808E+02, 8.308397E+02/

!c ... Temperature correction slopes (no units)
      data tcs_fy3e/                &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01/

!c ... Temperature correction intercepts (Kelvin)
      data tci_fy3e/                  &
       -0.475900E-00, -0.313900E-00,  &
       -0.266200E-00, -0.051300E-00,  &
       -0.073400E-00,  0.087500E-00/
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

! ... Set default return value
      bright = -1.0

! ... Check input parameters and return if they are bad
      if (rad .le. 0.0 .or.  &
          band .lt. 20 .or.  &
          band .gt. 25 ) return

! ... Get index into coefficient arrays
      index = band - 19

      
! ... Get the coefficients for fylat fy3/mersi_ii
      if (fylat_sensor_id == 1 .or. fylat_sensor_id == 2) then          ! modis to fy3/mersi_ii
         cwn = cwn_aqua(index)
         tcs = tcs_aqua(index)
         tci = tci_aqua(index)
  !      cwn = cwn_terra(index)
  !      tcs = tcs_terra(index)
  !      tci = tci_terra(index)
      else if (fylat_sensor_id == 3) then    ! npp/viirs to fy3/mersi_ii
         cwn = cwn_viirs(index)
         tcs = tcs_viirs(index)
         tci = tci_viirs(index)
      else if (fylat_sensor_id == 21) then   ! fy3d mersi_ii
         cwn = cwn_fy3d(index)
         tcs = tcs_fy3d(index)
         tci = tci_fy3d(index)
      else if (fylat_sensor_id == 22) then   ! fy3e mersi_ii
         cwn = cwn_fy3e(index)
         tcs = tcs_fy3e(index)
         tci = tci_fy3e(index)
      else
         print*,' sensor id is wrong !!!'
   !     call message('fylat_bright.f',          &
   !      'Platform name not recognized ' //     &
   !      '[OPERATOR ACTION: Contact SDST]', 0, 2)
      endif
     
! ... Compute brightness temperature
      if (units .eq. 1) then

! ...   Radiance units are
! ...   Watts per square meter per steradian per micron
        bright = (bright_m(1.0e+4 / cwn, rad) - tci) / tcs

      else

! ...   Radiance units are
! ...   milliWatts per square meter per steradian per wavenumber
        bright = (brite_m(cwn, rad) - tci) / tcs

      endif

end function fylat_planck_rad2tbb
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function bright_m(W, R) result(bright_value)

!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute brightness temperature given monochromatic Planck radiance
!    (Radiance units: Watts per square meter per steradian per micron)
!
!!INPUT parameterS:
!    W (real)           Wavelength (microns)
!    R (real)           Monochromatic Planck radiance (Watts per
!                       square meter per steradian per micron)
!
!!OUTPUT parameterS:
!    BRIGHT_M (real)    Brightness temperature (Kelvin)
!
!!REVISION HISTORY:
!
!!end
!-----------------------------------------------------------------------

      implicit none

! ... Include files
!      include 'fundamental_constants.inc'

! ... Planck constant (Joule second)
      real(kind=8), parameter :: h = 6.62606876d-34

! ... Speed of light in vacuum (meters per second)
      real(kind=8), parameter :: c = 2.99792458d+08

! ... Boltzmann constant (Joules per Kelvin)      
      real(kind=8), parameter :: k = 1.3806503d-23

! ... Derived constants      
      real(kind=8), parameter :: c1 = 2.0d+0 * h * c * c
      real(kind=8), parameter :: c2 = (h * c) / k
      
! ... Arguments
      real(kind=4) :: w, r, bright_value

! ... Local variables
      real(kind=8) :: ws

! ... Set default return value
      bright_value = -1.0
      
! ... Check input parameters and return if they are bad
      if (w .le. 0.0 .or. r .le. 0.0) return
                  
! ... Convert wavelength to meters
      ws = 1.0d-6 * dble(w)
      
! ... Compute brightness temperature
      bright_value = sngl(c2 / (ws * log(c1 / (1.0d+6 * dble(r) * ws**5) + 1.0d+0)))

end function bright_m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function brite_m(V, R) result(brite_value)

!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute brightness temperature given monochromatic Planck radiance
!    (Radiance units: milliWatts per square meter per steradian per
!    inverse centimeter)
!
!!INPUT parameterS:
!    V (real)          Wavenumber (inverse centimeters)
!    R (real)          Monochromatic Planck radiance (milliWatts per
!                      square meter per steradian per
!                      inverse centimeter)
!
!!OUTPUT parameterS:
!    BRITE_M (real)    Brightness temperature (Kelvin)
!
!!REVISION HISTORY:
!
!!end
!-----------------------------------------------------------------------

!      implicit none

! ... Include files
!      include 'fundamental_constants.inc'
! ... Planck constant (Joule second)
      real(kind=8), parameter :: h = 6.62606876d-34

! ... Speed of light in vacuum (meters per second)
      real(kind=8), parameter :: c = 2.99792458d+08

! ... Boltzmann constant (Joules per Kelvin)      
      real(kind=8), parameter :: k = 1.3806503d-23

! ... Derived constants      
      real(kind=8), parameter :: c1 = 2.0d+0 * h * c * c
      real(kind=8), parameter :: c2 = (h * c) / k
      
! ... Arguments
      real(kind=4) :: v, r, brite_value

! ... Local variables
      real(kind=8) :: vs

! ... Set default return value
      brite_value = -1.0
      
! ... Check input parameters and return if they are bad
      if (v .le. 0.0 .or. r .le. 0.0) return
                  
! ... Convert wavenumber to inverse meters
      vs = 1.0d+2 * dble(v)
      
! ... Compute brightness temperature
      brite_value = sngl(c2 * vs / log(c1 * vs**3 / (1.0d-5 * dble(r)) + 1.0d+0))
      
end function brite_m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
!~~~~~~~~~~~~~~~~ function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function fylat_planck_tbb2rad(temp, band, units) result(planck_value)

      real(kind=4)    :: temp
      integer(kind=4) :: band, units

!c ... Local variables
! 1 = modis to fy3 mersi_ii
      real(kind=4) :: cwn_terra(6), tcs_terra(6), tci_terra(6)
      real(kind=4) :: cwn_aqua(6), tcs_aqua(6), tci_aqua(6)
! 2 = viirs to fy3 mersi_ii
      real(kind=4) :: cwn_viirs(6), tcs_viirs(6), tci_viirs(6)
! 21 = real fy3d mersi_ii
      real(kind=4) :: cwn_fy3d(6), tcs_fy3d(6), tci_fy3d(6)      
! 22 = real fy3e mersi_ii
      real(kind=4) :: cwn_fy3e(6), tcs_fy3e(6), tci_fy3e(6)      
      
      real(kind=4) :: cwn, tcs, tci, planck_value
      integer index
      
!c ... External functions
      !real(kind=4) :: planck_m, planc_m
      !external planck_m, planc_m
            
!c ... Data statements

!-----------------------------------------------------------------------
!-----
! fylat_sensor_id = 1 / convert modis to mersi II
! ... Effective central wavenumbers (inverse centimeters)
! modis
!c      20v,  21,  22,   23v,
!c      24,   25,  27,   28v,
!c      29v,  30,  31v,  32v,
!c      33,   34,  35,   36,
      data cwn_terra/               &
       2.641767E+03, 2.465422E+03,  &
       1.362741E+03, 1.173198E+03,  &
       9.081998E+02, 8.315149E+02/

! ... Temperature correction slopes (no units)
      data tcs_terra/               &  
       9.993487E-01, 9.998701E-01,  &
       9.994937E-01, 9.995643E-01,  &
       9.995880E-01, 9.997388E-01/

! ... Temperature correction intercepts (Kelvin)
      data tci_terra/               & 
       4.744530E-01, 8.856134E-02,  & 
       2.037728E-01, 1.559624E-01,  &
       1.176660E-01, 6.856633E-02/


!-----------------------------------------------------------------------

!     AQUA MODIS DETECTOR-AVERAGED SPECTRAL RESPONSE
!     (LIAM GUMLEY 2003/06/05)

!     BAND 20 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K   v 20
!     BAND 21 TEMPERATURE RANGE WAS  180.00 K TO  400.00 K   
!     BAND 22 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
!     BAND 23 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K   v 21 
!     BAND 24 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 25 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 27 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 28 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K   v 22
!     BAND 29 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K   v 23
!     BAND 30 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K  
!     BAND 31 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K   v 24
!     BAND 32 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K   v 25
!     BAND 33 TEMPERATURE RANGE WAS  180.00 K TO  330.00 K
!     BAND 34 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
!     BAND 35 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K
!     BAND 36 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K

!     BANDS
! modis
!c      20v,  21,  22,   23v,
!c      24,   25,  27,   28v,
!c      29v,  30,  31v,  32v,
!c      33,   34,  35,   36,
! fylat
!       20,  21,  
!       22,  23,
!       24,  25


!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_aqua/                &
       2.647418E+03, 2.462446E+03,  &
       1.361638E+03, 1.169637E+03,  &
       9.076808E+02, 8.308397E+02/

!c ... Temperature correction slopes (no units)
      data tcs_aqua/                &
       9.993438E-01, 9.998729E-01,  &
       9.994894E-01, 9.995439E-01,  &
       9.995483E-01, 9.997404E-01/

!c ... Temperature correction intercepts (Kelvin)
      data tci_aqua/                &
       4.792821E-01, 8.659482E-02,  &
       2.053504E-01, 1.628724E-01,  &
       1.290129E-01, 6.810679E-02/
       
!-----
! fylat_sensor_id = 2 / convert npp/viirs to mersi II
!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_viirs/               &
       2.647418E+03, 2.462446E+03,  &
       1.361638E+03, 1.169637E+03,  &
       9.076808E+02, 8.308397E+02/

!c ... Temperature correction slopes (no units)
      data tcs_viirs/               &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01/

!c ... Temperature correction intercepts (Kelvin)
      data tci_viirs/               &
       0.000000E-00, 0.000000E-00,  &
       0.000000E-00, 0.000000E-00,  &
       0.000000E-00, 0.000000E-00/
       
!-----
! fylat_sensor_id = 21 / real fy3d/mersi II
!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_fy3d/                &
!       2.647418E+03, 2.462446E+03,  &
!       1.361638E+03, 1.169637E+03,  &
!       9.076808E+02, 8.308397E+02/
       2.6434359E+03, 2.471654E+03, & 
        1.382621E+03, 1.168182E+03, &
        9.333640E+02, 8.369410E+02/
!2.6434359E+03, 2.471654E+03, 1.382621E+03, 1.168182E+03, 9.33364E+02, 8.36941E+02
!c ... Temperature correction slopes (no units)
      data tcs_fy3d/                &
       0.9992917440, 0.9994814177,  &
       0.9989956900, 0.9997135336,  &
       0.9980397975, 0.9983777125/


!c ... Temperature correction intercepts (Kelvin)
      data tci_fy3d/                  &
       0.50718071650,  0.3493280160,  &
       0.40925130837,  0.1014073981,  &
       0.57633464244,  0.4317181810/    
        
!-----
! fylat_sensor_id = 22 / real fy3e/mersi II
!c ... Effective central wavenumbers (inverse centimeters)
      data cwn_fy3e/                &
       2.647418E+03, 2.462446E+03,  &
       1.361638E+03, 1.169637E+03,  &
       9.076808E+02, 8.308397E+02/

!c ... Temperature correction slopes (no units)
      data tcs_fy3e/                &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01,  &
       10.00000E-01, 10.00000E-01/

!c ... Temperature correction intercepts (Kelvin)
      data tci_fy3e/                  &
       -0.475900E-00, -0.313900E-00,  &
       -0.266200E-00, -0.051300E-00,  &
       -0.073400E-00,  0.087500E-00/
       
!-----------------------------------------------------------------------

! ... Check input parameters and return if they are bad

! ... Get index into coefficient arrays
      index = band - 19

! ... Check input parameters and return if they are bad
      if (temp .le. 0.0 .or.  &
          band .lt. 20 .or.   &
          band .gt. 25 ) return
      
! ... Get the coefficients for fylat fy3/mersi_ii
      if (fylat_sensor_id == 1 .or. fylat_sensor_id == 2) then          ! modis to fy3/mersi_ii
         cwn = cwn_aqua(index)
         tcs = tcs_aqua(index)
         tci = tci_aqua(index)
  !      cwn = cwn_terra(index)
  !      tcs = tcs_terra(index)
  !      tci = tci_terra(index)
      else if (fylat_sensor_id == 3) then    ! npp/viirs to fy3/mersi_ii
         cwn = cwn_viirs(index)
         tcs = tcs_viirs(index)
         tci = tci_viirs(index)
      else if (fylat_sensor_id == 21) then   ! fy3d mersi_ii
         cwn = cwn_fy3d(index)
         tcs = tcs_fy3d(index)
         tci = tci_fy3d(index)
      else if (fylat_sensor_id == 22) then   ! fy3e mersi_ii
         cwn = cwn_fy3e(index)
         tcs = tcs_fy3e(index)
         tci = tci_fy3e(index)
      else
         print*,' sensor id is wrong !!!'
   !     call message('fylat_planck.f',          &
   !      'Platform name not recognized ' //     &
   !      '[OPERATOR ACTION: Contact SDST]', 0, 2)
      endif
                
! ... Compute Planck radiance
      if (units .eq. 1) then

! ...   Radiance units are
! ...   Watts per square meter per steradian per micron
        planck_value = planck_m(1.0e+4 / cwn, temp * tcs + tci)

      else

!c ...   Radiance units are
!c ...   milliWatts per square meter per steradian per wavenumber
        planck_value = planc_m(cwn, temp * tcs + tci)

      endif

end function fylat_planck_tbb2rad
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function planc_m(V, T) result(planc_m_v)
      
!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute monochromatic Planck radiance given brightness temperature
!    (Radiance units: milliWatts per square meter per steradian per
!    inverse centimeter)
!
!!INPUT parameterS:
!    V (real)          Wavenumber (inverse centimeters)
!    T (real)          Brightness temperature (Kelvin)
!
!!OUTPUT parameterS:
!    PLANC_M (real)    Monochromatic Planck radiance (milliWatts per
!                      square meter per steradian per
!                      inverse centimeter)
!
!!REVISION HISTORY:
!
!
!!end
!-----------------------------------------------------------------------

!      implicit none

! ... Include files
!      include 'fundamental_constants.inc'
! ... Planck constant (Joule second)
      real(kind=8), parameter :: h = 6.62606876d-34

! ... Speed of light in vacuum (meters per second)
      real(kind=8), parameter :: c = 2.99792458d+08

! ... Boltzmann constant (Joules per Kelvin)      
      real(kind=8), parameter :: k = 1.3806503d-23

! ... Derived constants      
      real(kind=8), parameter :: c1 = 2.0d+0 * h * c * c
      real(kind=8), parameter :: c2 = (h * c) / k

! ... Arguments
      real(kind=4) :: v, t, planc_m_v

! ... Local variables
      double precision vs

! ... Set default return value
      planc_m_v = -1.0
      
! ... Check input parameters and return if they are bad
      if (v .le. 0.0 .or. t .le. 0.0) return
                  
! ... Convert wavenumber to inverse meters
      vs = 1.0d+2 * dble(v)
      
! ... Compute Planck radiance
      planc_m_v = sngl(1.0d+5 * (c1 * vs**3) / (exp(c2 * vs / dble(t)) - 1.0d+0))
            
end function planc_m

!~~~~~~~~~~~~~~~~ function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function planck_m(W, T) result(planck_m_v)

!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute monochromatic Planck radiance given brightness temperature
!    (Radiance units: Watts per square meter per steradian per micron)
!
!!INPUT parameterS:
!    W (real)           Wavelength (microns)
!    T (real)           Brightness temperature (Kelvin)
!
!!OUTPUT parameterS:
!    PLANCK_M (real)    Monochromatic Planck radiance (Watts per
!                       square meter per steradian per micron)
!
!!REVISION HISTORY:
!
!!end
!-----------------------------------------------------------------------

!      implicit none

! ... Include files
!      include 'fundamental_constants.inc'

! ... Planck constant (Joule second)
      real(kind=8), parameter :: h = 6.62606876d-34

! ... Speed of light in vacuum (meters per second)
      real(kind=8), parameter :: c = 2.99792458d+08

! ... Boltzmann constant (Joules per Kelvin)      
      real(kind=8), parameter :: k = 1.3806503d-23

! ... Derived constants      
      real(kind=8), parameter :: c1 = 2.0d+0 * h * c * c
      real(kind=8), parameter :: c2 = (h * c) / k
            
! ... Arguments
      real(kind=4) :: w, t, planck_m_v

! ... Local variables
      real(kind=8) :: ws

! ... Set default return value
      planck_m_v = -1.0
      
! ... Check input parameters and return if they are bad
      if (w .le. 0.0 .or. t .le. 0.0) return
                  
! ... Convert wavelength to meters
      ws = 1.0d-6 * dble(w)
      
! ... Compute Planck radiance
      planck_m_v = sngl(1.0d-6 * (c1 / ws**5) / (exp(c2 / (ws * dble(t))) - 1.0d+0))

end function  planck_m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ FUNCTION 1: planck_rad_fast ~~~~~~~~~~~~~~~~~~~~~~~
function planck_rad_fast(ichan, bt, dB_dT) RESULT(rad)

!-----------------------------------------------------------------------
! !F90 planck_rad_fast
!
! !Description:
!    This program is to look up planck radiance.
!
! !Input  parameters:
!  
! !Output parameters:
!
!----------------------------------------------------------------------- 

    integer (kind=4), intent(in) :: ichan
    real (kind=4), intent(in) :: bt
    real (kind=4), optional, intent(out) :: dB_dT
    integer :: l
    real (kind=4) :: rad, dB_dT_tmp
    
    l = (bt - planck_min_T)/planck_delta_T
    l = max(1,min(nplanck-1,l))
    
    dB_dT_tmp = (rutil%B_table(l+1,ichan)-rutil%B_table(l,ichan))/(rutil%T_planck(l+1)-rutil%T_planck(l))
    rad = rutil%B_table(l,ichan) + (bt-rutil%T_planck(l))*dB_dT_tmp
    
    if (present(dB_dT)) dB_dT = dB_dT_tmp
    
    return
    
! 3. end FUNCTION
end function planck_rad_fast

!------------------------------------------------------------------
! Subroutine to convert radiance to brightness temperature using a
! look-up table.
!------------------------------------------------------------------

FUNCTION planck_temp_fast(ichan, rad, dB_dT) result(bt)
    INTEGER (kind=4), intent(in) :: ichan
    REAL (kind=4), intent(in) :: rad
    REAL (kind=4), optional, intent(out) :: dB_dT
    INTEGER :: l
    REAL (kind=4) :: bt, dB_dT_tmp
    
    call locate(rutil%B_table(:,ichan),nplanck,rad,l)
    l = max(1,min(nplanck-1,l))
    
    dB_dT_tmp = (rutil%B_table(l+1,ichan)-rutil%B_table(l,ichan))/(rutil%T_planck(l+1)-rutil%T_planck(l))
    bt = rutil%T_planck(l) + (rad-rutil%B_table(l,ichan))/dB_dT_tmp
    
    if (present(dB_dT)) dB_dT = dB_dT_tmp
    
    return
END FUNCTION planck_temp_fast
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 1: planck main program ~~~~~~~~~~~~~~~~~
subroutine planck_main(chflg)

!-----------------------------------------------------------------------
! !F90 planck_main
!
! !Description:
!    This program is to get planck fast look up table.
!
! !Input  parameters
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
integer (kind=1), dimension(:), intent(in) :: chflg

! 2. begin program
  print*,'  ... load planck information  '

  !=== 2.1. load planck constant
  ! call load_calibration(Ialgo)

  !call load_planck_function_constants()

  !=== 2.2. load planck fast table
  call load_fast_planck(chflg)

! 3. end subroutine
end subroutine planck_main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 2: load planck FUNCTION constants ~~~~~
subroutine load_planck_function_constants()

!-----------------------------------------------------------------------
! !F90 planck_main
!
! !Description:
!    This program is to load planck constant.
!
! !Input  parameters
!    none
!
! !Output parameters:
!    none
!
!----------------------------------------------------------------------- 

! 1. define variables
integer(kind=4) :: idet, ichan
real(kind=4)    :: C_1_new_units, C_2_new_units
real(kind=4)    :: C_1_new_units_wlen, C_2_new_units_wlen
real(kind=4)    :: Wlen


! 2. begin program  
  !W to mW in final answer, m-1 to cm-1 in final answer, and convert nu**3 from (cm-1)**3 to (m-1)**3
  C_1_new_units = C_1*1.0e03*1.0e02*1.0e06
  
  !convert nu from cm-1 to m-1
  C_2_new_units = C_2*1.0e02
  
  !m to um in final answer
  C_1_new_units_wlen = C_1*1.0e-06
  
  !convert wlen from um to m
  C_2_new_units_wlen = C_2
  
!  do idet = 1, ndet_max
    do ichan = 20, 19+sat%nir
    
      if (sat%midwnum(ichan) /= missing_value_real4) then
        
         sat%planck_const1(ichan) = C_1_new_units*(sat%midwnum(ichan)**3)
         sat%planck_const2(ichan) = C_2_new_units*(sat%midwnum(ichan))
        
         !In meters
         Wlen = 1.0e-02/sat%midwnum(ichan)
        
         sat%planck_const1_wlen(ichan) = C_1_new_units_wlen/Wlen**5
         sat%planck_const2_wlen(ichan) = C_2_new_units_wlen/Wlen
        
      endif
      
    end do
!  end do

! 3. end subroutine  
end subroutine load_planck_function_constants
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 3: load fast_planck ~~~~~~~~~~~~~~~~~~~
subroutine load_fast_planck(chflg)

!-----------------------------------------------------------------------
! !F90 load_fast_planck
!
! !Description:
!    This program is to load fast planck table.
!
! !Input  parameters:
!    chflg           =
!    
! !Output parameters:
!
!----------------------------------------------------------------------- 

integer (kind=1), dimension(:), intent(in) :: chflg
integer (kind=4) ::  ichan, i , j
    
do i=1, nplanck
   rutil%T_planck(i) = planck_min_T + (i*planck_delta_T)
end do
    
do ichan = 20, 19+sat%nir
   if (chflg(ichan) > 0) then
      do i=1, nplanck
         call planck_rad_func(ichan, rutil%T_planck(i), rutil%B_table(i,ichan))
      end do
   endif
end do
  
! 3. end subroutine  
end subroutine load_fast_planck
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!------------------------------------------------------------------
! subroutine to convert brightness temperature to radiance for Aqua modis
!------------------------------------------------------------------
  
subroutine planck_rad_func(ichan, bt, rad)
    integer(kind=4), intent(in) :: ichan
    real(kind=4), intent(in) :: bt
    real(kind=4), intent(out) :: rad
    
    real(kind=4) :: Planck1
    real(kind=4) :: Planck2
    real(kind=4) :: a
    real(kind=4) :: b

    !a = sat%a(ichan,sat%idet) ! mid_wl (um)
    !b = sat%b(ichan,sat%idet)  
    !Planck1 = 3.7427E8*(1./(PI*a**5))
    !Planck2 = 14388.0

    if (bt > 0) then
      !rad = Planck1/(exp(Planck2/(bt*a + b))-1.0)
	  !rad = sat%planck_const1(ichan)/(exp((sat%planck_const2(ichan))/(sat%a(ichan)*bt+sat%b(ichan)))-1.0)
	  rad = fylat_planck_tbb2rad(bt, ichan, 1)
	 ! print*,'br',bt,rad
    else
      rad = 0.0
    endif
    
end subroutine planck_rad_func
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- end MODULE ---------------------------------
end module planck_module
