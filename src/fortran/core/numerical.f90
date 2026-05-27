module numerical

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is numerical module.
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
!C !END
!C----------------------------------------------------------------------

!use names_module
use constant

INTERFACE PACK_BYTES
     MODULE  PROCEDURE    &
         PACK_BYTES_I1,  &
         PACK_BYTES_I2
END INTERFACE 

INTERFACE median_filter
   MODULE  PROCEDURE      &
     median_filter_int8,  &
     median_filter_float32
END INTERFACE 

INTERFACE median_filter9
   MODULE  PROCEDURE       &
     median_filter9_int8,  &
     median_filter9_float32
END INTERFACE  

INTERFACE Get_Irregular_Lut_Index
  MODULE PROCEDURE &
    Get_Irregular_Lut_Index_Int32,   &
    Get_Irregular_Lut_Index_Float32, &
    Get_Irregular_Lut_Index_Float64
END INTERFACE  

INTERFACE Locate
  MODULE PROCEDURE  &
    Locate_Int32,   &
    Locate_Float32, &
    Locate_Float64
END INTERFACE      

CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~~~~ subroutine 1: ICNVRT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ICNVRT(WAY,NUM,STRING,LENGTH,IERR)

!-----------------------------------------------------------------------
! !F90 ICNVRT
!
! !Description:
!        This SUBROUTINE does an INTEGER-to-CHARACTER conversion
!        or a characater-to-INTEGER conversion depENDing on the
!        INTEGER WAY:
!                IF WAY = 0 THEN an INTEGER-to-CHARACTER conversion
!                is done. IF WAY .NE. 0 THEN a CHARACTER-to-INTEGER
!                conversion is done.
!
! !USAGE:
!U
!U        CALL ICNVRT(WAY,NUM,STRING)
!U             where WAY, NUM, STRING, and LENGTH are defined below.
!U
!U        Example: CALL ICNVRT(0,1000,STRING,LENGTH)
!U                 on RETURN STRING = '1000' and
!U                 LENGTH = 4.
!         
! !Input  parameters:
!    WAY - INTEGER::; Determines which way the conversion goes:
!              IF WAY = 0 THEN an INTEGER-to-CHARACTER conversion
!                         is performed;
!              IF WAY.NE.0 THEN a CHARACTER-to-INTEGER conversion
!                         is performed.
!
!    NUM - INTEGER::; an input only IF WAY = 0. NUM is the INTEGER
!                number to be converted to a CHARACTER expression.
!
!    STRING - CHARACTER; an input only IF WAY .NE. 0. STRING
!                is the CHARACTER expression to be converted to an
!                INTEGER value. It contain no decimal points or 
!                non-numeric characters other than possibly a
!                sign. IF STRING contains  a '+' sign, it will be
!                stripped of it on RETURN.
!
! !Output parameters:
!    NUM - INTEGER::; contains the INTEGER:: representation of 
!                STRING.
!
!    STRING - CHARACTER; contains the CHARACTER representation of NUM.
!
!    LENGTH - INTEGER::; The length of STRING to the first blank.
!                  The signIFicant part of STRING can be accessed with
!                  the declaration STRING(1:LENGTH).
!
!    IERR - INTEGER:: variable giving RETURN condition:
!                IERR = 0 for normal RETURN;
!                IERR = 1 IF NUM cannot be converted to STRING because
!                       STRING is too short or STRING cannot be
!                       converted to NUM because STRING is too long.
!                IERR = 2 IF STRING contained a non-numeric CHARACTER
!                       other than a leading sign or something went
!                       wrong with an INTEGER-to-CHARACTER conversion.
!
! !Other
!       ALGORITHM:
!A
!A         Nothing noteworthy, except that this SUBROUTINE will work
!A          for strange CHARACTER sets where the CHARACTER '1' doesn't
!A          follow '0', etc.
!A
!       MACHINE DEPENDENCIES: CM
!M          The parameter MAXINT (below) should be set to the
!M          number of digits that an INTEGER:: data type can have
!M          not including leading signs. For VAX FORTRAN V4.4-177
!M          MAXINT = 10.
!M
!M          NOTE: Under VAX FORTRAN V4.4-177, the
!M          error condition IERR = 1 will never occur for an
!M          INTEGER-to-CHARACTER conversion IF STRING
!M          is ALLOCATEd at least 11 bytes (CHARACTER*11).
!M
!       HISTORY:
!H
!H      written by:             bobby bodenheimer
!H      date:                   september 1986
!H      current version:        1.0
!H      modIFications:          none
!H
!       ROUTINES CALLED:
!C
!C          NONE.
!C
!----------------------------------------------------------------------
!       written for:    The CASCADE Project
!                       Oak Ridge National Laboratory
!                       U.S. Department of Energy
!                       contract number DE-AC05-840R21400
!                       subcontract number 37B-7685 S13
!                       organization:  The University of Tennessee
!----------------------------------------------------------------------
!       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
!       NO RESTRICTIONS ON ITS USE ARE IMPLIED
!----------------------------------------------------------------------


! Global Variables.
!
 INTEGER(KIND=4), INTENT(in):: WAY
 INTEGER(KIND=4), INTENT(out)::  LENGTH, IERR
 INTEGER(KIND=4), INTENT(inout):: NUM
 CHARACTER(LEN=*), INTENT(inout):: STRING
!
!
! Local Variables
!
      INTEGER(KIND=4)::       I
      INTEGER(KIND=4)::       MNUM
      INTEGER(KIND=4)::       M
      logical::       NEG
!
      INTEGER, parameter::MAXINT=10
!
      NEG = .FALSE.
      IERR = 0
!
!  INTEGER-to-CHARACTER conversion.
!
      IF (WAY == 0) THEN
         STRING = " "
         IF (NUM < 0) THEN
            NEG = .TRUE.
            MNUM = -NUM
            LENGTH = INT(LOG10(REAL(MNUM))) + 1
         ELSE IF (NUM == 0) THEN
            MNUM = NUM
            LENGTH = 1
         ELSE
            MNUM = NUM
            LENGTH = INT(LOG10(REAL(MNUM))) + 1
         END IF
         IF (LENGTH > LEN(STRING)) THEN
            IERR = 1
            RETURN
         END IF
ten:     DO I=LENGTH,1,-1    
            M=INT(REAL(MNUM)/10**(I-1))
            IF (M == 0) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "0"
            ELSE IF (M == 1) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "1"
            ELSE IF (M == 2) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "2"
            ELSE IF (M == 3) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "3"
            ELSE IF (M == 4) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "4"
            ELSE IF (M == 5) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "5"
            ELSE IF (M == 6) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "6"
            ELSE IF (M == 7) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "7"
            ELSE IF (M == 8) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "8"
            ELSE IF (M == 9) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "9"
            ELSE
               IERR = 2
               RETURN
            END IF
            MNUM = MNUM - M*10**(I-1)
         END DO ten

         IF (NEG .eqv. .true.) THEN
            STRING = "-"//STRING
            LENGTH = LENGTH + 1
         END IF
!
!  CHARACTER-to-INTEGER conversion.
!
      ELSE
         IF (STRING(1:1) == "-") THEN
            NEG = .TRUE.
            STRING = STRING(2:LEN(STRING))
         END IF
         IF (STRING(1:1) == "+") STRING = STRING(2:LEN(STRING))
         NUM = 0
         LENGTH = INDEX(STRING," ") - 1
         IF (LENGTH > MAXINT) THEN
            IERR = 1
            RETURN
         END IF
twenty:  DO I=LENGTH,1,-1
            IF (STRING(LENGTH-I+1:LENGTH-I+1) == "0") THEN
               M = 0
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "1") THEN
               M = 1
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "2") THEN
               M = 2
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "3") THEN
               M = 3
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "4") THEN
               M = 4
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "5") THEN
               M = 5
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "6") THEN
               M = 6
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "7") THEN
               M = 7
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "8") THEN
               M = 8
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "9") THEN
               M = 9
            ELSE
               IERR = 2
               RETURN
            END IF
            NUM = NUM + INT(10**(I-1))*M
         END DO twenty

         IF (NEG .eqv. .true.) THEN
            NUM = -NUM
            STRING = '-'//STRING
            LENGTH = LENGTH + 1
         END IF
      END IF
!
!  Last lines of ICNVRT
!
   RETURN
   
end subroutine ICNVRT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ FUNCTION 2: determine leap year ~~~~~~~~~~~~~~~~~~~
FUNCTION leap_year_fct(yyyy) RESULT(leap_flg)

!-----------------------------------------------------------------------
! !F90 leap_year_fct
!
! !Description:
!    This FUNCTION determine leap year
!
! !Input  parameters:
!    yyyy             = INTEGER year   
!   
! !Output parameters:
!    leap_flg         = leap year flag  [1 is leap year; 0 not]
!
!-----------------------------------------------------------------------

  INTEGER(KIND=4), INTENT(in) :: yyyy
  INTEGER(KIND=4) :: leap_flg
  
  leap_flg = 0
  
  IF ((modulo(yyyy,4) == 0 .and. modulo(yyyy,100) /= 0) .or. &
     modulo(yyyy,400) == 0) leap_flg = 1
  
  RETURN
  
END FUNCTION leap_year_fct
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ FUNCTION 3: compute day number in one year ~~~~~
FUNCTION compute_daynum(jmonth, jday, ileap) RESULT(day)

!-----------------------------------------------------------------------
! !F90 compute_daynum
!
! !Description:
!    This FUNCTION computes month.
!
! !Input  parameters:
!    jmonth            = INTEGER month  
!    jday              = INTEGER day   
!    ileap             = leap year flag
!
! !Output parameters:
!    day              = INTEGER day
!
!-----------------------------------------------------------------------

!  arguments
   INTEGER(KIND=4), INTENT(in) :: ileap
   INTEGER(KIND=4), INTENT(in) :: jmonth
   INTEGER(KIND=4), INTENT(in) :: jday
   INTEGER(KIND=4) :: day

!  Local variables

   IF (jmonth == 1) THEN
      day = jday
   ELSEIF (jmonth == 2) THEN
      day = 31 + jday
   ELSEIF (jmonth == 3) THEN
      day = 59 + jday + ileap
   ELSEIF (jmonth == 4) THEN
      day = 90 + jday + ileap
   ELSEIF (jmonth == 5) THEN
      day = 120 + jday + ileap
   ELSEIF (jmonth == 6) THEN
      day = 151 + jday + ileap
   ELSEIF (jmonth == 7) THEN
      day = 181 + jday + ileap
   ELSEIF (jmonth == 8) THEN
      day = 212 + jday + ileap
   ELSEIF (jmonth == 9) THEN
      day = 243 + jday + ileap
   ELSEIF (jmonth == 10) THEN
      day = 273 + jday + ileap
   ELSEIF (jmonth == 11) THEN
      day = 304 + jday + ileap
   ELSE
      day = 334 + jday + ileap
   ENDIF

END FUNCTION compute_daynum
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 1: Julian data calculation ~~~~~~~~~~~~~
SUBROUTINE  julian (IY, IM, ID, IH, MIT, JD)

!-----------------------------------------------------------------------
! !F90 julian
!
! !Description:
!    The main program computes julian day (1-365/366)
!
! !Input  parameters:
!    IY               = INTEGER year      
!    IM               = INTEGER month
!    ID               = INTEGER day
!    IH               = INTEGER hour
!    MIT              = INTEGER minute
!
!
! !Output parameters:
!    JD               = julian day
!
!-----------------------------------------------------------------------

! 1. define variables
INTEGER(kind=4), INTENT(in) :: IY, IM, ID, IH, MIT
REAL(kind=8), INTENT(out)   :: JD ! julian day
REAL(kind=4)    :: XI, XJ
INTEGER(kind=4) :: IY1, IM1

! 2. begin program
IF (IM <= 2) THEN   ! january & february
  IY1 = int(IY-1)
  IM1 = int(IM+12)
  JD = dble(int( 365.25*(IY1 + 4716.0)) + int( 30.6001*( IM1 + 1.0)) + 2.0 - &
       int( IY1/100.0 ) + int( int( IY1/100.0 )/4.0 ) + ID - 1524.5) + &
       dble((IH + MIT/60.+0./3600.)/24.)
     
ELSE

  JD = dble( int( 365.25*(IY + 4716.0)) + int( 30.6001*( IM + 1.0)) + 2.0 - &
       int( IY/100.0 ) + int( int( IY/100.0 )/4.0 ) + ID - 1524.5) + &
       dble((IH + MIT/60.+0./3600.)/24.)
     
ENDIF

! 3. END
END SUBROUTINE julian
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 2: Julian to date ~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE  julian_to_date (JD, year, month, day, hour, mint)

!-----------------------------------------------------------------------
! !F90 julian converter
!
! This function converts the Julian dates to Gregorian dates.
!
! Syntax:
! [day,month,year,hour,mint] = julian_to_date(JD)
!
! !Input  parameters:
!    JD               = julian day
!
! !Output parameters:
!    year             = INTEGER year      
!    month            = INTEGER month
!    day              = INTEGER day
!    hour             = INTEGER hour
!    mint             = INTEGER minute
!
!-----------------------------------------------------------------------

! 1. define variables
INTEGER(KIND=4), INTENT(out) :: year, month, day, hour, mint
REAL(KIND=8), INTENT(in)    :: JD ! julian day
REAL(KIND=4)     :: I, D, E, G   
REAL(KIND=4)    :: B, C, Fr
INTEGER(KIND=4)  :: A, a4

I = int(JD + 0.5)
Fr = abs( I - ( JD + 0.5) )

IF (I >= 2299160. ) THEN
     A = int( ( I- 1867216.25 ) / 36524.25 )
     a4 = int( A / 4 )
     B = I + 1. + float(A - a4)
ELSE
     B = I
ENDIF

C = B + 1524.
D = int( ( C - 122.1 ) / 365.25 )
E = int( 365.25 * D )
G = int( ( C - E ) / 30.6001 )
day = int( C - E + Fr - int( 30.6001 * G ) )

IF (G <= 13.5 ) THEN
    month = int(G - 1)
ELSE
    month = int(G - 13)
ENDIF

IF (month > 2.5) THEN
    year = int(D - 4716)
ELSE
    year = int(D - 4715)
ENDIF

hour = int( Fr * 24. )
mint = int( abs( hour -( Fr * 24. ) ) * 60. )

! 3. END
END SUBROUTINE julian_to_date
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ FUNCTION 12: pow x and y ~~~~~~~~~~~~~~~~~~~~~~~~
function pow(x,y) RESULT(res)

  real(kind=4) :: x,y,res
  
  res = x**y
  
end function pow
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ FUNCTION 5: compute sun to earth distance ~~~~~~~~~
FUNCTION compute_earth2sun(julday) RESULT(earth2sun)

!-----------------------------------------------------------------------
! !F90 compute_earth2sun
!
! !Description:
!    This FUNCTION computes the distance from sun to earth.
!
! !Input  parameters:
!    day              = number of day in one year
!
! !Output parameters:
!    earth2sun        = distance of earth to sun
!
!-----------------------------------------------------------------------
  INTEGER, INTENT(in) :: julday
  REAL(KIND=4) :: earth2sun


  earth2sun = 1.0 - 0.016729*cos(0.9856*(julday-4.0)*dtor)

  RETURN
 
END FUNCTION compute_earth2sun
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 7: compute_cos_zenith_angles ~~~~~~~~~~~
SUBROUTINE compute_cos_zenith_angles (satzen, solzen,  &
                                      cos_satzen, cos_solzen)
                                    
!-----------------------------------------------------------------------
! !F90 compute_cos_zenith_angles
!
! !Description:
!    This program is to compute_cos_zenith_angles.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
                                    
REAL(KIND=4), DIMENSION(:,:), INTENT(in)    :: satzen, solzen
!INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: space_mask
REAL(KIND=4), DIMENSION(:,:), INTENT(out)   :: cos_satzen, cos_solzen
INTEGER :: nx,ny,i,j

  nx = size(satzen,dim=1)
  ny = size(satzen,dim=2)
  
  cos_satzen = missing_value_4
  cos_solzen = missing_value_4
  
  DO j = 1, ny
  DO i = 1, nx
    
       cos_satzen(i,j) = cos(satzen(i,j)*dtor)
       cos_solzen(i,j) = cos(solzen(i,j)*dtor)
      
  END DO
  END DO

! 3. END SUBROUTINE    
END SUBROUTINE compute_cos_zenith_angles
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ SUBROUTINE 10: compute scat angle ~~~~~~~~~~~~~~
SUBROUTINE compute_scattering_angles(cos_satzen, cos_solzen, satzen, solzen, relaz, &
                                     scatzen)
                              
!-----------------------------------------------------------------------
! !F90 compute_scattering_angles
!
! !Description:
!    This program is to find p and z in the profile.
!
! !Input  parameters:
!    
!       
! !Output parameters:
!    
!
!-----------------------------------------------------------------------

REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: cos_satzen
REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: cos_solzen
REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: satzen
REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: solzen
REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: relaz
! INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: space_mask
! REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: glintzen
REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: scatzen
INTEGER(KIND=4) :: nx,ny, i, j
REAL(KIND=4) :: sin_solzen
REAL(KIND=4) :: sin_satzen
REAL(KIND=4) :: cos_relaz

nx = SIZE(satzen,dim=1)
ny = SIZE(satzen,dim=2)

scatzen = missing_value_4

DO j = 1, ny
DO i = 1, nx
   
   sin_solzen = sin(solzen(i,j)*dtor)
   sin_satzen = sin(satzen(i,j)*dtor)
   cos_relaz = cos(relaz(i,j)*dtor)
       
   scatzen(i,j) = compute_scat_zen(cos_solzen(i,j), cos_satzen(i,j), sin_solzen, sin_satzen, cos_relaz)
   !glintzen(i,j) = compute_glint_zen(cos_solzen(i,j), cos_satzen(i,j), sin_solzen, sin_satzen, cos_relaz)

END DO
END DO

END SUBROUTINE compute_scattering_angles
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!---------------------------------------------------------------------
! Compute the scattering angle.
!---------------------------------------------------------------------

FUNCTION compute_scat_zen(cos_solzen, cos_satzen, sin_solzen, sin_satzen, cos_relaz) RESULT(scatzen)

 REAL(KIND=4), INTENT(in) :: cos_satzen
 REAL(KIND=4), INTENT(in) :: cos_solzen
 REAL(KIND=4), INTENT(in) :: sin_satzen
 REAL(KIND=4), INTENT(in) :: sin_solzen
 REAL(KIND=4), INTENT(in) :: cos_relaz
 REAL(KIND=4) :: scatzen
 
  scatzen = missing_value_real4

  scatzen = -1.0 * (cos_solzen*cos_satzen - sin_solzen*sin_satzen*cos_relaz)
  
  IF (scatzen > 1.0) scatzen = 0.0
  scatzen = acos(scatzen)/dtor

  RETURN

END FUNCTION compute_scat_zen
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!-------------------------------------------------------------------------
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
SUBROUTINE LOCATE_FLOAT32(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real (kind=4),               intent(in)  :: x
    real (kind=4), dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

END SUBROUTINE LOCATE_FLOAT32

!-------------------------------------------------------------------------
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
SUBROUTINE LOCATE_FLOAT64(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real (kind=8),               intent(in)  :: x
    real (kind=8), dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

END SUBROUTINE LOCATE_FLOAT64
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!-------------------------------------------------------------------------
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
SUBROUTINE LOCATE_INT32(xx, n, x, j)

!   Arguments
    integer(kind=4),               intent(in)  :: n
    integer(kind=4),               intent(out) :: j
    integer(kind=4),               intent(in)  :: x
    integer(kind=4), dimension(:), intent(in)  :: xx

!   Local variables
    integer(kind=4) :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

END SUBROUTINE LOCATE_INT32
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-----------------------------------------------------------------
! Find a particular data bin in an irregularly spaced look-up table.
!-----------------------------------------------------------------

FUNCTION get_irregular_lut_index_int32(data, data_bin_start_lut, data_nbins_lut) result(index)

  INTEGER(KIND=4), INTENT(IN) :: data
  INTEGER(KIND=4), DIMENSION(:), INTENT(IN) :: data_bin_start_lut
  INTEGER(KIND=4), INTENT(IN) :: data_nbins_lut
  INTEGER(KIND=4) :: index
  
  CALL locate(data_bin_start_lut, data_nbins_lut, data, index)
  index = min(max(index, 1), data_nbins_lut)
  
  return

END FUNCTION get_irregular_lut_index_int32

!-----------------------------------------------------------------
! Find a particular data bin in an irregularly spaced look-up table.
!-----------------------------------------------------------------

FUNCTION get_irregular_lut_index_float32(data, data_bin_start_lut, data_nbins_lut) result(index)

  REAL(KIND=4), INTENT(IN) :: data
  REAL(KIND=4), DIMENSION(:), INTENT(IN) :: data_bin_start_lut
  INTEGER(KIND=4), INTENT(IN) :: data_nbins_lut
  INTEGER(KIND=4) :: index
  
  CALL locate(data_bin_start_lut, data_nbins_lut, data, index)
  index = min(max(index, 1), data_nbins_lut)
  
  return

END FUNCTION get_irregular_lut_index_float32

!-----------------------------------------------------------------
! Find a particular data bin in an irregularly spaced look-up table.
!-----------------------------------------------------------------

FUNCTION get_irregular_lut_index_float64(data, data_bin_start_lut, data_nbins_lut) result(index)

  REAL(KIND=8), INTENT(IN) :: data
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: data_bin_start_lut
  INTEGER(KIND=4), INTENT(IN) :: data_nbins_lut
  INTEGER(KIND=4) :: index
  
  CALL locate(data_bin_start_lut, data_nbins_lut, data, index)
  index = min(max(index, 1), data_nbins_lut)
  
  return

END FUNCTION get_irregular_lut_index_float64
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
SUBROUTINE median_filter_int8(image, width, mask)

  INTEGER(KIND=1), DIMENSION(:,:), INTENT(inout) :: image
  INTEGER(KIND=4), INTENT(in) :: width
  INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: mask
  
  SELECT CASE (width)
  CASE (3)
    CALL median_filter9(image, mask)
  CASE (5)
  !  call_median_filter25
  CASE DEFAULT
  !  call_median_filterx
  END SELECT  
  
END SUBROUTINE median_filter_int8

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------

SUBROUTINE median_filter_float32(image, width, mask)

  REAL(KIND=4), DIMENSION(:,:), INTENT(inout) :: image
  INTEGER(KIND=4), INTENT(in) :: width
  INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: mask
  
  SELECT CASE (width)
  CASE (3)
    CALL median_filter9(image, mask)
  CASE (5)
  !  call_median_filter25
  CASE DEFAULT
  !  call_median_filterx
  END SELECT 
  
END SUBROUTINE median_filter_float32
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
SUBROUTINE median_filter9_int8(image, mask)

  INTEGER(KIND=1), DIMENSION(:,:), INTENT(inout) :: image
  INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: mask
  
  INTEGER(KIND=4), parameter :: width = 3
  INTEGER(KIND=4) :: nx, ny, iline, ielem, astatus
  INTEGER(KIND=4) :: jm, jp, im, ip, dpix, arr_size
  INTEGER(KIND=1), DIMENSION(:,:), pointer :: arr
  
  dpix = (width - 1)/2
  arr_size = width*width
  
  nx = SIZE(image,1)
  ny = SIZE(image,2)
  
  ALLOCATE(arr(nx,ny),stat=astatus)
  IF (astatus /= 0) THEN
     PRINT*, "(a,'Not enough memory to ALLOCATE temporary array in median_filter.')"
     STOP
  ENDIF
  
  arr = image
  
  line_loop: DO iline=1, ny
    
    jm = iline - dpix
    jp = iline + dpix
    
    IF (jm < 1 .or. jp > ny) CYCLE
    
    element_loop: DO ielem=1, nx
    
      IF (mask(ielem,iline) == sym%YES) CYCLE
    
         im = ielem - dpix
         ip = ielem + dpix
      
      IF (im < 1 .or. ip > nx) CYCLE
      
      CALL opt_med9_int1(RESHAPE(image(im:ip,jm:jp),(/arr_size/)), arr(ielem,iline))
    
    END DO element_loop
  END DO line_loop
  
  image = arr
  
  DEALLOCATE(arr,stat=astatus)
  
  IF (astatus /= 0) THEN
     PRINT*, "(a,'Error deallocating temporary array in median_filter.')"
     STOP
  ENDIF
  
END SUBROUTINE median_filter9_int8

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------

SUBROUTINE median_filter9_float32(image, mask)

  REAL(KIND=4), DIMENSION(:,:), INTENT(inout) :: image
  INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: mask
  
  INTEGER(KIND=4), parameter :: width = 3
  INTEGER(KIND=4) :: nx, ny, iline, ielem, astatus
  INTEGER(KIND=4) :: jm, jp, im, ip, dpix, arr_size
  REAL(KIND=4), DIMENSION(:,:), pointer :: arr
  
  dpix = (width - 1)/2
  arr_size = width*width
  
  nx = SIZE(image,1)
  ny = SIZE(image,2)
  
  ALLOCATE(arr(nx,ny),stat=astatus)
  IF (astatus /= 0) THEN
     PRINT*, "(a,'Not enough memory to ALLOCATE temporary array in median_filter.')"
     STOP
  ENDIF
  
  arr = image
  
  line_loop: DO iline=1, ny
    
    jm = iline - dpix
    jp = iline + dpix
    
    IF (jm < 1 .or. jp > ny) CYCLE
    
    element_loop: DO ielem=1, nx
    
      IF (mask(ielem,iline) == sym%YES) CYCLE
    
      im = ielem - dpix
      ip = ielem + dpix
      
      IF (im < 1 .or. ip > nx) CYCLE
      
      CALL opt_med9(RESHAPE(image(im:ip,jm:jp),(/arr_size/)), arr(ielem,iline))
    
    END DO element_loop
  END DO line_loop
  
  image = arr
  
  DEALLOCATE(arr,stat=astatus)
  IF (astatus /= 0) THEN
     PRINT*, "(a,'Error deallocating temporary array in median_filter.')"
     STOP
  ENDIF
  
END SUBROUTINE median_filter9_float32
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE gradient2d(grid, nx, ny, mask, min_valid, max_valid, threshold_value, xmax, ymax, num_steps)

  REAL (kind=4), dimension(:,:), intent(in) :: grid
  INTEGER (kind=4), intent(in) :: nx, ny
  INTEGER (kind=1), dimension(:,:), intent(in) :: mask
  REAL (kind=4), intent(in) :: min_valid, max_valid, threshold_value
  INTEGER (kind=4), dimension(:,:), intent(inout) :: xmax, ymax
  INTEGER (kind=4), dimension(:,:), intent(inout), optional :: num_steps
  
  INTEGER (kind=4), parameter :: max_step = 150
  INTEGER (kind=4) :: ielem, iline, im, jm, ip, jp, dx_start, dy_start, index
  INTEGER (kind=4) :: direction, di, dj, i0, j0, i1, j1, ibad
  REAL (kind=4) :: min_grad, ref_value
  INTEGER (kind=4), dimension(8) :: icol, irow, di_default, dj_default
  REAL (kind=4), dimension(8) :: grad
  
  dx_start = 2
  dy_start = 2
  
  di_default = (/0,1,1,1,0,-1,-1,-1/)
  dj_default = (/-1,-1,0,1,1,1,0,-1/)
  
  xmax = missing_value_int4
  ymax = missing_value_int4
  
  if (present(num_steps)) then
    num_steps = missing_value_int4
  endif
  
  line_loop: do iline=1, ny
    
    jm = max(1,iline-dy_start)
    jp = min(ny,iline+dy_start)    
    
    element_loop: do ielem=1, nx
    
      if (mask(ielem,iline) == sym%YES .and. grid(ielem,iline) >= min_valid .and. grid(ielem,iline) <= max_valid) then
      
        im = max(1,ielem-dx_start)
        ip = min(nx,ielem+dx_start)
        
        icol = (/ielem,ip,ip,ip,ielem,im,im,im/)        
        irow = (/jm,jm,iline,jp,jp,jp,iline,jm/)
        
        direction = -999
        min_grad = 99999.0
        do index=1, 8 
          if (grid(icol(index),irow(index)) >= min_valid .and. grid(icol(index),irow(index)) <= max_valid) then
            grad(index) = grid(ielem,iline) - grid(icol(index),irow(index))
            if (grad(index) < min_grad) then
              min_grad = grad(index)
              direction = index
            endif
          endif
        end do
        
        if (direction >= 1 .and. direction <= 8) then
        
          di = di_default(direction)
          dj = dj_default(direction)
        
          do index = 1, max_step
            i0 = max(1,min(ielem + di*index,nx))
            j0 = max(1,min(iline + dj*index,ny))
            i1 = max(1,min(ielem + di*index + di,nx))
            j1 = max(1,min(iline + dj*index + dj,ny))
          
            ref_value = grid(i0,j0)
            ibad = 1
            if (grid(i1,j1) >= min_valid .and. grid(i1,j1) <= max_valid) then 
              ibad = 0
            endif
          
            if (grid(i1,j1) >= threshold_value .or. index == max_step .or. grid(i1,j1) < ref_value .or. ibad == 1) then
              xmax(ielem,iline) = i0
              ymax(ielem,iline) = j0
              if (present(num_steps)) then
                num_steps(ielem,iline) = index
              endif
              exit
            endif
          end do
          
        endif
        
      endif
    
    end do element_loop    
  end do line_loop

END SUBROUTINE gradient2d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Gradient2D_org(grid, nx, ny, mask, min_valid, max_valid, threshold_value, xmax, ymax)
            
  REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: grid
  INTEGER(KIND=4), INTENT(in) :: nx, ny
  INTEGER(KIND=1), DIMENSION(:,:), INTENT(in) :: mask
  REAL(KIND=4), INTENT(in) :: min_valid, max_valid, threshold_value
  INTEGER(KIND=4), DIMENSION(:,:), INTENT(inout) :: xmax, ymax
  
  INTEGER(KIND=4), parameter :: max_step = 150
  INTEGER(KIND=4) :: ielem, iline, im, jm, ip, jp, dx_start, dy_start, index
  INTEGER(KIND=4) :: direction, di, dj, i0, j0, i1, j1, ibad
  REAL(KIND=4) :: min_grad, ref_value
  INTEGER(KIND=4), DIMENSION(8) :: icol, irow, di_default, dj_default
  REAL(KIND=4), DIMENSION(8) :: grad
  
  dx_start = 2
  dy_start = 2
  
  di_default = (/0,1,1,1,0,-1,-1,-1/)
  dj_default = (/-1,-1,0,1,1,1,0,-1/)
  
  xmax = missing_value_int4
  ymax = missing_value_int4
  
  line_loop: DO iline=1, ny
    
    jm = max(1,iline-dy_start)
    jp = min(ny,iline+dy_start)    
    
    element_loop: DO ielem=1, nx
    
      IF (mask(ielem,iline) == sym%YES .and. grid(ielem,iline) >= min_valid .and. grid(ielem,iline) <= max_valid) THEN
      
        im = max(1,ielem-dx_start)
        ip = min(nx,ielem+dx_start)

        icol = (/ielem,ip,ip,ip,ielem,im,im,im/)
        irow = (/jm,jm,iline,jp,jp,jp,iline,jm/)

        min_grad = 99999.0
        DO index=1, 8 
           IF (grid(icol(index),irow(index)) >= min_valid .and. grid(icol(index),irow(index)) <= max_valid) THEN
              grad(index) = grid(ielem,iline) - grid(icol(index),irow(index))
              IF (grad(index) < min_grad) THEN
                 min_grad = grad(index)
                 direction = index
              ENDIF
           ENDIF
        END DO

        di = di_default(direction)
        dj = dj_default(direction)

        DO index = 1, max_step
           i0 = max(1,min(ielem + di*index,nx))
           j0 = max(1,min(iline + dj*index,ny))
           i1 = max(1,min(ielem + di*index + di,nx))
           j1 = max(1,min(iline + dj*index + dj,ny))
  
           ref_value = grid(i0,j0)
           ibad = 1
           IF (grid(i1,j1) >= min_valid .and. grid(i1,j1) <= max_valid) THEN 
              ibad = 0
           ENDIF
  
           IF (grid(i1,j1) >= threshold_value .or. index == max_step .or. grid(i1,j1) < ref_value .or. ibad == 1) THEN
              xmax(ielem,iline) = i0
              ymax(ielem,iline) = j0
              EXIT
           ENDIF
       END DO

      ENDIF
    
    END DO element_loop    
  END DO line_loop

END SUBROUTINE Gradient2D_org
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ SUBROUTINE 11:pack byte ~~~~~~~~~~~~~~~~~~~~~~~~~
!--- This Version packs into one byte words
SUBROUTINE PACK_BYTES_I1(input_bytes,bit_depth,output_bytes)

    INTEGER(KIND=1), DIMENSION(:), INTENT(in):: input_bytes
    INTEGER(KIND=4), DIMENSION(:), INTENT(in):: bit_depth
    INTEGER(KIND=1), DIMENSION(:), INTENT(out):: output_bytes
    INTEGER(KIND=1):: bit_start, bit_END, bit_offset
    INTEGER(KIND=1):: temp_byte
    INTEGER:: n_in,i_in,n_out,i_out
    INTEGER, parameter:: word_bit_depth = 8
  
!--- determine SIZE of vectors
    n_in = SIZE(input_bytes)
    n_out = SIZE(output_bytes)

!--- reset output byte
    output_bytes = 0

!--- initialize
    bit_offset = 0
    bit_start = 0
    bit_END = 0
    i_out = 1

!--- loop through input bytes
    DO i_in = 1, n_in

!--- determine starting and ENDing bit locations
       bit_start = bit_offset + 1
       bit_END = bit_start + bit_depth(i_in) - 1

!--- determine IF this input byte will fit on current output byte, IF not go to next
       IF (bit_END > word_bit_depth) THEN
          i_out = i_out + 1
          bit_offset = 0
          bit_start = bit_offset + 1
          bit_END = bit_start + bit_depth(i_in) - 1
       ENDIF

!--- check for exceeding the space allowed for the packed bytes
       IF (i_out > n_out) THEN
          PRINT *, "ERROR: Insufficient space for bit packing" 
          RETURN
       ENDIF

!--- place input byte into correct position
       temp_byte =0
       temp_byte = ishft(input_bytes(i_in),word_bit_depth-bit_depth(i_in))   !first ishft
       temp_byte = ishft(temp_byte,bit_END - word_bit_depth)                 !second ishft

!--- modIFy output byte
       output_bytes(i_out) = output_bytes(i_out) + temp_byte

!--- update bit offset
       bit_offset = bit_offset + bit_depth(i_in)

   END DO

END SUBROUTINE  PACK_BYTES_I1

!--- This Version packs into two byte words
SUBROUTINE PACK_BYTES_I2(input_bytes,bit_depth,output_bytes)

    INTEGER(KIND=1), DIMENSION(:), INTENT(in):: input_bytes
    INTEGER(KIND=4), DIMENSION(:), INTENT(in):: bit_depth
    INTEGER(KIND=2), DIMENSION(:), INTENT(out):: output_bytes
    INTEGER(KIND=1):: bit_start, bit_END, bit_offset
    INTEGER(KIND=2):: temp_byte                 
    INTEGER:: n_in,i_in,n_out,i_out
    INTEGER, parameter:: word_bit_depth = 16

!--- determine SIZE of vectors
    n_in = SIZE(input_bytes)
    n_out = SIZE(output_bytes)

!--- reset output byte
    output_bytes = 0

!--- initialize
    bit_offset = 0
    bit_start = 0
    bit_END = 0
    i_out = 1

!--- loop through input bytes
   DO i_in = 1, n_in

!--- determine starting and ENDing bit locations
      bit_start = bit_offset + 1
      bit_END = bit_start + bit_depth(i_in) - 1

!--- determine IF this input byte will fit on current output byte, IF not go to next
      IF (bit_END > word_bit_depth) THEN
         i_out = i_out + 1
         bit_offset = 0
         bit_start = bit_offset + 1
         bit_END = bit_start + bit_depth(i_in) - 1
      ENDIF

!--- check for exceeding the space allowed for the packed bytes
      IF (i_out > n_out) THEN
         PRINT *, "ERROR: Insufficient space for bit packing"
         RETURN
      ENDIF

!--- place input byte into correct position
      temp_byte =0
      temp_byte = ishft(input_bytes(i_in),word_bit_depth-bit_depth(i_in))   !first ishft
      temp_byte = ishft(temp_byte,bit_END - word_bit_depth)                 !second ishft

!--- modIFy output byte
      output_bytes(i_out) = output_bytes(i_out) + temp_byte

!--- update bit offset
      bit_offset = bit_offset + bit_depth(i_in)

   END DO

END SUBROUTINE  PACK_BYTES_I2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 4: compute spatial unIFormity ~~~~~~~~~~
SUBROUTINE compute_spatial_uniformity(dx, dy, space_mask, data1, data_mean, data_max, data_min, data_uni)

!-----------------------------------------------------------------------
! !F90 compute_spatial_unIFormity
!
! !Description:
!    This program is to output the message of reading file.
!
! !Input  parameters:
!    dx              = x
!    dy              = y
!    space_mask      = space mask
!    data1           = data
!
! !Output parameters:
!    data_mean       = mean value
!    data_max        = maximum value
!    data_min        = minimum value
!    data_uni        = unIFormity 
!
!-----------------------------------------------------------------------
                                                                                                          
  INTEGER(KIND=4), INTENT(in) :: dx, dy
  INTEGER(KIND=1), INTENT(in), DIMENSION(:,:) :: space_mask
  REAL(KIND=4), INTENT(in), DIMENSION(:,:) :: data1
  REAL(KIND=4), INTENT(out), DIMENSION(:,:), allocatable :: data_mean, data_max, data_min, data_uni
  INTEGER(KIND=4) :: nx, ny, nx_uni, ny_uni, nsub, astatus
  INTEGER(KIND=4) :: ielem, iline, ielem1, ielem2, iline1, iline2, n_good
  REAL(KIND=4), DIMENSION(:,:), allocatable :: temp
  INTEGER(KIND=4), DIMENSION(:,:), allocatable :: good
  INTEGER(KIND=1), DIMENSION(:,:), allocatable :: space_mask_temp
  
  nx = SIZE(data1,1)
  ny = SIZE(data1,2)
  
  nx_uni = 2*dx + 1
  ny_uni = 2*dy + 1
  nsub = nx_uni*ny_uni
  
  ALLOCATE(temp(nx_uni,ny_uni), good(nx_uni,ny_uni), &
          space_mask_temp(nx_uni,ny_uni), data_mean(nx,ny), data_max(nx,ny), &
          data_min(nx,ny),data_uni(nx,ny),stat=astatus)
          
  IF (astatus /= 0) THEN
     PRINT*, "(a,'Not enough memory to ALLOCATE spatial unIFormity arrays.')"
     STOP
  ENDIF

  data_mean = missing_value_real4
  data_max = missing_value_real4
  data_min = missing_value_real4
  data_uni = missing_value_real4
        
  line_loop: DO iline=1, ny
    
    iline1 = max(1,iline-dy)
    iline2 = min(ny,iline+dy)    
    ny_uni = (iline2 - iline1) + 1
    
    element_loop: DO ielem=1, nx
    
      !data_mean(ielem,iline) = missing_value_real4
      !data_max(ielem,iline) = missing_value_real4
      !data_min(ielem,iline) = missing_value_real4
      !data_uni(ielem,iline) = missing_value_real4
      
      IF (space_mask(ielem,iline) == sym%NO_SPACE) THEN
      
        ielem1 = max(1,ielem-dx)
        ielem2 = min(nx,ielem+dx)
        nx_uni = (ielem2 - ielem1) + 1
      
        space_mask_temp = sym%SPACE
        temp(1:nx_uni,1:ny_uni) = data1(ielem1:ielem2,iline1:iline2)
        space_mask_temp(1:nx_uni,1:ny_uni) = space_mask(ielem1:ielem2,iline1:iline2)
        n_good = nsub - sum(space_mask_temp)
      
        IF (n_good > 0) THEN
          temp = (1 - space_mask_temp)*temp
          data_mean(ielem,iline) =  sum(temp(1:nx_uni,1:ny_uni))/n_good
          data_uni(ielem,iline) = sqrt(max(0.0,(sum((temp(1:nx_uni,1:ny_uni))**2)/n_good - data_mean(ielem,iline)**2)))
          data_max(ielem,iline) = maxval(temp(1:nx_uni,1:ny_uni))
          data_min(ielem,iline) = minval(temp(1:nx_uni,1:ny_uni)) 
        ENDIF

      ENDIF
    
    END DO element_loop    
  END DO line_loop
  
  DEALLOCATE(temp, good, space_mask_temp,stat=astatus)
  IF (astatus /= 0) THEN
     PRINT*,"(a,'Error deallocating temporary spatial unIFormity arrays.')"
     STOP
  ENDIF

END SUBROUTINE compute_spatial_uniformity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 7: destory spatial unIFormity ~~~~~~~~~~
SUBROUTINE destroy_spatial_uniformity(data_mean, data_max, data_min, data_uni)

!-----------------------------------------------------------------------
! !F90 destroy_spatial_unIFormity
!
! !Description:
!    This program is to destory spatial unIFormity.
!
! !Input  parameters:
!
! !Output parameters:
!
!-----------------------------------------------------------------------
                                                                                                           
  REAL(KIND=4), INTENT(inout), DIMENSION(:,:), allocatable :: data_mean, data_max, data_min, data_uni
  INTEGER(KIND=4) :: astatus
  
  DEALLOCATE(data_mean, data_max, data_min, data_uni,stat=astatus)
  
  IF (astatus /= 0) THEN
     PRINT*,"(a,'Error deallocating spatial unIFormity arrays.')"
     STOP
  ENDIF
  
! 3. END SUBROUTINE
END SUBROUTINE  destroy_spatial_uniformity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ SUBROUTINE 8: profile computation t p z ~~~~~~~~~
SUBROUTINE prof_lookup_using_t(zlev, plev, tlev, z, p, t, ilev, a)
 
!-----------------------------------------------------------------------
! !F90 compute_levels
!
! !Description:
!    This program is to find p and z in the profile.
!
! !Input  parameters:
!    none
!       
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  REAL(kind=4), INTENT(in), DIMENSION(:) :: zlev, plev, tlev
  REAL(kind=4), INTENT(in) :: t
  REAL(kind=4), INTENT(out) :: p, z
  INTEGER(kind=4), INTENT(out) :: ilev
  REAL(kind=4), INTENT(out) :: a
  INTEGER(kind=4) :: nlev
  REAL(kind=4) :: dp, dt, dz
 
  nlev = SIZE(plev,1)
 
!   
  CALL LOCATE(tlev,nlev,t,ilev)
  ilev = max(1,min(nlev-1,ilev))
 
  dp = plev(ilev+1) - plev(ilev)
  dt = tlev(ilev+1) - tlev(ilev)
  dz = zlev(ilev+1) - zlev(ilev)
 
  IF (dp /= 0.0) THEN
    a = (t - tlev(ilev))/dt
    p = plev(ilev) + a*dp
    z = zlev(ilev) + a*dz
  ELSE
    a = 0.0
    p = plev(ilev)
    z = zlev(ilev)
  ENDIF
 
! 3. END SUBROUTINE 
END SUBROUTINE  prof_lookup_using_t


SUBROUTINE prof_lookup_using_p(zlev, plev, tlev, z, p, t, ilev, a)
   
  REAL(KIND=4), INTENT(in), DIMENSION(:) :: zlev, plev, tlev
  REAL(kind=4), INTENT(in) :: p
  REAL(kind=4), INTENT(out) :: z, t
  INTEGER(kind=4), INTENT(out) :: ilev
  REAL(kind=4), INTENT(out) :: a
  INTEGER(kind=4) :: nlev
  REAL(kind=4) :: dp, dt, dz

  nlev = SIZE(plev,1)
   
  CALL LOCATE(plev,nlev,p,ilev)
  ilev = max(1,min(nlev-1,ilev))

  dp = plev(ilev+1) - plev(ilev)
  dt = tlev(ilev+1) - tlev(ilev)
  dz = zlev(ilev+1) - zlev(ilev)

  IF (dp /= 0.0) THEN
    a = (p - plev(ilev))/dp
    t = tlev(ilev) + a*dt
    z = zlev(ilev) + a*dz
  ELSE
    a = 0.0
    t = tlev(ilev) 
    z = zlev(ilev)
  ENDIF
 
END SUBROUTINE prof_lookup_using_p


SUBROUTINE prof_lookup_using_z(zlev, plev, tlev, z, p, t, ilev, a)
   
  REAL(kind=4), INTENT(in), DIMENSION(:) :: zlev, plev, tlev
  REAL(kind=4), INTENT(in) :: z
  REAL(kind=4), INTENT(out) :: p, t
  INTEGER(kind=4), INTENT(out) :: ilev
  REAL(kind=4), INTENT(out) :: a
  INTEGER(kind=4) :: nlev
  REAL(kind=4) :: dp, dt, dz

  nlev = SIZE(plev,1)
   
  CALL LOCATE(zlev,nlev,z,ilev)
  ilev = max(1,min(nlev-1,ilev))

  dp = plev(ilev+1) - plev(ilev)
  dt = tlev(ilev+1) - tlev(ilev)
  dz = zlev(ilev+1) - zlev(ilev)

  IF (dp /= 0.0) THEN
    a = (z - zlev(ilev))/dz
    p = plev(ilev) + a*dp
    t = tlev(ilev) + a*dt
  ELSE
    a = 0.0
    p = plev(ilev) 
    t = tlev(ilev)
  ENDIF
 
END SUBROUTINE prof_lookup_using_z
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 5: INVERT 2x2 matrix ~~~~~~~~~~~~~~~~~~~
SUBROUTINE INVERT_2x2(A,A_inv,ierr)

!-----------------------------------------------------------------------
! !F90 INVERT 2x2 matrix
!
! !Description:
!    This program is to invert 2x2 matrixes.
!
! !Input  parameters:
!
! !Output parameters:
!
!-----------------------------------------------------------------------
  REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: A
  REAL(KIND=4), DIMENSION(:,:), INTENT(out):: A_inv
  REAL(KIND=4):: determinant
  INTEGER(KIND=4), INTENT(out):: ierr

!--- compute determinant
  ierr = 0
  determinant = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  
  IF (determinant == 0.0) THEN
!       PRINT *, "Singular Matrix in Invert 2x2"
        ierr = 1
  ENDIF

!--- compute inverse
  A_inv(1,1) = A(2,2)
  A_inv(1,2) = -A(1,2)
  A_inv(2,1) = -A(2,1)
  A_inv(2,2) = A(1,1)
  A_inv = A_inv / determinant

! 3. END SUBROUTINE
END SUBROUTINE INVERT_2x2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 6: INVERT 3x3 matrix ~~~~~~~~~~~~~~~~~~~
SUBROUTINE INVERT_3x3(A,A_inv,ierr)  

!-----------------------------------------------------------------------
! !F90 INVERT_3x3
!
! !Description:
!    This program is to invert 3x3 matrixes.
!
! !Input  parameters:
!
! !Output parameters:
!
!-----------------------------------------------------------------------
  REAL(KIND=4), DIMENSION(:,:), INTENT(in):: A
  REAL(KIND=4), DIMENSION(:,:), INTENT(out):: A_inv
  INTEGER(KIND=4), INTENT(out):: ierr
  REAL(KIND=4):: determinant

  ierr = 0
!--- compute determinant
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) - &
                A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3)) + &
                A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  IF (determinant == 0.0) THEN
!       PRINT *, "Singular Matrix in Invert 3x3"
        ierr = 1
  ENDIF

!--- compute inverse
  A_inv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
  A_inv(1,2) = A(1,3)*A(3,2) - A(3,3)*A(1,2)
  A_inv(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
  A_inv(2,1) = A(2,3)*A(3,1) - A(3,3)*A(2,1)
  A_inv(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
  A_inv(2,3) = A(1,3)*A(2,1) - A(2,3)*A(1,1)
  A_inv(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
  A_inv(3,2) = A(1,2)*A(3,1) - A(3,2)*A(1,1)
  A_inv(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
  A_inv = A_inv / determinant

! 3. END SUBROUTINE 
END SUBROUTINE INVERT_3x3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!-------------------------- END MODULE ---------------------------------
end module numerical
