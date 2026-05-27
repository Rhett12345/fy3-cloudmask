module frontend_module


!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is a frontend code for FY3/MERSI-II product code/
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

use names_module
use numerical
use data_arrays_module

implicit none


CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~FUNCTION 1: show system time ~~~~~~~~~~~~~~~~~~~~~~~~~~
function time_elapsed(clock_count_start, &
                      clock_count_end,   &
                      clock_rate,        &
                      option) RESULT(tdiff)

!-----------------------------------------------------------------------
! !F90 get_nwptime 
!
! !Description:
!    This program is to calculate and show total time for calculation.
!
! !Input  parameters:
!    clock_count_start = clock start
!    clock_count_end   = clock end
!    clock_rate        = clock rate
!    option            = option for print time 
!
! !Output parameters:
!
!-----------------------------------------------------------------------
                              
INTEGER, INTENT(in) :: clock_count_start, clock_count_END, clock_rate
INTEGER, INTENT(in) :: option
REAL(4) :: tdIFf
INTEGER :: dd, hh, mm, ss
    
    tdiff = (clock_count_end-clock_count_start)/real(clock_rate,kind=4)
    dd = int(tdIFf/86400)
    hh = int((int(modulo(tdIFf,86400.0)))/3600)
    mm = int((int(modulo(tdIFf,3600.0)))/60)
    ss = int(int(modulo(tdIFf,60.0)))
    
    IF (option == 1) THEN
    
      PRINT*," Total processing time: "
      PRINT*,'     Days    = ',dd
      PRINT*,'     Hours   = ',hh
      PRINT*,'     Minutes = ',mm
      PRINT*,'     Seconds = ',ss
      
    ELSE IF (option == 2) THEN
    
      PRINT*," Processing time for this L1b data: "
      PRINT*,'     Days    = ',dd
      PRINT*,'     Hours   = ',hh
      PRINT*,'     Minutes = ',mm
      PRINT*,'     Seconds = ',ss
      
    ENDIF
    
end function time_elapsed
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 1: extract satellite observation time ~~
subroutine extract_sattime (L1b_path,xname,year,month,day,hour,mint,dayn)

!-----------------------------------------------------------------------
! !F90 extract_sattime
!
! !Description:
!    This program is to extract satellite observation time.
!
! !Input  parameters:
!    L1b_path        = index number of senser option
!    xname           = satellite name
!
! !Output parameters:
!    year   sattime%year     = INTEGER year
!    month  sattime%month    = INTEGER month
!    day    sattime%day      = INTEGER day
!    hour   sattime%hour     = INTEGER hour 
!    mint   sattime%mint     = INTEGER minute
!    daynum                  = INTEGER number of day in one year
!
!-----------------------------------------------------------------------

!USE module

IMPLICIT NONE

! 1. define variables
!===== 1.1.input and out variables
!INTEGER(KIND=int1), INTENT(in):: index   ! index
CHARACTER*(*), INTENT(in)  :: L1b_path        ! path
CHARACTER*(*), INTENT(in)  :: xname     ! satname
INTEGER, INTENT(out) :: year, month, day, hour, mint, dayn

!===== 1.2.middle variables  
INTEGER(KIND=4), dimension(1:5) :: fp, ep      ! the first [fp] and END [ep] position of time string  

!===== 1.3.other variables
INTEGER length, ierr, leap_flg, day1, L1, L2
CHARACTER(len=200) :: fname, satname

!******* these parameters are set for Metesat-8 seviri *********
!data fp /21, 25, 29, 31/  ! year, day, hour, mint
!data ep /24, 27, 30, 32/
!***************************************************************
! 2. begin program
  PRINT*,'  ... fylat extract satellite observation time'
  
!xxxxxxxxxxxx
!Note: here, the satellite is FY3D/MERSI_II
!      IF we USE dIFferent satellite, we should change CHARACTER length.
!      Variable fp(x) and ep(x) of position of string should be changed. 
!xxxxxxxxxxxx
  fp(1:5) = (/20, 24, 26, 29, 31/) ! year, month, day, hour, mint
  ep(1:5) = (/23, 25, 27, 30, 32/)

  L1 = LEN(trim(L1b_path))
  L2 = LEN(trim(xname)) 
  satname = xname(L1+1:L2)   
  !print*,'  satellite L1b data name = ', satname
  
 !===== 2.1. read year
  fname = satname(fp(1):ep(1))
  CALL ICNVRT(1,year,fname,length,ierr)

  !===== 2.2. read month
  leap_flg = leap_year_fct(year) 
  fname    = satname(fp(2):ep(2))
  CALL ICNVRT(1,month,fname,length,ierr) 
  
  !===== 2.3. read day
  fname    = satname(fp(3):ep(3))
  CALL ICNVRT(1,day,fname,length,ierr) 
  dayn  = compute_daynum(month, day, leap_flg) 
 
  !===== 2.4. read hour
  fname = satname(fp(4):ep(4))
  CALL ICNVRT(1,hour,fname,length,ierr)

  !===== 2.5. read minute
  fname = satname(fp(5):ep(5))
  CALL ICNVRT(1,mint,fname,length,ierr)

  
! 3. END SUBROUTINE   
END SUBROUTINE extract_sattime
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module frontend_module