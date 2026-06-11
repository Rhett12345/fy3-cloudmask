module thresholds_read_module

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    cloud mask thresholds read module
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C    none
!C
!C !Author's information
!C    Author: Min Min/Wu Xiao/Zheng Zhaojun/Liu Ruixia/Zhang Miao/Yang Changjun/Qiu Hong
!C    E-mail: minmin@cma.gov.cn
!C    Tel   : 86-010-68406763
!C    National Satellite Meteorological Center, CMA 
!C  
!C !END
!C----------------------------------------------------------------------

use names_module
!use data_arrays_module
use constant

implicit none

!+++++++++++++++++++ step 1: Variables +++++++++++++++++++++++++++++++++
!     Variables
!     ver_string         String returned from thresholds file
!     ver_name           String searched for in thresholds file
!     thr_ver_id_Aqua    String containing current version id for Aqua
!     thr_ver_id_Terra   String containing current version id for Terra
!     thr_satid_terra    Identifier for Terra thresholds file 
!     thr_satid_aqua     Identifier for Aqua thresholds file 


CHARACTER(len=19),parameter :: ver_name='thresholds_file_ver'
CHARACTER*100 threshold_file_name
!character*6  thr_ver_id_Aqua, thr_ver_id_Terra
!character*5  thr_satid_terra, thr_satid_aqua

!parameter (ver_name='thresholds_file_ver')
!parameter (thr_ver_id_Aqua='v5.1.0')
!parameter (thr_ver_id_Terra='v5.1.0')
!parameter (thr_satid_terra='Terra')
!parameter (thr_satid_aqua='Aqua ')

integer(kind=4) :: lrn_thr_par 
! ANC THRESHOLDS
! --- Fortran OPEN
integer(kind=4), parameter   ::	ANC_THRESHOLD_UNIT = 41    
character(len= 3), parameter ::	ANC_THRESHOLD_STATUS = 'OLD'     
character(len= 9), parameter ::	ANC_THRESHOLD_FORM = 'FORMATTED' 
character(len=10), parameter ::	ANC_THRESHOLD_ACCESS = 'SEQUENTIAL' 

contains
!+++++++++++++++++++ step 2: Subroutines +++++++++++++++++++++++++++++++
subroutine thresholds_read( sensor_id )

!-------------------------------------------------------------------
! !DESCRIPTION:
!     Read all fylat cloud mask thresholds from parameter file.  
!     Also reads and verifies version and satellite id information 
!     from thresholds files.
!
! !INPUT PARAMETERS:
!     None
!
! !OUTPUT PARAMETERS:
!     sensor_id    Version string (RCS Id) for threshold parameter
!                  file.
!
!     Threshold values defined in the following include files:
!     LandDay_desert_c_thr.inc
!     PolarDay_desert_c_thr.inc
!     LandDay_desert_thr.inc
!     PolarDay_desert_thr.inc
!     LandDay_thr.inc
!     PolarDay_land_thr.inc
!     LandDay_coast_thr.inc
!     PolarDay_coast_thr.inc
!     LandNite_thr.inc
!     PolarNite_land_thr.inc
!     PolarDay_snow_thr.inc
!     Day_snow_thr.inc
!     PolarNite_snow_thr.inc
!     Nite_snow_thr.inc
!     ocean_day_thr.inc
!     PolarDay_ocean_thr.inc
!     ocean_nite_thr.inc
!     PolarNite_ocean_thr.inc
!     shadows_thr.inc
!     snglntr_thr.inc
!     spatial_var_thr.inc
!     noncld_obs_chk.inc
!     snow_mask.inc
!     Antarctic_day_thr.inc
!     swc_ndvi.inc
!     land_restoral.inc
!
!
! !END
!--------------------------------------------------------------------

!      implicit none
!      SAVE

! ... Arguments

      integer(kind=1) :: sensor_id, io_err
      character*80 version
      CHARACTER*80 ver_string    
        
! ... Local variables

      character*160 errmsg
      integer number

! ... MOD35 PCF number include file (defines LRN_THR_PAR)
!      include 'mod35.inc'
                  
! ... Version and satellite ID internal verification
!      include 'thresholds_ver.inc'
!      include 'platform_name.inc'

! ... Threshold include files

      include 'LandDay_desert_c_thr.inc'
      include 'PolarDay_desert_c_thr.inc'
      include 'LandDay_desert_thr.inc'
      include 'PolarDay_desert_thr.inc'
      include 'LandDay_thr.inc'
      include 'PolarDay_land_thr.inc'
      include 'LandDay_coast_thr.inc'
      include 'PolarDay_coast_thr.inc'
      include 'LandNite_thr.inc'
      include 'PolarNite_land_thr.inc'
      include 'PolarDay_snow_thr.inc'
      include 'Day_snow_thr.inc'
      include 'PolarNite_snow_thr.inc'
      include 'Nite_snow_thr.inc'
      include 'ocean_day_thr.inc'
      include 'PolarDay_ocean_thr.inc'
      include 'ocean_nite_thr.inc'
      include 'PolarNite_ocean_thr.inc'
      include 'shadows_thr.inc'
      include 'snglntr_thr.inc'
      include 'spatial_var_thr.inc'
      include 'noncld_obs_chk.inc'
      include 'snow_mask.inc'
      include 'Antarctic_day_thr.inc'
      include 'swc_ndvi.inc'
      include 'land_restoral.inc'
      include 'pfmft_nfmft_thr.inc'  ! added by minmin     

                 
! --- set the threshold unit number
      lrn_thr_par = ANC_THRESHOLD_UNIT
      
      if (sensor_id == 1 .or. sensor_id == 2) then
         threshold_file_name = trim(code_root_path)//'coeff/fylat_thresholds.mersi.aqua.v8'
      endif

      if (sensor_id == 21) then   ! fy3d-mersi-ii
         threshold_file_name = trim(code_root_path)//'coeff/fylat_thresholds.mersi.ii3d.v8'
      endif

!     open threshold file
      open( FILE=threshold_file_name,     &
            UNIT= ANC_THRESHOLD_UNIT,     &
            STATUS= ANC_THRESHOLD_STATUS, &
            FORM= ANC_THRESHOLD_FORM,     &
            ACCESS= ANC_THRESHOLD_ACCESS, &
            IOSTAT=io_err)

! ... Threshold parameter file version (for output to file)
      call param_string( lrn_thr_par, 'rcs_id', version )
      
! ... Get threshold ID info and perform internal verification.
      call param_string(lrn_thr_par, ver_name, ver_string)

! ... Daytime coastal desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'lds11_12hi_c', number, lds11_12hi_c )
      call param_real( lrn_thr_par, 'lds11_4hi_c',  number, lds11_4hi_c  )
      call param_real( lrn_thr_par, 'lds11_4lo_c',  number, lds11_4lo_c  )
      call param_real( lrn_thr_par, 'ldsco2_c',     number, ldsco2_c     )
      call param_real( lrn_thr_par, 'ldsh20_c',     number, ldsh20_c     )
      call param_real( lrn_thr_par, 'ldsref2_c',    number, ldsref2_c    )
      call param_real( lrn_thr_par, 'ldsref3_c',    number, ldsref3_c    )
      number = 2
      call param_real( lrn_thr_par, 'ldstci_c',    number, ldstci_c  )

! ... Daytime polar coastal desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pds11_12hi_c', number, pds11_12hi_c )
      call param_real( lrn_thr_par, 'pds11_4hi_c',  number, pds11_4hi_c  )
      call param_real( lrn_thr_par, 'pds11_4lo_c',  number, pds11_4lo_c  )
      call param_real( lrn_thr_par, 'pdsh20_c',     number, pdsh20_c     )
      call param_real( lrn_thr_par, 'pdsref2_c',    number, pdsref2_c    )
      call param_real( lrn_thr_par, 'pdsref3_c',    number, pdsref3_c    )
      number = 2
      call param_real( lrn_thr_par, 'pdstci_c',    number, pdstci_c  )

! ... Daytime desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'lds11_12hi', number, lds11_12hi )
      call param_real( lrn_thr_par, 'lds11_4hi',  number, lds11_4hi  )
      call param_real( lrn_thr_par, 'lds11_4lo',  number, lds11_4lo  )
      call param_real( lrn_thr_par, 'ldsco2',     number, ldsco2     )
      call param_real( lrn_thr_par, 'ldsh20',     number, ldsh20     )
      call param_real( lrn_thr_par, 'ldsref2',    number, ldsref2    )
      call param_real( lrn_thr_par, 'ldsref3',    number, ldsref3    )
      number = 2
      call param_real( lrn_thr_par, 'ldstci',    number, ldstci    )

! ... Daytime polar desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pds11_12hi', number, pds11_12hi )
      call param_real( lrn_thr_par, 'pds11_4hi',  number, pds11_4hi  )
      call param_real( lrn_thr_par, 'pds11_4lo',  number, pds11_4lo  )
      call param_real( lrn_thr_par, 'pdsh20',     number, pdsh20     )
      call param_real( lrn_thr_par, 'pdsref2',    number, pdsref2    )
      call param_real( lrn_thr_par, 'pdsref3',    number, pdsref3    )
      number = 2
      call param_real( lrn_thr_par, 'pdstci',    number, pdstci    )

! ... Daytime land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'dl11_12hi', number, dl11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'dl11_4lo',  number, dl11_4lo  )
      call param_real( lrn_thr_par, 'dlco2',     number, dlco2     )
      call param_real( lrn_thr_par, 'dlh20',     number, dlh20     )
      call param_real( lrn_thr_par, 'dlref1',    number, dlref1    )
      call param_real( lrn_thr_par, 'dlref3',    number, dlref3    )
      call param_real( lrn_thr_par, 'dlvrat',    number, dlvrat    )
      number = 2
      call param_real( lrn_thr_par, 'dltci',    number, dltci    )

! ... Daytime polar land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pdl11_12hi', number, pdl11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'pdl11_4lo',  number, pdl11_4lo  )
      call param_real( lrn_thr_par, 'pdlh20',     number, pdlh20     )
      call param_real( lrn_thr_par, 'pdlref1',    number, pdlref1    )
      call param_real( lrn_thr_par, 'pdlref3',    number, pdlref3    )
      call param_real( lrn_thr_par, 'pdlvrat',    number, pdlvrat    )
      number = 2
      call param_real( lrn_thr_par, 'pdltci',    number, pdltci    )

! ... Daytime coastal land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'dl11_12hi_t2', number, dl11_12hi_t2 )
      number = 4
      call param_real( lrn_thr_par, 'dl11_4lo_t2',  number, dl11_4lo_t2  )
      call param_real( lrn_thr_par, 'dlco2_t2',     number, dlco2_t2     )
      call param_real( lrn_thr_par, 'dlh20_t2',     number, dlh20_t2     )
      call param_real( lrn_thr_par, 'dlref1_t2',    number, dlref1_t2    )
      call param_real( lrn_thr_par, 'dlref3_t2',    number, dlref3_t2    )
      number = 2
      call param_real( lrn_thr_par, 'dltci_t2',    number, dltci_t2 )

! ... Daytime polar coastal land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pdl11_12hi_t2', number, pdl11_12hi_t2 )
      number = 4
      call param_real( lrn_thr_par, 'pdl11_4lo_t2',  number, pdl11_4lo_t2  )
      call param_real( lrn_thr_par, 'pdlh20_t2',     number, pdlh20_t2     )
      call param_real( lrn_thr_par, 'pdlref1_t2',    number, pdlref1_t2    )
      call param_real( lrn_thr_par, 'pdlref3_t2',    number, pdlref3_t2    )
      number = 2
      call param_real( lrn_thr_par, 'pdltci_t2',    number, pdltci_t2    )

! ... Nighttime land thresholds.

      number = 4
      call param_real( lrn_thr_par, 'nl4_12hi', number, nl4_12hi )
      call param_real( lrn_thr_par, 'nl4_12lo', number, nl4_12lo )
      call param_real( lrn_thr_par, 'nlco2',    number, nlco2    )
      call param_real( lrn_thr_par, 'nlh20',    number, nlh20    )
      call param_real( lrn_thr_par, 'nl7_11s',    number, nl7_11s  )
      call param_real( lrn_thr_par, 'nl_11_4l',    number, nl_11_4l  )
      call param_real( lrn_thr_par, 'nl_11_4h',    number, nl_11_4h  )
      call param_real( lrn_thr_par, 'nl_11_4m',    number, nl_11_4m  )
      number = 2
      call param_real( lrn_thr_par, 'bt_diff_bounds', number, bt_diff_bounds  )
      number = 1
      call param_real( lrn_thr_par, 'nl11_12hi', number, nl11_12hi)

! ... Nighttime polar land thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pnlh20',    number, pnlh20    )
      number = 1
      call param_real( lrn_thr_par, 'pnl11_12hi', number, pnl11_12hi)

! ... Daytime polar snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'dpsh20',  number, dpsh20  )
      call param_real( lrn_thr_par, 'dpsref1', number, dpsref1 )
      call param_real( lrn_thr_par, 'dpsref3', number, dpsref3 )
      call param_real( lrn_thr_par, 'dps4_11l', number, dps4_11l )
      call param_real( lrn_thr_par, 'dps4_11h', number, dps4_11h )
      call param_real( lrn_thr_par, 'dps4_11m1', number, dps4_11m1 )
      call param_real( lrn_thr_par, 'dps4_11m2', number, dps4_11m2 )
      call param_real( lrn_thr_par, 'dps4_11m3', number, dps4_11m3 )
      call param_real( lrn_thr_par, 'bt_11_bnds3', number, bt_11_bnds3 )
      number = 2
      call param_real( lrn_thr_par, 'dpstci',    number, dpstci    )
      number = 1
      call param_real( lrn_thr_par, 'dps11_12hi', number, dps11_12hi    )
      call param_real( lrn_thr_par, 'dps11_12adj', number, dps11_12adj    )

! ... Antarctic day thresholds.

      number = 4
      call param_real( lrn_thr_par, 'ant4_11l', number, ant4_11l )
      call param_real( lrn_thr_par, 'ant4_11h', number, ant4_11h )
      call param_real( lrn_thr_par, 'ant4_11m1', number, ant4_11m1 )
      call param_real( lrn_thr_par, 'ant4_11m2', number, ant4_11m2 )
      call param_real( lrn_thr_par, 'ant4_11m3', number, ant4_11m3 )
      call param_real( lrn_thr_par, 'bt_11_bnds4', number, bt_11_bnds4 )
      call param_real( lrn_thr_par, 'anth20',  number, anth20  )

! ... Daytime snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'ds4_11', number, ds4_11 )
      call param_real( lrn_thr_par, 'ds4_11hel', number, ds4_11hel )
      call param_real( lrn_thr_par, 'dsco2',  number, dsco2  )
      call param_real( lrn_thr_par, 'dsh20',  number, dsh20  )
      call param_real( lrn_thr_par, 'dsref3', number, dsref3 )
      number = 2
      call param_real( lrn_thr_par, 'dstci',    number, dstci    )
      number = 1
      call param_real( lrn_thr_par, 'ds11_12hi', number, ds11_12hi )
      call param_real( lrn_thr_par, 'ds11_12adj', number, ds11_12adj )

! ... Nighttime polar snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pn_4_12l',  number, pn_4_12l  )
      call param_real( lrn_thr_par, 'pn_4_12h',  number, pn_4_12h  )
      call param_real( lrn_thr_par, 'pn_4_12m1',  number, pn_4_12m1 )
      call param_real( lrn_thr_par, 'pn_4_12m2',  number, pn_4_12m2 )
      call param_real( lrn_thr_par, 'pn_4_12m3',  number, pn_4_12m3 )
      call param_real( lrn_thr_par, 'pn_7_11l',  number, pn_7_11l  )
      call param_real( lrn_thr_par, 'pn_7_11h',  number, pn_7_11h  )
      call param_real( lrn_thr_par, 'pn_7_11m1',  number, pn_7_11m1 )
      call param_real( lrn_thr_par, 'pn_7_11m2',  number, pn_7_11m2 )
      call param_real( lrn_thr_par, 'pn_7_11m3',  number, pn_7_11m3 )
      call param_real( lrn_thr_par, 'pn_7_11lw',  number, pn_7_11lw  )
      call param_real( lrn_thr_par, 'pn_7_11hw',  number, pn_7_11hw  )
      call param_real( lrn_thr_par, 'pn_7_11m1w',  number, pn_7_11m1w )
      call param_real( lrn_thr_par, 'pn_7_11m2w',  number, pn_7_11m2w )
      call param_real( lrn_thr_par, 'pn_7_11m3w',  number, pn_7_11m3w )
      call param_real( lrn_thr_par, 'pnsh20',    number, pnsh20    )
      call param_real( lrn_thr_par, 'pn_11_4l',  number, pn_11_4l  )
      call param_real( lrn_thr_par, 'pn_11_4h',  number, pn_11_4h  )
      call param_real( lrn_thr_par, 'pn_11_4m1',  number, pn_11_4m1 )
      call param_real( lrn_thr_par, 'pn_11_4m2',  number, pn_11_4m2 )
      call param_real( lrn_thr_par, 'pn_11_4m3',  number, pn_11_4m3 )
      call param_real( lrn_thr_par, 'bt_11_bounds', number, bt_11_bounds)
      call param_real( lrn_thr_par, 'bt_11_bnds2', number, bt_11_bnds2)
      number = 1
      call param_real( lrn_thr_par, 'pn65_11',    number, pn65_11   )
      call param_real( lrn_thr_par, 'pn13_11',    number, pn13_11   )
      call param_real( lrn_thr_par, 'pn7_11',    number, pn7_11   )
      call param_real( lrn_thr_par, 'pns11_12hi', number, pns11_12hi )
      call param_real( lrn_thr_par, 'pn11_12adj', number, pn11_12adj )

! ... Nighttime snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'ns11_4lo', number, ns11_4lo )
      call param_real( lrn_thr_par, 'ns4_12hi', number, ns4_12hi )
      call param_real( lrn_thr_par, 'nsco2',    number, nsco2    )
      call param_real( lrn_thr_par, 'nsh20',    number, nsh20    )
      number = 1
      call param_real( lrn_thr_par, 'n65_11',    number, n65_11   )
      call param_real( lrn_thr_par, 'ns11_12hi', number, ns11_12hi )
      call param_real( lrn_thr_par, 'ns11_12adj', number, ns11_12adj )

! ... Daytime ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'do11_12hi', number, do11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'do11_4lo',  number, do11_4lo  )
      call param_real( lrn_thr_par, 'dobt11',    number, dobt11    )
      call param_real( lrn_thr_par, 'doco2',     number, doco2     )
      call param_real( lrn_thr_par, 'doh20',     number, doh20     )
      call param_real( lrn_thr_par, 'doref2',    number, doref2    )
      call param_real( lrn_thr_par, 'doref3',    number, doref3    )
      call param_real( lrn_thr_par, 'dovrathi',  number, dovrathi  )
      call param_real( lrn_thr_par, 'dovratlo',  number, dovratlo  )
      number = 2
      call param_real( lrn_thr_par, 'dotci',    number, dotci    )

! ... Daytime polar ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pdo11_12hi', number, pdo11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'pdo11_4lo',  number, pdo11_4lo  )
      call param_real( lrn_thr_par, 'pdobt11',    number, pdobt11    )
      call param_real( lrn_thr_par, 'pdoh20',     number, pdoh20     )
      call param_real( lrn_thr_par, 'pdoref2',    number, pdoref2    )
      call param_real( lrn_thr_par, 'pdoref3',    number, pdoref3    )
      call param_real( lrn_thr_par, 'pdovrathi',  number, pdovrathi  )
      call param_real( lrn_thr_par, 'pdovratlo',  number, pdovratlo  )
      number = 2
      call param_real( lrn_thr_par, 'pdotci',    number, pdotci    )

! ... Nighttime ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'no11_12hi', number, no11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'no11_4lo',  number, no11_4lo  )
      call param_real( lrn_thr_par, 'nobt11',    number, nobt11    )
      call param_real( lrn_thr_par, 'noco2',     number, noco2     )
      call param_real( lrn_thr_par, 'noh20',     number, noh20     )
      call param_real( lrn_thr_par, 'no86_73',   number, no86_73   )
      call param_real( lrn_thr_par, 'no_11var',   number, no_11var  )

! ... Nighttime polar ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pno11_12hi', number, pno11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'pno11_4lo',  number, pno11_4lo  )
      call param_real( lrn_thr_par, 'pnobt11',    number, pnobt11    )
      call param_real( lrn_thr_par, 'pnoh20',     number, pnoh20     )
      call param_real( lrn_thr_par, 'pno86_73',   number, pno86_73   )
      call param_real( lrn_thr_par, 'pno_11var',   number, pno_11var  )

! ... Shadow Thresholds

      number = 2
      call param_real( lrn_thr_par, 'shadnir', number, shadnir )
      number = 1
      call param_real( lrn_thr_par, 'shavrat', number, shavrat )
      call param_real( lrn_thr_par, 'shad124', number, shad124 )

! ... Sun Glint Thresholds

      number = 2
      call param_real( lrn_thr_par, 'snglntv',   number, snglntv   )
      call param_real( lrn_thr_par, 'snglntvch', number, snglntvch )
      call param_real( lrn_thr_par, 'snglntvcl', number, snglntvcl )
      number = 1
      call param_real( lrn_thr_par, 'sg_tbdfl', number, sg_tbdfl )
      call param_real( lrn_thr_par, 'sg_tbdfh', number, sg_tbdfh )
      call param_real( lrn_thr_par, 'snglrat', number, snglrat )
      number = 4
      call param_real( lrn_thr_par, 'snglnt0', number, snglnt0 )
      call param_real( lrn_thr_par, 'snglnt10', number, snglnt10 )
      call param_real( lrn_thr_par, 'snglnt20', number, snglnt20 )
      call param_real( lrn_thr_par, 'snglnt_bounds', number, snglnt_bounds)

! ... Land Restoral Thresholds

      number = 1
      call param_real( lrn_thr_par, 'ldsr5_4_thr', number, ldsr5_4_thr )
      call param_real( lrn_thr_par, 'ldr5_4_thr', number, ldr5_4_thr )
      call param_real( lrn_thr_par, 'ld20m22', number, ld20m22 )
      call param_real( lrn_thr_par, 'ld22m31', number, ld22m31 )
      number = 3
      call param_real( lrn_thr_par, 'ldsbt11', number, ldsbt11 )
      call param_real( lrn_thr_par, 'ldsbt11bd', number, ldsbt11bd )
      call param_real( lrn_thr_par, 'lnbt11', number, lnbt11 )

! ... Day time ocean spatial variability threshold

      number = 1
      call param_real( lrn_thr_par, 'dovar11', number, dovar11 )

! ... Non-cloud Obstruction Thresholds

      number = 1
      call param_real( lrn_thr_par, 'nc_bt37',  number, nc_bt37 )
      call param_real( lrn_thr_par, 'nc37_11',  number, nc37_11 )
      call param_real( lrn_thr_par, 'nc21',     number, nc21 )
      call param_real( lrn_thr_par, 'nc11_12',  number, nc11_12 )
      call param_real( lrn_thr_par, 'ncrat',  number, ncrat )
      call param_real( lrn_thr_par, 'ncvrat',  number, ncvrat )
      call param_real( lrn_thr_par, 'ncsig',  number, ncsig )
     
! ... Snow Mask Thresholds

      number = 1
      call param_real( lrn_thr_par, 'sm_bt11',  number, sm_bt11 )
      call param_real( lrn_thr_par, 'sm_ndsi',  number, sm_ndsi )
      call param_real( lrn_thr_par, 'sm_ref2',  number, sm_ref2 )
      call param_real( lrn_thr_par, 'sm_ref3',  number, sm_ref3 )
      call param_real( lrn_thr_par, 'sm_co2',  number, sm_co2 )
      call param_real( lrn_thr_par, 'sm85_11',  number, sm85_11 )
      call param_real( lrn_thr_par, 'sm37_11',  number, sm37_11 )
      call param_real( lrn_thr_par, 'sm37_11hel',  number, sm37_11hel )
      call param_real( lrn_thr_par, 'sm_ndsi',  number, sm_ndsi )
      call param_real( lrn_thr_par, 'sm_mnir',  number, sm_mnir )

! ... Coast and shallow ocean ndvi thresholds.

      number = 2
      call param_real( lrn_thr_par, 'swc_ndvi', number, swc_ndvi )

! ... pfmft and nfmft thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pfmft_land',  number, pfmft_land  )
      call param_real( lrn_thr_par, 'pfmft_ocean', number, pfmft_ocean )
      call param_real( lrn_thr_par, 'pfmft_snow', number,  pfmft_snow )
      call param_real( lrn_thr_par, 'pfmft_cold', number, pfmft_cold )
      call param_real( lrn_thr_par, 'nfmft_land', number, nfmft_land )
      call param_real( lrn_thr_par, 'nfmft_ocean', number, nfmft_ocean )
      call param_real( lrn_thr_par, 'nfmft_snow', number, nfmft_snow )
      call param_real( lrn_thr_par, 'nfmft_desert', number, nfmft_desert)
      number = 1
      call param_real( lrn_thr_par, 'pfmft_11maxthre', number, pfmft_11maxthre )
      call param_real( lrn_thr_par, 'pfmft_btd_min', number, pfmft_btd_min )
      call param_real( lrn_thr_par, 'nfmft_maxthre', number, nfmft_maxthre )
      
      close(ANC_THRESHOLD_UNIT)

end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--------------------------------------------------------------------
subroutine param_string( PCF_NUM, NAME, STRING )

!-------------------------------------------------------------------
! !F77
!
! !DESCRIPTION:
!     Extract the string corresponding to a parameter
!     name in a parameter file. See function PARAM_READ_FILE for a
!     description of the file format.
!
! !INPUT PARAMETERS:
!     PCF_NUM    PCF number for parameter file
!     NAME       Name of parameter to extract
!
! !OUTPUT PARAMETERS:
!     STRING     String containing the extracted value
!
!
! Example:
!
!      implicit none
!
!      integer pcf_num
!      character*80 string
!      
!      pcf_num = 600300
!      call param_string( pcf_num, 'rcs_id', string )
!      write(*,'(a)') string
!      
!      end
!
! !END
!--------------------------------------------------------------------

!      implicit none
      
! ... Arguments

      integer pcf_num
      character*(*) name
      character*(*) string

! ... Local variables

      integer status
      integer param_max
      parameter( param_max = 10000 )
      character*255 param_list( param_max )
      integer param_num
      character*255 param_value
      integer ios
      character*12 pcf_text
      character*40 temp_string
            
! ... External variables
      
      integer param_read_file !, param_get_value
      external param_read_file !, param_get_value

! ... Read parameter file      

      status = param_read_file( pcf_num, param_max, param_num, param_list )
      !call param_read_file( pcf_num, param_max, param_num, param_list , status)
      if ( status .ne. 0 ) then
        write( pcf_text, '(i12)' ) pcf_num
!        call message( 'param_string',
!     &    'Error reading parameter file PCF#' // pcf_text //
!     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif

! ... Get parameter value

      !status = param_get_value( param_num, param_list, name, param_value )
      call  param_get_value( param_num, param_list, name, param_value, status )
      if ( status .ne. 0 ) then
        write( temp_string, '(a)' ) name( 1 : len( name ) )
!        call message( 'param_string',
!     &    'Error getting parameter value ' // temp_string //
!     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif
  
! ... Read parameter values into array

      read( param_value, '(a)', iostat = ios ) string
      if ( status .ne. 0 ) then
        write( temp_string, '(a)' ) name( 1 : len( name ) )
!        call message( 'param_string',
!     &    'Error reading string for parameter value ' // temp_string //
!     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif
      
      
end subroutine param_string
!--------------------------------------------------------------------
      
!--------------------------------------------------------------------
subroutine param_real( PCF_NUM, NAME, NUMBER, ARRAY )

!-------------------------------------------------------------------
! !DESCRIPTION:
!     Extract the real value(s) corresponding to a parameter
!     name in a parameter file. See function PARAM_READ_FILE for a
!     description of the file format.
!
! !INPUT PARAMETERS:
!     PCF_NUM    PCF number for parameter file
!     NAME       Name of parameter to extract
!     NUMBER     Number of values to extract for this parameter
!                (Maximum number of values in this version is 20).
!
! !OUTPUT PARAMETERS:
!     ARRAY      Array containing NUMBER extracted values
!
! !REVISION HISTORY:
!
! Example:
!
!      implicit none
!
!      integer pcf_num, number, i
!      real array(4)
!      
!      pcf_num = 600300
!      number = 4
!      call param_real( pcf_num, 'pds4_11', number, array )
!      write(*,*) ( array( i ), i = 1, number )
!      
!      end
!
! !END
!--------------------------------------------------------------------

!      implicit none
      
! ... Arguments

      integer pcf_num
      character*(*) name
      integer number
      real array( number )

! ... Local variables

      integer status
      integer, parameter ::  param_max = 10000 
      character*255 param_list( param_max )
      integer param_num
      character*255 param_value
      integer ios
      integer i                  
      character*12 pcf_text
      character*40 string
      
! ... External variables
      
      integer param_read_file !, param_get_value
      external param_read_file !, param_get_value

! ... Read parameter file      

      status = param_read_file( pcf_num, param_max, param_num, param_list )
      !call param_read_file( pcf_num, param_max, param_num, param_list, status )
      if ( status .ne. 0 ) then
        write( pcf_text, '(i12)' ) pcf_num
!        call message( 'param_real',
!     &    'Error reading parameter file PCF#' // pcf_text //
!     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif

! ... Get parameter value

      !status = param_get_value( param_num, param_list, name, param_value )
      call param_get_value( param_num, param_list, name, param_value, status )
      if ( status .ne. 0 ) then
        write( string, '(a)' ) name( 1 : len( name ) )
!        call message( 'param_real',
!     &    'Error getting parameter value ' // string //
!     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif

! ... Check that requested number of values does not exceed the
! ... number allowed by the FORMAT in the next READ statement

      if ( number .gt. 20 ) then
        write( string, '(a)' ) name( 1 : len( name ) )
!        call message( 'param_real',
!     &    'Too many values requested for parameter value ' // string //
!     &    ' [OPERATOR ACTION: Contact SDST]', ios, 2 )
      endif
      
! ... Read parameter values into array

      read( param_value, *, iostat=ios ) ( array( i ), i = 1, number )
      if ( ios .ne. 0 ) then
        write( string, '(a)' ) name( 1 : len( name ) )
!        call message( 'param_real',
!     &    'Error reading data for parameter value ' // string //
!     &    ' [OPERATOR ACTION: Contact SDST]', ios, 2 )
      endif
      
end subroutine param_real
!--------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine param_read_file2( PCF_NUM, PARAM_MAX, PARAM_NUM, PARAM_LIST, status ) 

!-------------------------------------------------------------------
! !DESCRIPTION:
!     Read a parameter file. A parameter file is an ASCII text file
!     containing 1 or more name/value pairs of the form
!
!     NAME : VALUE
!
!     A valid name/value pair must contain
!     - a name containing at least one character,
!     - a colon,
!     - at least one value. More than one value
!     may be defined by using commas to separate values, e.g.
!
!     ANGLES : 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0
!
!     Comments are identified by the '!' character, which may occur
!     at the beginning of a line, or after a name/value pair, thus
!
!     ! This is a comment
!     PI : 3.1415    ! This is also a comment
!
!     are both valid comments. Blank lines are ignored. 
!
! !INPUT PARAMETERS:
!     PCF_NUM       PCF number for parameter file
!     PARAM_MAX     Maximum number of parameters
!                   (dimension of output array PARAM_LIST)
!
! !OUTPUT PARAMETERS:
!     PARAM_NUM     Number of parameters read from FILE
!     PARAM_LIST    Array of parameter strings read from FILE
!
! !END
!--------------------------------------------------------------------
      
!      implicit none
!      save

! --- parameters
      integer 		param_max, status
      character*(*) 	param_list( param_max )
      integer 		pcf_num
      integer 		param_num

! --- internal variables
      character*255 	string
      integer           lun
      integer           param_len
      integer 		count

! --- Set number of parameters found
      param_num = 0

! --- file already opened by file_open 
      lun = pcf_num

! ... Get string length of parameter list
      param_len = len( param_list( 1 ) )

! ... Check that string length of parameter list does not exceed
! ... internal string length

      if ( param_len .gt. len( string ) ) then
        !param_read_file = -2
        status = -2
        return
      endif

! ... Read all lines from the input file, checking that maximum
! ... parameter element number is not exceeded

      count = 0
20    continue
        read( lun, '(a)', end = 40 ) string
        count = count + 1
        if ( count .gt. param_max ) then
          !param_read_file = -3
           status = -3
          return
        endif
        param_list( count ) = string( 1 : param_len )
      goto 20
40    continue

! --- rewind the file
     ! REWIND( lun )

! ... Set return values

      param_num = count
      !param_read_file = 0
      status = 0
      
end subroutine param_read_file2
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine param_get_value( PARAM_NUM, PARAM_LIST, PARAM_NAME, PARAM_VALUE, status )
           
!-------------------------------------------------------------------
! !DESCRIPTION:
!     Get the value string for a named parameter contained in an
!     array of parameter strings read by PARAM_READ_FILE.
!
! !INPUT PARAMETERS:
!     PARAM_NUM     Number of parameters in PARAM_LIST
!     PARAM_LIST    Array of parameter strings
!     PARAM_NAME    Name of parameter to extract
!
! !OUTPUT PARAMETERS:
!     PARAM_VALUE   String containing value for parameter PARAM_NAME
!
! 
! !END
!--------------------------------------------------------------------

      implicit none
      save
      
! ... Input arguments

      integer param_num, status
      character*(*) param_list( param_num )
      character*(*) param_name
      
! ... Output arguments

      character*(*) param_value
      
! ... Local variables

      character*255 name_string, curr_string, string
      integer name_length, curr_length
      integer i
      integer start_pos, end_pos, exc_pos

! ... External functions

      integer strpos, strlen
      EXTERNAL strpos, strlen
            
! ... Set return values

      !param_get_value = -1
      status = -1
      param_value = ' '

! ... If parameter name is empty, return

      if ( len( param_name ) .eq. 0 ) return

! ... Get lowercase compressed version of parameter name

      name_string( 1 : len( name_string ) ) = ' '
      name_string( 1 : len( param_name ) ) = param_name
      call strlower( name_string )
      call strcompress( name_string, .TRUE., name_length )

! ... Loop through parameter list until parameter name is found

      do i = 1, param_num

! ...   Get lowercase compressed version of current parameter      

        curr_string( 1 : len( curr_string ) ) = ' '
        curr_string( 1 : len( param_list( i ) ) ) = param_list( i )
        call strlower( curr_string )
        call strcompress( curr_string, .TRUE., curr_length )

! ...   Check that current parameter is valid

        if ( strpos( curr_string, ':' ) .ge. 1 .and.   &
             curr_string( 1 : 1 ) .ne. ':' .and.       &     
             curr_string( 1 : 1 ) .ne. '!' ) then

! ...     Check if current parameter matches parameter name

          if ( curr_string( 1 : name_length ) .eq.     &
               name_string( 1 : name_length ) ) then

! ...       Get start and end positions of parameter value

            start_pos = strpos( param_list( i ), ':' ) + 1
            end_pos = strlen( param_list( i ) )
            exc_pos = strpos( param_list( i ), '!' )
            if ( exc_pos .ge. 1 ) end_pos = exc_pos - 1

! ...       Extract parameter value and return

            string( 1 : len( string ) ) = ' '
            string = param_list( i )
            param_value( 1 : len( param_value ) ) = ' '
            param_value = string( start_pos : end_pos )
            !param_get_value = 0
            status = 0
            return

          endif
          
        endif

      end do
      
end subroutine param_get_value 
!--------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module thresholds_read_module
