module platform_module

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is to FY3/MERSI-II data.
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

use names_module
use data_arrays_module

implicit none

!---------- 1. fy3 mersi data public variables -------------------------
integer(kind=4), parameter,public       :: num_MERSI_II_Chan_fy3x       = 25   !通道数 total channel number 

!--- modis to mersi_ii [sensor_id = 1]
integer(kind=4), parameter,public       :: num_MERSI_II_Elem_1     = 1354 !列数 1=FY3D-MERSI_II data (convert modis to mersi_II)
integer(kind=4), parameter,public       :: num_MERSI_II_Line_1     = 2030 !行数 1=FY3D-MERSI_II data (convert modis to mersi_II)
integer(kind=4), parameter,public       :: num_vis_Chan_1          = 19   !可见通道数 1 and 2
integer(kind=4), parameter,public       :: num_ir_Chan_1           =  6   !红外通道数 1 and 2

!--- npp/viirs to mersi_ii [sensor_id = 3]
integer(kind=4), parameter,public       :: num_MERSI_II_Elem_2     = 2048 !列数 1=FY3D-MERSI_II data (convert npp/viirs to mersi_II)
integer(kind=4), parameter,public       :: num_MERSI_II_Line_2     = 2000 !行数 1=FY3D-MERSI_II data (convert npp/viirs to mersi_II)
integer(kind=4), parameter,public       :: num_vis_Chan_2          = 19   !可见通道数 1 and 2
integer(kind=4), parameter,public       :: num_ir_Chan_2           =  6   !红外通道数 1 and 2

!--- npp/viirs to mersi_ii [sensor_id = 3]
integer(kind=4), parameter,public       :: num_MERSI_II_Elem_3     = 3354 !列数 1=FY3D-MERSI_II data (convert npp/viirs to mersi_II)
integer(kind=4), parameter,public       :: num_MERSI_II_Line_3     = 3030 !行数 1=FY3D-MERSI_II data (convert npp/viirs to mersi_II)
integer(kind=4), parameter,public       :: num_vis_Chan_3          =  4   !可见通道数 1 and 2
integer(kind=4), parameter,public       :: num_ir_Chan_3           =  5   !红外通道数 1 and 2

!--- real fy3x mersi_ii [sensor_id = 21; 22; 23]
integer(kind=4), parameter,public       :: num_MERSI_II_Elem_fy3x  = 2048 !列数 2=fy3d/MERSI_II(real mersi_II data)
integer(kind=4), parameter,public       :: num_MERSI_II_Line_fy3x  = 2000 !行数 2=fy3d/MERSI_II(real mersi_II data)
integer(kind=4), parameter,public       :: num_vis_Chan_fy3x       = 19   !可见通道数 1 and 2
integer(kind=4), parameter,public       :: num_ir_Chan_fy3x        =  6   !红外通道数 1 and 2

CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_platform_info(ID)

integer(kind=1), intent(in) :: ID 
    
    print*,'**************'
    print*,' fylat Sensor '
    print*,'**************'
    
select case (ID)
  case (1)
    print*,'Sensor_id=1 / FY3D-MERSI_II data (convert modis to mersi_II)'
    sat%nLine = num_MERSI_II_Line_1
    sat%nElem = num_MERSI_II_Elem_1    
    sat%nChan = num_MERSI_II_Chan_fy3x
    sat%nvis  = num_vis_Chan_1
    sat%nir   = num_ir_Chan_1  
    sat%chan_flag = (/     1,     2,     3,     4,     5,     6,     7,      8,      9,     10,     11,     12,     13,     14,     15,     16,     17,     18,    19,           20,           21,           22,           23,           24,           25,   0,   0,   0,   0,   0/)
    sat%midwave   = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38,      3.77726,       4.0610,      7.34409,      8.54966,      11.0170,      12.0360,  0.,  0.,  0.,  0.,  0./)
    sat%midwnum   = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38, 2.647418E+03, 2.462446E+03, 1.361638E+03, 1.169637E+03, 9.076808E+02, 8.308397E+02,  0.,  0.,  0.,  0.,  0./)
!    sat%a         = (/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,   1.0, 9.993438E-01, 9.998729E-01, 9.994894E-01, 9.995439E-01, 9.995483E-01, 9.997404E-01, 1.0, 1.0, 1.0, 1.0, 1.0/)
!    sat%b         = (/   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0, 4.792821E-01, 8.659482E-02, 2.053504E-01, 1.628724E-01, 1.290129E-01, 6.810679E-02, 0.0, 0.0, 0.0, 0.0, 0.0/)
!standard_midwave = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38,         3.80,         4.05,         7.30,         8.55,         10.8,         12.0,   0,   0,   0,   0,   0/)

  case (2)
    print*,'Sensor_id=2 / FY3D-MERSI_II data (convert MODIS to mersi_II in mersi_II format)'
    sat%nLine = num_MERSI_II_Line_2
    sat%nElem = num_MERSI_II_Elem_2   
    sat%nChan = num_MERSI_II_Chan_fy3x 
    sat%nvis  = num_vis_Chan_2
    sat%nir   = num_ir_Chan_2  
    sat%chan_flag = (/     1,     2,     3,     4,     5,     6,     7,      8,      9,     10,     11,     12,     13,     14,     15,     16,     17,     18,    19,           20,           21,           22,           23,           24,           25,   0,   0,   0,   0,  0/)
    sat%midwave   = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38,      3.77726,       4.0610,      7.34409,      8.54966,      11.0170,      12.0360,  0.,  0.,  0.,  0.,  0./)
    sat%midwnum   = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38, 2.647418E+03, 2.462446E+03, 1.361638E+03, 1.169637E+03, 9.076808E+02, 8.308397E+02,  0.,  0.,  0.,  0.,  0./)
!    sat%a         = (/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,   1.0, 9.993438E-01, 9.998729E-01, 9.994894E-01, 9.995439E-01, 9.995483E-01, 9.997404E-01, 1.0, 1.0, 1.0, 1.0, 1.0/)
!    sat%b         = (/   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0, 4.792821E-01, 8.659482E-02, 2.053504E-01, 1.628724E-01, 1.290129E-01, 6.810679E-02, 0.0, 0.0, 0.0, 0.0, 0.0/)
!standard_midwave = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38,         3.80,         4.05,         7.30,         8.55,         10.8,         12.0,   0,   0,   0,   0,   0/)

  case (3)
    print*,'Sensor_id=3 / FY3D-MERSI_II data (convert NPP/VIIRS to mersi_II)'
    sat%nLine = num_MERSI_II_Line_3
    sat%nElem = num_MERSI_II_Elem_3   
    sat%nChan = num_MERSI_II_Chan_fy3x 
    sat%nvis  = num_vis_Chan_3
    sat%nir   = num_ir_Chan_3 
    sat%chan_flag = (/     1,     2,     3,     4,     5,     6,     7,      8,      9,     10,     11,     12,     13,     14,     15,     16,     17,     18,    19,           20,           21,           22,           23,           24,           25,   0,   0,   0,   0,   0/)
    sat%midwave   = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38,         3.80,         4.05,         7.20,         8.55,         10.8,         12.0,   0.,   0.,   0.,   0.,  0./)
!standard_midwave = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38,  3.80,  4.05,  7.20,  8.55,  10.8,  12.0,  0,  0,  0,  0,  0/)
    sat%midwnum   = (/  0.47,  0.55,  0.65,  0.86,  1.24,  1.64,  2.13,  0.412,  0.443,  0.490,  0.555,  0.670,  0.709,  0.746,  0.865,  0.905,  0.936,  0.940,  1.38, 2.647418E+03, 2.462446E+03, 1.361638E+03, 1.169637E+03, 9.076808E+02, 8.308397E+02,  0.,  0.,  0.,  0.,  0./)
!    sat%a         = (/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,   1.0,          1.0,          1.0,          1.0,          1.0,          1.0,          1.0,1.0,1.0,1.0,1.0,1.0/)
!    sat%b         = (/   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,0.0,0.0,0.0,0.0,0.0/) 

  case (21)
    print*,'Sensor_id=21 / FY3D-MERSI_II data (real FY-3D/mersi_II)'
    sat%nLine = num_MERSI_II_Line_fy3x
    sat%nElem = num_MERSI_II_Elem_fy3x  
    sat%nChan = num_MERSI_II_Chan_fy3x   
    sat%nvis  = num_vis_Chan_fy3x
    sat%nir   = num_ir_Chan_fy3x  
    sat%chan_flag = (/       1,       2,       3,       4,       5,       6,       7,        8,       9,       10,       11,       12,       13,       14,       15,       16,       17,       18,      19,            20,           21,           22,           23,          24,          25,  0,  0,  0,  0,  0/)
    sat%midwave   = (/  0.4712,  0.5548,  0.6536,  0.8687,  1.3814,  1.6451,  2.1255,  0.41130,  0.4442,  0.49095,  0.55602,  0.67032,  0.70948,  0.74651,  0.86568,  0.90583,  0.93696,  0.94085,  1.0301,       3.79599,      4.04587,      7.23264,      8.56031,     10.7139,    11.94827,   0.,   0.,   0.,   0.,  0./)
!standard_midwave = (/    0.47,    0.55,    0.65,    0.86,    1.24,    1.64,    2.13,    0.412,   0.443,    0.490,    0.555,    0.670,    0.709,    0.746,    0.865,    0.905,    0.936,    0.940,    1.38,          3.80,         4.05,         7.20,         8.55,        10.8,        12.0,  0,  0,  0,  0,  0/)
    sat%midwnum   = (/  0.4712,  0.5548,  0.6536,  0.8687,  1.3814,  1.6451,  2.1255,  0.41130,  0.4442,  0.49095,  0.55602,  0.67032,  0.70948,  0.74651,  0.86568,  0.90583,  0.93696,  0.94085,  1.0301, 2.6434359E+03, 2.471654E+03, 1.382621E+03, 1.168182E+03, 9.33364E+02, 8.36941E+02,   0.,   0.,   0.,   0.,  0./)
!    sat%a         = (/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,   1.0,          1.0,          1.0,          1.0,          1.0,          1.0,          1.0,1.0,1.0,1.0,1.0,1.0/)
!    sat%b         = (/   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,0.0,0.0,0.0,0.0,0.0/)

  case (22)
    print*,'Sensor_id=22 / FY3E-MERSI_II data (real FY-3E/mersi_II)'
    sat%nLine = num_MERSI_II_Line_fy3x
    sat%nElem = num_MERSI_II_Elem_fy3x  
    sat%nChan = num_MERSI_II_Chan_fy3x   
    sat%nvis  = num_vis_Chan_fy3x
    sat%nir   = num_ir_Chan_fy3x  
    sat%chan_flag = (/       1,       2,       3,       4,       5,       6,       7,        8,       9,       10,       11,       12,       13,       14,       15,       16,       17,       18,      19,            20,           21,           22,           23,          24,          25,  0,  0,  0,  0,  0/)
    sat%midwave   = (/  0.4712,  0.5548,  0.6536,  0.8687,  1.3814,  1.6451,  2.1255,  0.41130,  0.4442,  0.49095,  0.55602,  0.67032,  0.70948,  0.74651,  0.86568,  0.90583,  0.93696,  0.94085,  1.0301,       3.79599,      4.04587,      7.23264,      8.56031,     10.7139,    11.94827,   0.,   0.,   0.,   0.,  0./)
!standard_midwave = (/    0.47,    0.55,    0.65,    0.86,    1.24,    1.64,    2.13,    0.412,   0.443,    0.490,    0.555,    0.670,    0.709,    0.746,    0.865,    0.905,    0.936,    0.940,    1.38,          3.80,         4.05,         7.20,         8.55,        10.8,        12.0,  0,  0,  0,  0,  0/)
    sat%midwnum   = (/  0.4712,  0.5548,  0.6536,  0.8687,  1.3814,  1.6451,  2.1255,  0.41130,  0.4442,  0.49095,  0.55602,  0.67032,  0.70948,  0.74651,  0.86568,  0.90583,  0.93696,  0.94085,  1.0301, 2.6434359E+03, 2.471654E+03, 1.382621E+03, 1.168182E+03, 9.33364E+02, 8.36941E+02,   0.,   0.,   0.,   0.,  0./)
!    sat%a         = (/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,   1.0,          1.0,          1.0,          1.0,          1.0,          1.0,          1.0,1.0,1.0,1.0,1.0,1.0/)
!    sat%b         = (/   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,0.0,0.0,0.0,0.0,0.0/)


end select

print*,'Number of Element: ',sat%nElem
print*,'Number of Line   : ',sat%nLine
print*,'Number of Channel: ',sat%nChan
print*,'Number of VisChan: ',sat%nvis
print*,'Number of IR Chan: ',sat%nir
print*, '  '

end subroutine fylat_platform_info
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module platform_module
