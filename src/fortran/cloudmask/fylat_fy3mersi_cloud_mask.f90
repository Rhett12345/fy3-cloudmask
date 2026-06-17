module fylat_fy3mersi_cloud_mask

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    cloud mask algorithm for fy3/mersi sensor
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
use data_arrays_module
use constant
use thresholds_read_module
use cloudmask_data_arrays
!use io_module
! surface
use polar_module
use land_module
use water_module


implicit none

include 'global.inc'

!+++++++++++++++++++ step 1: Variables +++++++++++++++++++++++++++++++++


contains
!+++++++++++++++++++ step 2: Subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~~~~ subroutine 1: fy3mersi_cloud_mask ~~~~~~~~~~~~~~~~~
subroutine fy3mersi_cloud_mask()

!-----------------------------------------------------------------------
!                          PROCESSING SECTION
!
!        The Cloud Mask is processed pixel by pixel using a sliding
!        box approach. Regional context data is stored for use
!        in uniformity tests for uncertain pixels.  The edge of 
!        the area in which you are processing, (outline of region)
!        will be processed but will not include uniformity tests.
!----------------------------------------------------------------------



integer(kind=4) :: iline,ielem, i, j, k, ix_nwp, iy_nwp

real(kind=4)    :: pelev
integer(kind=1) :: eco_type, is_cold_sfc
!integer(kind=1) :: lsf
!real(kind=4)    :: lon,lat
real(kind=4)    :: satzen,tbadj
real(kind=4)    :: solzen
real(kind=4)    :: rela
real(kind=4)    :: glint
!real(kind=4)    :: u_wind
!real(kind=4)    :: v_wind
!real(kind=4)    :: tpw
integer(kind=4) :: xnwp,ynwp, i_sta,i_end,j_sta,j_end
real(kind=4),dimension(inband)     :: pxldat
real(kind=4),dimension(7)          :: btclr
real(kind=4),dimension(3,3,inband) :: indat
byte, dimension(:,:,:), pointer    :: qa
byte, dimension(:,:,:), pointer    :: cm

print*,'    ... fylat retrieve fy3/MERSI_II Cloud Mask start !!! '
!======================================================================
! read threshold file
!======================================================================
call thresholds_read( fylat_sensor_id )

!--- set local POINTERs to output structures
cm => cm_bitarray
qa => cm_qa_bitarray

allocate ( out_pwater (sat%nElem,sat%nLine) )            !lyj
allocate ( out_sfctmp (sat%nElem,sat%nLine) )            !lyj
allocate ( out_polar  (sat%nElem,sat%nLine) )            !lyj
allocate ( out_day    (sat%nElem,sat%nLine) )            !lyj
allocate ( out_night  (sat%nElem,sat%nLine) )            !lyj
allocate ( out_land   (sat%nElem,sat%nLine) )            !lyj
allocate ( out_water  (sat%nElem,sat%nLine) )            !lyj
allocate ( out_coast  (sat%nElem,sat%nLine) )            !lyj
allocate ( out_snglnt (sat%nElem,sat%nLine) )            !lyj
allocate ( out_snow   (sat%nElem,sat%nLine) )            !lyj
allocate ( out_ice    (sat%nElem,sat%nLine) )            !lyj
allocate ( out_desert (sat%nElem,sat%nLine) )            !lyj
allocate ( out_uniform(sat%nElem,sat%nLine) )            !lyj
allocate ( out_shadow (sat%nElem,sat%nLine) )            !lyj

out_pwater  = -999.0        !lyj
out_sfctmp  = -999.0        !lyj
out_polar   = -9            !lyj
out_day     = -9            !lyj
out_night   = -9            !lyj
out_land    = -9            !lyj
out_water   = -9            !lyj
out_coast   = -9            !lyj
out_snglnt  = -9            !lyj
out_snow    = -9            !lyj
out_ice     = -9            !lyj
out_desert  = -9            !lyj
out_uniform = -9            !lyj
out_shadow  = -9            !lyj

!======================================================================
! Loop over pixels in this segment
!======================================================================
!     open(1,file='test.txt')        !jincheng
line_loop_1: do iline= 1, sat%nLine
!line_loop_1: do iline= 200, 700                  !jincheng test
element_loop_1: do ielem= 1, sat%nElem
!element_loop_1: do ielem= 900, 1300               !jincheng test
!print*,'ij',iline,ielem         !jincheng
!open(1,file='test.txt')         !jincheng
!--- Initialize regional variables
    call pxinit(testbits,qa_bits,precip_water,vza,                 &
                sfctmp,pmsl,u_wind,v_wind,plat,plon,lsf,polar,     &
                day,night,land,water,coast,snglnt,visusd,vrused,   &
                snow,ice,desert,bad_value,bad_geo,uniform,         &
                shadow,smoke,cirrus_ir,cirrus_vis,nmtests,         &
                nbands,nbad_1km,nbad_250,hi_elev,antarctic,        &
                sh_ocean,sg_bad_data,map_ice,map_snow,sh_lake)
     
!--- define aliases
!   add clear sky bt
if (fylat_rtm_opt > 0 ) then
    ix_nwp   = sat%x_nwp(ielem,iline)
    iy_nwp   = sat%y_nwp(ielem,iline)
    btclr(1) = sat%bt_clr38(ix_nwp,iy_nwp)
    btclr(2) = sat%bt_clr40(ix_nwp,iy_nwp)
    btclr(3) = sat%bt_clr73(ix_nwp,iy_nwp)
    btclr(4) = sat%bt_clr86(ix_nwp,iy_nwp)
    btclr(5) = sat%bt_clr11(ix_nwp,iy_nwp)
    btclr(6) = sat%bt_clr12(ix_nwp,iy_nwp)
    btclr(7) = 0.0
endif
!   add clear sky bt
    pelev    = geo%dem(ielem,iline)
    lsf      = geo%lsm(ielem,iline)
    plon     = geo%lon(ielem,iline)
    plat     = geo%lat(ielem,iline)
    satzen   = geo%SensorZenith(ielem,iline)
    solzen   = geo%SolarZenith(ielem,iline)
    vza      = satzen 
    rela     = geo%RelAzimuth(ielem,iline)
    glint    = geo%GlintAngle(ielem,iline)
    refang   = glint
    eco_type = sat%eco(ielem,iline)
    ! data from nwp
    xnwp     = sat%x_nwp(ielem,iline)           !nwp longitude cell for org ncep
    ynwp     = sat%y_nwp(ielem,iline)           !nwp latitude cell
    if ( fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4 .or. fylat_nwp_opt == 5 ) then
       u_wind        = nwp26%u_sigma(xnwp,ynwp) 
       v_wind        = nwp26%v_sigma(xnwp,ynwp)
       precip_water  = nwp26%tpw(xnwp,ynwp)
       sfctmp        = nwp26%tsfc(xnwp,ynwp)
       pmsl          = nwp26%pmsl(xnwp,ynwp)
    endif
    if ( fylat_nwp_opt == 3 ) then
       u_wind        = nwp36%u_sigma(xnwp,ynwp) 
       v_wind        = nwp36%v_sigma(xnwp,ynwp)
       precip_water  = nwp36%tpw(xnwp,ynwp)
       sfctmp        = nwp36%tsfc(xnwp,ynwp) 
       pmsl          = nwp36%pmsl(xnwp,ynwp)
    endif

!--- define surface type              
    call get_pxldat(ielem,iline,pelev,pxldat,satzen,solzen,rela,glint,eco_type,u_wind,v_wind,precip_water,tbadj) ! revised by minmin 
 
!--  get water sfc temp
    if ((water .or. coast) .and. sat%sst(ielem,iline)>100.0) then
        sfctmp  = sat%sst(ielem,iline)
    endif
    
!--- check edge
    call chk_ele_lin_edge(ielem,iline,line_edge,ele_edge) ! revised by minmin 

!--- get indat
    indat = -999.0  
    call get_pxl3X3(ielem,iline,i_sta,i_end,j_sta,j_end)
    do k = 1, 19 
       indat(1:3,1:3,k) = sat%ref_vis(i_sta:i_end,j_sta:j_end,k)
    end do
    do k = 20, 25
       indat(1:3,1:3,k) = sat%tbb_ir(i_sta:i_end,j_sta:j_end,k-19)
    end do
    
!--- check uniformity
    call check_reg_uniformity(ielem,iline,line_edge,ele_edge, &
                              eco_type,day,land,water,        &
                              coast,snow,ice,uniform)   ! revised by minmin

!   Decision to process or not to process this pixel.
    if (sg_bad_data .and. water .and. snglnt) then
        process = .false.
    else
        process = .true.
    end if

!   Decision to cold surface to process this pixel.
    is_cold_sfc = 0
    if (sfctmp < 265.0 ) then
        is_cold_sfc = 1
    end if

    !if (iline ==1800 .and. ielem == 1100) then
    !   print*,'desert,hi_elev',desert,hi_elev
    !endif

    if(process) then
    !print*,'g',ielem,iline,polar,land,water,day

!    if (plat > 23 .and. plat < 25 .and. plon > 13 .and. plon < 15)  then !jincheng test
   
!    Decision tree for processing paths 
      if (polar) then

!        First polar processing

         if (day) then
!            Daytime processing
!print*,land,ice,snow,desert,coast,eco_type,uniform
             call polar_day(pxldat,vza,snglnt,                            &
                            visusd,refang,vrused,cirrus_vis,              &
                            land,ice,snow,desert,coast,                   &
                            eco_type,uniform,hi_elev,ielem,indat,         &
                            nmtests,testbits,tbadj,antarctic,             &
                            sh_ocean,sfctmp,qa_bits,confdnc,              &
                            btclr,is_cold_sfc)
                   
         else
!            Nighttime processing
             call polar_nite(pxldat,vza,land,ice,snow,desert,             &
                             hi_elev,sfctmp,eco_type,nmtests,testbits,    &
                             uniform,indat,ielem,antarctic,sh_ocean,      &
                             qa_bits,confdnc,                             &
                             btclr,is_cold_sfc)
         endif
                                      
      else if (land) then

!        Primarily land surface.

         if (day) then

!            Daytime processing.
             call land_day(pxldat,vza,visusd,vrused,           &
                           cirrus_vis,desert,coast,snow,       &
                           ice,hi_elev,tbadj,eco_type,         &
                           testbits,qa_bits,nmtests,confdnc,   &
                           btclr,is_cold_sfc)

         else
!            Nighttime processing.
             call land_nite(pxldat,plat,vza,ice,snow,coast,tbadj,          &
                            desert,hi_elev,sh_lake,sfctmp,eco_type,        &
                            nmtests,testbits,qa_bits,confdnc,precip_water, &
                            btclr,is_cold_sfc)
         endif

      else if (water) then

!          Primarily water surface.

         if (day) then
!            Daytime processing.
             call water_day(pxldat,vza,snglnt,visusd,refang,cirrus_vis,     &
                            sfctmp,hi_elev,uniform,ice,snow,ielem,iline,    &
                            line_edge,sh_ocean,indat,nmtests,               &
                            testbits,qa_bits,confdnc,btclr)
         else
!            Nighttime processing.
             call water_nite(pxldat,vza,uniform,ice,snow,indat,sfctmp,      &
                             sh_ocean,ielem,nmtests,testbits,qa_bits,       &
                             confdnc,btclr)
         end if

      end if

!     Test for shadows,if necessary. Set bit to indicate no shadow was found.
      if(.not.water .and. .not.coast .and. day .and.  &
         .not.polar .and. confdnc.ge.0.66) then
         call shadows(pxldat,shadow,visusd,qa_bits)
      end if

!     Test for possible non-cloud obstruction.
      if(land .and. day .and. .not. snow) then
         call noncld_obs_chk(indat,pxldat,confdnc,ielem,           &
                             line_edge,iline,qa_bits,              &
                             testbits,smoke)
      end if

!     Test of thin cirrus in the infrared
      if ( (.not. snow) .and. (.not. ice) ) then
         call thin_ci_chk_ir(pxldat,vza,cirrus_ir,qa_bits,testbits)
      end if

!     Set bits which indicate processing path through algorithm.
      call proc_path(water,land,day,ice,snow,snglnt,coast,   &
                     desert,smoke,shadow,testbits)

!     Set the bits which are not used in output array.
      call set_unused_bits(testbits)

!     Get cloud mask statistics
!      call get_stats(day,land,water,snglnt,snow,coast,desert,
!                     ice,shadow,smoke,cirrus_ir,cirrus_vis,
!                     confdnc,nmtests,nbands,bad_value,bad_geo,
!                     geo_flag,no,ns,nd,nt,ng,
!                     ni,nl,nw,nr,nv,nu,nn,nz,na,
!                     ne,npix,n1sm,n2sm,n3sm,n4sm,num_bad)
 
!     Set cloud mask quality bit flags
      call set_confdnc(confdnc,testbits)
     ! if (confdnc>1 .or. confdnc<0) then
     ! print*,'ff',ielem,iline,water,day,confdnc,testbits(1)
     ! endif

!     Set remaining QA bits 
      call set_quality_A(nmtests,nbands,lsf,qa_bits)

!     Store element in line array 
!               call fill_bit_line(nc,nmtests,nbands,bad_value,bad_geo,
!    +                             snglnt,desert,testbits,qa_bits,
!    +                             bitarray,qa_bitarray)
      call fill_bit_pixel(nmtests,nbands,bad_value,bad_geo,        &
                          snglnt,desert,testbits,qa_bits,          &
                          cm(ielem,iline,1:6),                     &
                          qa(ielem,iline,1:10))

!               call store_cm_scan(nc,confdnc,day,visusd,snglnt,snow,ice,
!    +                             coast,desert,land,process,
!    +                             scan_confdnc,scan_visusd,scan_snglnt,
!    +                             scan_snow,scan_ice,scan_coast,scan_desert,
!    +                             scan_land,scan_process,scan_day)
    

            
!    end if ! jincheng test
    end if !if(process) then

    out_pwater(ielem,iline) = precip_water                       !lyj
    out_sfctmp(ielem,iline) = sfctmp                             !lyj

    if ( polar   )  out_polar   (ielem,iline)  = 1               !lyj
    if ( day     )  out_day     (ielem,iline)  = 1               !lyj
    if ( night   )  out_night   (ielem,iline)  = 1               !lyj
    if ( land    )  out_land    (ielem,iline)  = 1               !lyj
    if ( water   )  out_water   (ielem,iline)  = 1               !lyj
    if ( coast   )  out_coast   (ielem,iline)  = 1               !lyj
    if ( snglnt  )  out_snglnt  (ielem,iline)  = 1               !lyj
    if ( snow    )  out_snow    (ielem,iline)  = 1               !lyj
    if ( ice     )  out_ice     (ielem,iline)  = 1               !lyj
    if ( desert  )  out_desert  (ielem,iline)  = 1               !lyj
    if ( uniform )  out_uniform (ielem,iline)  = 1               !lyj
    if ( shadow  )  out_shadow  (ielem,iline)  = 1               !lyj
    if ( .not.  polar   )  out_polar   (ielem,iline) = 0         !lyj
    if ( .not.  day     )  out_day     (ielem,iline) = 0         !lyj
    if ( .not.  night   )  out_night   (ielem,iline) = 0         !lyj
    if ( .not.  land    )  out_land    (ielem,iline) = 0         !lyj
    if ( .not.  water   )  out_water   (ielem,iline) = 0         !lyj
    if ( .not.  coast   )  out_coast   (ielem,iline) = 0         !lyj
    if ( .not.  snglnt  )  out_snglnt  (ielem,iline) = 0         !lyj
    if ( .not.  snow    )  out_snow    (ielem,iline) = 0         !lyj
    if ( .not.  ice     )  out_ice     (ielem,iline) = 0         !lyj
    if ( .not.  desert  )  out_desert  (ielem,iline) = 0         !lyj
    if ( .not.  uniform )  out_uniform (ielem,iline) = 0         !lyj
    if ( .not.  shadow  )  out_shadow  (ielem,iline) = 0         !lyj
    
!-----------------------------------------------------------------------
! END loop over pixels in segment
!-----------------------------------------------------------------------
end do element_loop_1
    
end do line_loop_1
!close(1)             !jincheng
cm => null()
qa => null()
   
end subroutine fy3mersi_cloud_mask
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 2: check surface type ~~~~~~~~~~~~~~~~~~
subroutine get_pxldat(i,j,pelev,pxldat0,satzen,solzen,rela,glint,eco_type,u_wind,v_wind,tpw,tbadj)

integer(kind=4) :: i,j,k
real(kind=4)    :: pelev
integer(kind=1) :: eco_type
!integer(kind=1) :: lsf
!real(kind=4)    :: plon,plat
real(kind=4)    :: satzen,tbadj
real(kind=4)    :: solzen,sza
real(kind=4)    :: rela
real(kind=4)    :: glint
real(kind=4)    :: u_wind
real(kind=4)    :: v_wind
real(kind=4)    :: tpw
real(kind=4)    :: ndvi
real(kind=4)    :: pxldat0(inband)
real(kind=4)    :: vrat
real(kind=4)    :: m05to19 ! added by minmin for correct real fy3d mersi-ii 1.38 channel
integer(kind=1) :: sim
integer(kind=1) :: sg_band(sg_bands_used)
!data sg_band /1,2,20,26,27,31,32,35/
data sg_band /3,4,20,24,25/
!data sg_band /3,4/
!pxldat02(1) = sat%ref_vis(i,j, 3)  !0.64
!pxldat02(2) = sat%ref_vis(i,j, 4)  !0.86
!pxldat02(3) = sat%ref_vis(i,j, 2)  !0.55
!pxldat02(4) = sat%ref_vis(i,j,19)  !1.38
!pxldat02(5) = sat%ref_vis(i,j, 7)  !2.25
!pxldat02(6) = sat%tbb_ir(i,j,1)    !3.80
!pxldat02(7) = sat%tbb_ir(i,j,4)    !8.50
!pxldat02(8) = sat%tbb_ir(i,j,5)    !11.0 
sim = sat%snow_mask(i,j)

!     First get 1km channels - check for bad data

      do  k = 1, sat%nChan

        if (k .le. 19) then
        
          if (sat%ref_vis(i,j,k) .gt. -99.0 .and. sat%ref_vis(i,j,k) .le. 2.3) then
            pxldat0(k) = sat%ref_vis(i,j,k) 
          else
            pxldat0(k) = bad_data
          end if

        else

          if (sat%tbb_ir(i,j,k-sat%nvis) .gt. 0.0 .and. sat%tbb_ir(i,j,k-sat%nvis) .lt. 1000.0) then
            pxldat0(k) = sat%tbb_ir(i,j,k-sat%nvis) 
          else
            pxldat0(k) = bad_data
          endif

        endif

      end do
      
      if (fylat_sensor_id > 20) then ! for real fy3d-mersi-ii 1.38 channel
         m05to19     = pxldat0(5) 
         pxldat0(5)  = pxldat0(19)
         pxldat0(19) = m05to19
      endif
      
      
! ... Loop to count number of good bands out of those used
      do k = 1 , sat%nChan
        if (nint(pxldat0(k)) .ne. nint(bad_data)) then
          nbands = nbands + 1
        else
          nbad_1km = nbad_1km + 1
          bad_value = .true.
        endif
      end do
      
!     get pixel latitude and longitude
      if (nint(plat) .eq. nint(bad_data)) then
         geo%flag(i,j,1) = 1
         bad_geo = .true.
      endif
      if (nint(plon) .eq. nint(bad_data)) then
         geo%flag(i,j,2) = 1
         bad_geo = .true.
      endif

! ... Check bands needed for sun-glint processing.  If any are missing or bad, set flag.
      do k = 1 , sg_bands_used
        if (nint(pxldat0(sg_band(k))) .eq. nint(bad_data)) then
          sg_bad_data = .true.
        endif
      end do

!     Get 11 um brightness temperature threshold adjustment for deserts.
      tbadj = (pelev / 1000.0) * 5.0
                         
!     Define "Antarctic" flags.
      if(plat .lt. -60.0) then
        antarctic = .true.
      end if
      
!     Determine whether or not current pixel will be processed as
!     desert.

!     Determine whether or not current pixel will be processed as desert.

!     Global.
      if(eco_type .eq.  8 .or. eco_type .eq. 46 .or.  &
         eco_type .eq. 50 .or. eco_type .eq. 51 .or.  &
         eco_type .eq. 59 .or. eco_type .eq. 71 .or.  &
         eco_type .eq. 11 .or. eco_type .eq. 9 .or.   &
         eco_type .eq. 52) then
        desert = .true.

!     High-elevation.
      else if( (pelev .gt. 2000.0) .and. (eco_type .eq. 42) ) then
        desert = .true.
!       Check for locations where there should be no high-elevation desert.
        if( (plat .le. 10.0 .and. plat .ge. -10.0) .and. (plon .ge. 90.0) ) then
          desert = .false.
        else if( (plat .ge. -30.0 .and. plat .le. -10.0) .and.   &
                 (plon .ge. 160.0 .and. plon .le. 180.0) ) then
          desert = .false.
        else if( (plat .ge. 10.0 .and. plat .le. 26.0) .and.     &
                 (plon .ge. 120.0 .and. plon .le. 180.0) ) then
          desert = .false.
        end if

!     Africa.
      else if(plat .le. 20.0 .and. plat .ge. -35.0 .and. plon .le. 60.0 .and. plon .ge. -20.0) then
        if(eco_type .eq.  7 .or. eco_type .eq. 41 .or.   &
           eco_type .eq. 43 .or. eco_type .eq. 58 .or.   &
           eco_type .eq. 36 .or. eco_type .eq. 91 .or.   &
           eco_type .eq. 32 .or. eco_type .eq. 29) then
           desert = .true.
        end if

!     Eurasia.
      else if(plat .le. 70.0 .and. plat .ge. -60.0 .and. plon .le. 180.0 .and. plon .ge. -20.0) then
        if(eco_type .eq. 11 .or. eco_type .eq. 2) then
           desert = .true.
        end if
!       Exclude New Zealand.
        if(plat .ge. -50.0 .and. plat .le. -30.0 .and. plon .ge. 160.0 .and. plon .le. 180.0) then
          desert = .false.
        end if

      end if
!     Add in Australia.
      if(plat .le. -11.0 .and. plat .gt. -40.0 .and. plon .le. 155.0 .and. plon .gt. 110.0) then
        if(eco_type .eq. 43 .or. eco_type .eq. 41 .or. eco_type .eq. 91) then
          desert = .true.
        end if
      end if

!     Determine whether or not visible ratio test may be used over
!     land surfaces.
      if(eco_type .eq. 2 .or. eco_type .eq. 8 .or.     &
         eco_type .eq. 11 .or. eco_type .eq. 40 .or.   &
         eco_type .eq. 41 .or. eco_type .eq. 46 .or.   &
         eco_type .eq. 51 .or. eco_type .eq. 52 .or.   &
         eco_type .eq. 59 .or. eco_type .eq. 71 .or.   &
         eco_type .eq. 50) then
        vrused = .false.
      endif   

!     Now we turn our attention to the ancillary data sets
!     get logical flag variables for given pixel
!     Day/night flag (current definition is > 85 degrees)
      sza = solzen
      if (sza .lt. 0.0) then
         bad_geo = .true.
         geo%flag(i,j,3) = 1
      else if (sza .gt. 85.0) then
         night = .true.
      else
         day = .true.
      endif

!     set polar flag (if lat is higher then +/-60)
      if (plat .gt. 60.0 .or. plat .lt. -60.0) then
         polar = .true.
      endif

!     set visusd flag (currently set to true if szen < 85 degrees)
      if (sza .gt. 85.0 ) then
         visusd = .false.
      endif

!     set the sunglnt flag
      if (nint(glint) .eq. 999.0) then
         geo%flag(i,j,5) = 1
      elseif (glint .le. 36.0) then
         snglnt = .true.
      endif   

!     set land/sea flag
! ... First make sure that it is not a missing value  [1]
      if (lsf .ne. -1) then
        if (lsf .eq. 1 .or. lsf .eq. 4) then
          land = .true.
          if(eco_type .eq. 14) then
!           Fix-up for missing ecosystem data in eastern Greenland and
!           north-eastern Siberia.  Without this, these regions become
!           completely "coast".
            if( (plat .lt. 64.0) .or. (plat .ge. 67.5 .and.      &
                ((plon .lt. -40.0 .and. plon .gt. -168.5) .or.   &
                                       plon .gt. -12.5)) .or.    &
                ((plat .ge. 64.0 .and. plat .lt. 67.5) .and.     &
                ((plon .lt. -40.0 .and. plon .gt. -168.5) .or.   &
                                        plon .gt. -30.0)) ) then
              coast = .true.
            end if
          end if
        elseif (lsf .eq. 2) then
          coast = .true.
          land = .true.
        elseif (lsf .eq. 3) then
          land = .true.
          sh_lake = .true.
!         Need shallow lakes to be processed as "coast" for day, but
!         not night.
          if(day) then
            coast = .true.
          end if
        else 
          water = .true.
          if(lsf .eq. 0) then
            sh_ocean = .true.
          end if
        endif        

!     If land/sea flag is missing, then calculate visible ratio to 
!     determine if land or water.

      elseif (nint(pxldat0(3)) .ne. nint(bad_data) .and.    &    ![2]
              nint(pxldat0(4)) .ne. nint(bad_data)) then 
        vrat = pxldat0(4) /pxldat0(3)
        if (vrat .gt. 0.9) then 
          land = .true.
        else
          water = .true.
        endif

! ... If all else fails, call it land.

      else
        land = .true.
        water = .false.
      endif    

!     Get pixel value of SST or model surface temp
!      if(land) then
!        sfctmp = contx_sfctmp(j,i)
!      else
!        sfctmp = contx_sst(j,i)
!      end if  

!     Set high elevation flag.
!     First, define "Greenland".
      Greenland = .false.
      if(land) then
        if(plat .ge. 60.0 .and. plat .lt. 67.0) then
          if(plon .ge. -60.0 .and. plon .lt. -30.0) then
            Greenland = .true.
          end if
        else if(plat .ge. 67.0 .and. plat .lt. 75.0) then
          if(plon .ge. -60.0 .and. plon .lt. -10.0) then
            Greenland = .true.
          end if
        else if(plat .ge. 75.0) then
          if(plon .ge. -70.0 .and. plon .lt. -10.0) then
            Greenland = .true.
          end if
        end if
      end if
      if( (pelev .gt. 2000.0) .or. (pelev .gt. 200.0 .and. Greenland  &
          .and. land) .or. (plat .ge. 75.7 .and. plat .le. 79.0       &
          .and. plon .ge. -73.0 .and. plon .le. -50.0 .and. land) ) then
        hi_elev = .true.
      end if

! ... Calculate raw NDVI value.
      if(nint(pxldat0(3)) .ne. nint(bad_data) .and. nint(pxldat0(4)) .ne. nint(bad_data)) then
        ndvi = (pxldat0(4)-pxldat0(3)) / (pxldat0(4)+pxldat0(3))
      end if

! ... Set the ancillary ice flag.
!      if ( sim .eq. sym%SEA_ICE .and. water ) then
!         map_ice = .true.
!      end if
!      
!      if ( sim .eq. sym%SNOW .and. water .and. plat .lt. -60.0) then
!         map_ice = .true.
!      end if
!      
! ... Set the ancillary snow flag.
!      if ( sim .eq. sym%SNOW ) then
!         map_snow = .true.
!      end if
!      if ( sim .eq. sym%SNOW .and. land) then
!         map_snow = .true.
!      end if
!      if( (sim .eq. sym%SNOW .and. land) .and. (plat .ge. 44.0 .or. plat .le. -40.0) ) then 
!         map_snow = .true.
!      endif

! ... Set the ancillary ice flag.
      if (sim .gt. 25 .and. sim .lt. 102 .and. water ) then
          map_ice = .true.
      endif
      if (sim .eq. 200 .and. water .and. plat .lt. -60.0) then
          map_ice = .true.
      endif

! ... Set the ancillary snow flag.
      if (sim .eq. 103 .or. sim .eq. 104 .or. sim .eq. 101) then
         map_snow = .true.
      endif
      if (sim.gt. 25 .and. sim .lt. 101 .and. land) then 
         map_snow = .true.
      endif
      if( (sim .eq. 200 .and. land) .and. (plat .ge. 44.0 .or. plat .le. -40.0) )  then
         map_snow = .true.
      endif

! ****************************************************************

      if(day) then

!       Run quick version of D.Hall's snow detection algorithm.
        call snow_mask(pxldat0,plat,land,snglnt,water,hi_elev,Greenland,ndsi_snow)

        if(water) then

          if(plat .ge. -60.0 .and. plat .le. 25.0) then
            if(map_ice .and. ndsi_snow) then
              ice = .true.
            end if
          else if( (sim .eq. 252 .or. sim .eq. 200 )  &
                 .and. (plat .ge. 44.0 .or. plat .le. -40.0) ) then
            if(ndsi_snow) then
              ice = .true.
            end if
          else if (plat .ge. 68.0 .or. plat .le. -60.0) then    ! added by minmin 20180424 for polar region
            if(ndsi_snow) then
              ice = .true.
            end if
          else if((lsf .eq. 3 .or. lsf .eq. 5) .and. ndsi_snow) then
            ice = .true.
          else if(map_ice .and. ndsi_snow) then
            ice = .true.
          endif

        else if(land) then

          if(plat .ge. -60.0 .and. plat .le. 25.0) then
!           Define New Zealand region which receives snow but snow
!           map does not show it.
            if( (plat .ge. -48.0 .and. plat .le. -34.0) .and.  &
                (plon .ge. 165.0 .and. plon .le. 180.0) ) then
              New_Zealand = .true.
            else
              New_Zealand = .false.
            end if
            if(map_snow .and. ndsi_snow) then
              snow = .true.
            else if(New_Zealand .and. ndsi_snow) then
              snow = .true.
            end if
          else if(plat .lt. -60.0) then
             snow = .true.
          else
            if(ndsi_snow) then
              snow = .true.
            end if
          endif

        endif

      else

        ice = map_ice
        snow = map_snow

      endif       
                  
end subroutine get_pxldat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine snow_mask(pxldat0,     &
                     plat,       &
                     land,       &
                     snglnt,     &
                     water,      &
                     hi_elev,    &
                     Greenland,  &
                     ndsi_snow)
      

!      save

! ... Common statement for debug purposes
!      common / bug / debug, *

!---------------------------------------------------------------------
!
!!Description:
!       Subroutine which implements a quick and dirty version of
!       the snow algorithm by Dorothy Hall, George Riggs and 
!       Vince Salmonson.
!
!!Input Parameters:
! pxldat0        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! plat          Latitude of current pixel
! land          Flag indicating a pixel to be processed as land
! snglnt        Logical variable flagging background as sunglint 
!               contaminated
! water         Logical variable indicating water surface for current 
!               pixel
! hi_elev       Logical flag indicating elevation > 2000 meters
! Greenland     Logical flag indicating location is "Greenland"
!
!!Output Parameters:
! ndsi-snow     logical variable indicating the prescence of snow in a
!               FOV according to the NDSI test
!
!!Revision History:
! 10/04  Collection 5b  R. Frey
! Removed definition of 'Greenland'.
! Removed tests on 1.38 um, 'sm_mnir', 3.7-11 um BTD in Greenland.
! Added 1.5K to 8.5-11 um BTD test threshold in Greenland.
!
!!Team-unique Header:
!
!!References and Credits:
! See snow product ATBD.
!
!!END
!---------------------------------------------------------------------

      include 'global.inc'
      include 'snow_mask.inc'
      !include 'platform_name.inc'
!
! ... scalar arguments ..
      logical snglnt,water,ndsi_snow,hi_elev,land,Greenland
      real plat
! ...
! ... array arguments ..
      real pxldat0(inband)
! ...
! ... local scalars ..
      real ndsi,masv88,masv55,masv188,masir13,masir11,     &
           masir85,masir37,diff,sth37_11,masnir,diff2,sth85_11
      integer debug,I_bad
!
! ... identification 
!      masv55 = pxldat0(4)     !0.55
!      masv88 = pxldat0(2)     !0.86
!      masv188 = pxldat0(26)   !1.38
!      masir37 = pxldat0(20)   !3.8
!      masir85 = pxldat0(29)   !8.5
!      masir11 = pxldat0(31)   !11.0
!      masir13 = pxldat0(35)   !13.9
      masv55 = pxldat0(2)      !0.55
      masv88 = pxldat0(4)      !0.86
      masv188 = pxldat0(19)    !1.38
      masir37 = pxldat0(20)    !3.8
      masir85 = pxldat0(23)    !8.5
      masir11 = pxldat0(24)    !11.0      
!     Select NIR channel for NDSI based on platform name.
!     Band 6 for Terra, band 7 for Aqua.
      !if(platform_name .eq. 'Terra') then
      !  masnir = pxldat0(6)
      !else if(platform_name .eq. 'Aqua ') then
      masnir = pxldat0(7)
      !masnir = pxldat0(6)
      !else
!        call message( 'snow_mask', 'Platform name not recognized' //
!     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      !endif

! ... debug statement ............................................
      !if (debug .gt. 0) then
      !  write(*,'(10x/,''Subroutine Snow_Mask '',/)')
      !endif
! ................................................................

! ... Initialize
      ndsi = 0.0
      ndsi_snow = .false.
      I_bad = nint(bad_data)

! ................................................................

! ... First, check the 11 micron brightness temperature.
      if (nint(masir11) .ne. I_bad .and. masir11 .le. sm_bt11(1)) then

! ...    perform the NDSI on the current pixel
         if (nint(masv55) .ne. I_bad .and. nint(masnir).ne.I_bad) then
            ndsi  = (masv55 - masnir) / (masv55 + masnir)

           if ( (ndsi .gt. sm_ndsi(1)) .and. (masv88 .gt. sm_ref2(1)) ) then

             ndsi_snow = .true.

!            Check for false snow detection.

!            Now, make sure NDSI is not flagging a thin cirrus
             if( .not. (Greenland .and. land) ) then
               !if(nint(masv188).ne.I_bad.and.nint(masir13).ne.I_bad) then
               !  if (masv188 .gt. sm_ref3(1) .and. masir13 .lt. sm_co2(1)) then  !   CHECK by Zhaojun  
               if(nint(masv188).ne.I_bad) then
                 if (masv188 .gt. sm_ref3(1)) then
                    ndsi_snow = .false.
                 endif
               endif
             endif
!            If in sunglint region, disregard if between -60 and 50 lat.
             if(water .and. snglnt .and. plat .le. 50.0 .and. plat .ge. -60.0) then
               ndsi_snow = .false.
             end if
!            Check for ice clouds mis-identified as snow. Note modified
!            BTD thresholds for Greenland.
             if(masir85 .ne. I_bad .and. masir11 .ne. I_bad) then
               diff = masir85 - masir11
               sth85_11 = sm85_11(1)
               
               if(masir37 .ne. I_bad .and. (.not. (Greenland .and. land)) ) then
                 diff2 = masir37 - masir11
                 sth37_11 = sm37_11(1)
                 if(diff .ge. sth85_11 .and. diff2 .ge. sth37_11) then
                   ndsi_snow = .false.
                 end if
               else
                 if(Greenland .and. land) then
                   sth85_11 = sm85_11(1) + 1.5
                 else
                   sth85_11 = sm85_11(1)
                 end if
                 if(diff .ge. sth85_11) then
                   ndsi_snow = .false.
                 end if
               end if
               
             end if
!            Check for water clouds mis-identified as snow.
             if(masir37 .ne. I_bad .and. masir11 .ne. I_bad) then
               if( .not. (land .and. Greenland) ) then
                 diff = masir37 - masir11
                 if(hi_elev) then
                   sth37_11 = sm37_11hel(1)
                 else
                   sth37_11 = sm37_11(1)
                 end if
                 if(diff .ge. sth37_11) then
                   ndsi_snow = .false.
                 end if
               end if
             end if

             if( (.not. hi_elev) .and. (.not. (Greenland .and. land)) ) then
               if(masnir .gt. sm_mnir(1)) then
                 ndsi_snow = .false.
               end if
             end if
             
           endif 
         endif
      endif

! ... debug statement ............................................
      !if (debug .gt. 0) then
      !  write(*,'(15x,'' snow mask variables '',L4)') ndsi_snow
      !  write(*,'(2x,''masv55 masnir ndsi masv188 masv88 masir13 snglnt hi_elev ndsi_snow '',/,6f8.4,f7.2,3L8)') masv55,masnir,ndsi, masv188,masv88,sth85_11,masir13,snglnt,hi_elev,ndsi_snow
      !  write(*,'(''IR vals '',3f10.2)') masir85,masir11,masir37
      !endif
! ................................................................

      return
      
end subroutine snow_mask
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine chk_ele_lin_edge(i,j,line_edge,ele_edge)

!      implicit none
!      save

!---------------------------------------------------------------------
!!F90
!
!!Description:
!     Routine for determining if the current element you are processing
!     is a border pixel.
!
!!Input parameters:
! nc            Counter for current pixel being processed
! s_pixels      Context of pixels in current scan
! pixels_in_edge Number of elements outside of processing region
! klin          Counts number of lines output to bit file
! nlin          Total number of lines to process
! lines_in_edge Number of lines outside of processing region
!
!!Output Parameters:
! ele_edge     Logical variable - true if processing border pixel
!
!---------------------------------------------------------------------
!      include 'global.inc'

!     scalar arguments
      integer i,j
      logical ele_edge
      logical line_edge
      
!     array arguments
!      integer s_pixels(nlcntx)

!     Initialize line_edge to false
      line_edge = .false.

!     Compare current line to border values
      if ((j .eq. 1) .or. (j .eq. sat%nLine)) then
          line_edge = .true.
      end if
      
!     Initialize ele_edge to false
      ele_edge = .false.

!     Compare current line to border values
      if ((i .eq. 1) .or. (i .eq. sat%nElem)) then
          ele_edge = .true.
      end if
  
      return
      
end subroutine chk_ele_lin_edge
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine check_reg_uniformity(ielem,iline,line_edge,ele_edge, &
                                eco_type,day,land,water,        &
                                coast,snow,ice,uniform)

!      implicit none
!      save


! ... Common statement for debug purposes
!      common / bug / debug, h_output

!--------------------------------------------------------------------
!!F90
!
!!DESCRIPTION:
! This entire routine checks is used to decide of the regional 
! processing box is consistent enough for running the uniformity check.
! The diffent backgroun variables are checked against the center
! pixel, and if they are not all the same, then the variable "uniform"
! is set to false.  All variables must be consistent for uniform to
! be set to true.
!
!!Input parameters:
! line_edge     Logical variable - true if processing border line
! ele_edge      Logical variable - true if processing border pixel
! contx_eco     Array containing context of lines of ecosystem values
! contx_topog   Array containing context of lines of land/sea values
! contx_snow    Array containing context of lines of snow values
! n_nadir       Nadir pixel number for lines in a context
! kele          Current processing element
! eco_type      Ecosystem type for current pixel
! day           Logical variable flagging day scenes
! land          Logical variable flagging land scenes
! water         Logical variable flagging water scenes
! coast         Logical variable flagging coast scenes
! snow          Logical variable flagging snow background scenes
! ice           Logical variable flagging ice background scenes
!
!!Output Parameters:
! uniform       Logical variable flagging uniform scenes
!              (Places where all pixels in context are similar)
!---------------------------------------------------------------------

!      include 'global.inc'
!
!     scalar arguments
      integer :: kele, iline, iElem
!       Number of lines in the processing box (context)
      integer,parameter :: nlcntx = 3

!       Number of elements in the processing box (context)
      integer,parameter :: necntx = 3
        
      logical day,uniform,line_edge,ele_edge,  &
              coast,land,water,snow,ice
      byte eco_type

!     array arguments
      !integer contx_topog(npixel,nlcntx),n_nadir(nlcntx), &
      !        contx_snow(npixel,nlcntx)
      !byte contx_eco(npixel,nlcntx)

!     local scalars
      integer i,j,imv,ide,itotal,nland,nwater,ncoast,h_output,debug, ii, jj,kk

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Subroutine check_reg_uniformity '')')
!        write(h_output,'(10x/,''Line Edge or Element Edge? '',2L5/)')
!     +        line_edge, ele_edge
!      endif
! ................................................................

! ... initialize variables
      imv = ((necntx - 1) / 2) + 1
      itotal = nlcntx * necntx
      nland = 0
      nwater = 0
      ncoast = 0
      kele = ielem


! ... Check all surface variables for consistency
!     If any of the checks fail, then set uniformity to zero

! ... First check if line or pixel is in a border region
      if (ele_edge .or. line_edge) then
         uniform = .false.

! ... Check if middle pixel has been flagged as snow or ice.
      else if(snow .or. ice) then
         uniform = .false.

      else

! ...    Check pixels in nlcntx by necntx region for consistency
         do 100 i = iline-1 , iline+1
            do 200 j = ielem-1 , ielem+1
               !ide = kele + (j - imv)
               ide = j
! ...          First check consistency of ecosystem
               if ((eco_type - sat%eco(ide,i)) .ne. 0) uniform = .false.
! ...          Next, check land/water consistency
               if (geo%lsm(ide,i) .eq. 1   .or.    &
                   geo%lsm(ide,i) .eq. 4) then
                  nland = nland + 1
               else if(geo%lsm(ide,i) .eq. 2 .or.  &
                       geo%lsm(ide,i) .eq. 3) then
                  ncoast = ncoast + 1
               else 
                  nwater = nwater + 1
               end if
! ...          Check snow consistency
               if (snow) then
                 !if (contx_snow(ide,i).ne.103    &
                 !    .and. contx_snow(ide,i).ne.104)
                 if (sat%snow_mask(ide,i) .ne. sym%SNOW) &
                  uniform = .false.
               elseif (.not. snow) then
                 !if (contx_snow(ide,i).eq.103    &
                 !    .or. contx_snow(ide,i).eq.104)
                 if (sat%snow_mask(ide,i) .eq. sym%SNOW) &
                  uniform = .false.
               endif
200         continue
! ...       One more check for a common nadir pixel.  Don't want to
! ...       invoke uniformity test if nadir pixels area not matched up.
            !if (n_nadir(i) .ne. n_nadir(nlcntx/2 + 1)) then
            !   uniform = .false.
            !endif
100      continue

!        At the end now, we want to decide if we have a coastline in our region.
         if(nwater + ncoast .eq. itotal ) then
           if (nwater .ne. 9) then
             uniform = .false.
!            Provide "double coastlines".
             coast = .true.
             land = .true.
             water = .false.
           endif
         else if(nland + ncoast .eq. itotal) then
           if (nland .ne. 9) then 
             uniform = .false.
           endif
         else 
           uniform = .false.
           coast = .true.
           land = .true.
           water = .false.
         end if

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(15x,'' Final Uniformity '',L4)') uniform
!        write(h_output,'(2x,''nwater nland ncoast itotal Day Land Water
!     +  Coast'',/,4I6,2x,4L5)') nwater,nland,ncoast,itotal,day,land,
!     + water,coast
!        write(h_output,'(15x,9i5)') ((contx_eco(jj,kk),jj=kele-1,
!     *         kele+1),kk=1,3)
!      endif
! ................................................................

      return
end subroutine check_reg_uniformity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine get_pxl3X3(i,j,i_sta,i_end,j_sta,j_end)

integer :: i, j
integer :: i_sta,i_end,j_sta,j_end

if (i == 1) then
    i_sta = 1
    i_end = 3
else if (i == sat%nElem) then
    i_sta = sat%nElem-2
    i_end = sat%nElem
else
    i_sta = i-1
    i_end = i+1
endif

if (j == 1) then
    j_sta = 1
    j_end = 3
else if (j == sat%nLine) then
    j_sta = sat%nLine-2
    j_end = sat%nLine
else
    j_sta = j-1
    j_end = j+1
endif

end subroutine get_pxl3X3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine convert_cloud_mask()

!-----------------------------------------------------------------------
!                          PROCESSING SECTION
!
!        The Cloud Mask is processed pixel by pixel using a sliding
!        box approach. Regional context data is stored for use
!        in uniformity tests for uncertain pixels.  The edge of 
!        the area in which you are processing, (outline of region)
!        will be processed but will not include uniformity tests.
!----------------------------------------------------------------------



integer(kind=4) :: iline,ielem, i, j, k, ix_nwp, iy_nwp
integer(kind=1) :: b0,b1,b2,b3

print*,'    ... convert Cloud Mask product !!! '

!--- set local POINTERs to output structures
!======================================================================
! Loop over pixels in this segment
!======================================================================
!b2b1
!00=云              0 
!01=可能晴空         1
!10=可信的晴空       2
!11=可信度高的晴空    3
! ++ convert cloud mask bit vales
line_loop_1: do iline= 1, sat%nLine
element_loop_1: do ielem= 1, sat%nElem

    cm_tmp(ielem,iline,1) = 5
    cm_tmp(ielem,iline,2) = 0
    
    b0 = ibits(cm_bitarray(ielem,iline,1),0,1)
    b1 = ibits(cm_bitarray(ielem,iline,1),1,1)
    b2 = ibits(cm_bitarray(ielem,iline,1),2,1)
    
    !print*,'ddd',cm_bitarray(ielem,iline,1),b0, b1, b2
    if (b0 == 1) then 
    
       cm_tmp(ielem,iline,2) = 1
       
       if (b2 ==0 .and. b1 == 0) then ! cloudy
          cm_tmp(ielem,iline,1) = 0
       endif
    
       if (b2 ==0 .and. b1 == 1) then ! pro cloudy
          cm_tmp(ielem,iline,1) = 1
       endif

       if (b2 ==1 .and. b1 == 0) then ! pro clear
          cm_tmp(ielem,iline,1) = 2
       endif

       if (b2 ==1 .and. b1 == 1) then ! clear
          cm_tmp(ielem,iline,1) = 3
       endif
                     
    endif
    
    if (b0 == 0) then 
    
       cm_tmp(ielem,iline,2) = 0
       cm_tmp(ielem,iline,1) = 0
       
    endif

!-----------------------------------------------------------------------
! END loop over pixels in segment
!-----------------------------------------------------------------------
end do element_loop_1
end do line_loop_1

! Apply spatial consistency filter to remove salt-and-pepper noise
call apply_spatial_filter()

end subroutine convert_cloud_mask
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!~~~~~~~~~~~~~~~~~~~ spatial consistency filter ~~~~~~~~~~~~~~~~~~~~~~~~
subroutine apply_spatial_filter()
    integer :: i, j, di, dj, center, same_count
    integer :: counts(0:3), max_count, max_class, c
    integer(kind=1), allocatable :: original(:,:)
    integer :: nx, ny

    nx = sat%nElem
    ny = sat%nLine
    allocate(original(nx, ny))

    do j = 1, ny
        do i = 1, nx
            original(i, j) = cm_tmp(i, j, 1)
        end do
    end do

    do j = 2, ny - 1
        do i = 2, nx - 1
            center = original(i, j)
            if (center == 5) cycle

            same_count = 0
            do c = 0, 3
                counts(c) = 0
            end do

            do dj = -1, 1
                do di = -1, 1
                    if (di == 0 .and. dj == 0) cycle
                    if (original(i+di, j+dj) == center) same_count = same_count + 1
                    c = original(i+di, j+dj)
                    if (c >= 0 .and. c <= 3) then
                        counts(c) = counts(c) + 1
                    end if
                end do
            end do

            if (same_count == 0) then
                max_count = 0
                max_class = center
                do c = 0, 3
                    if (counts(c) > max_count) then
                        max_count = counts(c)
                        max_class = c
                    end if
                end do
                cm_tmp(i, j, 1) = int(max_class, kind=1)
            end if
        end do
    end do

    deallocate(original)
end subroutine apply_spatial_filter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------- END MODULE ---------------------------------
end module fylat_fy3mersi_cloud_mask
