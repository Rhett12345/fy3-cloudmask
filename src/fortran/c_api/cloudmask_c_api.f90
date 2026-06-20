! =============================================================================
! cloudmask_c_api.f90
!
! ISO_C_BINDING wrapper for the FY-3D MERSI-II cloud mask algorithm.
! This file exposes a single C-callable function that processes the entire
! swath with OpenMP parallelization. It calls the existing Fortran subroutines
! from the original coeff codebase.
!
! Architecture: C++ (OpenMP pixel loop) -> this wrapper -> existing Fortran subroutines
! =============================================================================

module cloudmask_c_api_mod
    use, intrinsic :: iso_c_binding
    use cloudmask_data_arrays
    use names_module,       only: fylat_sensor_id, code_root_path
    use data_arrays_module, only: sat, geo, nwp26, nwp36
    use constant
    use thresholds_read_module
    use polar_module
    use land_module
    use water_module
    use fylat_fy3mersi_cloud_mask, only: snow_mask, chk_ele_lin_edge, &
        check_reg_uniformity
    implicit none

    include 'global.inc'

contains

    ! =========================================================================
    ! set_code_root_path_c -- Set the code root path for threshold file lookup.
    ! =========================================================================
    subroutine set_code_root_path_c(path_in, path_len) bind(C, name='set_code_root_path_c')
        character(c_char), intent(in) :: path_in(*)
        integer(c_int), value, intent(in) :: path_len
        integer :: i
        code_root_path = ' '
        do i = 1, min(path_len, 1000)
            code_root_path(i:i) = path_in(i)
        end do
    end subroutine set_code_root_path_c

    ! =========================================================================
    ! process_pixel_c -- C-callable per-pixel cloud mask algorithm.
    !
    ! This subroutine processes a single pixel given its data extracted from
    ! the swath arrays. All state is stored in the threadprivate module
    ! variables of cloudmask_data_arrays, making this thread-safe under OpenMP.
    ! =========================================================================
    subroutine process_pixel_c( &
        ! --- Input: pixel data (25 bands) ---
        pxldat_in,                  &
        ! --- Input: geometry ---
        lat_in, lon_in,             &
        satzen_in, solzen_in,       &
        relaz_in, glint_in,         &
        ! --- Input: NWP ---
        sfctmp_in, pmsl_in,         &
        uwind_in, vwind_in,         &
        tpw_in,                     &
        ! --- Input: ancillary ---
        pelev_in, eco_in,           &
        lsf_in, snow_mask_in,       &
        ! --- Input: clear-sky BT from RTM ---
        btclr_in,                   &
        ! --- Input: 3x3 context (3x3x25) ---
        indat_in,                   &
        ! --- Input: pixel indices ---
        ielem_in, iline_in,         &
        ! --- Output ---
        out_testbits, out_qa_bits,  &
        out_confidence, out_mask,   &
        out_nmtests, out_nbands,    &
        out_shadow_flag, out_smoke_flag &
    ) bind(C, name='process_pixel_c')


        ! --- Arguments ---
        real(c_float), intent(in)    :: pxldat_in(inband)
        real(c_float), value, intent(in) :: lat_in, lon_in
        real(c_float), value, intent(in) :: satzen_in, solzen_in
        real(c_float), value, intent(in) :: relaz_in, glint_in
        real(c_float), value, intent(in) :: sfctmp_in, pmsl_in
        real(c_float), value, intent(in) :: uwind_in, vwind_in
        real(c_float), value, intent(in) :: tpw_in
        real(c_float), value, intent(in) :: pelev_in
        integer(c_signed_char), value, intent(in) :: eco_in
        integer(c_signed_char), value, intent(in) :: lsf_in
        integer(c_signed_char), value, intent(in) :: snow_mask_in
        real(c_float), intent(in)    :: btclr_in(7)
        real(c_float), intent(in)    :: indat_in(necntx, nlcntx, inband)
        integer(c_int), value, intent(in) :: ielem_in, iline_in

        integer(c_signed_char), intent(out) :: out_testbits(6)
        integer(c_signed_char), intent(out) :: out_qa_bits(10)
        real(c_float), intent(out)   :: out_confidence
        integer(c_int), intent(out)  :: out_mask
        integer(c_int), intent(out)  :: out_nmtests, out_nbands
        integer(c_int), intent(out)  :: out_shadow_flag, out_smoke_flag

        ! --- Local variables ---
        real(c_float) :: pxldat_local(inband)
        real(c_float) :: tbadj_local
        integer(c_signed_char) :: eco_type_local
        logical :: process_local
        integer(c_signed_char) :: is_cold_sfc
        real(c_float) :: btclr_local(7)
        real(c_float) :: indat_local(necntx, nlcntx, inband)
        integer :: i_sta, i_end, j_sta, j_end

        ! Copy input data to local arrays
        pxldat_local = pxldat_in
        btclr_local  = btclr_in
        indat_local  = indat_in
        eco_type_local = eco_in

        ! Initialize pixel state (uses threadprivate module variables)
        call pxinit(testbits, qa_bits, precip_water, vza, sfctmp, pmsl, &
                     u_wind, v_wind, plat, plon, lsf,                   &
                     polar, day, night, land, water, coast, snglnt,      &
                     visusd, vrused, snow, ice, desert, bad_value,       &
                     bad_geo, uniform, shadow, smoke, cirrus_ir,         &
                     cirrus_vis, nmtests, nbands, nbad_1km, nbad_250,    &
                     hi_elev, antarctic, sh_ocean, sg_bad_data,          &
                     map_ice, map_snow, sh_lake)

        ! Set pixel geometry and ancillary data
        plat = lat_in
        plon = lon_in
        vza  = satzen_in
        refang = relaz_in  ! relative azimuth angle
        sfctmp = sfctmp_in
        pmsl   = pmsl_in
        u_wind = uwind_in
        v_wind = vwind_in
        precip_water = tpw_in

        ! Count valid bands (nbands)
        nbands = 0
        do k = 1, 19
            if (pxldat_local(k) > -99.0 .and. pxldat_local(k) <= 2.3) &
                nbands = nbands + 1
        end do
        do k = 20, 25
            if (pxldat_local(k) > 0.0 .and. pxldat_local(k) < 1000.0) &
                nbands = nbands + 1
        end do

        ! Use actual LSF from GEO data (passed from Python)
        lsf = lsf_in

        ! Compute processing flags from input data
        call compute_pixel_flags(pxldat_local, pelev_in, eco_in, &
                                  satzen_in, solzen_in, relaz_in,        &
                                  glint_in, snow_mask_in, tbadj_local)

        ! Check edge pixels
        call chk_ele_lin_edge(ielem_in, iline_in, line_edge, ele_edge)

        ! Check regional uniformity
        call check_reg_uniformity(ielem_in, iline_in, line_edge, ele_edge, &
                                   eco_type_local, day, land, water, coast, &
                                   snow, ice, uniform)

        ! Skip processing if bad data in sunglint over water
        process_local = .true.
        if (sg_bad_data .and. water .and. snglnt) process_local = .false.

        ! Cold surface flag
        is_cold_sfc = 0
        if (sfctmp < 265.0) is_cold_sfc = 1

        ! === Main decision tree ===
        if (process_local) then
            if (polar .and. day) then
                call polar_day(pxldat_local, vza, snglnt, visusd, refang, &
                    vrused, cirrus_vis, land, ice, snow, desert, coast,    &
                    eco_type_local, uniform, hi_elev, ielem_in, indat_local, &
                    nmtests, testbits, tbadj_local, antarctic, sh_ocean,   &
                    sfctmp, qa_bits, confdnc, btclr_local, is_cold_sfc)
            else if (polar .and. night) then
                call polar_nite(pxldat_local, vza, land, ice, snow,       &
                    desert, hi_elev, sfctmp, eco_type_local, nmtests, testbits, &
                    uniform, indat_local, ielem_in, antarctic, sh_ocean,   &
                    qa_bits, confdnc, btclr_local, is_cold_sfc)
            else if (land .and. day) then
                call land_day(pxldat_local, vza, visusd, vrused,           &
                    cirrus_vis, desert, coast, snow, ice, hi_elev,         &
                    tbadj_local, eco_type_local, testbits, qa_bits, nmtests, &
                    confdnc, btclr_local, is_cold_sfc,                     &
                    indat_local(:,:,24), ielem_in, iline_in)
            else if (land .and. night) then
                call land_nite(pxldat_local, plat, vza, ice, snow, coast,  &
                    tbadj_local, desert, hi_elev, sh_lake, sfctmp,         &
                    eco_type_local, nmtests, testbits, qa_bits, confdnc,   &
                    precip_water, btclr_local, is_cold_sfc,                &
                    indat_local(:,:,24), ielem_in, iline_in)
            else if (water .and. day) then
                call water_day(pxldat_local, vza, snglnt, visusd, refang,  &
                    cirrus_vis, sfctmp, hi_elev, uniform, ice, snow,       &
                    ielem_in, iline_in, line_edge, sh_ocean, indat_local,  &
                    nmtests, testbits, qa_bits, confdnc, btclr_local)
            else if (water .and. night) then
                call water_nite(pxldat_local, vza, uniform, ice, snow,     &
                    indat_local, sfctmp, sh_ocean, ielem_in, nmtests,      &
                    testbits, qa_bits, confdnc, btclr_local)
            end if

            ! === Post-processing restoral tests ===
            ! Shadow detection
            if (.not. water .and. .not. coast .and. day .and. &
                .not. polar .and. confdnc >= 0.66) then
                call shadows(pxldat_local, shadow, visusd, qa_bits)
            end if

            ! Non-cloud obstruction (smoke/dust)
            if (land .and. day .and. .not. snow) then
                call noncld_obs_chk(indat_local, pxldat_local, confdnc,    &
                    ielem_in, line_edge, iline_in, qa_bits, testbits, smoke)
            end if

            ! Thin cirrus IR check
            if (.not. snow .and. .not. ice) then
                call thin_ci_chk_ir(pxldat_local, vza, cirrus_ir,          &
                    qa_bits, testbits)
            end if

            ! Set processing path bits
            call proc_path(water, land, day, ice, snow, snglnt, coast,     &
                desert, smoke, shadow, testbits)

            ! Set unused bits and confidence
            call set_unused_bits(testbits)
            call set_confdnc(confdnc, testbits)
            call set_quality_A(nmtests, nbands, lsf, qa_bits)

            ! Fill output bit arrays
            call fill_bit_pixel(nmtests, nbands, bad_value, bad_geo,       &
                snglnt, desert, testbits, qa_bits,                         &
                out_testbits, out_qa_bits)
        end if

        ! === Convert cloud mask to integer value ===
        out_mask = convert_cloud_mask_value(out_testbits)

        ! === Copy results to output arguments ===
        out_confidence = confdnc
        out_nmtests    = nmtests
        out_nbands     = nbands
        out_shadow_flag = merge(1, 0, shadow)
        out_smoke_flag  = merge(1, 0, smoke)

    end subroutine process_pixel_c


    ! =========================================================================
    ! process_swath_c -- C-callable swath-level processing with OpenMP.
    !
    ! Processes all pixels in a 2048x2000 swath. This is the main entry point
    ! called from C++. Each pixel is independent and can be processed in parallel.
    ! =========================================================================
    subroutine process_swath_c( &
        ! --- Input arrays (nElem x nLine x ...) ---
        ref_vis, tbb_ir,           &  ! L1b: reflectance + BT
        lat_arr, lon_arr,          &  ! GEO
        satzen_arr, solzen_arr,    &  ! GEO angles
        relaz_arr, glint_arr,      &  ! GEO angles
        sfctmp_arr, pmsl_arr,      &  ! NWP surface
        uwind_arr, vwind_arr,      &  ! NWP wind
        tpw_arr,                   &  ! NWP TPW
        elev_arr, eco_arr,         &  ! Ancillary
        lsf_arr, snow_mask_arr,    &  ! Land-sea flag + Snow/ice mask
        btclr_arr,                 &  ! Clear-sky BT from RTM (nElem x nLine x 7)
        ! --- Dimensions ---
        nElem, nLine,              &
        ! --- Output arrays (nElem x nLine x ...) ---
        out_cm_bitarray,           &  ! 6-byte cloud mask
        out_qa_bitarray,           &  ! 10-byte QA
        out_cloud_mask,            &  ! Integer cloud mask (0-3)
        out_confidence,            &  ! Confidence (0-1)
        out_nmtests_arr,           &
        out_nbands_arr,            &
        out_shadow_arr,            &
        out_smoke_arr              &
    ) bind(C, name='process_swath_c')

        ! --- Arguments ---
        integer(c_int), value, intent(in) :: nElem, nLine
        real(c_float), intent(in)     :: ref_vis(nElem, nLine, 19)
        real(c_float), intent(in)     :: tbb_ir(nElem, nLine, 6)
        real(c_float), intent(in)     :: lat_arr(nElem, nLine)
        real(c_float), intent(in)     :: lon_arr(nElem, nLine)
        real(c_float), intent(in)     :: satzen_arr(nElem, nLine)
        real(c_float), intent(in)     :: solzen_arr(nElem, nLine)
        real(c_float), intent(in)     :: relaz_arr(nElem, nLine)
        real(c_float), intent(in)     :: glint_arr(nElem, nLine)
        real(c_float), intent(in)     :: sfctmp_arr(nElem, nLine)
        real(c_float), intent(in)     :: pmsl_arr(nElem, nLine)
        real(c_float), intent(in)     :: uwind_arr(nElem, nLine)
        real(c_float), intent(in)     :: vwind_arr(nElem, nLine)
        real(c_float), intent(in)     :: tpw_arr(nElem, nLine)
        real(c_float), intent(in)     :: elev_arr(nElem, nLine)
        integer(c_signed_char), intent(in) :: eco_arr(nElem, nLine)
        integer(c_signed_char), intent(in) :: lsf_arr(nElem, nLine)
        integer(c_signed_char), intent(in) :: snow_mask_arr(nElem, nLine)
        real(c_float), intent(in)     :: btclr_arr(nElem, nLine, 7)

        integer(c_signed_char), intent(out) :: out_cm_bitarray(nElem, nLine, 6)
        integer(c_signed_char), intent(out) :: out_qa_bitarray(nElem, nLine, 10)
        integer(c_int), intent(out)   :: out_cloud_mask(nElem, nLine)
        real(c_float), intent(out)    :: out_confidence(nElem, nLine)
        integer(c_int), intent(out)   :: out_nmtests_arr(nElem, nLine)
        integer(c_int), intent(out)   :: out_nbands_arr(nElem, nLine)
        integer(c_int), intent(out)   :: out_shadow_arr(nElem, nLine)
        integer(c_int), intent(out)   :: out_smoke_arr(nElem, nLine)

        ! --- Local variables ---
        integer :: iline, ielem
        real(c_float) :: pxldat_local(inband)
        real(c_float) :: btclr_local(7)
        real(c_float) :: indat_local(necntx, nlcntx, inband)
        real(c_float) :: tbadj_local
        integer(c_signed_char) :: eco_type_local
        logical :: process_local
        integer(c_signed_char) :: is_cold_sfc
        integer :: i_sta, i_end, j_sta, j_end, ii, jj, k
        integer(c_signed_char) :: out_tb(6), out_qa(10)
        real(c_float) :: out_conf
        integer(c_int) :: out_mask, out_nm, out_nb, out_sh, out_sm

        ! Set sensor ID (default: FY-3D = 21)
        ! This should be passed from Python, but we use a default for now
        if (fylat_sensor_id == 0) fylat_sensor_id = 21

        ! Set code root path for threshold file lookup
        ! The threshold file is at: {code_root_path}/coeff/fylat_thresholds.mersi.ii3d.v8
        if (len_trim(code_root_path) == 0) then
            call get_environment_variable('FY3_CODE_ROOT', code_root_path)
            if (len_trim(code_root_path) == 0) then
                code_root_path = '../coeff/'
            end if
        end if

        ! Load thresholds (must be done before parallel region)
        call thresholds_read(fylat_sensor_id)

        ! Set sat structure dimensions (needed by chk_ele_lin_edge, check_reg_uniformity)
        sat%nElem = nElem
        sat%nLine = nLine

        ! Allocate sat pointer arrays used by check_reg_uniformity
        if (associated(sat%eco)) deallocate(sat%eco)
        if (associated(sat%snow_mask)) deallocate(sat%snow_mask)
        allocate(sat%eco(nElem, nLine))
        allocate(sat%snow_mask(nElem, nLine))
        sat%eco = eco_arr
        sat%snow_mask = snow_mask_arr

        ! Allocate geo%lsm used by check_reg_uniformity
        if (associated(geo%lsm)) deallocate(geo%lsm)
        allocate(geo%lsm(nElem, nLine))
        ! Use actual LSF from GEO data
        geo%lsm = lsf_arr

        ! Initialize output arrays
        out_cm_bitarray  = 0
        out_qa_bitarray  = 0
        out_cloud_mask   = 5  ! fill value
        out_confidence   = 0.0
        out_nmtests_arr  = 0
        out_nbands_arr   = 0
        out_shadow_arr   = 0
        out_smoke_arr    = 0

        ! === Main pixel loop with OpenMP parallelization ===
        !$omp parallel do schedule(dynamic, 8) &
        !$omp   private(iline, ielem, pxldat_local, btclr_local, indat_local, &
        !$omp            tbadj_local, eco_type_local, process_local,          &
        !$omp            is_cold_sfc, i_sta, i_end, j_sta, j_end, ii, jj, k, &
        !$omp            out_tb, out_qa, out_conf, out_mask, out_nm, out_nb,  &
        !$omp            out_sh, out_sm) &
        !$omp   shared(ref_vis, tbb_ir, lat_arr, lon_arr, satzen_arr,        &
        !$omp          solzen_arr, relaz_arr, glint_arr, sfctmp_arr,          &
        !$omp          pmsl_arr, uwind_arr, vwind_arr, tpw_arr, elev_arr,    &
        !$omp          eco_arr, snow_mask_arr, btclr_arr,                     &
        !$omp          out_cm_bitarray, out_qa_bitarray, out_cloud_mask,      &
        !$omp          out_confidence, out_nmtests_arr, out_nbands_arr,       &
        !$omp          out_shadow_arr, out_smoke_arr, nElem, nLine)

        do iline = 1, nLine
            do ielem = 1, nElem

                ! --- Extract 25-band pixel data ---
                pxldat_local(1:19) = ref_vis(ielem, iline, 1:19)
                pxldat_local(20:25) = tbb_ir(ielem, iline, 1:6)

                ! FY-3D band 5/19 swap
                if (fylat_sensor_id > 20) then
                    pxldat_local(5) = ref_vis(ielem, iline, 19)
                    pxldat_local(19) = ref_vis(ielem, iline, 5)
                end if

                ! Range check VIS channels
                do k = 1, 19
                    if (pxldat_local(k) <= -99.0 .or. pxldat_local(k) > 2.3) &
                        pxldat_local(k) = bad_data
                end do

                ! Range check IR channels
                do k = 20, 25
                    if (pxldat_local(k) <= 0.0 .or. pxldat_local(k) >= 1000.0) &
                        pxldat_local(k) = bad_data
                end do

                ! --- Clear-sky BT ---
                btclr_local(1:7) = btclr_arr(ielem, iline, 1:7)

                ! --- Extract 3x3 context ---
                i_sta = max(1, ielem - 1)
                i_end = min(nElem, ielem + 1)
                j_sta = max(1, iline - 1)
                j_end = min(nLine, iline + 1)

                indat_local = bad_data
                do jj = j_sta, j_end
                    do ii = i_sta, i_end
                        do k = 1, 19
                            indat_local(ii - ielem + 2, jj - iline + 2, k) = &
                                ref_vis(ii, jj, k)
                        end do
                        do k = 20, 25
                            indat_local(ii - ielem + 2, jj - iline + 2, k) = &
                                tbb_ir(ii, jj, k - 19)
                        end do
                    end do
                end do

                ! --- Process this pixel ---
                call process_pixel_c( &
                    pxldat_local,               &
                    lat_arr(ielem, iline),      &
                    lon_arr(ielem, iline),      &
                    satzen_arr(ielem, iline),   &
                    solzen_arr(ielem, iline),   &
                    relaz_arr(ielem, iline),    &
                    glint_arr(ielem, iline),    &
                    sfctmp_arr(ielem, iline),   &
                    pmsl_arr(ielem, iline),     &
                    uwind_arr(ielem, iline),    &
                    vwind_arr(ielem, iline),    &
                    tpw_arr(ielem, iline),      &
                    elev_arr(ielem, iline),     &
                    eco_arr(ielem, iline),      &
                    lsf_arr(ielem, iline),      &
                    snow_mask_arr(ielem, iline), &
                    btclr_local,                &
                    indat_local,                &
                    ielem, iline,               &
                    out_tb, out_qa,             &
                    out_conf, out_mask,         &
                    out_nm, out_nb, out_sh, out_sm)

                ! --- Store results ---
                out_cm_bitarray(ielem, iline, 1:6)  = out_tb(1:6)
                out_qa_bitarray(ielem, iline, 1:10) = out_qa(1:10)
                out_cloud_mask(ielem, iline)  = out_mask
                out_confidence(ielem, iline)  = out_conf
                out_nmtests_arr(ielem, iline) = out_nm
                out_nbands_arr(ielem, iline)  = out_nb
                out_shadow_arr(ielem, iline)  = out_sh
                out_smoke_arr(ielem, iline)   = out_sm

            end do
        end do
        !$omp end parallel do

        ! Apply 3x3 median filter on confidence to remove salt-and-pepper.
        ! Smooths the continuous confidence field then reclassifies cloud mask.
        call smooth_conf_reclassify(out_confidence, out_cloud_mask, nElem, nLine)

        ! Apply 7/8 majority filter as secondary cleanup.
        call apply_spatial_consistency(out_cloud_mask, nElem, nLine)

    end subroutine process_swath_c


    ! =========================================================================
    ! Helper: compute_pixel_flags -- set surface type flags from input data
    ! =========================================================================
    subroutine compute_pixel_flags(pxldat0, pelev, eco_type_in, &
                                    satzen, solzen, rela, glint, &
                                    snow_mask_val, tbadj)
        real(c_float), intent(in)       :: pxldat0(inband)
        real(c_float), intent(in)       :: pelev
        integer(c_signed_char), intent(in) :: eco_type_in
        real(c_float), intent(in)       :: satzen, solzen, rela, glint
        integer(c_signed_char), intent(in) :: snow_mask_val
        real(c_float), intent(out)      :: tbadj

        include 'snow_mask.inc'

        real :: ndvi_val, ndsi_val
        integer :: sg_band(5)
        logical :: ndsi_snow_local
        integer :: eco_int

        ! Ecosystem type
        eco_int = int(eco_type_in)

        ! Threshold adjustment for elevation
        tbadj = (pelev / 1000.0) * 5.0

        ! Antarctic flag
        antarctic = (plat < -60.0)

        ! Desert determination (simplified from get_pxldat)
        desert = .false.
        if (eco_int >= 7 .and. eco_int <= 10) desert = .true.
        if (eco_int == 16) desert = .true.

        ! Visible ratio test disable for certain ecosystems
        vrused = .true.
        select case (eco_int)
            case (2, 8, 11, 40, 41, 46, 50, 51, 52, 59, 71)
                vrused = .false.
        end select

        ! Day/night
        if (solzen < 0.0) then
            bad_geo = .true.
            day = .false.
            night = .false.
        else if (solzen > 85.0) then
            night = .true.
            day = .false.
        else
            day = .true.
            night = .false.
        end if

        ! Polar
        polar = (abs(plat) > 60.0)

        ! Visible data usable
        visusd = .not. night

        ! Sunglint
        snglnt = (glint <= 36.0)

        ! Land/water from land-sea flag
        ! (Assumes lsf is already set from ancillary data)
        ! Default: if lsf not set externally, use simple heuristic
        if (lsf == 1 .or. lsf == 4) then
            land = .true.
            water = .false.
            coast = .false.
        else if (lsf == 2) then
            coast = .true.
            land = .true.
            water = .false.
        else if (lsf == 3) then
            land = .true.
            water = .false.
            coast = .false.
            sh_lake = .true.
        else
            water = .true.
            land = .false.
            coast = .false.
            sh_ocean = .true.
        end if

        ! High elevation
        hi_elev = (pelev > 2000.0)
        if (plat > 60.0 .and. pelev > 500.0) then
            if ((plon > -80.0 .and. plon < -20.0) .or. &
                (plon > 20.0 .and. plon < 180.0)) then
                hi_elev = .true.
                Greenland = .true.
            end if
        end if

        ! Snow/ice from ancillary map
        map_snow = .false.
        map_ice  = .false.
        select case (int(snow_mask_val))
            case (1, 3)
                map_snow = .true.
            case (2, 4)
                map_ice = .true.
            case (5)
                map_snow = .true.
                map_ice = .true.
        end select

        ! Daytime snow detection via NDSI
        if (day .and. land .and. .not. snglnt .and. .not. water) then
            call snow_mask(pxldat0, plat, land, snglnt, water, hi_elev, &
                            Greenland, ndsi_snow_local)
            ndsi_snow = ndsi_snow_local

            if (ndsi_snow .and. .not. map_ice) then
                snow = .true.
            else if (ndsi_snow .and. map_ice) then
                ice = .true.
            else if (map_snow) then
                snow = .true.
            else if (map_ice) then
                ice = .true.
            end if
        else
            ! Night: use ancillary map only
            ice  = map_ice
            snow = map_snow
        end if

        ! Sun-glint band check
        sg_band = (/3, 4, 20, 24, 25/)
        sg_bad_data = .false.
        if (snglnt .and. water) then
            if (any(pxldat0(sg_band) == bad_data)) sg_bad_data = .true.
        end if

    end subroutine compute_pixel_flags


    ! =========================================================================
    ! Helper: convert_cloud_mask_value -- extract 2-bit cloud mask from testbits
    ! =========================================================================
    function convert_cloud_mask_value(tb) result(cm_val)
        integer(c_signed_char), intent(in) :: tb(6)
        integer(c_int) :: cm_val

        integer :: b0, b1, b2
        integer :: byte_idx, bit_offset

        ! Extract bits 0, 1, 2 from the first byte
        ! Bit 0: quality flag (processed or not)
        ! Bit 1: confidence LSB
        ! Bit 2: confidence MSB
        b0 = 0; b1 = 0; b2 = 0

        if (btest(tb(1), 0)) b0 = 1
        if (btest(tb(1), 1)) b1 = 1
        if (btest(tb(1), 2)) b2 = 1

        if (b0 == 0) then
            cm_val = 5  ! Not processed -> fill value (was 0=cloudy, bug)
        else
            if (b2 == 0 .and. b1 == 0) cm_val = 0  ! Cloudy
            if (b2 == 0 .and. b1 == 1) cm_val = 1  ! Probably cloudy
            if (b2 == 1 .and. b1 == 0) cm_val = 2  ! Probably clear
            if (b2 == 1 .and. b1 == 1) cm_val = 3  ! Confident clear
        end if
    end function convert_cloud_mask_value

    ! =========================================================================
    ! smooth_conf_reclassify -- 3x3 median filter on confidence, reclassify mask
    ! =========================================================================
    subroutine smooth_conf_reclassify(confidence, cloud_mask, nElem, nLine)
        integer, intent(in) :: nElem, nLine
        real(c_float), intent(inout) :: confidence(nElem, nLine)
        integer, intent(inout) :: cloud_mask(nElem, nLine)

        real(c_float) :: window(9), tmp, conf_val
        real(c_float), allocatable :: smoothed(:,:)
        integer :: i, j, di, dj, n, p, q

        allocate(smoothed(nElem, nLine))
        smoothed = confidence

        ! 3x3 median filter
        do j = 2, nLine - 1
            do i = 2, nElem - 1
                n = 0
                do dj = -1, 1
                    do di = -1, 1
                        conf_val = confidence(i + di, j + dj)
                        if (conf_val >= 0.0 .and. cloud_mask(i + di, j + dj) /= 5) then
                            n = n + 1
                            window(n) = conf_val
                        end if
                    end do
                end do
                if (n > 0) then
                    do p = 1, n - 1
                        do q = p + 1, n
                            if (window(p) > window(q)) then
                                tmp = window(p)
                                window(p) = window(q)
                                window(q) = tmp
                            end if
                        end do
                    end do
                    if (mod(n, 2) == 1) then
                        smoothed(i, j) = window((n + 1) / 2)
                    else
                        smoothed(i, j) = (window(n / 2) + window(n / 2 + 1)) * 0.5
                    end if
                end if
            end do
        end do

        ! Reclassify from smoothed confidence (preserve fill=5 for unprocessed)
        do j = 1, nLine
            do i = 1, nElem
                if (cloud_mask(i, j) == 5) cycle
                confidence(i, j) = smoothed(i, j)
                conf_val = smoothed(i, j)
                if (conf_val > 0.99) then
                    cloud_mask(i, j) = 3
                else if (conf_val > 0.95) then
                    cloud_mask(i, j) = 2
                else if (conf_val > 0.66) then
                    cloud_mask(i, j) = 1
                else
                    cloud_mask(i, j) = 0
                end if
            end do
        end do

        deallocate(smoothed)
    end subroutine smooth_conf_reclassify

    ! =========================================================================
    ! apply_spatial_consistency -- 7/8 majority filter for isolated clusters
    ! =========================================================================
    subroutine apply_spatial_consistency(cloud_mask, nElem, nLine)

        integer, intent(in) :: nElem, nLine
        integer, intent(inout) :: cloud_mask(nElem, nLine)

        integer :: i, j, di, dj, center, counts(0:3), max_count, max_class, c
        integer, allocatable :: original(:,:)

        allocate(original(nElem, nLine))
        original = cloud_mask

        do j = 2, nLine - 1
            do i = 2, nElem - 1
                center = original(i, j)
                if (center == 5) cycle

                do c = 0, 3
                    counts(c) = 0
                end do
                do dj = -1, 1
                    do di = -1, 1
                        if (di == 0 .and. dj == 0) cycle
                        c = original(i + di, j + dj)
                        if (c >= 0 .and. c <= 3) then
                            counts(c) = counts(c) + 1
                        end if
                    end do
                end do

                max_count = 0
                max_class = center
                do c = 0, 3
                    if (counts(c) > max_count) then
                        max_count = counts(c)
                        max_class = c
                    end if
                end do

                ! Flip if >= 7 of 8 neighbors agree on a different class
                if (max_class /= center .and. max_count >= 7) then
                    cloud_mask(i, j) = max_class
                end if
            end do
        end do

        deallocate(original)
    end subroutine apply_spatial_consistency

end module cloudmask_c_api_mod
