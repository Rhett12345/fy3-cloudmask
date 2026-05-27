! =============================================================================
! cloudmask_data_arrays.f90 -- Thread-safe version with OpenMP threadprivate
!
! This is a modified version of the original cloudmask_data_arrays module.
! All per-pixel state variables are declared with the threadprivate directive
! so that each OpenMP thread has its own private copy.
!
! Original: retrieval_system_V3.1_cldmask/src/cloudmask/cloudmask_data_arrays.f90
! =============================================================================

module cloudmask_data_arrays
    implicit none

    ! --- Per-pixel logical flags (threadprivate for OpenMP safety) ---
    logical :: line_edge, ele_edge, polar, land, day, night, ice, snglnt, visusd, water, bad_value, &
               coast, desert, vrused, snow, bad_geo, map_ice, map_snow, ndsi_snow,                  &
               hi_elev, antarctic, sh_ocean, sg_bad_data, sh_lake,                                  &
               New_Zealand, Greenland, process, cirrus_ir, cirrus_vis, no_250,                       &
               uniform, shadow, smoke

    !$omp threadprivate(line_edge, ele_edge, polar, land, day, night, ice, snglnt, visusd, water, &
    !$omp               bad_value, coast, desert, vrused, snow, bad_geo, map_ice, map_snow,       &
    !$omp               ndsi_snow, hi_elev, antarctic, sh_ocean, sg_bad_data, sh_lake,            &
    !$omp               New_Zealand, Greenland, process, cirrus_ir, cirrus_vis, no_250,           &
    !$omp               uniform, shadow, smoke)

    integer         :: lsf, nmtests, nbands, nbad_1km, nbad_250
    real(kind=4)    :: confdnc, precip_water, vza, plat, plon, sfctmp, pmsl, u_wind, v_wind, refang

    !$omp threadprivate(lsf, nmtests, nbands, nbad_1km, nbad_250)
    !$omp threadprivate(confdnc, precip_water, vza, plat, plon, sfctmp, pmsl, u_wind, v_wind, refang)

    ! --- Output pointer arrays (allocated at runtime, NOT threadprivate) ---
    real(kind=4), dimension(:,:), pointer :: out_pwater
    real(kind=4), dimension(:,:), pointer :: out_sfctmp
    integer, dimension(:,:), pointer      :: out_polar, out_day, out_night, out_land, out_water,      &
                                             out_coast, out_snglnt, out_snow, out_ice, out_desert,    &
                                             out_uniform, out_shadow

    ! --- Per-pixel bit arrays (threadprivate for OpenMP safety) ---
    byte :: testbits(6)    ! 48-bit cloud mask
    byte :: qa_bits(10)    ! 80-bit QA flags

    !$omp threadprivate(testbits, qa_bits)

end module cloudmask_data_arrays
