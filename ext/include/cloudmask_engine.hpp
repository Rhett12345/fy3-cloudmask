/**
 * cloudmask_engine.hpp
 *
 * C++ header for the FY-3D MERSI-II Cloud Mask Engine.
 * Declares the Fortran ISO_C_BINDING interface and provides a C++ wrapper
 * with OpenMP parallelization.
 */

#pragma once

#include <cstdint>
#include <vector>
#include <stdexcept>

// Fortran array dimension constants (must match global.inc)
constexpr int INBAND = 25;
constexpr int NLCNTX = 3;
constexpr int NECNTX = 3;
constexpr int IR_BAND = 6;
constexpr int VIS_BAND = 19;

/**
 * Result structure for a single pixel.
 */
struct PixelResult {
    int8_t testbits[6];      // 48-bit cloud mask test results
    int8_t qa_bits[10];      // 80-bit QA flags
    float  confidence;       // Unobstructed FOV confidence [0, 1]
    int    cloud_mask;       // Integer cloud mask: 0=cloudy, 1=prob_cloudy, 2=prob_clear, 3=clear
    int    n_tests;          // Number of spectral tests applied
    int    n_bands;          // Number of good bands
    int    shadow;           // Shadow flag (0 or 1)
    int    smoke;            // Smoke/dust flag (0 or 1)
};

/**
 * Result structure for a full swath.
 */
struct SwathResult {
    std::vector<int8_t> cm_bitarray;    // (nElem * nLine * 6)
    std::vector<int8_t> qa_bitarray;    // (nElem * nLine * 10)
    std::vector<int>    cloud_mask;     // (nElem * nLine)
    std::vector<float>  confidence;     // (nElem * nLine)
    std::vector<int>    nmtests;        // (nElem * nLine)
    std::vector<int>    nbands;         // (nElem * nLine)
    std::vector<int>    shadow;         // (nElem * nLine)
    std::vector<int>    smoke;          // (nElem * nLine)

    SwathResult(int nElem, int nLine)
        : cm_bitarray(nElem * nLine * 6, 0)
        , qa_bitarray(nElem * nLine * 10, 0)
        , cloud_mask(nElem * nLine, 5)      // 5 = fill value
        , confidence(nElem * nLine, 0.0f)
        , nmtests(nElem * nLine, 0)
        , nbands(nElem * nLine, 0)
        , shadow(nElem * nLine, 0)
        , smoke(nElem * nLine, 0)
    {}
};

// ============================================================================
// Fortran ISO_C_BINDING interface
// ============================================================================

extern "C" {
    /**
     * Process a single pixel through the cloud mask algorithm.
     * Thread-safe: each call uses threadprivate Fortran module state.
     */
    void process_pixel_c(
        const float* pxldat_in,           // [25] band data
        float lat_in, float lon_in,
        float satzen_in, float solzen_in,
        float relaz_in, float glint_in,
        float sfctmp_in, float pmsl_in,
        float uwind_in, float vwind_in,
        float tpw_in,
        float pelev_in,
        int8_t eco_in,
        int8_t snow_mask_in,
        const float* btclr_in,            // [7] clear-sky BT
        const float* indat_in,            // [3x3x25] context
        int ielem_in, int iline_in,
        int8_t* out_testbits,             // [6]
        int8_t* out_qa_bits,              // [10]
        float*  out_confidence,
        int*    out_mask,
        int*    out_nmtests,
        int*    out_nbands,
        int*    out_shadow,
        int*    out_smoke
    );

    /**
     * Process an entire swath with OpenMP parallelization.
     * This is the high-performance entry point.
     */
    void process_swath_c(
        const float* ref_vis,             // [nElem * nLine * 19]
        const float* tbb_ir,              // [nElem * nLine * 6]
        const float* lat_arr,             // [nElem * nLine]
        const float* lon_arr,
        const float* satzen_arr,
        const float* solzen_arr,
        const float* relaz_arr,
        const float* glint_arr,
        const float* sfctmp_arr,
        const float* pmsl_arr,
        const float* uwind_arr,
        const float* vwind_arr,
        const float* tpw_arr,
        const float* elev_arr,
        const int8_t* eco_arr,
        const int8_t* snow_mask_arr,
        const float* btclr_arr,           // [nElem * nLine * 7]
        int nElem, int nLine,
        int8_t* out_cm_bitarray,          // [nElem * nLine * 6]
        int8_t* out_qa_bitarray,          // [nElem * nLine * 10]
        int*    out_cloud_mask,           // [nElem * nLine]
        float*  out_confidence,           // [nElem * nLine]
        int*    out_nmtests_arr,
        int*    out_nbands_arr,
        int*    out_shadow_arr,
        int*    out_smoke_arr
    );
}

// ============================================================================
// C++ wrapper class
// ============================================================================

class CloudMaskEngine {
public:
    /**
     * Process a full swath and return results.
     *
     * @param ref_vis    Visible reflectance [nElem * nLine * 19], float32
     * @param tbb_ir     IR brightness temperature [nElem * nLine * 6], float32
     * @param lat        Latitude [nElem * nLine], float32
     * @param lon        Longitude [nElem * nLine], float32
     * @param satzen     Satellite zenith angle [nElem * nLine], float32
     * @param solzen     Solar zenith angle [nElem * nLine], float32
     * @param relaz      Relative azimuth angle [nElem * nLine], float32
     * @param glint      Glint angle [nElem * nLine], float32
     * @param sfctmp     Surface temperature [nElem * nLine], float32
     * @param pmsl       Mean sea level pressure [nElem * nLine], float32
     * @param uwind      U-wind component [nElem * nLine], float32
     * @param vwind      V-wind component [nElem * nLine], float32
     * @param tpw        Total precipitable water [nElem * nLine], float32
     * @param elev       Elevation [nElem * nLine], float32
     * @param eco        Ecosystem type [nElem * nLine], int8
     * @param snow_mask  Snow/ice mask [nElem * nLine], int8
     * @param btclr      Clear-sky BT [nElem * nLine * 7], float32
     * @param nElem      Number of elements (columns)
     * @param nLine      Number of lines (rows)
     * @return           SwathResult with all output arrays
     */
    static SwathResult process_swath(
        const float* ref_vis,
        const float* tbb_ir,
        const float* lat,
        const float* lon,
        const float* satzen,
        const float* solzen,
        const float* relaz,
        const float* glint,
        const float* sfctmp,
        const float* pmsl,
        const float* uwind,
        const float* vwind,
        const float* tpw,
        const float* elev,
        const int8_t* eco,
        const int8_t* snow_mask,
        const float* btclr,
        int nElem,
        int nLine
    ) {
        SwathResult result(nElem, nLine);

        process_swath_c(
            ref_vis, tbb_ir,
            lat, lon, satzen, solzen, relaz, glint,
            sfctmp, pmsl, uwind, vwind, tpw,
            elev, eco, snow_mask, btclr,
            nElem, nLine,
            result.cm_bitarray.data(),
            result.qa_bitarray.data(),
            result.cloud_mask.data(),
            result.confidence.data(),
            result.nmtests.data(),
            result.nbands.data(),
            result.shadow.data(),
            result.smoke.data()
        );

        return result;
    }
};
