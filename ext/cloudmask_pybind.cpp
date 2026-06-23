/**
 * cloudmask_pybind.cpp
 *
 * pybind11 bindings for the FY-3D MERSI-II Cloud Mask Engine.
 * Exposes the Fortran/C++ hybrid engine to Python with zero-copy numpy arrays.
 *
 * Memory layout handling:
 *   - Input: Python/numpy arrays are C-order (row-major). We transpose them
 *     to Fortran column-major order before passing to the Fortran engine.
 *   - Output: Fortran writes results in column-major order. We reshape the
 *     flat output vectors using numpy's order='F' to get correct (nElem, nLine)
 *     arrays for Python.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <cstring>

#include "include/cloudmask_engine.hpp"

namespace py = pybind11;

// ============================================================================
// C-order to Fortran-order transpose helpers (cache-blocked)
// ============================================================================

static constexpr int BLOCK = 32;

/**
 * Transpose 2D C-order (nElem, nLine) -> Fortran-order memory layout.
 * Uses cache-blocked tiling for better memory access patterns.
 */
template<typename T>
py::array_t<T> transpose_2d(py::array_t<T> arr, int nElem, int nLine) {
    auto buf = arr.request();
    auto out = py::array_t<T>({nElem, nLine});
    auto out_buf = out.request();
    const T* src = static_cast<const T*>(buf.ptr);
    T* dst = static_cast<T*>(out_buf.ptr);

    for (int jj = 0; jj < nLine; jj += BLOCK) {
        int jEnd = std::min(jj + BLOCK, nLine);
        for (int ii = 0; ii < nElem; ii += BLOCK) {
            int iEnd = std::min(ii + BLOCK, nElem);
            for (int j = jj; j < jEnd; j++) {
                for (int i = ii; i < iEnd; i++) {
                    dst[j * nElem + i] = src[i * nLine + j];
                }
            }
        }
    }
    return out;
}

/**
 * Transpose 3D C-order (nElem, nLine, K) -> Fortran-order memory layout.
 * For small K (<=32), processes all K slices in inner loop for cache locality.
 */
template<typename T>
py::array_t<T> transpose_3d(py::array_t<T> arr, int nElem, int nLine, int K) {
    auto buf = arr.request();
    auto out = py::array_t<T>({nElem, nLine, K});
    auto out_buf = out.request();
    const T* src = static_cast<const T*>(buf.ptr);
    T* dst = static_cast<T*>(out_buf.ptr);

    for (int jj = 0; jj < nLine; jj += BLOCK) {
        int jEnd = std::min(jj + BLOCK, nLine);
        for (int ii = 0; ii < nElem; ii += BLOCK) {
            int iEnd = std::min(ii + BLOCK, nElem);
            for (int j = jj; j < jEnd; j++) {
                for (int i = ii; i < iEnd; i++) {
                    const T* src_px = src + (i * nLine + j) * K;
                    for (int k = 0; k < K; k++) {
                        dst[(k * nLine + j) * nElem + i] = src_px[k];
                    }
                }
            }
        }
    }
    return out;
}

// ============================================================================
// Main processing function
// ============================================================================

py::dict process_swath_py(
    py::array_t<float, py::array::c_style | py::array::forcecast> ref_vis,
    py::array_t<float, py::array::c_style | py::array::forcecast> tbb_ir,
    py::array_t<float, py::array::c_style | py::array::forcecast> lat,
    py::array_t<float, py::array::c_style | py::array::forcecast> lon,
    py::array_t<float, py::array::c_style | py::array::forcecast> satzen,
    py::array_t<float, py::array::c_style | py::array::forcecast> solzen,
    py::array_t<float, py::array::c_style | py::array::forcecast> relaz,
    py::array_t<float, py::array::c_style | py::array::forcecast> glint,
    py::array_t<float, py::array::c_style | py::array::forcecast> sfctmp,
    py::array_t<float, py::array::c_style | py::array::forcecast> pmsl,
    py::array_t<float, py::array::c_style | py::array::forcecast> uwind,
    py::array_t<float, py::array::c_style | py::array::forcecast> vwind,
    py::array_t<float, py::array::c_style | py::array::forcecast> tpw,
    py::array_t<float, py::array::c_style | py::array::forcecast> elev,
    py::array_t<int8_t, py::array::c_style | py::array::forcecast> eco,
    py::array_t<int8_t, py::array::c_style | py::array::forcecast> lsf,
    py::array_t<int8_t, py::array::c_style | py::array::forcecast> snow_mask,
    py::array_t<float, py::array::c_style | py::array::forcecast> btclr,
    int nElem, int nLine,
    std::string code_root_path
) {
    // Validate input shapes
    auto rv_buf = ref_vis.request();
    if (rv_buf.ndim != 3 || rv_buf.shape[0] != nElem || rv_buf.shape[1] != nLine || rv_buf.shape[2] != 19)
        throw std::runtime_error("ref_vis must be (nElem, nLine, 19)");

    auto ir_buf = tbb_ir.request();
    if (ir_buf.ndim != 3 || ir_buf.shape[0] != nElem || ir_buf.shape[1] != nLine || ir_buf.shape[2] != 6)
        throw std::runtime_error("tbb_ir must be (nElem, nLine, 6)");

    auto bt_buf = btclr.request();
    if (bt_buf.ndim != 3 || bt_buf.shape[0] != nElem || bt_buf.shape[1] != nLine || bt_buf.shape[2] != 7)
        throw std::runtime_error("btclr must be (nElem, nLine, 7)");

    // --- Set code root path for threshold file lookup ---
    set_code_root_path_c(code_root_path.c_str(), static_cast<int>(code_root_path.size()));

    // --- Transpose inputs: C-order -> Fortran column-major order ---
    auto ref_vis_f = transpose_3d<float>(ref_vis, nElem, nLine, 19);
    auto tbb_ir_f  = transpose_3d<float>(tbb_ir,  nElem, nLine, 6);
    auto btclr_f   = transpose_3d<float>(btclr,   nElem, nLine, 7);

    auto lat_f     = transpose_2d<float>(lat,     nElem, nLine);
    auto lon_f     = transpose_2d<float>(lon,     nElem, nLine);
    auto satzen_f  = transpose_2d<float>(satzen,  nElem, nLine);
    auto solzen_f  = transpose_2d<float>(solzen,  nElem, nLine);
    auto relaz_f   = transpose_2d<float>(relaz,   nElem, nLine);
    auto glint_f   = transpose_2d<float>(glint,   nElem, nLine);
    auto sfctmp_f  = transpose_2d<float>(sfctmp,  nElem, nLine);
    auto pmsl_f    = transpose_2d<float>(pmsl,    nElem, nLine);
    auto uwind_f   = transpose_2d<float>(uwind,   nElem, nLine);
    auto vwind_f   = transpose_2d<float>(vwind,   nElem, nLine);
    auto tpw_f     = transpose_2d<float>(tpw,     nElem, nLine);
    auto elev_f    = transpose_2d<float>(elev,    nElem, nLine);

    auto eco_f       = transpose_2d<int8_t>(eco,       nElem, nLine);
    auto lsf_f       = transpose_2d<int8_t>(lsf,       nElem, nLine);
    auto snow_mask_f = transpose_2d<int8_t>(snow_mask, nElem, nLine);

    // --- Call the Fortran engine via C++ wrapper ---
    SwathResult result = CloudMaskEngine::process_swath(
        ref_vis_f.data(),
        tbb_ir_f.data(),
        lat_f.data(),
        lon_f.data(),
        satzen_f.data(),
        solzen_f.data(),
        relaz_f.data(),
        glint_f.data(),
        sfctmp_f.data(),
        pmsl_f.data(),
        uwind_f.data(),
        vwind_f.data(),
        tpw_f.data(),
        elev_f.data(),
        eco_f.data(),
        lsf_f.data(),
        snow_mask_f.data(),
        btclr_f.data(),
        nElem, nLine
    );

    // --- Convert results: Fortran flat vectors -> numpy arrays ---
    // Fortran stores (nElem, nLine) as column-major: [i0_j0, i1_j0, ..., iN_j0, i0_j1, ...]
    // We copy into a flat numpy array and reshape with order='F' to get correct (nElem, nLine).

    auto np = py::module_::import("numpy");

    auto reshape_2d_int = [&](const std::vector<int>& vec) {
        auto flat = py::array_t<int>(vec.size());
        std::memcpy(flat.mutable_data(), vec.data(), vec.size() * sizeof(int));
        return np.attr("reshape")(flat, py::make_tuple(nElem, nLine), py::arg("order") = "F");
    };
    auto reshape_2d_float = [&](const std::vector<float>& vec) {
        auto flat = py::array_t<float>(vec.size());
        std::memcpy(flat.mutable_data(), vec.data(), vec.size() * sizeof(float));
        return np.attr("reshape")(flat, py::make_tuple(nElem, nLine), py::arg("order") = "F");
    };
    auto reshape_3d_int8 = [&](const std::vector<int8_t>& vec, int K) {
        auto flat = py::array_t<int8_t>(vec.size());
        std::memcpy(flat.mutable_data(), vec.data(), vec.size() * sizeof(int8_t));
        return np.attr("reshape")(flat, py::make_tuple(nElem, nLine, K), py::arg("order") = "F");
    };

    py::dict out;
    out["cm_bitarray"] = reshape_3d_int8(result.cm_bitarray, 6);
    out["qa_bitarray"] = reshape_3d_int8(result.qa_bitarray, 10);
    out["cloud_mask"]  = reshape_2d_int(result.cloud_mask);
    out["confidence"]  = reshape_2d_float(result.confidence);
    out["nmtests"]     = reshape_2d_int(result.nmtests);
    out["nbands"]      = reshape_2d_int(result.nbands);
    out["shadow"]      = reshape_2d_int(result.shadow);
    out["smoke"]       = reshape_2d_int(result.smoke);
    return out;
}


// ============================================================================
// Transpose verification (test-only, returns pass/fail dict)
// ============================================================================

py::dict verify_transpose_py(int nElem, int nLine) {
    auto np = py::module_::import("numpy");
    py::dict result;

    // --- 2D verification ---
    // Create C-order (nElem, nLine) with arr[i,j] = 10000*i + j
    auto arr2d = py::array_t<int>({nElem, nLine});
    auto arr2d_buf = arr2d.request();
    int* arr2d_ptr = static_cast<int*>(arr2d_buf.ptr);
    for (int i = 0; i < nElem; i++)
        for (int j = 0; j < nLine; j++)
            arr2d_ptr[i * nLine + j] = 10000 * i + j;

    auto arr2d_f = transpose_2d<int>(arr2d, nElem, nLine);
    auto arr2d_f_buf = arr2d_f.request();
    int* arr2d_f_ptr = static_cast<int*>(arr2d_f_buf.ptr);

    bool pass_2d = true;
    for (int i = 0; i < nElem && pass_2d; i++)
        for (int j = 0; j < nLine && pass_2d; j++)
            if (arr2d_f_ptr[j * nElem + i] != 10000 * i + j)
                pass_2d = false;

    result["pass_2d"] = pass_2d;

    // --- 2D round-trip: Fortran order -> C-order (via reshape with order='F') ---
    auto flat_2d = py::array_t<int>(nElem * nLine);
    std::memcpy(flat_2d.mutable_data(), arr2d_f_ptr, nElem * nLine * sizeof(int));
    auto arr2d_rt_f = np.attr("reshape")(flat_2d, py::make_tuple(nElem, nLine), py::arg("order") = "F");
    auto arr2d_rt = py::cast<py::array_t<int>>(
        np.attr("ascontiguousarray")(arr2d_rt_f)
    );
    auto arr2d_rt_buf = arr2d_rt.request();
    int* arr2d_rt_ptr = static_cast<int*>(arr2d_rt_buf.ptr);

    bool pass_2d_rt = true;
    for (int i = 0; i < nElem && pass_2d_rt; i++)
        for (int j = 0; j < nLine && pass_2d_rt; j++)
            if (arr2d_rt_ptr[i * nLine + j] != 10000 * i + j)
                pass_2d_rt = false;

    result["pass_2d_roundtrip"] = pass_2d_rt;

    // --- 3D verification ---
    int K = 7;  // use K=7 like btclr
    auto arr3d = py::array_t<int>({nElem, nLine, K});
    auto arr3d_buf = arr3d.request();
    int* arr3d_ptr = static_cast<int*>(arr3d_buf.ptr);
    for (int i = 0; i < nElem; i++)
        for (int j = 0; j < nLine; j++)
            for (int k = 0; k < K; k++)
                arr3d_ptr[(i * nLine + j) * K + k] = 10000 * i + 10 * j + k;

    auto arr3d_f = transpose_3d<int>(arr3d, nElem, nLine, K);
    auto arr3d_f_buf = arr3d_f.request();
    int* arr3d_f_ptr = static_cast<int*>(arr3d_f_buf.ptr);

    bool pass_3d = true;
    for (int i = 0; i < nElem && pass_3d; i++)
        for (int j = 0; j < nLine && pass_3d; j++)
            for (int k = 0; k < K && pass_3d; k++)
                if (arr3d_f_ptr[(k * nLine + j) * nElem + i] != 10000 * i + 10 * j + k)
                    pass_3d = false;

    result["pass_3d"] = pass_3d;

    // --- 3D round-trip: Fortran order -> C-order ---
    auto flat_3d = py::array_t<int>(nElem * nLine * K);
    std::memcpy(flat_3d.mutable_data(), arr3d_f_ptr, nElem * nLine * K * sizeof(int));
    auto arr3d_rt_f = np.attr("reshape")(flat_3d, py::make_tuple(nElem, nLine, K), py::arg("order") = "F");
    auto arr3d_rt = py::cast<py::array_t<int>>(
        np.attr("ascontiguousarray")(arr3d_rt_f)
    );
    auto arr3d_rt_buf = arr3d_rt.request();
    int* arr3d_rt_ptr = static_cast<int*>(arr3d_rt_buf.ptr);

    bool pass_3d_rt = true;
    for (int i = 0; i < nElem && pass_3d_rt; i++)
        for (int j = 0; j < nLine && pass_3d_rt; j++)
            for (int k = 0; k < K && pass_3d_rt; k++)
                if (arr3d_rt_ptr[(i * nLine + j) * K + k] != 10000 * i + 10 * j + k)
                    pass_3d_rt = false;

    result["pass_3d_roundtrip"] = pass_3d_rt;

    return result;
}


PYBIND11_MODULE(_cloudmask_native, m) {
    m.doc() = R"doc(
        FY-3D MERSI-II Cloud Mask Engine -- Native C++/Fortran backend.

        This module provides the high-performance cloud mask algorithm
        implemented in Fortran with C++ OpenMP parallelization.

        Architecture:
            Python (config/CLI/IO) -> C++ (OpenMP pixel loop) -> Fortran (core algorithm)
    )doc";

    m.def("process_swath", &process_swath_py,
        py::arg("ref_vis"),
        py::arg("tbb_ir"),
        py::arg("lat"),
        py::arg("lon"),
        py::arg("satzen"),
        py::arg("solzen"),
        py::arg("relaz"),
        py::arg("glint"),
        py::arg("sfctmp"),
        py::arg("pmsl"),
        py::arg("uwind"),
        py::arg("vwind"),
        py::arg("tpw"),
        py::arg("elev"),
        py::arg("eco"),
        py::arg("lsf"),
        py::arg("snow_mask"),
        py::arg("btclr"),
        py::arg("nElem"),
        py::arg("nLine"),
        py::arg("code_root_path") = "",
        R"doc(
            Process a full swath through the cloud mask algorithm.

            Parameters
            ----------
            ref_vis : ndarray, shape (nElem, nLine, 19), float32
                Visible/NIR reflectance for 19 channels.
            tbb_ir : ndarray, shape (nElem, nLine, 6), float32
                IR brightness temperature for 6 channels.
            lat, lon : ndarray, shape (nElem, nLine), float32
                Geolocation.
            satzen, solzen, relaz, glint : ndarray, shape (nElem, nLine), float32
                Geometry angles (degrees).
            sfctmp, pmsl : ndarray, shape (nElem, nLine), float32
                NWP surface temperature (K) and pressure (hPa).
            uwind, vwind : ndarray, shape (nElem, nLine), float32
                NWP wind components (m/s).
            tpw : ndarray, shape (nElem, nLine), float32
                Total precipitable water (cm).
            elev : ndarray, shape (nElem, nLine), float32
                Elevation (m).
            eco : ndarray, shape (nElem, nLine), int8
                IGBP ecosystem type.
            lsf : ndarray, shape (nElem, nLine), int8
                Land-sea flag (0=water, 1=land, 2=coast, 3=shallow_lake, 4=land).
            snow_mask : ndarray, shape (nElem, nLine), int8
                NISE snow/ice mask.
            btclr : ndarray, shape (nElem, nLine, 7), float32
                Clear-sky brightness temperatures from RTM.
            nElem, nLine : int
                Swath dimensions (columns, rows).

            Returns
            -------
            dict with keys:
                cm_bitarray : ndarray, shape (nElem, nLine, 6), int8
                    48-bit cloud mask test results.
                qa_bitarray : ndarray, shape (nElem, nLine, 10), int8
                    80-bit QA flags.
                cloud_mask : ndarray, shape (nElem, nLine), int32
                    Integer cloud mask (0=cloudy, 3=confident clear).
                confidence : ndarray, shape (nElem, nLine), float32
                    Confidence value [0, 1].
                nmtests : ndarray, shape (nElem, nLine), int32
                    Number of spectral tests applied.
                nbands : ndarray, shape (nElem, nLine), int32
                    Number of good bands.
                shadow : ndarray, shape (nElem, nLine), int32
                    Shadow detection flag.
                smoke : ndarray, shape (nElem, nLine), int32
                    Smoke/dust obstruction flag.
        )doc"
    );

    m.attr("__version__") = "3.6.0";
    m.attr("__backend__") = "C++/Fortran (OpenMP)";

    m.def("verify_transpose", &verify_transpose_py,
        py::arg("nElem"), py::arg("nLine"),
        "Verify C-to-Fortran transpose correctness (2D + 3D round-trip). Returns dict of pass/fail.");
}
