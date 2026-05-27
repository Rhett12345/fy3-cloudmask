# =============================================================================
# CMakeLists.txt -- FY-3D MERSI-II Cloud Mask Native Engine
#
# Builds the C++/Fortran hybrid library with OpenMP parallelization.
# Architecture: Python -> pybind11 -> C++ (OpenMP) -> Fortran (core algorithm)
# =============================================================================

cmake_minimum_required(VERSION 3.16)
project(cloudmask_native LANGUAGES C Fortran)

# Find pybind11
find_package(pybind11 REQUIRED)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find required libraries
find_package(HDF5 REQUIRED COMPONENTS C Fortran)

# ============================================================================
# Source directories
# ============================================================================
set(FY3_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src/fortran" CACHE PATH
    "Fortran source directory (core modules)")

set(CLOUDMASK_SRC_DIR "${FY3_SRC_DIR}/cloudmask" CACHE PATH
    "Cloud mask algorithm source directory")

set(CAPI_DIR "${FY3_SRC_DIR}/c_api" CACHE PATH
    "Modified ISO_C_BINDING wrappers (with OpenMP)")

set(UTILS_DIR "${FY3_SRC_DIR}/utils" CACHE PATH
    "Utility sources (strings + C)")

# ============================================================================
# Include directories for Fortran compiler
# ============================================================================
# c_api/ takes priority (has modified cloudmask_data_arrays.f90 with threadprivate)
set(FORTRAN_INCLUDE_DIRS
    ${CAPI_DIR}                                    # Modified wrappers
    ${CLOUDMASK_SRC_DIR}                           # Algorithm + .inc files
    ${FY3_SRC_DIR}/core                            # Foundation modules
    ${UTILS_DIR}                                   # String utils
)

# ============================================================================
# Collect Fortran source files
# ============================================================================

# --- Core modules ---
set(CORE_MODULES
    ${FY3_SRC_DIR}/core/names_module.f90
    ${FY3_SRC_DIR}/core/data_arrays_module.f90
    ${FY3_SRC_DIR}/core/platform_module.f90
    ${FY3_SRC_DIR}/core/constant.f90
    ${FY3_SRC_DIR}/core/numerical.f90
    ${FY3_SRC_DIR}/core/planck_module.f90
    ${FY3_SRC_DIR}/core/frontend_module.f90
)

# --- Cloud mask data arrays (our modified version with threadprivate) ---
set(CLOUDMASK_DATA
    ${CAPI_DIR}/cloudmask_data_arrays.f90
)

# --- Threshold reader ---
set(THRESHOLD_READER
    ${CLOUDMASK_SRC_DIR}/thresholds_read_module.f90
    ${CLOUDMASK_SRC_DIR}/param_read_file.f
)

# --- Bit manipulation ---
set(BIT_OPS
    ${CLOUDMASK_SRC_DIR}/set_bit.f
    ${CLOUDMASK_SRC_DIR}/clear_bit.f
    ${CLOUDMASK_SRC_DIR}/check_bits.f
    ${CLOUDMASK_SRC_DIR}/check_qa_bits.f
    ${CLOUDMASK_SRC_DIR}/set_qa_bit.f
    ${CLOUDMASK_SRC_DIR}/set_confdnc.f
    ${CLOUDMASK_SRC_DIR}/set_quality_A.f
    ${CLOUDMASK_SRC_DIR}/set_unused_bits.f
    ${CLOUDMASK_SRC_DIR}/proc_path.f
    ${CLOUDMASK_SRC_DIR}/fill_bit_pixel.f90
)

# --- Core algorithm functions ---
set(ALGORITHM_CORE
    ${CLOUDMASK_SRC_DIR}/conf_test.f
    ${CLOUDMASK_SRC_DIR}/conf_test_2val.f
    ${CLOUDMASK_SRC_DIR}/pxinit.f
    ${CLOUDMASK_SRC_DIR}/tview.f
    ${CLOUDMASK_SRC_DIR}/trispc.f
    ${CLOUDMASK_SRC_DIR}/get_regdif.f
    ${CLOUDMASK_SRC_DIR}/get_regstd.f
    ${CLOUDMASK_SRC_DIR}/spatial_var.f
    ${CLOUDMASK_SRC_DIR}/chk_spatial_var.f
    ${CLOUDMASK_SRC_DIR}/chk_spatial2.f
    ${CLOUDMASK_SRC_DIR}/get_sg_thresholds.f90
    ${CLOUDMASK_SRC_DIR}/get_pn_thresholds.f
    ${CLOUDMASK_SRC_DIR}/get_nl_thresholds.f90
)

# --- Restoral tests ---
set(RESTORAL_TESTS
    ${CLOUDMASK_SRC_DIR}/chk_land.f90
    ${CLOUDMASK_SRC_DIR}/chk_land_nite.f90
    ${CLOUDMASK_SRC_DIR}/chk_coast.f90
    ${CLOUDMASK_SRC_DIR}/chk_sunglint.f90
    ${CLOUDMASK_SRC_DIR}/chk_shallow_water.f
    ${CLOUDMASK_SRC_DIR}/shadows.f90
    ${CLOUDMASK_SRC_DIR}/noncld_obs_chk.f90
    ${CLOUDMASK_SRC_DIR}/thin_ci_chk_ir.f90
)

# --- Spectral test modules ---
set(SPECTRAL_TESTS
    ${CLOUDMASK_SRC_DIR}/LandDay.f90
    ${CLOUDMASK_SRC_DIR}/LandNite.f90
    ${CLOUDMASK_SRC_DIR}/LandDay_desert.f90
    ${CLOUDMASK_SRC_DIR}/LandDay_desert_c.f90
    ${CLOUDMASK_SRC_DIR}/LandDay_coast.f90
    ${CLOUDMASK_SRC_DIR}/ocean_day.f90
    ${CLOUDMASK_SRC_DIR}/ocean_nite.f90
    ${CLOUDMASK_SRC_DIR}/Day_snow.f90
    ${CLOUDMASK_SRC_DIR}/Nite_snow.f90
    ${CLOUDMASK_SRC_DIR}/PolarDay_land.f90
    ${CLOUDMASK_SRC_DIR}/PolarDay_desert.f90
    ${CLOUDMASK_SRC_DIR}/PolarDay_desert_c.f90
    ${CLOUDMASK_SRC_DIR}/PolarDay_coast.f90
    ${CLOUDMASK_SRC_DIR}/PolarDay_snow.f90
    ${CLOUDMASK_SRC_DIR}/PolarDay_ocean.f90
    ${CLOUDMASK_SRC_DIR}/PolarNite_land.f90
    ${CLOUDMASK_SRC_DIR}/PolarNite_snow.f90
    ${CLOUDMASK_SRC_DIR}/PolarNite_ocean.f90
    ${CLOUDMASK_SRC_DIR}/Antarctic_day.f90
)

# --- Decision tree modules ---
set(DECISION_MODULES
    ${CLOUDMASK_SRC_DIR}/polar_module.f90
    ${CLOUDMASK_SRC_DIR}/land_module.f90
    ${CLOUDMASK_SRC_DIR}/water_module.f90
    ${CLOUDMASK_SRC_DIR}/fylat_fy3mersi_cloud_mask.f90
)

# --- Our ISO_C_BINDING wrapper ---
set(C_API_WRAPPER
    ${CAPI_DIR}/cloudmask_c_api.f90
)

# ============================================================================
# All Fortran sources
# ============================================================================
set(ALL_FORTRAN_SOURCES
    ${CORE_MODULES}
    ${CLOUDMASK_DATA}
    ${THRESHOLD_READER}
    ${BIT_OPS}
    ${ALGORITHM_CORE}
    ${RESTORAL_TESTS}
    ${SPECTRAL_TESTS}
    ${DECISION_MODULES}
    ${C_API_WRAPPER}
)

# ============================================================================
# C source files
# ============================================================================
set(C_SOURCES
    ${UTILS_DIR}/optmed.c
)

# ============================================================================
# Build the Fortran/C static library
# ============================================================================
add_library(cloudmask_fortran STATIC ${ALL_FORTRAN_SOURCES} ${C_SOURCES})

target_include_directories(cloudmask_fortran PRIVATE ${FORTRAN_INCLUDE_DIRS})

# Compiler flags
set_target_properties(cloudmask_fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules
)

# OpenMP for Fortran
if(OpenMP_Fortran_FOUND)
    target_link_libraries(cloudmask_fortran PRIVATE OpenMP::OpenMP_Fortran)
endif()

# HDF5
if(HDF5_FOUND)
    target_include_directories(cloudmask_fortran PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cloudmask_fortran PRIVATE ${HDF5_C_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
endif()

# Position-independent code for shared library
set_target_properties(cloudmask_fortran PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ============================================================================
# Build the pybind11 module
# ============================================================================
pybind11_add_module(_cloudmask_native
    ${CMAKE_CURRENT_SOURCE_DIR}/cloudmask_pybind.cpp
)

target_include_directories(_cloudmask_native PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules
)

target_link_libraries(_cloudmask_native PRIVATE
    cloudmask_fortran
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(_cloudmask_native PRIVATE OpenMP::OpenMP_CXX)
endif()

# ============================================================================
# Install
# ============================================================================
install(TARGETS _cloudmask_native
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/fy3_cloudmask
)
