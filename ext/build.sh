#!/bin/bash
# =============================================================================
# build.sh -- Build the FY-3D Cloud Mask Native Engine (GNU toolchain)
#
# Usage: ./build.sh [--debug] [--clean] [--install]
#
# All Fortran sources live under: ../src/fortran/
#   core/       Foundation modules (names, constant, planck, etc.)
#   cloudmask/  Cloud mask algorithm + threshold .inc files
#   utils/      String utilities + C sources
#   c_api/      Modified ISO_C_BINDING + OpenMP wrappers
#
# Builds: _cloudmask_native.cpython-*.so
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

# Parse arguments
DEBUG=0
CLEAN=0
INSTALL=0
for arg in "$@"; do
    case $arg in
        --debug)   DEBUG=1 ;;
        --clean)   CLEAN=1 ;;
        --install) INSTALL=1 ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

if [ $CLEAN -eq 1 ]; then
    echo "Cleaning build directory..."
    rm -rf "$BUILD_DIR"
    echo "Done."
    exit 0
fi

# =============================================================================
# Toolchain
# =============================================================================
echo "=== Checking toolchain ==="

# Find compilers: prefer system, fallback to conda
find_compiler() {
    local name="$1"
    local conda_name="$2"
    # Try conda prefix first
    if [ -n "$CONDA_PREFIX" ] && [ -x "$CONDA_PREFIX/bin/$conda_name" ]; then
        echo "$CONDA_PREFIX/bin/$conda_name"
        return
    fi
    # Try system
    if command -v "$name" &>/dev/null; then
        echo "$(command -v "$name")"
        return
    fi
    # Try conda default prefix
    local default_conda="$(dirname $(dirname $(which python3 2>/dev/null || echo /usr/bin/python3)))"
    if [ -x "$default_conda/bin/$conda_name" ]; then
        echo "$default_conda/bin/$conda_name"
        return
    fi
    echo ""
}

GFORTRAN=$(find_compiler gfortran gfortran)
GXX=$(find_compiler g++ x86_64-conda-linux-gnu-g++)
GCC=$(find_compiler gcc x86_64-conda-linux-gnu-gcc)

# Fallback: plain g++/gcc if conda wrappers not found
[ -z "$GXX" ] && GXX=$(find_compiler g++ g++)
[ -z "$GCC" ] && GCC=$(find_compiler gcc gcc)

for cmd in "$GFORTRAN" "$GXX" "$GCC"; do
    if [ -z "$cmd" ] || [ ! -x "$cmd" ]; then
        echo "ERROR: Compiler not found: $cmd"
        echo "Install: conda install gfortran_linux-64 gxx_linux-64 gcc_linux-64"
        exit 1
    fi
done

echo "  gfortran: $($GFORTRAN --version | head -1)"
echo "  g++:      $($GXX --version | head -1)"
echo "  gcc:      $($GCC --version | head -1)"

# Python / pybind11
PYTHON=${PYTHON:-python3}
if ! $PYTHON -c "import pybind11" 2>/dev/null; then
    echo "  Installing pybind11..."
    $PYTHON -m pip install pybind11 --quiet
fi
PYBIND11_INCLUDES=$($PYTHON -m pybind11 --includes)
PYTHON_EXT=$($PYTHON -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
PYTHON_VERSION=$($PYTHON -c "import sysconfig; print(sysconfig.get_config_var('SOABI'))")
echo "  Python: $($PYTHON --version) ($PYTHON_VERSION)"

# HDF5 - find dynamically
find_hdf5_lib() {
    local search_dirs="${CONDA_PREFIX:-/usr}/lib /usr/lib/x86_64-linux-gnu /usr/lib64"
    for dir in $search_dirs; do
        if [ -f "$dir/libhdf5.so" ] || [ -f "$dir/libhdf5.so.310" ] || [ -f "$dir/libhdf5.so.200" ]; then
            echo "$dir"
            return
        fi
    done
    echo ""
}

HDF5_LIB_DIR=$(find_hdf5_lib)
if [ -z "$HDF5_LIB_DIR" ]; then
    echo "ERROR: HDF5 not found. Install: conda install hdf5"
    exit 1
fi
echo "  HDF5: $HDF5_LIB_DIR"

# =============================================================================
# Source paths (all Fortran now lives under src/fortran/)
# =============================================================================
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FORTRAN_SRC="$PROJECT_ROOT/src/fortran"
CORE_DIR="$FORTRAN_SRC/core"
CLOUDMASK_DIR="$FORTRAN_SRC/cloudmask"
UTILS_DIR="$FORTRAN_SRC/utils"
CAPI_DIR="$FORTRAN_SRC/c_api"
INCLUDE_DIR="$SCRIPT_DIR/include"

echo "  Fortran source root: $FORTRAN_SRC"

# =============================================================================
# Compiler flags
# =============================================================================
if [ $DEBUG -eq 1 ]; then
    FFLAGS="-g -O0 -fopenmp -fPIC -ffree-line-length-none -fno-range-check -Wall -Wextra -fcheck=all"
    CXXFLAGS="-g -O0 -fPIC -std=c++17 -Wall -Wextra"
    CFLAGS="-g -O0 -fPIC -Wall"
else
    FFLAGS="-O3 -fopenmp -fPIC -ffree-line-length-none -fno-range-check -DNDEBUG"
    CXXFLAGS="-O2 -fPIC -std=c++17 -DNDEBUG"
    CFLAGS="-O2 -fPIC -DNDEBUG"
fi

# =============================================================================
# Create build directories + clean stale .mod
# =============================================================================
mkdir -p "$BUILD_DIR/fortran_modules"
mkdir -p "$BUILD_DIR/obj"

echo ""
echo "=== Cleaning stale .mod files ==="
find "$FORTRAN_SRC" -name "*.mod" -delete 2>/dev/null || true
echo "  Cleaned."

# =============================================================================
# Fortran compilation helpers
# =============================================================================
FINCLUDES="-I$BUILD_DIR/fortran_modules -I$CAPI_DIR -I$CLOUDMASK_DIR -I$CORE_DIR -I$UTILS_DIR"
MODOUT="-J$BUILD_DIR/fortran_modules"

echo ""
echo "=== Compiling Fortran ==="

FCOMPILE() {
    local src="$1"
    local obj="$BUILD_DIR/obj/$(basename "${src%.*}").o"
    if [ ! -f "$obj" ] || [ "$src" -nt "$obj" ]; then
        echo "  FC  $(basename $src)"
        $GFORTRAN $FFLAGS $FINCLUDES $MODOUT -c "$src" -o "$obj"
    fi
    echo "$obj"
}

FCOMPILE_F77() {
    local src="$1"
    local obj="$BUILD_DIR/obj/$(basename "${src%.*}").o"
    if [ ! -f "$obj" ] || [ "$src" -nt "$obj" ]; then
        echo "  FC  $(basename $src) [f77]"
        $GFORTRAN $FFLAGS -ffixed-form $FINCLUDES $MODOUT -c "$src" -o "$obj"
    fi
    echo "$obj"
}

# --- Core modules (order matters for dependencies) ---
FCOMPILE "$CORE_DIR/names_module.f90"
FCOMPILE "$CORE_DIR/data_arrays_module.f90"
FCOMPILE "$CORE_DIR/platform_module.f90"
FCOMPILE "$CORE_DIR/constant.f90"
FCOMPILE "$CORE_DIR/numerical.f90"
FCOMPILE "$CORE_DIR/planck_module.f90"
FCOMPILE "$CORE_DIR/frontend_module.f90"

# --- Modified cloudmask_data_arrays (threadprivate) from c_api/ ---
FCOMPILE "$CAPI_DIR/cloudmask_data_arrays.f90"

# --- Threshold reader ---
FCOMPILE "$CLOUDMASK_DIR/thresholds_read_module.f90"
FCOMPILE_F77 "$CLOUDMASK_DIR/param_read_file.f"

# --- String utilities ---
for f in strlower.f strcompress.f strlen.f strpos.f; do
    FCOMPILE_F77 "$UTILS_DIR/$f"
done

# --- Bit operations ---
for f in set_bit.f clear_bit.f check_bits.f check_qa_bits.f set_qa_bit.f \
         set_confdnc.f set_quality_A.f set_unused_bits.f proc_path.f \
         fill_bit_pixel.f90; do
    case "$f" in
        *.f90) FCOMPILE "$CLOUDMASK_DIR/$f" ;;
        *.f)   FCOMPILE_F77 "$CLOUDMASK_DIR/$f" ;;
    esac
done

# --- Core algorithm ---
for f in conf_test.f conf_test_2val.f pxinit.f tview.f trispc.f \
         get_regdif.f get_regstd.f spatial_var.f chk_spatial_var.f \
         chk_spatial2.f get_sg_thresholds.f90 get_pn_thresholds.f \
         get_nl_thresholds.f90; do
    case "$f" in
        *.f90) FCOMPILE "$CLOUDMASK_DIR/$f" ;;
        *.f)   FCOMPILE_F77 "$CLOUDMASK_DIR/$f" ;;
    esac
done

# --- Restoral tests ---
for f in chk_land.f90 chk_land_nite.f90 chk_coast.f90 chk_sunglint.f90 \
         chk_shallow_water.f shadows.f90 noncld_obs_chk.f90 thin_ci_chk_ir.f90; do
    case "$f" in
        *.f90) FCOMPILE "$CLOUDMASK_DIR/$f" ;;
        *.f)   FCOMPILE_F77 "$CLOUDMASK_DIR/$f" ;;
    esac
done

# --- Spectral tests ---
for f in LandDay.f90 LandNite.f90 LandDay_desert.f90 LandDay_desert_c.f90 \
         LandDay_coast.f90 ocean_day.f90 ocean_nite.f90 Day_snow.f90 \
         Nite_snow.f90 PolarDay_land.f90 PolarDay_desert.f90 \
         PolarDay_desert_c.f90 PolarDay_coast.f90 PolarDay_snow.f90 \
         PolarDay_ocean.f90 PolarNite_land.f90 PolarNite_snow.f90 \
         PolarNite_ocean.f90 Antarctic_day.f90; do
    FCOMPILE "$CLOUDMASK_DIR/$f"
done

# --- Decision tree modules ---
for f in polar_module.f90 land_module.f90 water_module.f90 \
         fylat_fy3mersi_cloud_mask.f90; do
    FCOMPILE "$CLOUDMASK_DIR/$f"
done

# --- C API wrapper (ISO_C_BINDING + OpenMP) ---
FCOMPILE "$CAPI_DIR/cloudmask_c_api.f90"

# =============================================================================
# C compilation
# =============================================================================
echo ""
echo "=== Compiling C ==="

for f in optmed.c optmed_int1.c; do
    echo "  CC  $f"
    $GCC $CFLAGS -c "$UTILS_DIR/$f" -o "$BUILD_DIR/obj/$(basename ${f%.*}).o"
done

# =============================================================================
# C++ pybind11
# =============================================================================
echo ""
echo "=== Compiling C++ pybind11 ==="

echo "  CXX cloudmask_pybind.cpp"
$GXX $CXXFLAGS $PYBIND11_INCLUDES \
    -I"$INCLUDE_DIR" -I"$BUILD_DIR/fortran_modules" \
    -c "$SCRIPT_DIR/cloudmask_pybind.cpp" -o "$BUILD_DIR/obj/cloudmask_pybind.o"

# =============================================================================
# Link
# =============================================================================
echo ""
echo "=== Linking ==="

OUTPUT="$BUILD_DIR/_cloudmask_native${PYTHON_EXT}"
FOBJS=$(find "$BUILD_DIR/obj" -name "*.o" | sort)

# Use $ORIGIN-relative RPATH for portability
$GXX -shared -fPIC $CXXFLAGS -fopenmp \
    $FOBJS \
    -o "$OUTPUT" \
    -L"$HDF5_LIB_DIR" \
    -lhdf5 -lhdf5_fortran \
    -lgfortran -lquadmath \
    -Wl,-rpath,"\$ORIGIN/../lib" \
    -Wl,-rpath,"$HDF5_LIB_DIR"

echo ""
echo "=== Build successful ==="
echo "  Output: $OUTPUT"
echo "  Size:   $(du -h "$OUTPUT" | cut -f1)"

# =============================================================================
# Verify
# =============================================================================
echo ""
echo "=== Verification ==="
echo -n "  Intel symbols: "
nm "$OUTPUT" 2>/dev/null | grep -c "for_cpystr\|for_trim\|for_write" && echo "" || echo "0 (OK)"
echo -n "  libifcore in ldd: "
if ldd "$OUTPUT" 2>/dev/null | grep -qi "libifcore\|libifport"; then
    echo "FOUND (BAD!)"
else
    echo "not found (OK)"
fi
echo -n "  RPATH: "
readelf -d "$OUTPUT" 2>/dev/null | grep RPATH || echo "(none)"

# =============================================================================
# Install
# =============================================================================
if [ $INSTALL -eq 1 ]; then
    echo ""
    echo "=== Installing ==="
    INSTALL_DIR="$PROJECT_ROOT/src/fy3_cloudmask"
    cp "$OUTPUT" "$INSTALL_DIR/"
    echo "  Installed to: $INSTALL_DIR/$(basename $OUTPUT)"

    echo ""
    echo "=== Quick test ==="
    $PYTHON -c "
import os
os.environ.setdefault('FY3_CODE_ROOT', '$PROJECT_ROOT/coeff/')
from fy3_cloudmask.algorithm.native_backend import is_native_available, get_backend_info
print('  Native available:', is_native_available())
info = get_backend_info()
print('  Backend:', info['backend'])
print('  Version:', info['version'])
" 2>&1
fi

echo ""
echo "Done."
