# =============================================================================
# Makefile -- FY-3D MERSI-II Cloud Mask Native Engine
#
# Builds the C++/Fortran hybrid library with OpenMP parallelization.
# Usage: make [all|lib|clean|install]
#
# Prerequisites: ifort/gfortran, g++, OpenMP, HDF5, pybind11, numpy
# =============================================================================

# --- Compiler configuration ---
FC       ?= gfortran
CC       ?= gcc
CXX      ?= g++
PYTHON   ?= python3

# --- Paths ---
PROJECT_ROOT = ..
FORTRAN_SRC  = $(PROJECT_ROOT)/src/fortran
CORE_DIR     = $(FORTRAN_SRC)/core
CLOUDMASK_DIR = $(FORTRAN_SRC)/cloudmask
UTILS_DIR    = $(FORTRAN_SRC)/utils
CAPI_DIR     = $(FORTRAN_SRC)/c_api
EXT_DIR      = .
BUILD_DIR    = build

# --- Include paths (c_api/ has priority for modified modules) ---
FINCLUDES = -I$(CAPI_DIR) -I$(CLOUDMASK_DIR) -I$(CORE_DIR) -I$(UTILS_DIR) -J$(BUILD_DIR)

# --- Compiler flags ---
ifeq ($(FC),ifort)
    FFLAGS   = -O3 -qopenmp -fPIC -module $(BUILD_DIR)
    OMP_FLAG = -qopenmp
else
    FFLAGS   = -O3 -fopenmp -fPIC -ffree-line-length-none -fdefault-integer-8
    OMP_FLAG = -fopenmp
endif

CFLAGS   = -O2 -fPIC
CXXFLAGS = -O2 -fPIC -std=c++17

# pybind11 flags
PYBIND11_FLAGS = $(shell $(PYTHON) -m pybind11 --includes 2>/dev/null || \
                   python3 -c "import pybind11; print('-I' + pybind11.get_include())")
PYTHON_EXT     = $(shell $(PYTHON)-config --extension-suffix 2>/dev/null || echo ".so")

# HDF5
HDF5_CFLAGS  = $(shell pkg-config --cflags hdf5 2>/dev/null || echo "-I/usr/include/hdf5/serial")
HDF5_LDFLAGS = $(shell pkg-config --libs hdf5 2>/dev/null || echo "-lhdf5_serial")

# ============================================================================
# Source files
# ============================================================================

# Core modules
CORE_MODS = \
    $(CORE_DIR)/names_module.f90 \
    $(CORE_DIR)/data_arrays_module.f90 \
    $(CORE_DIR)/platform_module.f90 \
    $(CORE_DIR)/constant.f90 \
    $(CORE_DIR)/numerical.f90 \
    $(CORE_DIR)/planck_module.f90 \
    $(CORE_DIR)/frontend_module.f90

# Our modified cloudmask_data_arrays (with threadprivate)
CLOUDMASK_DATA = $(CAPI_DIR)/cloudmask_data_arrays.f90

# Threshold reader
THRESHOLD_READER = \
    $(CLOUDMASK_DIR)/thresholds_read_module.f90 \
    $(CLOUDMASK_DIR)/param_read_file.f

# Bit operations
BIT_OPS = \
    $(CLOUDMASK_DIR)/set_bit.f \
    $(CLOUDMASK_DIR)/clear_bit.f \
    $(CLOUDMASK_DIR)/check_bits.f \
    $(CLOUDMASK_DIR)/check_qa_bits.f \
    $(CLOUDMASK_DIR)/set_qa_bit.f \
    $(CLOUDMASK_DIR)/set_confdnc.f \
    $(CLOUDMASK_DIR)/set_quality_A.f \
    $(CLOUDMASK_DIR)/set_unused_bits.f \
    $(CLOUDMASK_DIR)/proc_path.f \
    $(CLOUDMASK_DIR)/fill_bit_pixel.f90

# Core algorithm
ALG_CORE = \
    $(CLOUDMASK_DIR)/conf_test.f \
    $(CLOUDMASK_DIR)/conf_test_2val.f \
    $(CLOUDMASK_DIR)/pxinit.f \
    $(CLOUDMASK_DIR)/tview.f \
    $(CLOUDMASK_DIR)/trispc.f \
    $(CLOUDMASK_DIR)/get_regdif.f \
    $(CLOUDMASK_DIR)/get_regstd.f \
    $(CLOUDMASK_DIR)/spatial_var.f \
    $(CLOUDMASK_DIR)/chk_spatial_var.f \
    $(CLOUDMASK_DIR)/chk_spatial2.f \
    $(CLOUDMASK_DIR)/get_sg_thresholds.f90 \
    $(CLOUDMASK_DIR)/get_pn_thresholds.f \
    $(CLOUDMASK_DIR)/get_nl_thresholds.f90

# Restoral tests
RESTORAL = \
    $(CLOUDMASK_DIR)/chk_land.f90 \
    $(CLOUDMASK_DIR)/chk_land_nite.f90 \
    $(CLOUDMASK_DIR)/chk_coast.f90 \
    $(CLOUDMASK_DIR)/chk_sunglint.f90 \
    $(CLOUDMASK_DIR)/chk_shallow_water.f \
    $(CLOUDMASK_DIR)/shadows.f90 \
    $(CLOUDMASK_DIR)/noncld_obs_chk.f90 \
    $(CLOUDMASK_DIR)/thin_ci_chk_ir.f90

# Spectral tests
SPECTRAL = \
    $(CLOUDMASK_DIR)/LandDay.f90 \
    $(CLOUDMASK_DIR)/LandNite.f90 \
    $(CLOUDMASK_DIR)/LandDay_desert.f90 \
    $(CLOUDMASK_DIR)/LandDay_desert_c.f90 \
    $(CLOUDMASK_DIR)/LandDay_coast.f90 \
    $(CLOUDMASK_DIR)/ocean_day.f90 \
    $(CLOUDMASK_DIR)/ocean_nite.f90 \
    $(CLOUDMASK_DIR)/Day_snow.f90 \
    $(CLOUDMASK_DIR)/Nite_snow.f90 \
    $(CLOUDMASK_DIR)/PolarDay_land.f90 \
    $(CLOUDMASK_DIR)/PolarDay_desert.f90 \
    $(CLOUDMASK_DIR)/PolarDay_desert_c.f90 \
    $(CLOUDMASK_DIR)/PolarDay_coast.f90 \
    $(CLOUDMASK_DIR)/PolarDay_snow.f90 \
    $(CLOUDMASK_DIR)/PolarDay_ocean.f90 \
    $(CLOUDMASK_DIR)/PolarNite_land.f90 \
    $(CLOUDMASK_DIR)/PolarNite_snow.f90 \
    $(CLOUDMASK_DIR)/PolarNite_ocean.f90 \
    $(CLOUDMASK_DIR)/Antarctic_day.f90

# Decision tree modules
DECISION = \
    $(CLOUDMASK_DIR)/polar_module.f90 \
    $(CLOUDMASK_DIR)/land_module.f90 \
    $(CLOUDMASK_DIR)/water_module.f90 \
    $(CLOUDMASK_DIR)/fylat_fy3mersi_cloud_mask.f90

# C API wrapper
C_API = $(CAPI_DIR)/cloudmask_c_api.f90

# C source
C_SRC = $(UTILS_DIR)/optmed.c

# All Fortran sources
ALL_FSRC = $(CORE_MODS) $(CLOUDMASK_DATA) $(THRESHOLD_READER) $(BIT_OPS) \
           $(ALG_CORE) $(RESTORAL) $(SPECTRAL) $(DECISION) $(C_API)

# ============================================================================
# Object files
# ============================================================================
FOBJS  = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(filter %.f90,$(ALL_FSRC)))
FOBJS += $(patsubst %.f,$(BUILD_DIR)/%.o,$(filter %.f,$(ALL_FSRC)))
COBJS  = $(patsubst %.c,$(BUILD_DIR)/%.o,$(C_SRC))
CXXOBJ = $(BUILD_DIR)/cloudmask_pybind.o

# ============================================================================
# Targets
# ============================================================================

.PHONY: all lib clean install

all: lib

lib: $(BUILD_DIR)/_cloudmask_native$(PYTHON_EXT)

# Create build directory structure
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)/$(CORE_DIR) $(BUILD_DIR)/$(CLOUDMASK_DIR) $(BUILD_DIR)/$(CAPI_DIR)

# Fortran compilation rules
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(FINCLUDES) -c $< -o $@

$(BUILD_DIR)/%.o: %.f | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(FINCLUDES) -c $< -o $@

# C compilation
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# C++ pybind11 compilation
$(BUILD_DIR)/cloudmask_pybind.o: cloudmask_pybind.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(PYBIND11_FLAGS) -Iinclude -I$(BUILD_DIR) -c $< -o $@

# Link shared library
$(BUILD_DIR)/_cloudmask_native$(PYTHON_EXT): $(FOBJS) $(COBJS) $(CXXOBJ)
	$(CXX) -shared -fPIC $(CXXFLAGS) $(OMP_FLAG) \
		$^ \
		-o $@ \
		$(HDF5_LDFLAGS) \
		-lgfortran -lquadmath 2>/dev/null || \
	$(CXX) -shared -fPIC $(CXXFLAGS) $(OMP_FLAG) \
		$(FOBJS) $(COBJS) $(CXXOBJ) \
		-o $@ \
		-L$(shell $(FC) -print-file-name=.) \
		$(HDF5_LDFLAGS)

# Install
install: lib
	cp $(BUILD_DIR)/_cloudmask_native$(PYTHON_EXT) ../src/fy3_cloudmask/

# Clean
clean:
	rm -rf $(BUILD_DIR)
