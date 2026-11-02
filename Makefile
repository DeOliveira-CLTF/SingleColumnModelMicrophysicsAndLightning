# Makefile for mp_offline
#
# Builds the two drivers of this project:
#   mp_offline_driver  -- GOAMAZON single-column case (Thompson or WSM6, MYNN optional,
#                          including the unified-EDMF plume core inside MYNN)
#   ensemble_driver     -- ERA5 ensemble x saturation-perturbation case
#
# Requires: gfortran, and a NetCDF-Fortran install findable by nf-config (on the
# PATH). Everything else (Thompson, WSM6, MYNN...) is extracted physics already
# checked into src/ -- nothing is downloaded by this Makefile.
#
# Usage:
#   make                  # builds both drivers into build/
#   make mp_offline_driver
#   make ensemble_driver
#   make clean            # removes build/
#   make DEBUG=1          # -O0 -g -fbacktrace instead of -O2 (faster to rebuild, slower to run)
#
# Then run from build/, e.g.:
#   cd build && ./mp_offline_driver     (reads ./namelist.mp_offline)

FC      := gfortran
BUILD   := build
SRC     := src
MYNN    := $(SRC)/mynn
WSM6    := $(SRC)/wsm6
UEDMF   := $(SRC)/unified_edmf

NETCDF_DIR = $(shell nf-config --prefix)
NETCDF_INC = $(shell nf-config  --includedir)
NETCDF_LIB = $(NETCDF_DIR)/lib
NETCDF_LIB_FLAGS = $(NETCDF_LIB) -lnetcdff -lnetcdf

#NF_FFLAGS := $(shell nf-config --fflags 2>/dev/null)
#NF_FLIBS  := $(shell nf-config --flibs 2>/dev/null)
NF_FFLAGS := -L$(NETCDF_LIB_FLAGS) 
NF_FLIBS  := -L$(NETCDF_LIB) -lnetcdff -lnetcdf
ifeq ($(strip $(NF_FFLAGS)),)
$(warning nf-config not found on PATH -- NetCDF-Fortran flags will be empty; the link step will fail. Install NetCDF-Fortran and make sure nf-config is on PATH.)
endif

ifdef DEBUG
OPT := -O0 -g -fbacktrace
else
OPT := -O2 -fbacktrace
endif

FFLAGS := -cpp -Dmpas -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 $(OPT) $(NF_FFLAGS) -I$(NETCDF_INC) -I$(BUILD) -J$(BUILD)

# ---- object files, in dependency order (see the comment block at the end for why) ----
COMMON_OBJS := \
  $(BUILD)/mpas_stubs_offline.o \
  $(BUILD)/mpas_constants.o \
  $(BUILD)/mpas_atmphys_constants.o \
  $(BUILD)/ccpp_kind_types.o \
  $(BUILD)/mpas_atmphys_functions.o \
  $(BUILD)/module_mp_radar.o \
  $(BUILD)/module_mp_thompson.o \
  $(BUILD)/module_libmassv.o \
  $(BUILD)/mp_wsm6.o \
  $(BUILD)/mp_wsm6_effectRad.o \
  $(BUILD)/module_mp_wsm6.o \
  $(BUILD)/mpas_atmphys_lightning.o \
  $(BUILD)/output_mp_netcdf.o

MYNN_OBJS := \
  $(BUILD)/mynn_shared.o \
  $(BUILD)/bl_mynn_subroutines.o \
  $(BUILD)/bl_mynn.o \
  $(BUILD)/sf_mynn.o \
  $(BUILD)/plume_microphysics_thompson.o \
  $(BUILD)/mynn_offline_coupling.o

MP_DRIVER_OBJS := $(COMMON_OBJS) $(MYNN_OBJS) \
  $(BUILD)/goamazon_bridge.o \
  $(BUILD)/mp_offline_driver.o

ENSEMBLE_DRIVER_OBJS := $(COMMON_OBJS) \
  $(BUILD)/era5_bridge.o \
  $(BUILD)/saturation_perturbation.o \
  $(BUILD)/updraft_forcing.o \
  $(BUILD)/output_ensemble_netcdf.o \
  $(BUILD)/ensemble_driver.o

.PHONY: all clean
all: $(BUILD)/mp_offline_driver $(BUILD)/ensemble_driver

$(BUILD):
	mkdir -p $(BUILD)

# ---- link ----
$(BUILD)/mp_offline_driver: $(MP_DRIVER_OBJS)
	$(FC) $(FFLAGS) -o $@ $(MP_DRIVER_OBJS) $(NF_FLIBS)

$(BUILD)/ensemble_driver: $(ENSEMBLE_DRIVER_OBJS)
	$(FC) $(FFLAGS) -o $@ $(ENSEMBLE_DRIVER_OBJS) $(NF_FLIBS)

# ---- compile (one pattern rule per source directory; order-only deps on $(BUILD)
#      just make sure the directory exists, not the compile order -- that's fixed
#      by the object lists above, which `make` compiles depth-first as needed) ----
$(BUILD)/%.o: $(SRC)/%.F90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@
$(BUILD)/%.o: $(SRC)/%.F | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@
$(BUILD)/%.o: $(MYNN)/%.F90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@
$(BUILD)/%.o: $(MYNN)/%.F | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@
$(BUILD)/%.o: $(WSM6)/%.F90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@
$(BUILD)/%.o: $(WSM6)/%.F | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@
$(BUILD)/%.o: $(UEDMF)/%.F90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

# ---- explicit dependency edges (Fortran `use` statements) so a from-scratch
#      `make -j` builds every file only after the modules it uses exist ----
$(BUILD)/mpas_constants.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/mpas_atmphys_constants.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/mpas_atmphys_functions.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/module_mp_radar.o: $(BUILD)/ccpp_kind_types.o
$(BUILD)/module_mp_thompson.o: $(BUILD)/mpas_stubs_offline.o $(BUILD)/mpas_atmphys_functions.o $(BUILD)/module_mp_radar.o
$(BUILD)/mp_wsm6.o: $(BUILD)/ccpp_kind_types.o $(BUILD)/module_libmassv.o $(BUILD)/module_mp_radar.o
$(BUILD)/mp_wsm6_effectRad.o: $(BUILD)/ccpp_kind_types.o $(BUILD)/mp_wsm6.o
$(BUILD)/module_mp_wsm6.o: $(BUILD)/mpas_stubs_offline.o $(BUILD)/mp_wsm6.o $(BUILD)/mp_wsm6_effectRad.o
$(BUILD)/mpas_atmphys_lightning.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/output_mp_netcdf.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/mynn_shared.o: $(BUILD)/ccpp_kind_types.o
$(BUILD)/bl_mynn_subroutines.o: $(BUILD)/ccpp_kind_types.o
$(BUILD)/bl_mynn.o: $(BUILD)/ccpp_kind_types.o $(BUILD)/bl_mynn_subroutines.o
$(BUILD)/sf_mynn.o: $(BUILD)/ccpp_kind_types.o $(BUILD)/mynn_shared.o
$(BUILD)/plume_microphysics_thompson.o: $(BUILD)/mpas_stubs_offline.o $(BUILD)/module_mp_thompson.o
$(BUILD)/mynn_offline_coupling.o: $(BUILD)/mpas_stubs_offline.o $(BUILD)/bl_mynn_subroutines.o $(BUILD)/bl_mynn.o $(BUILD)/sf_mynn.o
$(BUILD)/goamazon_bridge.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/mp_offline_driver.o: $(filter-out $(BUILD)/mp_offline_driver.o,$(MP_DRIVER_OBJS))
$(BUILD)/era5_bridge.o: $(BUILD)/mpas_stubs_offline.o $(BUILD)/mpas_atmphys_constants.o
$(BUILD)/saturation_perturbation.o: $(BUILD)/mpas_stubs_offline.o $(BUILD)/mpas_atmphys_functions.o
$(BUILD)/updraft_forcing.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/output_ensemble_netcdf.o: $(BUILD)/mpas_stubs_offline.o
$(BUILD)/ensemble_driver.o: $(filter-out $(BUILD)/ensemble_driver.o,$(ENSEMBLE_DRIVER_OBJS))

clean:
	rm -rf $(BUILD)
