# Source lists, compilation rules and targets for FastEarth3D.
# Included by config/Makefile after the flag sets are assembled.

# --- Library object list (in module-dependency order) ------------------------
obj_fastearth = \
	$(objdir)/fe_precision.o \
	$(objdir)/fe_constants.o \
	$(objdir)/fe_params.o \
	$(objdir)/fe_sht.o \
	$(objdir)/fe_field.o \
	$(objdir)/fe_earth_structure.o \
	$(objdir)/fe_radial_integrals.o \
	$(objdir)/fe_band.o \
	$(objdir)/fe_radial_fe.o \
	$(objdir)/fe_viscoelastic.o \
	$(objdir)/fe_response.o \
	$(objdir)/fe_gravity.o \
	$(objdir)/fe_sle.o \
	$(objdir)/fe_timestep.o \
	$(objdir)/fe_rotation.o \
	$(objdir)/fe_coupling.o \
	$(objdir)/fe_io.o \
	$(objdir)/fastearth.o

# --- Inter-module dependencies (so `make -j` stays correct) ------------------
$(objdir)/fe_constants.o:        $(objdir)/fe_precision.o
$(objdir)/fe_params.o:           $(objdir)/fe_precision.o $(objdir)/fe_constants.o
$(objdir)/fe_sht.o:              $(objdir)/fe_precision.o
$(objdir)/fe_field.o:            $(objdir)/fe_precision.o $(objdir)/fe_sht.o
$(objdir)/fe_earth_structure.o:  $(objdir)/fe_precision.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_params.o
$(objdir)/fe_radial_integrals.o: $(objdir)/fe_precision.o
$(objdir)/fe_band.o:             $(objdir)/fe_precision.o
$(objdir)/fe_radial_fe.o:        $(objdir)/fe_constants.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_radial_integrals.o $(objdir)/fe_band.o
$(objdir)/fe_viscoelastic.o:     $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o
$(objdir)/fe_response.o:         $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_viscoelastic.o
$(objdir)/fe_gravity.o:          $(objdir)/fe_earth_structure.o
$(objdir)/fe_sle.o:              $(objdir)/fe_sht.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_response.o
$(objdir)/fe_timestep.o:         $(objdir)/fe_response.o $(objdir)/fe_sle.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_viscoelastic.o \
                                 $(objdir)/fe_precision.o
$(objdir)/fe_rotation.o:         $(objdir)/fe_sht.o $(objdir)/fe_constants.o
$(objdir)/fe_coupling.o:         $(objdir)/fe_response.o $(objdir)/fe_sle.o \
                                 $(objdir)/fe_rotation.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_params.o \
                                 $(objdir)/fe_timestep.o $(objdir)/fe_viscoelastic.o
$(objdir)/fe_io.o:               $(objdir)/fe_coupling.o $(objdir)/fe_response.o \
                                 $(objdir)/fe_viscoelastic.o $(objdir)/fe_sht.o \
                                 $(objdir)/fe_constants.o
# The umbrella module re-exports every component, so it compiles last.
$(objdir)/fastearth.o:           $(objdir)/fe_coupling.o $(objdir)/fe_io.o

# --- Pattern rule ------------------------------------------------------------
$(objdir)/%.o: $(srcdir)/%.f90 | $(objdir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

$(objdir):
	mkdir -p $(objdir)

$(bindir):
	mkdir -p $(bindir)

# --- Library -----------------------------------------------------------------
fastearth-static: $(obj_fastearth)
	ar rcs $(objdir)/libfastearth.a $(obj_fastearth)
	@echo ""
	@echo "    $(objdir)/libfastearth.a is ready."
	@echo ""

# --- Tests -------------------------------------------------------------------
test_params: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_params.f90 \
		-o $(bindir)/test_params.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_params.x is ready."

test_band: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_band.f90 \
		-o $(bindir)/test_band.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_band.x is ready."

test_sht: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sht.f90 \
		-o $(bindir)/test_sht.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sht.x is ready."

test_earth: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_earth.f90 \
		-o $(bindir)/test_earth.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_earth.x is ready."

test_mesh: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_mesh.f90 \
		-o $(bindir)/test_mesh.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_mesh.x is ready."

test_integrals: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_integrals.f90 \
		-o $(bindir)/test_integrals.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_integrals.x is ready."

test_assembly: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_assembly.f90 \
		-o $(bindir)/test_assembly.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_assembly.x is ready."

test_love: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_love.f90 \
		-o $(bindir)/test_love.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_love.x is ready."

test_relax: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_relax.f90 \
		-o $(bindir)/test_relax.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_relax.x is ready."

test_etd1: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_etd1.f90 \
		-o $(bindir)/test_etd1.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_etd1.x is ready."

# Coupling-order characterization (§3c): measures the convergence order of the
# strain<->memory coupling. Standalone like test_etd1 -- a dt-sweep diagnostic,
# intentionally NOT in TESTS / `make check`. Build + run directly.
test_couple_order: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_couple_order.f90 \
		-o $(bindir)/test_couple_order.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_couple_order.x is ready."

# SLE<->memory coupling-order characterization (§3c 3b): drives a fast-evolving load
# through the full sea-level driver and measures the order restored by co-converging
# the ocean load σ and the end-of-step memory τ. Standalone dt-sweep diagnostic, like
# test_couple_order -- intentionally NOT in TESTS / `make check`. Build + run directly.
test_sle_couple_order: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_couple_order.f90 \
		-o $(bindir)/test_sle_couple_order.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_couple_order.x is ready."

# Adaptive-Δt controller (§3c): field step-doubling estimate order + the adaptive
# stepper converging to a fine reference with far fewer steps. Standalone diagnostic,
# NOT in `make check`.
test_timestep: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_timestep.f90 \
		-o $(bindir)/test_timestep.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_timestep.x is ready."

test_response: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_response.f90 \
		-o $(bindir)/test_response.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_response.x is ready."

test_sle: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle.f90 \
		-o $(bindir)/test_sle.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle.x is ready."

test_flotation: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_flotation.f90 \
		-o $(bindir)/test_flotation.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_flotation.x is ready."

test_ve_response: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_ve_response.f90 \
		-o $(bindir)/test_ve_response.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_ve_response.x is ready."

test_sle_ve: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_ve.f90 \
		-o $(bindir)/test_sle_ve.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_ve.x is ready."

test_benchmark_love: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_love.f90 \
		-o $(bindir)/test_benchmark_love.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_love.x is ready."

test_coupling: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_coupling.f90 \
		-o $(bindir)/test_coupling.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_coupling.x is ready."

test_restart: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_restart.f90 \
		-o $(bindir)/test_restart.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_restart.x is ready."

test_benchmark_disc: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_disc.f90 \
		-o $(bindir)/test_benchmark_disc.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_disc.x is ready."

test_benchmark_martinec: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_martinec.f90 \
		-o $(bindir)/test_benchmark_martinec.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_martinec.x is ready."

test_field: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_field.f90 \
		-o $(bindir)/test_field.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_field.x is ready."

test_flotation_load: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_flotation_load.f90 \
		-o $(bindir)/test_flotation_load.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_flotation_load.x is ready."

# Standalone SLE benchmark (Martinec-2018 case E2): ~750 steps at lmax=128, runs
# in minutes -- intentionally NOT in TESTS / `make check`. Build with `make
# openmp=1 test_benchmark_sle` and run $(bindir)/test_benchmark_sle.x directly.
test_benchmark_sle: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_sle.f90 \
		-o $(bindir)/test_benchmark_sle.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_sle.x is ready."

test_sle_eustatic: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_eustatic.f90 \
		-o $(bindir)/test_sle_eustatic.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_eustatic.x is ready."

test_sle_subgrid: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_subgrid.f90 \
		-o $(bindir)/test_sle_subgrid.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_subgrid.x is ready."

TESTS = test_params test_band test_sht test_earth test_mesh test_integrals test_assembly test_love test_relax test_response test_sle test_flotation test_flotation_load test_ve_response test_sle_ve test_benchmark_love test_coupling test_restart test_benchmark_disc test_benchmark_martinec test_field test_sle_subgrid

check: $(TESTS)
	@echo ""
	@echo "=== Running FastEarth3D test suite ==="
	@for t in $(TESTS); do \
		echo "--- $$t ---"; \
		$(bindir)/$$t.x || exit 1; \
	done
	@echo ""
	@echo "=== All tests passed ==="

# --- Housekeeping ------------------------------------------------------------
.PHONY: usage check clean showconfig

usage:
	@echo ""
	@echo "    * FastEarth3D build *"
	@echo ""
	@echo " make fastearth-static : build libfastearth.a"
	@echo " make check            : build + run the test suite"
	@echo " make test_sht         : build the SHT round-trip test"
	@echo " make clean            : remove objects and binaries"
	@echo " make showconfig       : show the active build configuration"
	@echo ""
	@echo "   switches:  debug=0|1|2   openmp=0|1"
	@echo ""

showconfig:
	@echo "----------------------"
	@echo "FastEarth3D build configuration"
	@echo "----------------------"
	@echo "compiler  : $(FC)"
	@echo "host      : $(shell hostname)"
	@echo "openmp    : $(openmp)"
	@echo "debug     : $(debug)"
	@echo "FFLAGS    : $(FFLAGS)"
	@echo "LFLAGS    : $(LFLAGS)"

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a
	rm -f $(bindir)/*.x
	rm -rf $(bindir)/*.dSYM
