# Source lists, compilation rules and targets for FastEarth3D.
# Included by config/Makefile after the flag sets are assembled.

# --- Library object list (in module-dependency order) ------------------------
obj_fastearth = \
	$(objdir)/fe_precision.o \
	$(objdir)/fe_constants.o \
	$(objdir)/fe_sht.o \
	$(objdir)/fe_earth_structure.o \
	$(objdir)/fe_radial_integrals.o \
	$(objdir)/fe_lis.o \
	$(objdir)/fe_radial_fe.o \
	$(objdir)/fe_viscoelastic.o \
	$(objdir)/fe_response.o \
	$(objdir)/fe_gravity.o \
	$(objdir)/fe_sle.o \
	$(objdir)/fe_rotation.o \
	$(objdir)/fe_coupling.o \
	$(objdir)/fastearth.o

# --- Inter-module dependencies (so `make -j` stays correct) ------------------
$(objdir)/fe_constants.o:        $(objdir)/fe_precision.o
$(objdir)/fe_sht.o:              $(objdir)/fe_precision.o
$(objdir)/fe_earth_structure.o:  $(objdir)/fe_precision.o $(objdir)/fe_constants.o
$(objdir)/fe_radial_integrals.o: $(objdir)/fe_precision.o
$(objdir)/fe_lis.o:              $(objdir)/fe_precision.o
$(objdir)/fe_radial_fe.o:        $(objdir)/fe_constants.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_radial_integrals.o $(objdir)/fe_lis.o
$(objdir)/fe_viscoelastic.o:     $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o
$(objdir)/fe_response.o:         $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_constants.o
$(objdir)/fe_gravity.o:          $(objdir)/fe_earth_structure.o
$(objdir)/fe_sle.o:              $(objdir)/fe_sht.o $(objdir)/fe_constants.o
$(objdir)/fe_rotation.o:         $(objdir)/fe_sht.o $(objdir)/fe_constants.o
$(objdir)/fe_coupling.o:         $(objdir)/fe_viscoelastic.o $(objdir)/fe_gravity.o \
                                 $(objdir)/fe_sle.o $(objdir)/fe_rotation.o
# The umbrella module re-exports every component, so it compiles last.
$(objdir)/fastearth.o:           $(objdir)/fe_coupling.o

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

test_response: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_response.f90 \
		-o $(bindir)/test_response.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_response.x is ready."

TESTS = test_sht test_earth test_mesh test_integrals test_assembly test_love test_relax test_response

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
