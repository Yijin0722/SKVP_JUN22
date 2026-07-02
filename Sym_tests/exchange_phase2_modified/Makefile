default: sym_phase2 baseline

FC = gfortran

FFLAGS = -fopenmp -O2 -fimplicit-none -fbacktrace
F77FLAGS = -fopenmp -O2 -fbacktrace

common_objs = mod_gen_parameters.o mod_bsplines.o module_skvp_AtomDiatom.o \
              h4bmkp.o sub_potential_BMKP.o Potential_Interface_3.o \
              sub_basic_aux_mat_calcul.o sub_potential_aux_mat_cacul.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.f
	$(FC) $(F77FLAGS) -c $<

h4bmkp.o: h4bmkp.f
	$(FC) $(F77FLAGS) -c $<

baseline: $(common_objs) skvp_AtomDiatom_baseline.o
	$(FC) $(FFLAGS) -o $@ $^ -L /usr/lib -llapack -lblas

sym_phase1: $(common_objs) skvp_AtomDiatom_sym_phase1.o
	$(FC) $(FFLAGS) -o $@ $^ -L /usr/lib -llapack -lblas

sym_phase2: $(common_objs) skvp_AtomDiatom_sym_phase2.o
	$(FC) $(FFLAGS) -o $@ $^ -L /usr/lib -llapack -lblas

clean:
	rm -f *.o *.mod baseline sym_phase1 sym_phase2
