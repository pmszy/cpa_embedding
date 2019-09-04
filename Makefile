RM = rm -f

#FLAGS= -CB -fast -O2 -132 -r8 -funroll-loops -inline all -g -traceback -fpconstant -fpe3 -fp-stack-check -check bounds
FLAGS= -O2 -132 -r8 -funroll-loops -inline all 
LINK=mpif90 
F90=mpif90
EXE=AVEDISDCA_intel
LFLAGS =
LFLAGS = -r8 -132
#FFLAGS = -openmp -fpp -i4
FFLAGS =
#MKLPATH=/usr/local/packages/lapack/3.4.2/INTEL-140-MVAPICH2-2.0/lib/
#MKLPATH=/opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64
#MKLINCLUDE=/opt/intel/composer_xe_2013_sp1.1.106/mkl/include
SPRNGLIB=~/lib -lsprng -lstdc++
SPRNGINC=~/include


OBJS = mod_tprof.o module_Global.o dateregupg.o geod.o ginit.o \
	main.o meas.o put.o  readin.o sumup.o tables.o 
EXETMDCA=AVEDISDCA_intel

AVEDISDCA: $(OBJS)
	#$(LINK) $(FLAGS) $(OBJS) -L $(SPRNGLIB) -L $(MKLPATH) -I $(SPRNGINC) -llapack -lrefblas $(FFLAGS) $(LFLAGS) -o $(EXE) 	
	$(LINK) $(FLAGS) $(OBJS) -L $(SPRNGLIB) -L $(MKLPATH) -I $(SPRNGINC) $(FFLAGS) $(LFLAGS) -o $(EXE) 	

mod_tprof.o: mod_tprof.F
	$(F90) -c mod_tprof.F

module_Global.o: module_Global.for
	$(F90) $(FLAGS) $(FFLAGS) -c module_Global.for

dateregupg.o: module_Global.for dateregupg.F
	$(F90) $(FLAGS) $(FFLAGS) -L $(SPRNGLIB) -I $(SPRNGINC) -c dateregupg.F

geod.o: module_Global.for geod.F
	$(F90) $(FLAGS) $(FFLAGS) -c geod.F

ginit.o: module_Global.for ginit.F
	$(F90) $(FLAGS) -L $(SPRNGLIB) -I $(SPRNGINC) -c ginit.F
	
	
main.o: module_Global.for main.F
	$(F90) $(FLAGS) $(FFLAGS) -c main.F

meas.o: module_Global.for meas.F
	$(F90) $(FLAGS) $(FFLAGS) -c meas.F

put.o: module_Global.for put.F
	$(F90) $(FLAGS) $(FFLAGS) -c put.F


readin.o: module_Global.for readin.F
	$(F90) $(FLAGS) $(FFLAGS) -L $(SPRNGLIB) -I $(SPRNGINC) -c readin.F

sumup.o: module_Global.for sumup.F
	$(F90) $(FLAGS) $(FFLAGS) -c sumup.F

tables.o: module_Global.for tables.F
	$(F90) $(FLAGS) $(FFLAGS) -c tables.F
	
	

clean:
	$(RM) $(OBJS) $(EXE) $(SPEC)
realclean:
	 $(RM) *.o core dump.F global.mod module_Global.d work.pc *.F *~ *.f 
	

