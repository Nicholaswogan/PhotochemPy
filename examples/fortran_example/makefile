
FC=gfortran
FCFLAGS=-g -freal-4-real-8 -g -fopenmp -O3
SRC=../../src/photochem_data.f90 \
		../../src/photochem_vars.f90 \
		../../src/photochem_wrk.f90 \
		../../src/photochem_clima.f90 \
		../../src/photochem.f90 \
		../../src/cvode_funcs.f90	\
		../../src/lin_alg.f \
		example.f90
		
OBJ=photochem_data.o \
		photochem_vars.o \
		photochem_wrk.o \
		photochem_clima.o \
		photochem.o \
		cvode_funcs.o	\
		lin_alg.o \
		example.o
	
IDIR=../../src/dependencies/include
MODDIR=../../src/dependencies/modules
LDIR=../../src/dependencies/lib
FFLAGS=-I$(IDIR) -I$(MODDIR) -L$(LDIR)
LIBS=-lgomp -lm -lsundials_fcvode -lsundials_cvode \
 		 -lsundials_fnvecserial -lsundials_nvecserial -lyaml -lminpack
		 
%.o: ../../src/%.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)
	
%.o: ../../src/%.f
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

photo.run: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(FFLAGS) $(LIBS)

clean:
	@rm -f *.mod *.o photo.run
	@touch photo.run.dSYM
	@rm -r photo.run.dSYM
