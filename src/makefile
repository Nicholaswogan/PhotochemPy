
FC=gfortran
FCFLAGS=-g -freal-4-real-8# -Wunused-variable# -fbounds-check # I think fbounds-check is finding a false-positive
SRC=photochem_data.f90 \
		photochem_vars.f90 \
		photochem_wrk.f90 \
		photochem_lightning.f90 \
		photochem_clima.f90 \
		photochem.f90 \
		lin_alg.f \
		main.f90 \
		test_all.f90
# OBJ=$(SRC:.f90=.o)
OBJ=photochem_data.o \
		photochem_vars.o \
		photochem_wrk.o \
		photochem_lightning.o \
		photochem_clima.o \
		photochem.o \
		lin_alg.o \
		main.o 
OBJJ=photochem_data.o \
		photochem_vars.o \
		photochem_wrk.o \
		photochem_lightning.o \
		photochem_clima.o \
		photochem.o \
		lin_alg.o \
		test_all.o
	
IDIR=dependencies/include
MODDIR=dependencies/modules
LDIR=dependencies/lib
FFLAGS=-I$(IDIR) -I$(MODDIR) -L$(LDIR) -Idependencies/fortran
LIBS=-lgomp -lm -lsundials_fcvode_mod -lsundials_cvode \
 		  -lyaml -lminpack

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

%.o: %.f
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

photo.run: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(FFLAGS) $(LIBS)

test.all: $(OBJJ)
	gfortran $(FCFLAGS) -o $@ $(OBJJ) $(FFLAGS) $(LIBS)

clean:
	@rm -f *.mod *.o photo.run test.all
	@touch photo.run.dSYM test.all.dSYM
	@rm -r photo.run.dSYM test.all.dSYM
