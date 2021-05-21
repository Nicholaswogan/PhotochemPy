
FC=gfortran
FCFLAGS=-g -freal-4-real-8 -Wunused-variable# -fbounds-check # I think fbounds-check is finding a false-positive
SRC=photochem_data.f90 \
		photochem_vars.f90 \
		photochem_wrk.f90 \
		photochem_clima.f90 \
		photochem.f90 \
		cvode_funcs.f90	\
		lin_alg.f \
		main.f90
OBJ=$(SRC:.f90=.o)
	
IDIR=dependencies/include
MODDIR=dependencies/modules
LDIR=dependencies/lib
FFLAGS=-I$(IDIR) -I$(MODDIR) -L$(LDIR)
LIBS=-lgomp -lm -lsundials_fcvode -lsundials_cvode \
 		 -lsundials_fnvecserial -lsundials_nvecserial -lyaml -lminpack

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

photo.run: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(FFLAGS) $(LIBS)

clean:
	@rm -f *.mod *.o photo.run 
	@touch photo.run.dSYM
	@rm -r photo.run.dSYM
