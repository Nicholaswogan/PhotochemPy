
FC=gfortran
FCFLAGS=-g -freal-4-real-8 -O #-fbounds-check # I think fbounds-check is finding a false-positive
SRC=photochem_data.f90 \
		photochem_vars.f90 \
		photochem_wrk.f90 \
		photochem.f90 \
		cvode_funcs.f90	\
		lin_alg.f \
		main.f90
OBJ=$(SRC:.f90=.o)

IDIR=cvode-5.7.0/install/include
LDIR=cvode-5.7.0/install/lib
FFLAGS=-I$(IDIR) -L$(LDIR)
LIBS=-lgomp -lm -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

photo.run: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(FFLAGS) $(LIBS)

clean:
	@rm -f *.mod *.o photo.run
