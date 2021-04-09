FCFLAGS="-O3 -freal-4-real-8"
OTHERFLAGS="-I../../src/cvode-5.7.0/install/include -L../../src/cvode-5.7.0/install/lib -lm -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial"

# Provide path to Fortran source code.
gfortran -c ../../src/modules/Rainout_vars.f90 $FCFLAGS
gfortran -c ../../src/modules/reading_vars.f90 $FCFLAGS
gfortran -c ../../src/lin_alg.f $FCFLAGS
gfortran -c ../../src/Photochem.f90 $FCFLAGS
gfortran -c ../../src/cvode_funcs.f90 $FCFLAGS

# compile the main program
gfortran -c Hadean.f90 $FCFLAGS
gfortran -c Archean.f90 $FCFLAGS
gfortran -c Modern.f90 $FCFLAGS
gfortran -c Archean2Proterozoic.f90 $FCFLAGS

# link everything
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o cvode_funcs.o Hadean.o -o Hadean.run $FCFLAGS $OTHERFLAGS
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o cvode_funcs.o Archean.o -o Archean.run $FCFLAGS $OTHERFLAGS
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o cvode_funcs.o Modern.o -o Modern.run $FCFLAGS $OTHERFLAGS
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o cvode_funcs.o Archean2Proterozoic.o -o Archean2Proterozoic.run $FCFLAGS $OTHERFLAGS


# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
rm lin_alg.o
rm photochem.mod
rm photochem.o
rm Modern.o
rm Archean.o
rm Hadean.o
