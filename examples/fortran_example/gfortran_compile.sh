FCFLAGS="-O3 -fopenmp -freal-4-real-8"

# Provide path to Fortran source code.
gfortran -c ../../src/modules/Rainout_vars.f90 $FCFLAGS
gfortran -c ../../src/modules/reading_vars.f90 $FCFLAGS
gfortran -c ../../src/lin_alg.f $FCFLAGS
gfortran -c ../../src/Photochem.f90 $FCFLAGS

# compile the main program
gfortran -c Hadean.f90 $FCFLAGS
gfortran -c Archean.f90 $FCFLAGS
gfortran -c Modern.f90 $FCFLAGS

# link everything
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o Hadean.o -o Hadean.run $FCFLAGS
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o Archean.o -o Archean.run $FCFLAGS
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o Modern.o -o Modern.run $FCFLAGS

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
rm lin_alg.o
rm photochem.mod
rm photochem.o
rm Modern.o
rm Archean.o
rm Hadean.o
