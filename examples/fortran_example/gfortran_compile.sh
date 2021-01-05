FCFLAGS="-O3 -fopenmp -freal-4-real-8"

# Provide path to Fortran source code.
gfortran -c ../../src/modules/Rainout_vars.f90 $FCFLAGS
gfortran -c ../../src/modules/reading_vars.f90 $FCFLAGS
gfortran -c ../../src/lin_alg.f $FCFLAGS
gfortran -c ../../src/Photochem.f90 $FCFLAGS

# compile the main program
gfortran -c PhotoMain.f90 $FCFLAGS

# link everything
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o PhotoMain.o -o Photo.run $FCFLAGS

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
rm lin_alg.o
rm photochem.mod
rm photochem.o
rm PhotoMain.o
