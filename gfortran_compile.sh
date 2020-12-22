
# compile sub-modules
gfortran -c src/modules/Rainout_vars.f90
gfortran -c src/modules/reading_vars.f90

# compile main module
gfortran -c src/Photochem.f90 src/lin_alg.f


# gfortran PhotoMain.o -o testrun

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
rm lin_alg.o
rm photochem.mod
rm photochem.o
