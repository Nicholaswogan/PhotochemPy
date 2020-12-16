
# compile sub-modules
gfortran -c modules/Rainout_vars.f90
gfortran -c modules/reading_vars.f90

# compile main module
gfortran -c Photochem.f90 lin_alg.f
gfortran -c PhotoMain.f90

# gfortran PhotoMain.o -o testrun

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
rm lin_alg.o
