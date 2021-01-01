gfortran -c src/modules/Rainout_vars.f90 -O3 -fopenmp
gfortran -c src/modules/reading_vars.f90 -O3 -fopenmp
gfortran -c src/lin_alg.f -O3 -fopenmp
gfortran -c src/Photochem.f90 -O3 -fopenmp -freal-4-real-8
gfortran -c PhotoMain.f90 -O3 -fopenmp
gfortran Rainout_vars.o reading_vars.o Photochem.o lin_alg.o PhotoMain.o -o Photo.run -O3 -fopenmp -freal-4-real-8

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
rm lin_alg.o
rm photochem.mod
rm photochem.o
rm PhotoMain.o
