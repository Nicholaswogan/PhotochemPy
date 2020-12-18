
# compile sub-modules
gfortran -c src/modules/Rainout_vars.f90
gfortran -c src/modules/reading_vars.f90

# compile main module
f2py -c src/Photochem.f90 src/lin_alg.f Rainout_vars.o reading_vars.o \
-m Photochem \
--opt="-O3" \
--f90flags='-fopenmp -freal-4-real-8' \
--f77flags='-freal-4-real-8' \
-lgomp \
only: allocate_memory right_hand_side jacobian read_species read_reactions \
read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet \
read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt \
dochem photo setup integrate

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
