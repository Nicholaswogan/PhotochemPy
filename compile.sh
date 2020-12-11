
# compile sub-modules
gfortran -c modules/Rainout_vars.f90
gfortran -c modules/reading_vars.f90

# compile main module
f2py -c Photochem.f90 lin_alg.f Rainout_vars.o reading_vars.o \
-m Photochem \
only: allocate_memory right_hand_side jacobian read_species read_reactions \
read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet \
read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt \
dochem photo

# clean stuff up
rm rainout_vars.mod Rainout_vars.o
rm reading_vars.mod reading_vars.o
