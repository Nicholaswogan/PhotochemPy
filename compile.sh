
python -m numpy.f2py -c src/modules/Rainout_vars.f90 src/modules/reading_vars.f90 src/Photochem.f90 src/lin_alg.f \
-m Photochem \
--opt="-O3" \
--f90flags='-fopenmp -freal-4-real-8' \
--f77flags='-fopenmp -freal-4-real-8' \
-lgomp \
only: allocate_memory right_hand_side jacobian read_species read_reactions \
read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet \
read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt \
dochem photo setup integrate

# f2py -c src/modules/Rainout_vars.f90 src/modules/reading_vars.f90 src/Photochem.f90 src/lin_alg.f \
# -m Photochem \
# --opt="-O3" \
# --f90flags='-freal-4-real-8 ' \
# --f77flags='-freal-4-real-8' \
# only: allocate_memory right_hand_side jacobian read_species read_reactions \
# read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet \
# read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt \
# dochem photo setup integrate
