# gfortran -c src/photochem_clima.f90 -O3 -lminpack

python -m numpy.f2py -c src/photochem_data.f90 \
                    		src/photochem_vars.f90 \
                    		src/photochem_wrk.f90 \
                        src/photochem_lightning.f90 \
                    		src/photochem_clima.f90 \
                        src/photochem.f90 \
                    		src/cvode_funcs.f90	\
                    		src/lin_alg.f \
-m Photochem \
--opt="-O3" \
--f90flags='-freal-4-real-8 -fopenmp' \
--f77flags='-freal-4-real-8 -fopenmp' \
-lgomp \
-Isrc/dependencies/include \
-Isrc/dependencies/modules \
-Lsrc/dependencies/lib \
-lm \
-lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial -lyaml -lminpack \
only: setup right_hand_side jacobian \
integrate cvode cvode_save cvode_equilibrium \
steam2photochem

# rm photochem_clima.mod photochem_clima.o
# python -m numpy.f2py -c src/modules/Rainout_vars.f90 src/modules/reading_vars.f90 src/Photochem.f90 src/lin_alg.f src/cvode_funcs.f90 \
# -m Photochem \
# --opt="-O3" \
# --f90flags='-freal-4-real-8 ' \
# --f77flags='-freal-4-real-8' \
# -Isrc/cvode-5.7.0/install/include \
# -Lsrc/cvode-5.7.0/install/lib \
# -lm \
# -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial \
# only: allocate_memory right_hand_side jacobian read_species read_reactions \
# read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet \
# read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt \
# dochem photo setup integrate cvode
