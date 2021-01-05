
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_data_dir(('PhotochemPy/DATA','PhotochemPy/DATA'))
    return config



from numpy.distutils.core import setup, Extension
import subprocess


only = '''only: allocate_memory right_hand_side jacobian read_species read_reactions
 read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet
 read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt
 dochem photo setup integrate :'''

option = 1 # Parallel and faster with openmp.
# option = 2 # Serial and slower.

if option == 1:
    extensions = [
    Extension(name="Photochem",
              sources=['src/modules/Rainout_vars.f90', 'src/modules/reading_vars.f90','src/Photochem.f90','src/lin_alg.f'],
              extra_f90_compile_args = ['-O3','-fopenmp', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              libraries=['gomp'],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          packages=['PhotochemPy'],
          version='0.0.1',
          ext_modules=extensions,
          configuration=configuration)

if option == 2:
    extensions = [
    Extension(name="Photochem",
              sources=['src/modules/Rainout_vars.f90', 'src/modules/reading_vars.f90','src/Photochem.f90','src/lin_alg.f'],
              extra_f90_compile_args = ['-O3', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          packages=['PhotochemPy'],
          version='0.0.1',
          ext_modules=extensions,
          configuration=configuration)
