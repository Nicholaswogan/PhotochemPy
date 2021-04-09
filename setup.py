
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_data_dir(('PhotochemPy/DATA','PhotochemPy/DATA'))
    return config

from numpy.distutils.core import setup, Extension
import numpy as np
import subprocess
import os
import sys

version = '0.1.0'

only = '''only: allocate_memory right_hand_side jacobian read_species read_reactions
 read_atmosphere photgrid rates gridw readflux initphoto initmie read_planet
 read_photochem rainout ltning aertab densty aercon photsatrat difco sedmnt
 dochem photo setup integrate cvode :'''

option = 1 # (Default) Parallel version (fast)
# option = 2 # Serial version (slow)

if option == 1: # installing with parallel computation
    extensions = [
    Extension(name="Photochem",
              sources=['src/modules/Rainout_vars.f90', 'src/modules/reading_vars.f90','src/Photochem.f90','src/lin_alg.f','src/cvode_funcs.f90'],
              extra_f90_compile_args = ['-O3','-fopenmp', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              libraries=['gomp','m','sundials_fcvode','sundials_cvode','sundials_fnvecserial','sundials_nvecserial'],
              library_dirs =["src/cvode-5.7.0/install/lib"],
              include_dirs = ["src/cvode-5.7.0/install/include"],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          python_requires='>3.6.0',
          packages=['PhotochemPy'],
          version=version,
          ext_modules=extensions,
          configuration=configuration)

if option == 2: # istalling with serial compuation
    extensions = [
    Extension(name="Photochem",
              sources=['src/modules/Rainout_vars.f90', 'src/modules/reading_vars.f90','src/Photochem.f90','src/lin_alg.f'],
              extra_f90_compile_args = ['-O3', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              libraries=['m','sundials_fcvode','sundials_cvode','sundials_fnvecserial','sundials_nvecserial'],
              library_dirs =["src/cvode-5.7.0/install/lib"],
              include_dirs = ["src/cvode-5.7.0/install/include"],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          python_requires='>3.6.0',
          packages=['PhotochemPy'],
          version=version,
          ext_modules=extensions,
          configuration=configuration)
