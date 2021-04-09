
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
# option = 3 # Parallel with Spike solver (super duper fast)

# For option 3 follow these instructions: https://photochempy.readthedocs.io/en/latest/install.html

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
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          python_requires='>3.6.0',
          packages=['PhotochemPy'],
          version=version,
          ext_modules=extensions,
          configuration=configuration)


if option == 3: # installing with spike (fastest! But hard to install)
    # spike directory
    spike_lib = 'src/spike-1.0/lib/x64/'

    # install libraries
    install1 = "python -m pip install mkl"
    install2 = "python -m pip install mkl-include"
    install3 = "python -m pip install intel-openmp"
    subprocess.run(install1.split())
    subprocess.run(install2.split())
    subprocess.run(install3.split())

    # This should find where the mkl libraries are installed (hopefully!)
    cmd = 'python -m pip show mkl'
    out = subprocess.check_output(cmd.split()).decode("utf-8")
    ind = [i for i,a in enumerate(out.split()) if a=='Location:'][0]+1
    mkl_location = out.split()[ind]
    root_dir = mkl_location[:mkl_location.find('/lib/python')]
    mkl_lib = root_dir+'/lib'
    mkl_include = root_dir+'/include'

    # must have numpy >1.16.0
    if float('.'.join(np.__version__.split('.')[:2]))<1.16:
        sys.exit("Must have numpy > 1.16.0")

    extensions = [
    Extension(name="Photochem",
              sources=['src/modules/Rainout_vars.f90', 'src/modules/reading_vars.f90','src/Photochem_spike.f90','src/lin_alg.f'],
              extra_f90_compile_args = ['-O3','-fopenmp', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              libraries=['spike','mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'iomp5', 'pthread', 'm', 'dl'],
              runtime_library_dirs=[mkl_lib],
              library_dirs =[mkl_lib,spike_lib],
              include_dirs = [mkl_include],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          python_requires='>3.6.0',
          packages=['PhotochemPy'],
          version=version,
          ext_modules=extensions,
          configuration=configuration)
