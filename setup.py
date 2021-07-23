
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_data_dir(('PhotochemPy/DATA','PhotochemPy/DATA'))
    return config

from numpy.distutils.core import setup, Extension
import numpy as np
# from subprocess import call
# 
# comile dependencies
# err = call('sh dependencies.sh'.split())
# if err:
#     raise Exception('The dependencies failed to compile.')

version = '0.2.0'

only = 'only: setup right_hand_side jacobian'+ \
        ' integrate cvode cvode_save cvode_equilibrium steam2photochem :'
        
sources = ['src/photochem_data.f90', \
           'src/photochem_vars.f90', \
           'src/photochem_wrk.f90', \
           'src/photochem_lightning.f90', \
           'src/photochem_clima.f90', \
           'src/photochem.f90', \
           'src/cvode_funcs.f90', \
           'src/lin_alg.f']

option = 1 # (Default) Parallel version (fast)
# option = 2 # Serial version (slow)

if option == 1: # installing with parallel computation
    extensions = [
    Extension(name="Photochem",
              sources=sources,
              extra_f90_compile_args = ['-O3','-fopenmp', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              libraries=['gomp','m','sundials_fcvode','sundials_cvode','sundials_fnvecserial','sundials_nvecserial','yaml','minpack'],
              library_dirs =["src/dependencies/lib"],
              include_dirs = ["src/dependencies/include","src/dependencies/modules"],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          python_requires='>3.6.0',
          packages=['PhotochemPy'],
          version=version,
          ext_modules=extensions,
          install_requires=['pathos'],
          configuration=configuration)

if option == 2: # istalling with serial compuation
    extensions = [
    Extension(name="Photochem",
              sources=sources,
              extra_f90_compile_args = ['-O3', '-freal-4-real-8'],
              extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
              libraries=['m','sundials_fcvode','sundials_cvode','sundials_fnvecserial','sundials_nvecserial','yaml'],
              library_dirs =["src/dependencies/lib"],
              include_dirs = ["src/dependencies/include","src/dependencies/modules"],
              f2py_options=only.split())
              ]

    setup(name = 'PhotochemPy',
          python_requires='>3.6.0',
          packages=['PhotochemPy'],
          version=version,
          ext_modules=extensions,
          install_requires=['pathos'],
          configuration=configuration)
