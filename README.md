
# PhotochemPy [![Documentation Status](https://readthedocs.org/projects/photochempy/badge/?version=latest)](https://photochempy.readthedocs.io/en/latest/?badge=latest) [![Memcheck](https://img.shields.io/badge/memcheck-clean-green.svg?style=flat)]() [![Build Status](https://travis-ci.com/Nicholaswogan/PhotochemPy.svg?branch=main)](https://travis-ci.com/Nicholaswogan/PhotochemPy)
PhotochemPy is a photochemical model of rocky planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere.

PhotochemPy is a Python wrapper to Fortran source code. This makes the code very speedy, but also user-friendly.

## Installation

**Requirements:**
To install PhotochemPy, you must have the following installed on your system.
- `Python` with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

**Python Module:** After satisfying the requirements, the following command will install the PhotochemPy package to your Python installation.

`python -m pip install git+git://github.com/Nicholaswogan/PhotochemPy.git`

Alternatively, if you want maximum performance, then install PhotochemPy with the [Spike](http://www.ecs.umass.edu/~polizzi/spike/index.htm) solver following [these instructions](https://photochempy.readthedocs.io/en/latest/install.html#python).

**Fortran source:** If you prefer to use the code exclusively in Fortran, that is OK too. An example is provided in the folder `examples/fortran_example`.

## Documentation
Read the [documentation here](https://photochempy.readthedocs.io/en/latest/). The best way to get started is [with this Tutorial in the documentation](https://photochempy.readthedocs.io/en/latest/Tutorial.html).

## History
PhotochemPy is an updated version of the [Atmos](https://github.com/VirtualPlanetaryLaboratory/atmos) photochemical model, originally developed by Jim Kasting and Kevin Zahnle and further developed by many of their students and colleges.

## Code integrity
I regularly check for memory issues (uninitialized variables, accessing an array beyond its bounds) with the [valgrind memcheck](http://valgrind.org) tool, and make sure all errors are dealt with.
