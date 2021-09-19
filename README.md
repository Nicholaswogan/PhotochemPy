
# PhotochemPy [![Memcheck](https://img.shields.io/badge/memcheck-clean-green.svg?style=flat)]() [![Build Status](https://travis-ci.com/Nicholaswogan/PhotochemPy.svg?branch=main)](https://travis-ci.com/Nicholaswogan/PhotochemPy)
PhotochemPy is a photochemical model of planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere, or evolve atmospheres through time.

<!-- [![Documentation Status](https://readthedocs.org/projects/photochempy/badge/?version=latest)](https://photochempy.readthedocs.io/en/latest/?badge=latest) -->

PhotochemPy is a Python wrapper to Fortran source code. This makes the code very speedy, but also user-friendly.

## Installation

**Requirements:**
- MacOS or any Linux OS.
- `Python` (3.7, 3.8, or 3.9)
- Fortran and C compiler. I suggest the GNU compiler collection (includes `gfortran`, `gcc`, etc.). If you are using a Mac, install it with Homebrew: `brew install gcc.`

**Python Module:** 

- Clone or download the github repository.
- Navigate to the root directory of PhotochemPy, then install with `python -m pip install .`

**Fortran library:** If you prefer to use the code exclusively in Fortran, that is OK too. You can build `libphotochem` with CMake. Download or clone this respository, then from the root directory of the repository run

```sh
mkdir build
cd build
cmake ..
make -j
```

## Examples/Tutorial

See the `examples` directory. Also check out [this tutorial](https://github.com/Nicholaswogan/PhotochemPy/blob/main/docs/source/Tutorial.ipynb)

## History
PhotochemPy is a distant fork of the [Atmos](https://github.com/VirtualPlanetaryLaboratory/atmos) photochemical model, originally developed by Jim Kasting and Kevin Zahnle and further developed by many of their students and colleges.

## Contact + Publications
If you have questions email me: wogan@uw.edu

Also, if you plan on using PhotochemPy for a publication please email me before you submit anything to a journal, just so we can confirm your planned application of the model is reasonable.
