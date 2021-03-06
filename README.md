
# PhotochemPy [![Memcheck](https://img.shields.io/badge/memcheck-clean-green.svg?style=flat)]() [![Build Status](https://travis-ci.com/Nicholaswogan/PhotochemPy.svg?branch=main)](https://travis-ci.com/Nicholaswogan/PhotochemPy)
PhotochemPy is a photochemical model of rocky planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere, or evolve atmospheres through time.

<!-- [![Documentation Status](https://readthedocs.org/projects/photochempy/badge/?version=latest)](https://photochempy.readthedocs.io/en/latest/?badge=latest) -->

PhotochemPy is a Python wrapper to Fortran source code. This makes the code very speedy, but also user-friendly.

## Installation

**Requirements:**
To install PhotochemPy, you must have the following installed on your system.
- `Python` (>3.6.0) with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection, version >4.9.4 (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

**Python Module:** After satisfying the requirements, then follow these setups to install PhotochemPy

- Clone or download the github repository [`https://github.com/Nicholaswogan/PhotochemPy`](https://github.com/Nicholaswogan/PhotochemPy)
- In a terminal, navigate to the folder `src/dependencies`, and run the shell script `compile.sh`.
- Navigate to the root directory of PhotochemPy, then install with `python -m pip install .`

**Fortran source:** If you prefer to use the code exclusively in Fortran, that is OK too. An example is provided in the folder `examples/fortran_example`.

## Examples/Tutorial

See the `examples` directory. Also check out [this tutorial](https://github.com/Nicholaswogan/PhotochemPy/blob/main/docs/source/Tutorial.ipynb)

<!-- ## Documentation
Read the [documentation here](https://photochempy.readthedocs.io/en/latest/). The best way to get started is [with this Tutorial in the documentation](https://photochempy.readthedocs.io/en/latest/Tutorial.html). -->

## History
PhotochemPy is a distant fork of the [Atmos](https://github.com/VirtualPlanetaryLaboratory/atmos) photochemical model, originally developed by Jim Kasting and Kevin Zahnle and further developed by many of their students and colleges.

## Contact + Publications
If you have questions email me: wogan@uw.edu

Also, if you plan on using PhotochemPy for a publication please email me before you submit anything to a journal, just so we can confirm your planned application of the model is reasonable.
