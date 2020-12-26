# PhotochemPy
PhotochemPy is a photochemical model of rocky planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere.

PhotochemPy is a Python wrapper to the Fortran source code.

## Installation
**Python Module:** The following command will install the PhotochemPy Python package.

`pip install git+git://github.com/Nicholaswogan/PhotochemPy.git`

**Fortran source:** If you prefer to use the code exclusively in Fortran, that is OK too. An example is provided in the file `PhotoMain.f90`, which can be compiled with `./gfortran_compile.sh`.

## Usage
See `Tutorial.ipynb`.

## History
PhotochemPy is an updated version of the [Atmos](https://github.com/VirtualPlanetaryLaboratory/atmos) photochemical model, originally developed by Jim Kasting and Kevin Zahnle and further developed by many of their students and colleges.
