# Installation

**Requirements:**
To install PhotochemPy, you must have the following installed on your system.
- `Python` with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

## Python
After satisfying the requirements, the following command will install the PhotochemPy package to your Python installation.

`python -m pip install git+git://github.com/Nicholaswogan/PhotochemPy.git`

Alternatively, you can instead compile PhotochemPy in this directory by running `./compile.sh`. This is useful for if you are adding a new feature to the source code (`src` directory) and need to compile a bunch of times to test it.

## Fortran
If you prefer to use the code exclusively in Fortran, that is OK too. Examples are provided in the folder [examples/fortran_example](https://github.com/Nicholaswogan/PhotochemPy/tree/main/examples/fortran_example).

## Note on parallel processing
PhotochemPy is a parallel program. To set the number of threads used by PhotochemPy set the `OMP_NUM_THREADS` environment variable: `export OMP_NUM_THREADS=6`. This command will force PhotochemPy to use 6 threads unless the computer has less than 6 threads, then it will use the maximum available threads. Generally speaking, using more threads will make the program faster.
