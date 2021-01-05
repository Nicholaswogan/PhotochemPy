# Installation

**Requirements:**
To install PhotochemPy, you must have the following installed on your system.
- `Python` with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

**Python Module:** After satisfying the requirements, the following command will install the PhotochemPy package to your Python installation.

`Python -m pip install git+git://github.com/Nicholaswogan/PhotochemPy.git`

Alternatively, you can instead compile PhotochemPy in this directory by running `./compile.sh`. This is useful for if you are adding a new feature to the source code (`src` directory) and need to compile a bunch of times to test it.

**Fortran source:** If you prefer to use the code exclusively in Fortran, that is OK too. An example is provided in the folder `examples/fortran_example`. If you only use the Fortran source, you do not need `Python` with `numpy` installed on your system (you only need the GNU compiler collection).

## Parallel programing and maximizing performance
PhotochemPy is a parallel program. To set the number of threads used by PhotochemPy export the `OMP_NUM_THREADS` environment variable: `export OMP_NUM_THREADS=6`. This command will force PhotochemPy to use 6 threads (unless the computer has less than 6 threads, then it will use the maximum available threads).

Generally speaking, you will get the best performance when `OMP_NUM_THREADS` is set to the number of available cores. The number of cores is often less than the number of available threads!
