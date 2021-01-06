# Installation

**Requirements:**
To install PhotochemPy, you must have the following installed on your system.
- `Python` with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

## Python
After satisfying the above requirements, the following command will install the PhotochemPy package to your Python installation with **default** options.

`Python -m pip install git+git://github.com/Nicholaswogan/PhotochemPy.git`

**Other install options:**

There are 3 total install options. Each option has more or less parallel processing.
- Option 1 (default): Parallel jacobian calculation. (pretty fast)
- Option 2: Parallel jacobian calculation and parallel radiative transfer calculations. (fastest)
- Option 3: No parallel processing. (slowest)

To install one of the non-default options:
1. Download the PhotochemPy respository
2. Open the setup.py file located in the main directory and edit the [option](https://github.com/Nicholaswogan/PhotochemPy/blob/main/setup.py#L19) variable to the option you desire.
3. Open a terminal in the main PhotochemPy directory and install with `python -m pip install .`

Note! *Option 2 will not have perfect reproducibility from run-to-run*. Calculations will differ by a factor close to machine precision (i.e. something like a factor of ~1.000000000000001). This is not a bug in the code. It is a fundamental aspect to some kinds of parallel processing.

Alternatively, if you don't want to add PhotochemPy to you Python installation you can instead compile it in this directory by running `./compile.sh`. This is useful for if you are adding a new feature to the source code (`src directory`) and need to compile a bunch of times to test it.

## Fortran
If you prefer to use the code exclusively in Fortran, that is OK too. Examples are provided in the folder [examples/fortran_example](https://github.com/Nicholaswogan/PhotochemPy/tree/main/examples/fortran_example).

## A note on parallel processing
PhotochemPy is a parallel program (for install options 1 and 2). To set the number of threads used by PhotochemPy export the `OMP_NUM_THREADS` environment variable: `export OMP_NUM_THREADS=6`. This command will force PhotochemPy to use 6 threads unless the computer has less than 6 threads, then it will use the maximum available threads. Generally speaking, using more threads will make the program faster.
