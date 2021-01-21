# Installation

## Requirements
To install PhotochemPy, you must have the following installed on your system.
- `Python` (>3.6.0) with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection, version >4.9.4 (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

## Python
### The easy way
After satisfying the requirements, the following command will install the standard PhotochemPy package to your Python installation.

`python -m pip install git+git://github.com/Nicholaswogan/PhotochemPy.git`

### The hard way (if you desperately need more speed)
Another option is to install PhotochemPy with [Spike](http://www.ecs.umass.edu/~polizzi/spike/index.htm). Spike is a method for solving banded linear systems in parallel, and it is faster than the LINPACK routine used in the standard PhotochemPy installation. This is not the default installation option because installing it is complicated **and might fail or get messed up!** Installing with Spike is probably only worth it if you are considering really large problems (e.g. >100 species, and >200 atmospheric layers), and are using a computer with more than about 4 cores.

To install PhotochemPy with Spike follow these instructions (These instructions will only work for MacOS)
1. Install gfortran 6 with the Homebrew command `brew install gcc@6`. Other version of gfortran might not work.
2. Symlink gfortran-6 with gfortran: `sudo ln -s /usr/local/bin/gfortran-6 /usr/local/bin/gfortran`
3. Download or clone PhotochemPy from Github: `git clone https://github.com/Nicholaswogan/PhotochemPy.git`
4. Navigate to the directory `src/spike-1.0/src` (relative to the root-directory of PhotochemPy) with a terminal. Then build spike with the command `make -j all`
5. Navigate up one directory to `src/spike-1.0` and create a new anaconda environment called "photochem" with the `photochem.yml` file: `conda env create -f photochem.yml`. After that, activate it with `conda activate photochem`. **Note!** Installing lots of additional packages to the `photochem` environment might cause PhotochemPy to not run.
6. Open the `setup.py` (in the PhotochemPy root directory) in a text editor, and change `option = 3` at [this line](https://github.com/Nicholaswogan/PhotochemPy/blob/main/setup.py#L19).
7. Navigate the terminal to the root directory of PhotochemPy, and install with `python -m pip install .` (make sure the `photochem` anaconda environment is active)
9. Run the test to make sure the installation worked: `python test_PhotochemPy.py`

**Two things to Note!** First, Spike is fastest when the number of threads used is the same as the number of cores. My computer has 6 cores (and 12 total threads). To guarantee 6 threads are used I add the line `export OMP_NUM_THREADS=6` to the file `~/.bash_profile`. Also, Spike will cause PhotochemPy to have thread-dependent reproducibility. In other words, the answer will be a very tiny amount different depending on the number of threads used. This isn't a bug. It is just an aspect of how the Spike algorithm works.

### For development
You can also compile PhotochemPy in PhotochemPy root directory by running `./compile.sh`. This is useful for if you are adding a new feature to the source code (`src` directory) and need to compile a bunch of times to test it.

## Fortran
If you prefer to use the code exclusively in Fortran, that is OK too. Examples are provided in the folder [examples/fortran_example](https://github.com/Nicholaswogan/PhotochemPy/tree/main/examples/fortran_example).

## Note on parallel processing
PhotochemPy is a parallel program. To set the number of threads used by PhotochemPy set the `OMP_NUM_THREADS` environment variable: `export OMP_NUM_THREADS=6`. This command will force PhotochemPy to use 6 threads unless the computer has less than 6 threads, then it will use the maximum available threads.
