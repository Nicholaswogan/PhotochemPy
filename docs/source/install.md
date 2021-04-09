# Installation

## Requirements
To install PhotochemPy, you must have the following installed on your system.
- `Python` (>3.6.0) with the `numpy` package. I suggest using [anaconda](https://www.anaconda.com/) to install these regardless of your operating system.
- The GNU compiler collection, version >4.9.4 (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

## Python
After satisfying the requirements, then follow these setups to install PhotochemPy

- Clone or download the github repository [`https://github.com/Nicholaswogan/PhotochemPy`](https://github.com/Nicholaswogan/PhotochemPy)
- In a terminal, navigate to the folder `src/cvode-5.7.0/build_dir`, and run the shell script `compile.sh`. This compiles CVODE.
- Navigate to the root directory of PhotochemPy, then install with `python -m pip install .`

### For development
You can also compile PhotochemPy in PhotochemPy root directory by running `./compile.sh`. This is useful for if you are adding a new feature to the source code (`src` directory) and need to compile a bunch of times to test it. For this to work,  you must first compile CVODE by running the shell script `src/cvode-5.7.0/build_dir/compile.sh`.

<!-- ## Fortran
If you prefer to use the code exclusively in Fortran, that is OK too. Examples are provided in the folder [examples/fortran_example](https://github.com/Nicholaswogan/PhotochemPy/tree/main/examples/fortran_example). -->

## Note on parallel processing
PhotochemPy is a parallel program. To set the number of threads used by PhotochemPy set the `OMP_NUM_THREADS` environment variable: `export OMP_NUM_THREADS=6`. This command will force PhotochemPy to use 6 threads unless the computer has less than 6 threads, then it will use the maximum available threads.
