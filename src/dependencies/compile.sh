
# first clean
./clean.sh

cd ../cvode-5.7.0
mkdir build_dir
cd build_dir
cmake \
-DEXAMPLES_ENABLE_CXX=OFF \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
-DCMAKE_BUILD_TYPE=Release \
-DEXAMPLES_ENABLE_F77=OFF \
-DEXAMPLES_ENABLE_F90=OFF \
-DBUILD_FORTRAN77_INTERFACE=ON \
-DBUILD_SHARED_LIBS=OFF \
-DCMAKE_Fortran_COMPILER=gfortran \
../
make install

cd ../../fortran-yaml
mkdir build_dir
cd build_dir
cmake \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
../
make install

cp -r modules ../../dependencies

cd ../../minpack
gfortran -c *.f -O3
ar rcs libminpack.a *.o
rm *.o
mv libminpack.a ../dependencies/lib




