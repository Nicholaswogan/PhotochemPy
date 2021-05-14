

cd ../cvode-5.7.0
mkdir build_dir
cd build_dir
cmake \
-DBUILD_EXAMPLES=OFF \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
-DEXAMPLES_ENABLE_F77=ON \
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




