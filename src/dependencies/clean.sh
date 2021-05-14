# cvode
cd ../cvode-5.7.0
if [ -d "build_dir" ] 
then
  rm -r build_dir
fi

# yaml 
cd ../fortran-yaml
if [ -d "build_dir" ] 
then
  rm -r build_dir
fi

# clean
cd ../dependencies
if [ -d "lib" ] 
then
  find . -type f ! -name '*.sh' -delete
  rm -r -- ./*/
fi
