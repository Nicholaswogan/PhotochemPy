############################################################################
#  Compilers
############################################################################

#####################################################################
########### Two possible options (choose one)
########### Option 1- (any) Fortran Compiler is used alone
########### Option 2- Intel Fortran + (any) C Compiler to remove the Fortran runtime library dependencies
#######################################################################

OPTION = 1
#OPTION = 2

################################################################
########## Select C compiler and flags (mandatory for option #2)
################################################################
#CC     = /usr/local/bin/gcc-10
#CFLAGS = -O3 -c -fopenmp -fPIC

##########################################################
########## Select F90 compiler and flags  (Mandatory)
##########################################################

#>>>>> Remark: three examples F1-F2-F3 (choose one or create your own)

#F1- Intel Fortran Compiler (Mandatory for Option #2 above)
#F90 = /opt/intel/bin/ifort
#F90FLAGS = -c -O3 -cpp  -qopenmp -fPIC

##F2- Gnu Fortran Compiler
F90 = gfortran
F90FLAGS = -c -O3 -fopenmp -ffree-line-length-none -ffixed-line-length-none -cpp

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !Warning: Spike has not been tested with this compiler!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##F3-Portland Group Fortran Compiler
#F90 = pgf90
#F90FLAGS = -O1 -mp
