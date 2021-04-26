

program main
  use photochem, only: allocate_memory, setup, integrate, cvode_equilibrium
  use photochem_vars, only: rootdir
  implicit none
  integer :: nnz, nnq, nnp, nnsp, nnr, kks, kkj
  integer converged
  logical success

  rootdir = '../PhotochemPy/'

  nnz = 200 ! number of vertical layers in the atmosphere
  nnq = 71 ! number of long-lived species (does not include inert and short-lived species)
  nnp = 4 ! number of particles
  nnsp = 86 ! total number of species
  nnr = 477 ! number of reactions (as given in reactions file)
  kks = 37 ! number of photolysis species
  kkj = 64 ! number of photolysis reactions
  call allocate_memory(nnz,nnq,nnp,nnsp,nnr,kks,kkj) ! allocates memory

  ! Now load all input files. here I use the Archean+haze template
  call setup('../input/templates/Hadean+HCN/species.dat', &
             '../input/templates/Hadean+HCN/reactions.rx', &
             '../input/templates/Hadean+HCN/planet.dat', &
             '../input/templates/Hadean+HCN/input_photchem.dat', &
             '../input/templates/Hadean+HCN/atmosphere.txt', &
             '../input/templates/Hadean+HCN/Sun_4.0Ga.txt')
  call integrate(100,converged)
  call cvode_equilibrium(1.d-3, 1.d-30, .true., success)
end program
