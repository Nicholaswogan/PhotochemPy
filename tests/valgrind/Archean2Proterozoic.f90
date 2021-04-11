Program PhotoMain
  use Photochem
  implicit none
  integer :: nnw, nnz, nnq, nnp, nnsp, nnr, kks, kkj
  integer i, converged

  double precision, dimension(1) :: teval
  double precision, dimension(1,71,200) :: solution
  logical :: use_fast_jacobian
  logical :: success


  ! Point the Photochem module to the directory where the DATA folder is located
  rootdir = '../../PhotochemPy/'

  ! Establish the dimensions of your problem
  ! and then allocate memory to variables.
  nnz = 200 ! number of vertical layers in the atmosphere
  nnq = 71 ! number of long-lived species (does not include inert and short-lived species)
  nnp = 4 ! number of particles
  nnsp = 86 ! total number of species
  nnr = 477 ! number of reactions (as given in reactions file)
  kks = 37 ! number of photolysis species
  kkj = 64 ! number of photolysis reactions
  call allocate_memory(nnz,nnq,nnp,nnsp,nnr,kks,kkj) ! allocates memory

  ! Now load all input files. here I use the Archean+haze template
  call setup('../../input/templates/Archean2Proterozoic/species.dat', &
             '../../input/templates/Archean2Proterozoic/reactions.rx', &
             '../../input/templates/Archean2Proterozoic/PLANET.dat', &
             '../../input/templates/Archean2Proterozoic/input_photchem.dat', &
             '../../input/templates/Archean2Proterozoic/atmosphere.txt', &
             '../../input/templates/Archean2Proterozoic/Sun_2.7Ga.txt')

  teval(1) = 1.d17
  use_fast_jacobian = .true. ! use fast jacobian

  ! goe.txt is the output file
  call cvode(0.0d0,usol_init,nnq,nnz,teval,1,1.d-4,1.d-22, use_fast_jacobian ,solution,success )

end Program
