  Program PhotoMain
    use Photochem
    implicit none
    integer :: nnw, nnz, nnq, nnp, nnsp, nnr, kks, kkj
    integer i, converged

    ! Point the Photochem module to the directory where the DATA folder is located
    rootdir = '../../PhotochemPy/'

    ! Establish the dimensions of your problem
    ! and then allocate memory to variables.
    nnz = 200 ! number of vertical layers in the atmosphere
    nnq = 69 ! number of long-lived species (does not include inert and short-lived species)
    nnp = 4 ! number of particles
    nnsp = 86 ! total number of species
    nnr = 477 ! number of reactions (as given in reactions file)
    kks = 37 ! number of photolysis species
    kkj = 64 ! number of photolysis reactions
    call allocate_memory(nnz,nnq,nnp,nnsp,nnr,kks,kkj) ! allocates memory

    ! Now load all input files. here I use the Archean+haze template
    call setup('../../input/templates/Hadean+HCN/species.dat', &
               '../../input/templates/Hadean+HCN/reactions.rx', &
               '../../input/templates/Hadean+HCN/planet.dat', &
               '../../input/templates/Hadean+HCN/input_photchem.dat', &
               '../../input/templates/Hadean+HCN/atmosphere.txt', &
               '../../input/templates/Hadean+HCN/Sun_4.0Ga.txt')

    ! integrate to photochemical equilbirium
    call integrate(1,converged)
    ! 1000 = number of steps the integrator takes
    ! converged = 1 or 0. If 1 then integegrator found photochemical equilbrium
    ! or 0 if it didn't find it.


    ! Now you could save some output.
    ! for example, the atmosphere at equilbrium is stored in the variable
    ! usol_out


  end Program
