  Program PhotoMain
    use Photochem
    implicit none
    integer :: nnw, nnz, nnq, nnp, nnsp, nnr, kks, kkj
    integer i,converged

    nnz = 200
    nnq = 61
    nnp = 4
    nnsp = 74
    nnr = 392
    kks = 33
    kkj = 60

    call allocate_memory(nnz,nnq,nnp,nnsp,nnr,kks,kkj)
    call setup('input/templates/Archean+haze/species.dat', &
               'input/templates/Archean+haze/reactions.rx', &
               'input/templates/Archean+haze/planet.dat', &
               'input/templates/Archean+haze/input_photchem.dat', &
               'input/templates/Archean+haze/atmosphere.txt', &
               'input/templates/Archean+haze/Sun_2.7Ga.txt')


    call integrate(1000,converged)


  end Program
