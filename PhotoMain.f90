  Program PhotoMain
    use Photochem
    implicit none
    integer :: nnw, nnz, nnq, nnp, nnsp, nnr, kks, kkj
    integer i

    nnz = 200
    nnq = 61
    nnp = 4
    nnsp = 74
    nnr = 392
    kks = 33
    kkj = 60

    call allocate_memory(nnz,nnq,nnp,nnsp,nnr,kks,kkj)
    call setup('input/speciesOG.dat', &
               'input/reactions.rx', &
               'input/planet.dat', &
               'input/input_photchem.dat', &
               'input/atmosphere.txt', &
               'input/Sun_2.7Ga.txt')


    call integrate

  end Program
