  Program PhotoMain
    use Photochem
    implicit none
    integer :: nnw, nnz, nnq, nnp, nnsp, nnr, kks, kkj

    nnw = 118
    nnz = 200
    nnq = 63
    nnp = 4
    nnsp = 74
    nnr = 392
    kks = 33
    kkj = 60

    call allocate_memory(nnw,nnz,nnq,nnp,nnsp,nnr,kks,kkj)
    call read_species("input/species.dat")
    call read_reactions("input/reactions.rx")
    call read_atmosphere("input/atmosphere.txt")
    call photogrid(100.D0)
    call rates
    ! etc....

  end Program
