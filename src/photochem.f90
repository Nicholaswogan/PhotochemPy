module photochem
  implicit none

  public
  ! public integrate, allocate_memory, setup, jacobian, right_hand_side, cvode_equilibrium
  integer, parameter, private :: err_len = 1000

contains
  ! Module subroutines go here.
  ! include "fortran_subroutine_name.f90"
  ! e.g. include "Photo.f90"
  ! etc...

  ! ALL THESE WORK!!!
  include "determine_dimensions.f90"
  include "read_species.f90" ! reads species.dat
  include "read_reactions.f90" ! reads reactions.rx
  include "read_atmosphere.f90" ! reads atmosphere.txt
  include "read_planet.f90" ! reads planet.dat
  include "read_photochem.f90" ! reads input_photochem.dat
  include "photgrid.f90" ! step up grid for photolysis calculations
  include "Initphoto.f90"
  include "Initmie.f90"
  ! above is data read-in stuff
  include "Rates.f90" ! calculates reaction rates
  include "Xsections.f90"
  include "Rainout.f90"
  include "Aqueous.f90"
  include "Ltning.f90" ! Needs work for time dependent model
  include "Aertab.f90"
  include "Densty.f90"
  include "Aercon.f90"
  include "PhotSatrat.f90"
  include "Difco.f90"
  include "Sedmnt.f90"
  include "Dochem.f90"
  include "Chempl.f90"
  include "Photo.f90" ! need to deal with precision problem
  include "Rayleigh.f90"
  include "Twostr.f90"
  include "setup.f90"
  include "integrate.f90"
  include "right_hand_side.f90"
  include "jacobian.f90"
  include "cvode.f90"

  include "redox_conservation.f90"

  subroutine allocate_memory(nnz, nnq, nnp, nnsp,&
     nnr, kks, kkj)
     use photochem_data
     use photochem_vars
     use photochem_wrk
    implicit none
    integer, intent(in) ::  nnz, nnq, nnp, nnsp, nnr, kks, kkj

    ! The dimensions.
    nz = nnz
    nz1 = nz-1
    nq  = nnq
    nq1 = nq
    np = nnp
    nsp = nnsp
    nsp2 = nnsp+2
    ks = kks
    kj = kkj
    nr = nnr
    LDA=3*NQ+1
    NEQ=NQ*NZ

    ! if allocated, then deallocate
    if (allocated(ISPEC)) then

      deallocate(ISPEC)
      deallocate(LBOUND)
      deallocate(vdep0)
      deallocate(vdep)
      deallocate(fixedmr)
      deallocate(distflux)
      deallocate(sgflux)
      deallocate(distheight)
      deallocate(MBOUND)
      deallocate(SMFLUX)
      deallocate(VEFF0)
      deallocate(VEFF)
      deallocate(atomsO)
      deallocate(atomsH)
      deallocate(atomsC)
      deallocate(atomsS)
      deallocate(atomsN)
      deallocate(atomsCl)
      deallocate(redoxstate)
      deallocate(mass)
      deallocate(chemj)
      deallocate(jchem)
      deallocate(reactype)
      deallocate(rateparams)
      deallocate(iloss)
      deallocate(iprod)
      deallocate(photoreac)
      deallocate(photospec)
      deallocate(photonums)
      deallocate(usol_init)
      deallocate(den)
      deallocate(T)
      deallocate(EDD)
      deallocate(aersol)
      deallocate(wfall)
      deallocate(rpar)
      deallocate(aersol_init)
      deallocate(wfall_init)
      deallocate(rpar_init)
      deallocate(numl)
      deallocate(nump)
      deallocate(Press)
      deallocate(P)
      deallocate(z)
      deallocate(dz)
      deallocate(sq)
      deallocate(VH2O)
      deallocate(VH2SO4)
      deallocate(FSULF)
      deallocate(H2SO4S)
      deallocate(S8S)
      deallocate(QEXTT)
      deallocate(W0T)
      deallocate(GFT)
      deallocate(SL)
      deallocate(A)
      deallocate(H)
      deallocate(RAINGC)
      deallocate(RAIN)
      deallocate(XSAVE)
      deallocate(tauedd)
      deallocate(hscale)
      deallocate(H_ATM)
      deallocate(DK)
      deallocate(SCALE_H)
      deallocate(BHN2)
      deallocate(BH2N2)
      deallocate(h2osat)
      deallocate(DD)
      deallocate(DL)
      deallocate(DU)
      deallocate(ADL)
      deallocate(ADU)
      deallocate(ADD)
      deallocate(usol_out)
      deallocate(flow)
      deallocate(fluxo)
      deallocate(yp)
      deallocate(yl)

    endif

    ! allocate memory

    ! Defined in species.dat
    allocate(ISPEC(nsp2)) ! issue with this one
    allocate(LBOUND(nq))
    allocate(vdep0(nq))
    allocate(vdep(nq))
    allocate(fixedmr(nq))
    allocate(distflux(nq))
    allocate(sgflux(nq))
    allocate(distheight(nq))
    allocate(MBOUND(nq))
    allocate(SMFLUX(nq))
    allocate(VEFF0(nq))
    allocate(VEFF(nq))
    allocate(atomsO(nsp2))
    allocate(atomsH(nsp2))
    allocate(atomsC(nsp2))
    allocate(atomsS(nsp2))
    allocate(atomsN(nsp2))
    allocate(atomsCl(nsp2))
    allocate(redoxstate(nsp2))
    allocate(mass(nq))
    ispec = ''
    lbound = 0
    vdep0 = 0.d0
    vdep = 0.d0
    fixedmr = 0.d0
    distflux = 0.d0
    sgflux = 0.d0
    distheight = 0.d0
    mbound = 0
    smflux = 0.d0
    veff0 = 0.d0
    veff = 0.d0
    atomso = 0
    atomsh = 0
    atomsc = 0
    atomss = 0
    atomsn = 0
    atomscl = 0
    redoxstate = 0.d0
    redox_factor = 0.d0
    mass = 0.d0

    ! definined in reactions.rx
    allocate(chemj(5,nr))
    allocate(jchem(5,nr))
    allocate(reactype(nr))
    allocate(rateparams(4,nr))
    allocate(iloss(2,nsp,nmax))
    allocate(iprod(nsp,nmax))
    allocate(photoreac(kj))
    allocate(photospec(ks))
    allocate(photonums(kj))
    jchem = 0
    rateparams = 0.d0
    iloss = 0
    iprod = 0
    photoreac = 0
    photospec = 0
    photonums = 0

    ! needed in atmosphere.txt
    allocate(usol_init(nq,nz))
    allocate(den(nz))
    allocate(T(nz))
    allocate(EDD(nz))
    allocate(aersol(nz,np))
    allocate(wfall(nz,np))
    allocate(rpar(nz,np))
    allocate(aersol_init(nz,np))
    allocate(wfall_init(nz,np))
    allocate(rpar_init(nz,np))
    allocate(numl(nsp))
    allocate(nump(nsp))
    usol_init = 0.d0
    den = 0.d0
    t = 0.d0
    edd = 0.d0
    aersol = 0.d0
    wfall = 0.d0
    rpar = 0.d0
    aersol_init = 0.d0
    wfall_init = 0.d0
    rpar_init = 0.d0
    numl = 0
    nump = 0


    ! needed in Densty.f90
    allocate(Press(nz))
    allocate(P(nz))
    press = 0.d0
    p = 0.d0

    ! needed in photogrid.f90
    allocate(z(nz))
    allocate(dz(nz))
    z = 0.d0
    dz = 0.d0

    ! needed in initphoto.f90.
    allocate(sq(kj,nz,kw))
    sq = 0.d0

    ! needed in Aertab.f90
    allocate(VH2O(Nf,nz))
    allocate(VH2SO4(Nf,nz))
    VH2O = 0.0d0
    VH2SO4 = 0.0d0

    ! needed in Aercon.f90
    allocate(FSULF(nz))
    allocate(H2SO4S(nz))
    allocate(S8S(nz))
    fsulf = 0.0d0
    h2so4s = 0.0d0
    s8s = 0.0d0

    ! needed in Photo.f90
    allocate(QEXTT(kw,nz,np))
    allocate(W0T(kw,nz,np))
    allocate(GFT(kw,nz,np))
    qextt = 0.0d0
    w0t = 0.0d0
    gft = 0.0d0

    ! needed in dochem.f90
    allocate(SL(NSP,NZ))
    sl = 0.0d0


    ! needed in rates.f90
    allocate(A(NR,NZ))
    ! zero out
    A = 0.0d0

    ! needed in rainout.f90
    allocate(H(NQ,NZ))
    allocate(RAINGC(NQ,NZ))
    allocate(RAIN(NZ))
    allocate(XSAVE(naq,nz))
    h = 0.0d0
    raingc = 0.0d0
    rain = 0.0d0
    xsave = 0.0d0

    ! needed in Difco.f90
    allocate(tauedd(nz))
    allocate(hscale(nz))
    allocate(H_ATM(nz))
    allocate(DK(nz))
    allocate(SCALE_H(nq,nz))
    allocate(BHN2(nz))
    allocate(BH2N2(nz))
    tauedd = 0.0d0
    hscale = 0.0d0
    h_atm = 0.0d0
    dk = 0.0d0
    scale_h = 0.0d0
    bhn2 = 0.0d0
    bh2n2 = 0.0d0

    ! needed in PhotSatrat.f90
    allocate(h2osat(nz))
    h2osat = 0.0d0

    ! needed in setup.f90
    allocate(DD(NQ1,NZ))
    allocate(DL(NQ1,NZ))
    allocate(DU(NQ1,NZ))
    allocate(ADL(NQ,NZ))
    allocate(ADU(NQ,NZ))
    allocate(ADD(NQ,NZ))
    dd = 0.0d0
    dl = 0.0d0
    du = 0.0d0
    adl = 0.0d0
    adu = 0.0d0
    add = 0.0d0

    ! integrate.f90
    allocate(usol_out(nq,nz))
    allocate(flow(nq))
    allocate(fluxo(nq,nz))
    allocate(yp(nq,nz))
    allocate(yl(nq,nz))
    usol_out = 0.0d0
    flow = 0.0d0
    fluxo = 0.0d0
    yp = 0.0d0
    yl = 0.0d0

  end subroutine allocate_memory

  ! for julia wrapper
  subroutine get_usol_init(usol)
    use photochem_data, only: nq, nz
    use photochem_vars, only: usol_init
    implicit none
    double precision, dimension(nq,nz),intent(out) ::usol
    usol = usol_init
  end subroutine

end module
