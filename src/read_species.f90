
subroutine read_species(species_dat, err)
  use photochem_data, only: iSL, ISPEC, mass, CO2_inert, N2_inert, redoxstate, &
                            atomsO, atomsH, atomsC, atomsS, atomsN, atomsCl, &
                            LSO2, LH2CO, lh2so4, lso4aer, lh2s, &
                            LCO, LH2O, LH2, LCH4, LO2, LH, &
                            Ls8aer, Lhcaer, Lhcaer2, ls2, ls3, ls4, &
                            lno, lo, LCO2, ls
  use photochem_vars, only: LBOUND, VDEP0, veff0, FIXEDMR, distflux, sgflux, &
                            distheight, mbound, smflux, vdep, veff, FCO2, FN2
  implicit none

  ! input
  character(len=*), intent(in) :: species_dat
  character(len=err_len), intent(out) :: err
  integer :: iLL, iSP, iIN
  integer :: io, io1
  character(len=8) :: species, spectype
  character(len=8) :: AX
  integer :: LX
  integer :: LA,LB,LD,LE,LF,LM
  integer :: LBC, LG
  real*8 :: XX,YY,ZZ,XXX,YYY,ZZZ
  real*8, dimension(6) :: atom_mass
  err = ''
  ! open species.dat
  open(4, file=trim(species_dat),status='OLD')

  iLL=0
  ! counter for short lived species
  iSL=0
  ! counter for number of lines in species.dat file
  iSP=0
  ! counter for inert species
  iIN=0

  FCO2=0.
  FN2=0.
  CO2_inert = 0
  N2_inert = 0
  io = 0
  do while (io == 0)
    read(4,*,iostat=io) SPECIES,SPECTYPE
    if ((scan(species,'*').LT.1) .and. (io == 0)) then  
      iSP=iSP+1
      ISPEC(iSP)=species
      ! This loads the "Lnumbers" for ease of use later in the code
      ! call LNUM(ISPEC(isP),iSP)
      if(trim(species).eq.'SO2') LSO2=iSP
      if(trim(species).eq.'H2CO') LH2CO=iSP
      if(trim(species).eq.'H2SO4') LH2SO4=iSP
      if(trim(species).eq.'SO4AER') LSO4AER=iSP
      if(trim(species).eq.'H2S') LH2S=iSP
      if(trim(species).eq.'CO') LCO=iSP
      if(trim(species).eq.'H2O') LH2O=iSP
      if(trim(species).eq.'H2') LH2=iSP
      if(trim(species).eq.'CH4') LCH4=iSP
      if(trim(species).eq.'O2') LO2=iSP
      if(trim(species).eq.'S8AER') Ls8aer=iSP
      if(trim(species).eq.'HCAER') Lhcaer=iSP
      if(trim(species).eq.'HCAER2') lhcaer2=iSP
      if(trim(species).eq.'S2') ls2=iSP
      if(trim(species).eq.'S3') ls3=iSP
      if(trim(species).eq.'S4') ls4=iSP
      if(trim(species).eq.'S') ls=iSP
      if(trim(species).eq.'NO') lno=iSP
      if(trim(species).eq.'O') lo=iSP
      if(trim(species).eq.'H') lh=iSP
      if(trim(species).eq.'CO2') lco2=iSP
      ! Return to previous line in species.dat file
      backspace 4

      ! read in atmoic number data, NEVER use LC,LH,LN,LO,LS as placeholders
      ! as they mean something else...
      read(4,*,iostat = io1) AX,AX,LA,LB,LD,LE,LF,LM
      if (io1 /= 0) then
        err = 'Problem reading number of atoms for '//trim(species)//' in '// &
              trim(species_dat)
        return
      endif

      if (SPECTYPE.EQ.'LL') then
        iLL=iLL+1
        ! Return to previous line in species.dat file
        backspace 4

        ! This section reads in the boundary conditions from species.dat.
        read(4,*,iostat=io1) AX,AX,LX,LX,LX,LX,LX,LX,LBC,XX,YY,ZZ,XXX,LG,YYY,ZZZ
        if (io1 /= 0) then
          err = 'Problem reading boundary conditions for '//trim(species)//' in '// &
                trim(species_dat)
          return
        endif
        if ((lbc < 0) .or. (lbc > 3)) then
          write(ax,'(i8)') lbc
          err = 'lbound = '//trim(adjustl(AX))//' for '//trim(species)//' in '// &
                trim(species_dat)//' is not a valid boundary condition'
          return
        endif
        if ((LG < 0) .or. (LG > 2)) then
          write(ax,'(i8)') lbc
          err = 'mbound = '//trim(adjustl(AX))//' for '//trim(species)//' in '// &
                trim(species_dat)//' is not a valid boundary condition'
          return
        endif
        if (YY < -1.d-100) then
          err = 'fixedmr for '//trim(species)//' in '// &
                trim(species_dat)//' can not be negative.'
          return
        endif

        LBOUND(iLL)=LBC
        VDEP0(iLL)=XX
        FIXEDMR(iLL)=YY
        if (LBOUND(iLL).eq.3) then
          ! distributed flux
          distflux(iLL)=ZZ
        else
          ! lower boundary flux
          SGFLUX(iLL)=ZZ
        endif
        distheight(iLL)=XXX
        MBOUND(iLL)=LG
        SMFLUX(iLL)=YYY
        VEFF0(iLL)=ZZZ
      elseif (SPECTYPE.EQ.'IN') then
        iIN=iIN+1
        backspace 4
        read(4,*,iostat=io1) AX,AX,LX,LX,LX,LX,LX,LX,XX
        if (io1 /= 0) then
          err = 'Problem reading boundary conditions for '//trim(species)//' in '// &
                trim(species_dat)
          return
        endif
        if (species.EQ.'CO2') then
          FCO2=XX
          CO2_inert = 1
        endif
        if (species.EQ.'N2') then
          FN2=XX
          N2_inert = 1
        endif
      endif
      if (SPECTYPE.EQ.'SL') iSL=iSL+1
      if (SPECTYPE.EQ.'HV') iIN=iIN+1
      if (SPECTYPE.EQ.'M') iIN=iIN+1

      atomsO(iLL+iSL+iIN)=LA
      atomsH(iLL+iSL+iIN)=LB
      atomsC(iLL+iSL+iIN)=LD
      atomsS(iLL+iSL+iIN)=LE
      atomsN(iLL+iSL+iIN)=LF
      atomsCL(iLL+iSL+iIN)=LM

      ! molar mass of each atom [O, H, C, S, N, Cl]
      atom_mass = (/15.999, 1.00784, 12.0107, 32.065, 14.0067, 35.453/)
      mass(iLL+iSL+iIN) = la*atom_mass(1) + lb*atom_mass(2) + &
                          ld*atom_mass(3) + le*atom_mass(4) + &
                          lf*atom_mass(5) + lm*atom_mass(6)

    endif
  enddo
  ! close the file
  close(4)
  
  redoxstate = atomsO*1.0 + atomsH*(-0.5) + atomsS*(-2.) + &
               atomsCL*(-1.0) + atomsC*(-2)

  VDEP = VDEP0
  VEFF = VEFF0

end subroutine
