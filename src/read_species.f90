
      subroutine read_species(species_dat)
        use reading_vars
        implicit none

        ! module variables
        ! integer :: iSL
        ! character(len=8), allocatable, dimension(:) :: ISPEC
        ! integer, allocatable, dimension(:) :: LBOUND
        ! real*8, allocatable, dimension(:) :: VDEP0
        ! real*8, allocatable, dimension(:) :: FIXEDMR
        ! real*8, allocatable, dimension(:) :: distflux
        ! real*8, allocatable, dimension(:) :: SGFLUX
        ! real*8, allocatable, dimension(:) :: distheight
        ! integer, allocatable, dimension(:) :: MBOUND
        ! real*8, allocatable, dimension(:) :: SMFLUX
        ! real*8, allocatable, dimension(:) :: VEFF0
        ! integer, allocatable, dimension(:) :: atomsO
        ! integer, allocatable, dimension(:) :: atomsH
        ! integer, allocatable, dimension(:) :: atomsC
        ! integer, allocatable, dimension(:) :: atomsS
        ! integer, allocatable, dimension(:) :: atomsN
        ! integer, allocatable, dimension(:) :: atomsCL
        ! integer, allocatable, dimension(:) :: mass
        ! lnums

        ! local variables
        character(len=*) :: species_dat
!f2py   intent(in) :: species_dat
        integer :: iLL, iSP, iIN
        integer :: i
        character(len=8) :: species, spectype
        character(len=8) :: AX
        integer :: LX
        integer :: LA,LB,LD,LE,LF,LM
        integer :: LBC, LG
        real*8 :: XX,YY,ZZ,XXX,YYY,ZZZ
        real*8, dimension(6) :: atom_mass

        ! open species.dat
        open(4, file=trim(species_dat),status='OLD')

        iLL=0
        ! counter for short lived species
        iSL=0
        ! counter for number of lines in species.dat file
        iSP=0
        ! counter for inert species
        iIN=0

        i = 0
        FCO2=0.
        FN2=0.
        CO2_inert = 0
        N2_inert = 0
        do while (I.LT.1000)
        ! Note: Below will crash if species.dat is longer than 1000 lines.
          read(4,*, end=96) SPECIES,SPECTYPE
          if (scan(species,'*').LT.1) then  ! else ignore comments in species.dat file (lines that start with *)
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
            read(4,*) AX,AX,LA,LB,LD,LE,LF,LM

            if (SPECTYPE.EQ.'LL') then
              iLL=iLL+1
              ! Return to previous line in species.dat file
              backspace 4

              ! This section reads in the boundary conditions from species.dat.
          read(4,*) AX,AX,LX,LX,LX,LX,LX,LX,LBC,XX,YY,ZZ,XXX,LG,YYY,ZZZ

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
            endif
            if (SPECTYPE.EQ.'IN') then
              iIN=iIN+1
              backspace 4
              read(4,*) AX,AX,LX,LX,LX,LX,LX,LX,XX
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
          I=I+1
        enddo
   96   CONTINUE
        ! close the file
        close(4)

        redoxstate = atomsO*1.0 + atomsH*(-0.5) + atomsS*(-2.) + &
                     atomsCL*(-1.0) + atomsC*(-2)

        do i=1,NQ
          VDEP(i) = VDEP0(i)
          VEFF(i) = VEFF0(i)
        enddo

      end subroutine
