  subroutine setup(species_dat,reactions_rx,planet_dat,&
                 & photochem_dat, atmosphere_txt, flux_txt)

    use photochem_data, only: nz, nr, nz1 , nq, nq1, nw, isl, kw, frak, ihztype, jtrop, &
                              lh, lh2, planet, wavl, dz
    use photochem_vars, only: usol_init, veff, mbound, den, T
    use photochem_wrk, only: scale_h, bhn2, bh2n2, hscale, du, dl, dk, dd, &
                             adu, adl, add, h_atm, A
    implicit none

    ! local variables
    character(len=*),intent(in) :: species_dat
    character(len=*),intent(in) :: reactions_rx
    character(len=*),intent(in) :: planet_dat
    character(len=*),intent(in) :: photochem_dat
    character(len=*),intent(in) :: atmosphere_txt
    character(len=*),intent(in) :: flux_txt
    integer i,j
    real*8, dimension(nq,nz) :: Fval
    real*8, dimension(nz) :: H2O


    ! this subroutine will load all the data into memory
    ! (e.g. cross sections, rate date, etc.)

    call read_species(species_dat)
    call read_reactions(reactions_rx)
    call read_planet(planet_dat)
    call read_photochem(photochem_dat)
    call read_atmosphere(atmosphere_txt)
    call photgrid(100.0D5)
    call densty
    call rates(nz, nr, T, den, A)
    call difco
    call photsatrat(jtrop,nz,h2o)
    call dochem(Fval,-1,jtrop,isl,usol_init,nq,nz)

    if (mbound(LH2) .gt. 0) then
      do i=1,nz
!        !don't use molecular diffusion
        bHN2(i) = 0.0d0
        bH2N2(i) = 0.0d0

      enddo
    else
!      !use effusion velocity formulation of diffusion limited flux
      Veff(LH) = 1.0*bhN2(nz)/DEN(NZ)*(1./Hscale(nz) &
       - 1./scale_H(LH,nz))
!      !diff lim flux
      Veff(LH2) = 1.0*bH2N2(nz)/DEN(NZ)*(1./Hscale(nz) &
      - 1./scale_H(LH2,nz))

    endif

    if (planet .eq. 'EARTH') call ltning(usol_init,nq,nz)
    call aertab
    call initphoto(flux_txt)
    call initmie(nw,wavl,kw,frak,ihztype)

    ! Compute the jacobian coefficents
    do i=1,nq1
      DU(i,1) = DK(1)/DEN(1)/DZ(1)**2
      DL(i,NZ) = DK(NZ1)/DEN(NZ)/DZ(NZ)**2
      DD(i,1) = DU(i,1)
      DD(i,NZ) = DL(i,NZ)
      do J=2,NZ1
        DU(i,J) = DK(J)/DEN(J)/DZ(J)**2
        DL(i,J) = DK(J-1)/DEN(J)/DZ(J)**2
        DD(i,J) = DU(i,J) + DL(i,J)
      enddo
    enddo
    do i=1,nq
      do j=1,nz
        ADU(i,j) = 0.0D0
        ADL(i,j) = 0.0D0
        ADD(i,j) = 0.0D0
      enddo
    enddo

    if (mbound(LH2).eq.0) then
      ! diff limited flux implemented as effusion velocity
      DU(LH,1) = DU(LH,1) + bHN2(1)/Den(1)/DZ(1)**2
      ADU(LH,1) = bHN2(1)/Den(1)/DZ(1)/2.* &
          (1./scale_H(LH,1)-1./H_atm(1))
      DU(LH2,1) = DU(LH2,1) + bH2N2(1)/Den(1)/DZ(1)**2
      ADU(LH2,1) = bH2N2(1)/Den(1)/DZ(1)/2.* &
          (1./scale_H(LH2,1)-1./H_atm(1))
      ! upper boundary condition
      DL(LH,NZ) = DL(LH,NZ) + bHN2(nz1)/Den(nz)/DZ(NZ)**2
      ADL(LH,NZ) = -bHN2(nz1)/Den(nz)/DZ(nz)/2.* &
          (1./scale_H(LH,nz1)-1./H_atm(nz1))
      DL(LH2,NZ) = DL(LH2,NZ) + bH2N2(nz1)/Den(nz)/DZ(NZ)**2
      ADL(LH2,NZ) = -bH2N2(nz1)/Den(nz)/DZ(NZ)/2.* &
          (1./scale_H(LH2,nz1)-1./H_atm(nz1))
      !  unused....
      DD(LH,1) = DU(LH,1)
      ADD(LH,1) = -ADU(LH,1)
      DD(LH2,1) = DU(LH2,1)
      ADD(LH2,1) = -ADU(LH2,1)

      ! interior grid points   fixed 8-13-05
      do j=2,nz1
        DU(LH,j) = DU(LH,j) + bHN2(j)/Den(j)/DZ(j)**2
        ADU(LH,j) = bHN2(j)/Den(j)/DZ(j)/2.* &
                (1./scale_H(LH,j)-1./H_atm(j))
        DL(LH,j) = DL(LH,j) + bHN2(j-1)/Den(j)/DZ(j)**2
        ADL(LH,j) = -bHN2(j-1)/Den(j)/DZ(j)/2.* &
                (1./scale_H(LH,j-1)-1./H_atm(j-1))
        DU(LH2,j) = DU(LH2,j) + bH2N2(j)/Den(j)/DZ(j)**2
        ADU(LH2,j) = bH2N2(j)/Den(j)/DZ(j)/2.* &
                 (1./scale_H(LH2,j)-1./H_atm(j))
        DL(LH2,j) = DL(LH2,j) + bH2N2(j-1)/Den(j)/DZ(j)**2
        ADL(LH2,j) = -bH2N2(j-1)/Den(j)/DZ(j)/2.* &
                 (1./scale_H(LH2,j-1)-1./H_atm(j-1))
        DD(LH,j) = DU(LH,j) + DL(LH,j)
        ADD(LH,j) = -ADU(LH,j) - ADL(LH,j)
        DD(LH2,j) = DU(LH2,j) + DL(LH2,j)
        ADD(LH2,j) = -ADU(LH2,j) - ADL(LH2,j)
       enddo
    endif  !end molecular diffusion for H and H2 loop

    call rainout(Jtrop,0,Usol_init,nq,nz) ! help rainout not break on first step




  end subroutine
