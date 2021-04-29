  subroutine setup(species_dat,reactions_rx,planet_dat,&
                 & photochem_dat, atmosphere_txt, flux_txt, err)

    use photochem_data, only: nz, nr, nz1 , nq, nq1, nw, isl, &
                              kw, frak, ihztype, jtrop, &
                              lh, lh2, planet, wavl, dz, z, ztrop
    use photochem_vars, only: usol_init, veff, mbound, den, T
    use photochem_wrk, only: scale_h, bhn2, bh2n2, hscale, du, dl, dk, dd, &
                             adu, adl, add, h_atm, A, rain, raingc
    implicit none

    ! local variables
    character(len=*),intent(in) :: species_dat
    character(len=*),intent(in) :: reactions_rx
    character(len=*),intent(in) :: planet_dat
    character(len=*),intent(in) :: photochem_dat
    character(len=*),intent(in) :: atmosphere_txt
    character(len=*),intent(in) :: flux_txt
    character(len=err_len), intent(out) :: err
    integer i,j
    real*8, dimension(:,:), allocatable :: Fval
    real*8, dimension(:), allocatable :: H2O
    integer :: nnq, nnsp, nnp, nnr, kks, kkj, nnw, nnz

    ! this subroutine will load all the data into memory
    ! (e.g. cross sections, rate date, etc.)
    err = ''
    call determine_dimensions(species_dat,reactions_rx,planet_dat, &
                              photochem_dat, atmosphere_txt, flux_txt, &
                              nnq, nnsp, nnp, nnr, kks, kkj, nnw, nnz, err)
    if (len_trim(err) /= 0) return
    call allocate_memory(nnz,nnq,nnp,nnsp,nnr,kks,kkj)
    allocate(Fval(nq,nz))
    allocate(H2O(nz))
    call read_species(species_dat,err)
    if (len_trim(err) /= 0) return
    call read_reactions(reactions_rx, err)
    if (len_trim(err) /= 0) return
    call read_planet(planet_dat,err)
    if (len_trim(err) /= 0) return
    call read_photochem(photochem_dat,err) 
    if (len_trim(err) /= 0) return
    call read_atmosphere(atmosphere_txt, err)
    if (len_trim(err) /= 0) return
    ! all above DO NOT depend on anything like T, den etc.
    ! only need to be read in once at the beginning.
    
    ! all below depend on something (T, density, eddy, etc.)
    ! so all of it should always be run before a calculation
    
    ! next step is to improve arguments in this functions, and to identify functions which depend on usol
    call photgrid(100.0D5, nz, ztrop, z, dz, jtrop)
    call densty ! DEPENDS ON USOL
    call rates(nz, nr, T, den, A) ! DEPENDS ON USOL
    call difco ! DEPENDS ON USOL
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
    call initphoto(flux_txt,err)
    if (len_trim(err) /= 0) return
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

    call rainout(.true.,Jtrop,Usol_init,nq,nz, rain, raingc) ! help rainout not break on first step
    



  end subroutine
