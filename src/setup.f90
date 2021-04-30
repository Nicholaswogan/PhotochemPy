  subroutine setup(species_dat,reactions_rx,planet_dat,&
                 & photochem_dat, atmosphere_txt, flux_txt, err)

    use photochem_data, only: nz, nr, nq, nw, &
                              kw, frak, ihztype, jtrop, &
                              lh, lh2, planet, wavl, wav, wavu, dz, z, ztrop, &
                              lgrid, flux, background_spec
    use photochem_vars, only: usol_init, veff, mbound, den, T, P, press, edd, H2Osat
    use photochem_wrk, only: scale_h, bhn2, bh2n2, tauedd, hscale, du, dl, dk, dd, &
                             adu, adl, add, h_atm, A, rain, raingc, &
                             zapNO, zapO2, proNOP, zapCO, zapH2, zapO
    implicit none

    ! local variables
    character(len=*),intent(in) :: species_dat
    character(len=*),intent(in) :: reactions_rx
    character(len=*),intent(in) :: planet_dat
    character(len=*),intent(in) :: photochem_dat
    character(len=*),intent(in) :: atmosphere_txt
    character(len=*),intent(in) :: flux_txt
    character(len=err_len), intent(out) :: err
    ! integer i,j
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
    call gridw(nw,wavl,wav,wavu,lgrid, err) ! makes grid (depends on nothing)
    if (len_trim(err) /= 0) return
    call readflux(flux_txt,nw,wavl,flux,err) ! reads flux (depnds on gridw)
    if (len_trim(err) /= 0) return
    call initmie(nw,wavl,kw,frak,ihztype)
    call photgrid(100.0D5, nz, z, dz) 
    ! end depends on nothing
    ! This stuff depends on T
    JTROP=minloc(Z,1, Z .ge. ztrop)-1 ! depends on ztrop (sorta like T)
    call initphoto(err) ! depends on T
    if (len_trim(err) /= 0) return
    call aertab
    ! end Stuff that depends on T
    
    ! begin stuff that needs to be inizialized
    call densty(nq, nz, usol_init, T, den, P, press) ! DEPENDS ON USOL
    call rainout(.true.,Jtrop,Usol_init,nq,nz, T,den, rain, raingc) ! initial conditions for rainout solve
    ! end stuff that needs to be inizialized
    
    ! begin stuff that should be in RHS (depends on usol)
    call densty(nq, nz, usol_init, T, den, P, press) ! DEPENDS ON USOL
    call rates(nz, nr, T, den, A) ! DEPENDS ON USOL (via den)
    call difco(nq,nz,usol_init, T, den, edd, &
              hscale, tauedd, DK, H_atm, bhn2, bh2n2, scale_H)
    call photsatrat(nz, T, P, den, Jtrop, H2Osat, H2O) ! depends on usol
    if (planet .eq. 'EARTH') call ltning(nq, nz, usol_init, &
                                  zapNO, zapO2, proNOP, zapCO, zapH2, zapO)
    call diffusion_coeffs(nq, nz, den, dz, DK, bhN2, bh2N2, scale_H, H_atm, &
                         DU, DL, DD, ADU, ADL, ADD)
                                  
    ! below is H escape
    if (background_spec /= 'H2') then ! then we can consider its upper and lower boundary
      if (mbound(LH2) == 0) then
  !       do i=1,nz
  ! !        !don't use molecular diffusion
  !         bHN2(i) = 0.0d0
  !         bH2N2(i) = 0.0d0
  !       enddo
  !     else
  !      !use effusion velocity formulation of diffusion limited flux
        Veff(LH) = 1.0*bhN2(nz)/DEN(NZ)*(1./Hscale(nz) &
         - 1./scale_H(LH,nz))
  !      !diff lim flux
        Veff(LH2) = 1.0*bH2N2(nz)/DEN(NZ)*(1./Hscale(nz) &
        - 1./scale_H(LH2,nz))

      endif
    endif    
  end subroutine
  
  subroutine diffusion_coeffs(nq, nz, den, dz, DK, bhN2, bh2N2, scale_H, H_atm, &
                              DU, DL, DD, ADU, ADL, ADD)
    use photochem_data, only: lh2, lh, nq1, nz1, background_spec
    implicit none
    
    integer, intent(in) :: nq, nz
    real(8), intent(in) :: den(nz), dz(nz), dk(nz), bhn2(nz), bh2n2(nz)
    real(8), intent(in) :: scale_H(nq,nz), H_atm(nz)
    
    real(8), intent(out) :: DU(nq,nz), DL(nq,nz), DD(nq,nz)
    real(8), intent(out) :: ADU(nq,nz), ADL(nq,nz), ADD(nq,nz)
    
    integer :: i, j
  
    ! eddy diffusion stuff
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

    if (background_spec /= 'H2') then ! there is H2 molecular diffusion
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
  
  
  
  end subroutine
