  subroutine setup(species_dat,reactions_rx,set_file,&
                 & atmosphere_txt, flux_txt, err)

    use photochem_data, only: np, nq, nz, nw, &
                              frak, ihztype, jtrop, &
                              wavl, wav, wavu, dz, z, ztrop, &
                              lgrid, flux, &
                              top_atmos, bottom_atmos, background_mu, mass
                              
    use photochem_vars, only: den, P, press, T, edd, usol_init, rpar_init, &
                              wfall_init, aersol_init
    use photochem_wrk, only: rain, raingc, rpar, wfall, aersol
    implicit none

    ! local variables
    character(len=*),intent(in) :: species_dat
    character(len=*),intent(in) :: reactions_rx
    character(len=*),intent(in) :: set_file
    character(len=*),intent(in) :: atmosphere_txt
    character(len=*),intent(in) :: flux_txt
    character(len=1000), intent(out) :: err
    integer :: nnq, nnsp, nnp, nnr, kks, kkj, nnw
    real(8), allocatable :: mubar_z(:)
    integer :: i
    
    ! this subroutine will load all the data into memory
    ! (e.g. cross sections, rate date, etc.)
    err = ''
    call determine_dimensions(species_dat,reactions_rx, set_file, &
                              atmosphere_txt, flux_txt, &
                              nnq, nnsp, nnp, nnr, kks, kkj, nnw, err)
    if (len_trim(err) /= 0) return
    call allocate_memory(nnq,nnp,nnsp,nnr,kks,kkj,nnw)
    call read_species(species_dat,err)
    if (len_trim(err) /= 0) return
    call read_reactions(reactions_rx, err)
    if (len_trim(err) /= 0) return
    call read_settings(set_file,err)
    if (len_trim(err) /= 0) return
    call read_atmosphere_file(atmosphere_txt, err)
    if (len_trim(err) /= 0) return
    call gridw(nw,wavl,wav,wavu,lgrid, err) ! makes grid (depends on nothing)
    if (len_trim(err) /= 0) return
    call readflux(flux_txt,nw,wavl,flux,err) ! reads flux (depnds on gridw)
    if (len_trim(err) /= 0) return
    call initmie(nw,wavl,frak,ihztype)
    ! end depends on nothing
    
    ! this stuff depends on z
    call allocate_memory_z(nz,err)
    if (len_trim(err) /= 0) return
    call photgrid(top_atmos, bottom_atmos, nz, z, dz) 
    call interp2atmosfile(nz, nq, np, z, T, edd, &
                          usol_init, rpar_init, wfall_init, &
                          aersol_init, err)
    if (len_trim(err) /= 0) return
    ! end stuff dthat depends on just z
    
    ! This stuff depends on T an z
    JTROP=minloc(Z,1, Z .ge. ztrop)-1 ! depends on ztrop (sorta like T)
    call initphoto(err) ! depends on T
    if (len_trim(err) /= 0) return
    call aertab
    ! end Stuff that depends on T
    
    ! begin stuff that needs to be inizialized
    rpar = rpar_init
    wfall = wfall_init
    aersol = aersol_init
    allocate(mubar_z(nz))
    do i = 1,nz
      call mean_molecular_weight(nq, usol_init(:,i), mass, background_mu, mubar_z(i))
    enddo
    call densty(nz, mubar_z, T, den, P, press)
    call rainout(.true.,Jtrop,usol_init,nq,nz, T,den, rain, raingc,err)
    if (len_trim(err) /= 0) return 
    ! end stuff that needs to be inizialized
  end subroutine

  
  subroutine diffusion_coeffs(nq, nz, den, dz, DK, bx1x2, scale_H, H_atm, &
                              DU, DL, DD, ADU, ADL, ADD)
    use photochem_data, only: nq1, nz1, np
    implicit none
    
    integer, intent(in) :: nq, nz
    real(8), intent(in) :: den(nz), dz(nz), dk(nz), bx1x2(nq,nz)
    real(8), intent(in) :: scale_H(nq,nz), H_atm(nz)
    
    real(8), intent(out) :: DU(nq,nz), DL(nq,nz), DD(nq,nz)
    real(8), intent(out) :: ADU(nq,nz), ADL(nq,nz), ADD(nq,nz)
    
    integer :: i, j, k
  
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
    
    do k=1,(NQ-NP)
      !       if(mbound(k).eq.0) then
      DU(k,1) = DU(k,1) + bX1X2(k,1)/Den(1)/DZ(1)**2
      ADU(k,1) = bX1X2(k,1)/Den(1)/DZ(1)/2.* &
         (1./scale_H(k,1)-1./H_atm(1))
      ! upper boundary condition
      DL(k,NZ) = DL(k,NZ) + bX1X2(k,nz1)/Den(nz)/DZ(NZ)**2
      ADL(k,NZ) = -bX1X2(k,nz1)/Den(nz)/DZ(nz)/2.* &
        (1./scale_H(k,nz1)-1./H_atm(nz1))
      !  unused...
      DD(k,1) = DU(k,1)
      ADD(k,1) = -ADU(k,1)

      ! interior grid points   ?fixed 8-13-05
      do j=2,nz1
         DU(k,j) = DU(k,j) + bX1X2(k,j)/Den(j)/DZ(j)**2
         ADU(k,j) = bX1X2(k,j)/Den(j)/DZ(j)/2.* &
               (1./scale_H(k,j)-1./H_atm(j))
         DL(k,j) = DL(k,j) + bX1X2(k,j-1)/Den(j)/DZ(j)**2
         ADL(k,j) = -bX1X2(k,j-1)/Den(j)/DZ(j)/2.* &
               (1./scale_H(k,j-1)-1./H_atm(j-1))
         DD(k,j) = DU(k,j) + DL(k,j)
         ADD(k,j) = -ADU(k,j) - ADL(k,j)
      enddo
    enddo

  end subroutine
