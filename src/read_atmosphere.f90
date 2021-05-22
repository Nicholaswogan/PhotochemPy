
subroutine read_atmosphere_file(atmosphere_txt, err)
  use photochem_data, only: ispec, np, nq, &
                            nzf, z_file, &
                            T_file, edd_file, usol_file, rpar_file, &
                            wfall_file, aersol_file
  implicit none

  ! local variables
  character(len=*), intent(in) :: atmosphere_txt
  character(len=err_len), intent(out) :: err
  
  
  character(len=10000) :: line
  character(len=8), dimension(1000) :: arr1
  character(len=24), dimension(1000) :: arr11
  character(len=24),allocatable, dimension(:) :: labels
  integer :: ind(1)
  real*8, allocatable :: temp(:)
  integer :: i, n, nn, io, j, k, ii, iii
  
  err = ''
  
  open(4, file=trim(atmosphere_txt),status='old',iostat=io)
  if (io /= 0) then
    err = 'Can not open file '//trim(atmosphere_txt)
    return
  endif
  read(4,'(A)') line
  
  nzf = -1
  io = 0
  do while (io == 0)
    read(4,*,iostat=io)
    nzf = nzf + 1
  enddo

  if (allocated(z_file)) then
    deallocate(z_file, T_file)
    deallocate(edd_file, usol_file)
    deallocate(rpar_file, wfall_file, aersol_file)
  endif
  allocate(z_file(nzf), T_file(nzf))
  allocate(edd_file(nzf), usol_file(nq,nzf))
  allocate(rpar_file(nzf,np), wfall_file(nzf,np), aersol_file(nzf,np))
  z_file = 0.d0
  T_file = 0.d0
  edd_file = 0.d0
  usol_file = 0.d0
  rpar_file = 0.d0
  wfall_file = 0.d0
  aersol_file = 0.d0
  
  rewind(4)
  read(4,'(A)') line
  n = 0
  nn = 0
  do i=1,1000
    read(line,*,iostat=io) arr1(1:i)
    if (io==-1) exit
    n = n+1
  enddo
  read(4,'(A)') line
  do i=1,1000
    read(line,*,iostat=io) arr11(1:i)
    if (io==-1) exit
    nn = nn+1
  enddo
  if (n /= nn) then
    err = 'There is a missing column label in the file '//trim(atmosphere_txt)
    return
  endif
  
  allocate(labels(n))
  allocate(temp(n))
  rewind(4)
  read(4,'(A)') line
  read(line,*) (labels(i),i=1,n)

  ! reads in mixing ratios
  iii = 0
  do i=1,nq
    do j=1,n
      if (labels(j).eq.ispec(i)) then
        iii = iii+1
        do k = 1,nzf
          read(4,*,iostat=io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in initial atmosphere in '//trim(atmosphere_txt)
            return
          endif
          usol_file(i,k) = temp(j)
        enddo
        rewind(4) ! rewind!
        read(4,*) ! skip first line
        exit
      endif
    enddo
  enddo

  if (iii.ne.nq) then
    print*,'Warning: Did not find initial data for some species in '// &
            trim(atmosphere_txt)//' . The program will assume initial mixing ratios of zero.'
  endif

  rewind(4)
  read(4,*)
  ! reads in temperature
  ind = findloc(labels,'temp')
  if (ind(1) /= 0) then
    do k=1,nzf
      read(4,*,iostat = io) (temp(ii),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in temperature in '//trim(atmosphere_txt)
        return
      endif
      T_file(k) = temp(ind(1))
    enddo
  else
    err = 'temp was not found in input file '//trim(atmosphere_txt)
    return
  endif
  
  rewind(4)
  read(4,*)
  ! reads in alt
  ind = findloc(labels,'alt')
  if (ind(1) /= 0) then
    do k=1,nzf
      read(4,*,iostat = io) (temp(ii),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in altitude in '//trim(atmosphere_txt)
        return
      endif
      z_file(k) = temp(ind(1))*1.d5 ! conver to cm
    enddo
  else
    err = '"alt" was not found in input file '//trim(atmosphere_txt)
    return
  endif

  rewind(4)
  read(4,*)
  ! reads in eddy diffusion?
  ind = findloc(labels,'eddy')
  if (ind(1) /= 0) then
    do k=1,nzf
      read(4,*,iostat = io) (temp(ii),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in eddy diffusion in '//trim(atmosphere_txt)
        return
      endif
      edd_file(k) = temp(ind(1))
    enddo
  else
    err = 'eddy was not found in input file '//trim(atmosphere_txt)
    return
  endif

  rewind(4)
  read(4,*)
  iii = 0
  ! reads in aersol parameters
  if (np.gt.0) then
  do j=1,n
    do i=1,np
      if (trim(labels(j)).eq.trim(ISPEC(nq-np+i))//'_AERSOL') then
        iii = iii + 1
        do k=1,nzf
          read(4,*,iostat = io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in aersol parameters in '//trim(atmosphere_txt)
            return
          endif
          aersol_file(k,i) = temp(j)
        enddo
        rewind(4)
        read(4,*)
        exit
      else if (trim(labels(j)).eq.trim(ISPEC(nq-np+i))//'_WFALL') then
        iii = iii + 1
        do k=1,nzf
          read(4,*,iostat = io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in aersol parameters in '//trim(atmosphere_txt)
            return
          endif
          wfall_file(k,i) = temp(j)
        enddo
        rewind(4)
        read(4,*)
        exit
      else if (trim(labels(j)).eq.trim(ISPEC(nq-np+i))//'_RPAR') then
        iii = iii + 1
        do k=1,nzf
          read(4,*,iostat = io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in aersol parameters in '//trim(atmosphere_txt)
            return
          endif
          rpar_file(k,i) = temp(j)
        enddo
        rewind(4)
        read(4,*)
        exit
      endif
    enddo
  enddo
  endif
  close(4)

  if(iii /= np*3) then
    err = 'Missing aersol parameters in '//trim(atmosphere_txt)
    return
  endif
  
  

end subroutine

subroutine interp2atmosfile(nz, nq, np, z, T, edd, usol, rpar, wfall, aersol, err)
  use photochem_data, only: nzf, z_file, &
                            T_file, edd_file, usol_file, rpar_file, &
                            wfall_file, aersol_file
  implicit none
  integer, intent(in) :: nz, nq, np
  real(8), intent(in) :: z(nz)
  real(8), intent(out) :: T(nz), edd(nz), usol(nq,nz)
  real(8), intent(out) :: rpar(nz,np), wfall(nz,np), aersol(nz,np)
  character(len=err_len), intent(out) :: err
  
  integer :: i
  
  err = ''
  
  call interp(nz, nzf, z, z_file, T_file, T, err)
  if (len_trim(err) /= 0) return
  
  call interp(nz, nzf, z, z_file, dlog10(dabs(edd_file)), edd, err)
  if (len_trim(err) /= 0) return
  edd = 10.d0**edd
  
  do i = 1,nq
    call interp(nz, nzf, z, z_file, dlog10(dabs(usol_file(i,:))), usol(i,:), err)
    if (len_trim(err) /= 0) return
  enddo
  usol = 10.d0**usol
  
  do i = 1,np
    call interp(nz, nzf, z, z_file, dlog10(dabs(rpar_file(:,i))), rpar(:,i), err)
    if (len_trim(err) /= 0) return
    call interp(nz, nzf, z, z_file, dlog10(dabs(wfall_file(:,i))), wfall(:,i), err)
    if (len_trim(err) /= 0) return
    call interp(nz, nzf, z, z_file, dlog10(dabs(aersol_file(:,i))), aersol(:,i), err)
    if (len_trim(err) /= 0) return
  enddo
  rpar = 10.d0**rpar
  wfall = 10.d0**wfall
  aersol = 10.d0**aersol
  
  
  if (z(1) < z_file(1)) then
    print*,'Warning: vertical grid is being extrapolated below where there is input data.'
  endif
  
  if (z(nz) > z_file(nzf)) then
    print*,'Warning: vertical grid is being extrapolated above where there is input data.'
  endif
  
  
end subroutine


subroutine steam2photochem(nq, np, nz_in, P_surf, P_top, usol_layer, rpar_layer, err)
  use photochem_clima, only: pahlevan_H2_clima
  use photochem_data, only: g, background_mu, background_spec, mass, nz
  use photochem_vars, only: P, edd
  
  implicit none
  integer, intent(in) :: nq, np, nz_in
  real(8), intent(in) :: P_surf, P_top
  real(8), intent(in) :: usol_layer(nq), rpar_layer(np)
  character(len=1000), intent(out) :: err
  
  real(8) :: z_in(nz_in), P_in(nz_in), T_in(nz_in), ztrop_in
  real(8) :: edd_in(nz_in), mubar
  
  real(8) :: P_in_rev(nz_in), edd_in_rev(nz_in)
  real(8) :: P_rev(nz), edd_rev(nz)
  
  integer :: i, j
  err = ""
  
  ! background must be H2
  if (background_spec /= "H2") then
    err = "steam2photochem can not be called unless the background atmosphere is H2"
    return
  endif
  if (P_top > P_surf) then
    err = "P_top must be less than P_surf"
    return
  endif
  if (sum(usol_layer) > 0.5d0) then
    err = "usol_layer sums to >50% of the atmosphere. "// &
          "This means assuming background H2 atmosphere is not OK."
    return
  endif

  ! compute mubar
  call mean_molecular_weight(nq, usol_layer, mass, background_mu, mubar)
  ! Use simple climate model to compute T(z)
  ! z_in, T_in, ztrop_in
  call pahlevan_H2_clima(P_surf, mubar, g, P_top, nz_in, &
                         z_in, P_in, T_in, ztrop_in, err)
  if (len_trim(err) /= 0) return
  
  ! compute eddy diffusion by interpolating via log10(P). Note, P_in is bars
  ! reverse P and eddy then interpolate to them
  do i = 1,nz
    j = nz-i+1
    P_rev(i) = P(j)
    edd_rev(i) = edd(j)
  enddo
  do i = 1,nz_in
    j = nz_in-i+1
    P_in_rev(i) = P_in(j)
  enddo

  call interp(nz_in, nz, dlog10(P_in_rev), dlog10(P_rev), dlog10(edd_rev), edd_in_rev, err)
  if (len_trim(err) /= 0) return
  edd_in_rev = 10.d0**edd_in_rev
  do i = 1,nz_in
    j = nz_in-i+1
    edd_in(i) = edd_in_rev(j)
  enddo  
  ! Set up new grid (no outputs except err)
  call new_z_grid(nq, np, nz_in, ztrop_in, P_surf, z_in, &
                  T_in, edd_in, usol_layer, rpar_layer, err)
  if (len_trim(err) /= 0) return


end subroutine


subroutine new_z_grid(nq, np, nz_in, ztrop_in, P_surf_in, z_in, T_in, edd_in, usol_layer, rpar_layer, err)
  use photochem_data, only: nz, z, dz, jtrop, ztrop, P0, LH2O, mass, background_mu
  use photochem_vars, only: T, edd, usol_init, rpar_init, &
                            wfall_init, aersol_init
  implicit none
  
  integer, intent(in) :: nq, np, nz_in
  real(8), intent(in) :: ztrop_in, P_surf_in, z_in(nz_in), T_in(nz_in), edd_in(nz_in)
  real(8), intent(in) :: usol_layer(nq), rpar_layer(np)
  character(len=err_len), intent(out) :: err
  
  real(8) :: den_tmp(nz_in), P_tmp(nz_in), press_tmp(nz_in), H2O_tmp(nz_in), h2osat_tmp(nz_in)
  real(8), allocatable :: mubar_z(:)
  integer :: i
  
  P0 = P_surf_in
  nz = nz_in
  ! this stuff depends on z
  call allocate_memory_z(nz,err)
  if (len_trim(err) /= 0) return
  z = z_in
  dz = z_in(2)-z_in(1)
  T = T_in
  edd = edd_in
  do i = 1,nq
    usol_init(i,:) = usol_layer(i)
  enddo
  ztrop = ztrop_in
  jtrop=minloc(z,1, z .ge. ztrop)-1 ! depends on ztrop (sorta like T)
  ! treat water vapor special
  allocate(mubar_z(nz))
  do i = 1,nz
    call mean_molecular_weight(nq, usol_init(:,i), mass, background_mu, mubar_z(i))
  enddo
  call densty(nz, mubar_z, T, den_tmp, P_tmp, press_tmp) 
  call photsatrat(nz, T, P_tmp, den_tmp, Jtrop, H2Osat_tmp, H2O_tmp)
  do i=1,jtrop
    usol_init(LH2O,i) = H2O_tmp(i) 
  enddo
  do i = jtrop+1,nz
    usol_init(lH2O,i) = H2O_tmp(jtrop)
  enddo
  do i = 1,np
    rpar_init(:,i) = rpar_layer(i)
  enddo
  wfall_init = 0.d0
  aersol_init = 0.d0
  ! end stuff dthat depends on just z
  
  ! This stuff depends on T an z
  call initphoto(err) ! depends on T
  if (len_trim(err) /= 0) return
  call aertab
  ! end Stuff that depends on T
  
end subroutine






