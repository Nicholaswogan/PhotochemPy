
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



