
subroutine read_atmosphere(atmosphere_txt, err)
  use photochem_data, only: np, nq, nsp2, nz, ispec
  use photochem_vars, only: T, edd, usol_init, &
                            rpar_init, wfall_init, aersol_init
  use photochem_wrk, only: rpar, wfall, aersol
                           
  implicit none

  ! local variables
  character(len=*), intent(in) :: atmosphere_txt
  character(len=err_len), intent(out) :: err
  character(len=10000) :: line
  character(len=8), dimension(1000) :: arr1
  character(len=24), dimension(1000) :: arr11
  character(len=24),allocatable, dimension(:) :: labels
  integer :: ind(1)
  real*8, allocatable :: temp(:), file_mix(:), file_z(:)
  integer :: i, n, nn, io, j, k, ii, iii
  integer :: nzf
  
  
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
    nzf = nzf + 1
    read(4,*,iostat=io)
  enddo
  ! allocate this guy
  allocate(file_mix(nzf), &
           file_z(nzf))

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
          usol_init(i,k) = temp(j)
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
    do k=1,nz
      read(4,*,iostat = io) (temp(ii),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in temperature in '//trim(atmosphere_txt)
        return
      endif
      T(k) = temp(ind(1))
    enddo
  else
    err = 'temp was not found in input file '//trim(atmosphere_txt)
    return
  endif

  ! do not read in den. It is established later
  ! rewind(4)
  ! read(4,*)
  ! reads in density
  ! ind = findloc(labels,'density')
  ! if (ind(1) /= 0) then
    ! do k=1,nz
      ! read(4,*,iostat=io) (temp(ii),ii=1,n)
      ! DEN(k) = temp(ind(1))
    ! enddo
  ! endif


  rewind(4)
  read(4,*)
  ! reads in eddy diffusion?
  ind = findloc(labels,'eddy')
  if (ind(1) /= 0) then
    do k=1,nz
      read(4,*,iostat = io) (temp(ii),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in eddy diffusion in '//trim(atmosphere_txt)
        return
      endif
      edd(k) = temp(ind(1))
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
      if (trim(labels(j)).eq.trim(ISPEC(size(ISPEC)-(nsp2-nq)-np+i))//'_AERSOL') then
        iii = iii + 1
        do k=1,nz
          read(4,*,iostat = io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in aersol parameters in '//trim(atmosphere_txt)
            return
          endif
          aersol(k,i) = temp(j)
          aersol_init(k,i) = temp(j)
        enddo
        rewind(4)
        read(4,*)
        exit
      else if (trim(labels(j)).eq.trim(ISPEC(size(ISPEC)-(nsp2-nq)-np+i))//'_WFALL') then
        iii = iii + 1
        do k=1,nz
          read(4,*,iostat = io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in aersol parameters in '//trim(atmosphere_txt)
            return
          endif
          wfall(k,i) = temp(j)
          wfall_init(k,i) = temp(j)
        enddo
        rewind(4)
        read(4,*)
        exit
      else if (trim(labels(j)).eq.trim(ISPEC(size(ISPEC)-(nsp2-nq)-np+i))//'_RPAR') then
        iii = iii + 1
        do k=1,nz
          read(4,*,iostat = io) (temp(ii),ii=1,n)
          if (io /= 0) then
            err = 'Problem reading in aersol parameters in '//trim(atmosphere_txt)
            return
          endif
          rpar(k,i) = temp(j)
          rpar_init(k,i) = temp(j)
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
