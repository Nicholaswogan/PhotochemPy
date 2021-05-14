

subroutine determine_dimensions(species_dat,reactions_rx,planet_dat, &
                                photochem_dat, atmosphere_txt, flux_txt, &
                                nq, nsp, np, nr, ks, kj, nw, nz, nzf, err)
  use photochem_vars, only: rootdir
  implicit none
  integer, parameter :: str_length  = 1000
  ! input
  character(len=*),intent(in) :: species_dat
  character(len=*),intent(in) :: reactions_rx
  character(len=*),intent(in) :: planet_dat
  character(len=*),intent(in) :: photochem_dat
  character(len=*),intent(in) :: atmosphere_txt
  character(len=*),intent(in) :: flux_txt
  ! output
  integer, intent(out) :: nq, nsp, np, nr, ks, kj, nw, nz, nzf
  character(len=err_len), intent(out) :: err
  ! local
  integer :: io
  character(len=str_length) :: line, spname, sptype
  character(len=str_length), allocatable :: photosp(:) 
  integer :: i, j
  logical :: unique
  
  nz = 200 ! for now we will hard-code number of layers
  
  nq = 0
  nsp = 0
  np = 0
  ! determine nq, nsp, np
  open(100, file=trim(species_dat),status='old',iostat=io)
  if (io /= 0) then
    err = 'The input file '//trim(species_dat)//'  does not exist'
    return
  endif
  do while (io == 0)
    read(100,'(A)',iostat=io) line
    if (line(1:1) /= '*' .and. io == 0) then
      read(line,*) spname, sptype
      if (sptype == 'LL') then
        nq = nq + 1
        nsp = nsp + 1
        if (index(spname,'AER') /= 0) then
          np = np + 1
        endif
      else if (sptype == 'SL') then
        nsp = nsp + 1
      else if (sptype == 'IN') then
        nsp = nsp + 1
      else if (sptype == "M" .or. sptype == "HV") then
        ! nothing
      else
        err = "Species type "//trim(sptype)//' in '//trim(species_dat)//' is not a valid species type.'
        return
      endif
    endif
  enddo  
  close(100)
  
  ! count reactions
  nr = 0
  kj = 0
  open(100, file=trim(reactions_rx),status='old',iostat=io)
  if (io /= 0) then
    err = 'The input file '//trim(reactions_rx)//'  does not exist'
    return
  endif
  do while (io == 0)
    read(100,'(A)',iostat=io) line
    if (io == 0) then
      if (len_trim(line) == 0) then
        write(spname,'(i5)') nr+1
        err = 'Line '//trim(adjustl(spname))//' in '//trim(reactions_rx)//' is blank. This is not allowed.'
        return
      endif
      nr = nr + 1
      if (trim(line(49:55)) == 'PHOTO') then
        kj = kj + 1
      endif
    endif
  enddo
  ! find number of unique photolysis species
  ks = 0
  allocate(photosp(kj))
  kj = 0
  rewind(100)
  do i = 1, nr
    read(100,'(A)') line
    if (trim(line(49:55)) == 'PHOTO') then
      kj = kj + 1
      photosp(kj) = line(1:10)
      ! compare to all previous
      unique = .true.
      do j = 1, kj-1
        if (photosp(kj) == photosp(j)) then
          unique = .false.
          exit
        endif
      enddo
      if (unique) then
        ks = ks + 1
      endif
    endif
  enddo
  close(100)
  
  ! find number of wavelengths
  open(100, file=trim(photochem_dat),status='old',iostat=io)
  if (io /= 0) then
    err = 'The input file '//trim(photochem_dat)//'  does not exist'
    return
  endif
  do while (io == 0)
    read(100,'(A)',iostat=io) line
    if (line(1:1) == '*' .or. line(1:1) == 'C') then
      cycle
    else
      do i = 1,str_length
        if (line(1:5) == "LGRID") then
          read(line,'(a,i2)') sptype, j
        endif
      enddo
    endif
  enddo
  close(100)
  if (j == 0) then
    open(100,file=trim(rootdir)//'DATA/GRIDS/wogan.grid',status = 'old',iostat = io)
    if (io /= 0) then
      err = 'Was not able to find wogan.grid'
      return
    endif
  else
    err = 'LGRID = 0, in input_photochem.dat, is the only option.'
    return 
  endif
  nw = -2 ! skip the header
  do while(io == 0)
    read(100,'(A)',iostat=io) line
    if (io == 0) nw = nw + 1
  enddo
  close(100)
  
  ! check that the other files exist
  open(100, file=trim(planet_dat),status='old',iostat=io)
  if (io /= 0) then
    err = 'The input file '//trim(planet_dat)//'  does not exist'
    return
  endif
  close(100)
  open(100, file=trim(atmosphere_txt),status='old',iostat=io)
  if (io /= 0) then
    err = 'The input file '//trim(atmosphere_txt)//'  does not exist'
    return
  endif
  nzf = -2
  do while(io == 0)
    read(100,*,iostat=io)
    nzf = nzf + 1
  enddo

  close(100)
  open(100, file=trim(flux_txt),status='old',iostat=io)
  if (io /= 0) then
    err = 'The input file '//trim(flux_txt)//'  does not exist'
    return
  endif
  close(100)

end subroutine