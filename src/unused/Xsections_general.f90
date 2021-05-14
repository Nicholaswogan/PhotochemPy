subroutine Xsections_general(species,reactions_rx,kj,nz,nw,kw,wavl,T,jn,sq)
  implicit none

  ! module variables

  ! inputs
  character(len=8), intent(in) :: species
  character(len=*), intent(in) :: reactions_rx
  integer, intent(in) :: kj, nz, nw, kw
  double precision, dimension(nw+1), intent(in) :: wavl
  double precision, dimension(nz), intent(in) :: T

  ! in/out
  integer, intent(inout) :: jn
  double precision, dimension(kj,nz,kw), intent(inout) :: sq

  ! other vars
  character(len=50) :: fmt
  character(len=10) :: reac1, reac2, prod1, prod2, prod3
  character(len=25) :: reacs, prods
  character(len=5) :: label
  character(len=500) :: xsname, qyname
  character(len=1000) :: temperatureline, line
  character(len=8), dimension(100) :: temperatures
  character(len=8) :: temp_dumb
  integer :: io, stat, numtempcols, kdata, i, ierr, iw
  double precision, allocatable, dimension(:,:) :: yy, yy_grid
  double precision, allocatable, dimension(:) :: temperature_nums
  double precision, allocatable, dimension(:) :: x1, y1, x2, qy
  double precision, allocatable, dimension(:) :: dumby
  double precision :: deltax = 1.d-4, biggest=1.d+36, zero=0.0d0
  double precision, dimension(nw) :: yg1, yq1
  double precision, dimension(nw,nz) :: yg2
  double precision, dimension(nz) :: p_interp

  character(len=500) :: rootdir
  rootdir = '../PhotochemPy/'

  ! initialize
  ierr = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! Read in cross sections !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open xs file
  xsname = trim(rootdir)//'DATA/XSECTIONS/'//trim(species)//'/'// &
            trim(species)//'.XS.dat'
  open(11, file=trim(xsname), status='OLD')

  ! determine number of temperature columns
  read(11,*)
  read(11,*)
  read(11,'(A)') temperatureline
  do i=1,100
    read(temperatureline,*,iostat=io) temperatures(1:i)
    if (io == -1) exit
  enddo
  if (i == 101) then
    print*,'More temperature data than allowed for ',species
  endif
  ! now i-2 is number of temperature columns
  numtempcols = max(i-2,1)

  ! determine how many xsections points there are
  read(11,*)
  read(11,*)
  io = 0
  kdata = 0
  do while(io == 0)
    kdata = kdata + 1
    READ(11, *, IOSTAT=io)
  enddo

  ! allocate memory
  allocate(yy(kdata+4,numtempcols))
  allocate(yy_grid(nw,numtempcols+2))
  allocate(temperature_nums(numtempcols+2))
  allocate(y1(kdata+4))
  allocate(x1(kdata+4))
  allocate(dumby(numtempcols+1))

  ! now read in data
  rewind(11)
  do i=1,4
    read(11,*)
  enddo
  do i=1,kdata
    read(11,'(A)') line
    read(line,*) dumby(1:numtempcols+1)
    x1(i) = dumby(1)
    yy(i,:) = dumby(2:)
  enddo
  close(11)

  ! re-bin data to the grid
  y1 = 0.d0
  y1 = yy(:,1)
  call addpnt(x1,y1,kdata+4,kdata,x1(1)*(1.d0-deltax),zero)
  call addpnt(x1,y1,kdata+4,kdata,               zero,zero)
  call addpnt(x1,y1,kdata+4,kdata,x1(kdata)*(1.d0+deltax),zero)
  call addpnt(x1,y1,kdata+4,kdata,            biggest,zero)
  call inter2(nw+1,wavl,yg1,kdata,x1,y1,ierr)
  kdata = kdata - 4
  yy_grid(:,1) = yg1
  do i=1,numtempcols
    x1(:kdata) = x1(3:kdata+2)
    x1(kdata+1:kdata+4) = 0.d0
    y1 = 0.d0
    y1 = yy(:,i)
    call addpnt(x1,y1,kdata+4,kdata,x1(1)*(1.d0-deltax),zero)
    call addpnt(x1,y1,kdata+4,kdata,               zero,zero)
    call addpnt(x1,y1,kdata+4,kdata,x1(kdata)*(1.d0+deltax),zero)
    call addpnt(x1,y1,kdata+4,kdata,            biggest,zero)
    call inter2(nw+1,wavl,yg1,kdata,x1,y1,ierr)
    kdata = kdata-4
    yy_grid(:,i+1) = yg1
  enddo
  x1(:kdata) = x1(3:kdata+2)
  x1(kdata+1:kdata+4) = 0.d0
  y1 = 0.d0
  y1 = yy(:,numtempcols)
  call addpnt(x1,y1,kdata+4,kdata,x1(1)*(1.d0-deltax),zero)
  call addpnt(x1,y1,kdata+4,kdata,               zero,zero)
  call addpnt(x1,y1,kdata+4,kdata,x1(kdata)*(1.d0+deltax),zero)
  call addpnt(x1,y1,kdata+4,kdata,            biggest,zero)
  call inter2(nw+1,wavl,yg1,kdata,x1,y1,ierr)
  kdata = kdata-4
  yy_grid(:,numtempcols+2) = yg1
  if (ierr .NE. 0) then
    print*, "failed while binning cross sections for ", species
  endif

  ! find temperatures of each column
  if (numtempcols == 1) then
    temperature_nums(2) = 300.d0 ! if 1 temp column then assume its 300K
  else
    do i=2,numtempcols+1
      temp_dumb = temperatures(i)
      read(temp_dumb(:index(temp_dumb,'K')-1),*) temperature_nums(i)
    enddo
  endif
  temperature_nums(1) = 0.d0
  temperature_nums(numtempcols+2) = 3000.d0
  do i=1,nz
    if ((T(i) > 3000.d0) .or. (T(i) < 0.d0)) then
      print*, 'Atmospheric temperature can not be above 3000 K or below 0 K'
      stop
    endif
  enddo

  ! now interpolate to atmospheric temperature for each bin
  ! and at each temperature
  do i = 1, nw
    call interp_linear(1, numtempcols+2, temperature_nums, yy_grid(i,:), nz, &
                       T, p_interp)
    yg2(i,:) = p_interp
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! Quantum Yields !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! should use this intead of reading reactions again
  ! j is the index of sq(j,i,iw)
  ! spec_ind = index(ispec,species)
  ! do j=1,kj
  !   if (photoreac(j) == spec_ind) then
  !     reac1 = chemj(1,photonums(j))
  !     reac2 = chemj(2,photonums(j))
  !     prod1 = chemj(3,photonums(j))
  !     prod2 = chemj(4,photonums(j))
  !     prod3 = chemj(5,photonums(j))
  !     print*,reac1,species
  !   endif
  ! enddo

  fmt = '(A10,A10,A10,A10,A8,A5)'
  open(12, file=trim(reactions_rx),status='OLD')
  io = 0
  do while(io == 0)
    read(12,trim(fmt),IOSTAT=io) reac1,reac2,prod1,prod2,prod3,label

    if ((label.eq.'PHOTO').and.(trim(species).eq.trim(reac1))) then

      reacs = trim(reac1) // " " // trim(reac2)
      prods = trim(prod1) // " " // trim(prod2) // " " &
             // trim(prod3)

      do i=1, len(trim(reacs))
        if(reacs(i:i) == ' ') reacs(i:i) = '_'
      end do

      do i=1, len(trim(prods))
        if(prods(i:i) == ' ') prods(i:i) = '_'
      enddo

      ! open QY file
      qyname = trim(rootdir)//'DATA/XSECTIONS/'//trim(species)//'/'// &
                trim(reacs)//"_"//trim(prods)//".QY.dat"

      ! open qy
      open(13, file=trim(qyname), status='OLD')
      ! count lines
      kdata = 0
      stat = 0
      do i=1,5
        read(13,*)
      enddo
      do while(stat == 0)
        read(13, *, IOSTAT=stat)
        kdata = kdata + 1
      end do

      allocate(qy(kdata+4))
      allocate(x2(kdata+4))

      ! rewind and read in data
      rewind(13)
      do i=1,4
        read(13,*)
      enddo
      do i=1,kdata
        read(13,*) x2(i), qy(i)
      enddo
      close(13)

      ! re-bin data to grid
      call addpnt(x2,qy,kdata+4,kdata,x2(1)*(1.d0-deltax),zero)
      call addpnt(x2,qy,kdata+4,kdata,               zero,zero)
      call addpnt(x2,qy,kdata+4,kdata,x2(kdata)*(1.d0+deltax),zero)
      call addpnt(x2,qy,kdata+4,kdata,            biggest,zero)
      call inter2(nw+1,wavl,yq1,kdata,x2,qy,ierr)
      if (ierr .ne. 0) THEN
        print*, "failed while reading in quantum yields for ", species
      endif

      do iw = 1, nw
        do i = 1, nz
          sq(jn,i,iw) = yg2(iw,i)*yq1(iw)
        enddo
      enddo

      jn = jn+1

      deallocate(qy)
      deallocate(x2)

    endif
  enddo
  close(12)

  ! deallocate
  ! all should auto-deallocate after going out of scope.
  deallocate(yy,yy_grid)
  deallocate(temperature_nums)
  deallocate(x1, y1)
  deallocate(dumby)

end subroutine
