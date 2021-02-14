subroutine XS(species,reactions_rx,kj,nz,nw,wavl,wav,T,sq)
  implicit none

  ! module variables

  ! inputs
  character(len=8), intent(in) :: species
  character(len=*), intent(in) :: reactions_rx
  integer, intent(in) :: kj, nz, nw
  double precision, dimension(nw), intent(in) :: wavl
  double precision, dimension(nw+1), intent(in) :: wav
  double precision, dimension(nz), intent(in) :: T

  ! in/out
  integer, dimension(kj,nz,nw), intent(inout) :: sq

  ! other vars
  character(len=50) :: fmt
  character(len=10) :: reac1, reac2, prod1, prod2, prod3
  character(len=5) :: label
  character(len=500) :: xsname
  character(len=1000) :: temperatureline, line
  character(len=8), dimension(100) :: temperatures
  integer io, stat, numtempcols, kdata
  double precision, allocatable, dimension(:,:) :: yy
  double precision, allocatable, dimension(:) :: x1
  double precision, allocatable, dimension(:) :: dumby

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
    if (io==-1) exit
  enddo
  ! now i-2 is number of temperature columns
  numtempcols = max(i-2,1)

  ! determine how many xsections points there are
  read(11,*)
  io = 0
  kdata = 0
  do while(io == 0)
    kdata = kdata + 1
    READ(11, *, IOSTAT=io)
  enddo

  ! allocate memory
  allocate(yy(kdata,numtempcols))
  allocate(xx(kdata))
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
    yy(i) = dumby(2:)
  enddo

  close(11)

  ! for now no T dependence. Only use T = 300 K.
  if (numtempcols.eq.1) then
    y1 = yy(:,1)
  else
    y1 = yy(:,)
  endif

  ! interpolate XS data
  CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
  CALL addpnt(x1,y1,kdata,n1,               zero,zero)
  CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
  CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
  CALL inter2(nw+1,wavl,yg1,n1,x1,y1,ierr)

  ! Quantum Yields
  fmt = '(A10,A10,A10,A10,A10,A5)'
  open(12, file=trim(reactions_rx),status='OLD')

  read(12,trim(fmt),IOSTAT=stat) reac1,reac2,prod1,prod2,prod3,label
  do while(stat == 0)
    read(12,trim(fmt),IOSTAT=stat) reac1,reac2,prod1,prod2,prod3,label

    if ((label.eq.'PHOTO').and.(species.eq.reac1)) then

      labnum = labnum + 1
      reacs = trim(reac1) // " " // trim(reac2)
      prods = trim(prod1) // " " // trim(prod2) // " " &
             // trim(prod3)

      do i=1, len(trim(reacs))
        if(reacs(i:i) == ' ') then
          reacs(i:i) = '_'
        endif
      end do

      do i=1, len(trim(prods))
        if(prods(i:i) == ' ') then
          prods(i:i) = '_'
        endif
      enddo

      ! open QY file

      ! fill sq for that photolysis reaction

    endif
  enddo

end subroutine XS


! subroutine read_xsection_file(species,nw,wavl,wav,xsection_data)
!   implicit none
! end subroutine
