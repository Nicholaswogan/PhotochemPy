subroutine Xsections_general(species,reactions_rx,kj,nz,nw,kw,wavl,jn,sq)
  implicit none

  ! module variables

  ! inputs
  character(len=8), intent(in) :: species
  character(len=*), intent(in) :: reactions_rx
  integer, intent(in) :: kj, nz, nw, kw
  double precision, dimension(nw+1), intent(in) :: wavl
  ! double precision, dimension(nz), intent(in) :: T

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
  integer io, stat, numtempcols, kdata, i, ierr, iw
  double precision, allocatable, dimension(:,:) :: yy
  double precision, allocatable, dimension(:) :: x1, y1, x2, qy
  double precision, allocatable, dimension(:) :: dumby
  double precision :: deltax = 1.E-4, biggest=1.E+36, zero=0.0
  double precision, dimension(nw) :: yg1, yq1

  character(len=500) :: rootdir
  rootdir = '../PhotochemPy/'

  ! initialize
  ierr = 0

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

  ! for now no T dependence. Only use T = 300 K.
  if (numtempcols.eq.1) then
    y1 = yy(:,1)
  else
    do i=1,numtempcols
      if (trim(temperatures(i)).eq.'300K') exit
    enddo
    y1 = yy(:,i)
  endif

  ! if (numtempcols.eq.-1) then
  !   ! make an array of cross sections
  !   allocate(temperature_nums(numtempcols+2))
  !   allocate(yyy(numtempcols+2))
  !
  !   do i=2,numtempcols
  !     read(temperatures(:index(temperatures,'K')-1),*) temperature_nums(i)
  !   enddo
  !   temperature_nums(1) = 0.d0
  !   temperature_nums(numtempcols+1) = 3000.d0
  !
  !   do i=1,kdata
  !     yyy(2:numtempcols) = yy(i,:)
  !     yyy(1) = yy(i,1)
  !     yyy(numtempcols+1) = yy(i,numtempcols)
  !     do j=1,nz
  !       call inter2(nw+1,T(j),yg11,kdata,x1,y1,ierr)
  !       yg1(j,i) =  yg11
  !     enddo
  !   enddo
  ! endif

  ! interpolate XS data
  call addpnt(x1,y1,kdata+4,kdata,x1(1)*(1.-deltax),zero)
  call addpnt(x1,y1,kdata+4,kdata,               zero,zero)
  call addpnt(x1,y1,kdata+4,kdata,x1(kdata)*(1.+deltax),zero)
  call addpnt(x1,y1,kdata+4,kdata,            biggest,zero)
  call inter2(nw+1,wavl,yg1,kdata,x1,y1,ierr)
  IF (ierr .NE. 0) THEN
    print*, "failed while reading in cross sections for ", species
  ENDIF

  ! Quantum Yields
  fmt = '(A10,A10,A10,A10,A8,A5)'
  open(12, file=trim(reactions_rx),status='OLD')
  io = 0
  do while(io == 0)
    read(12,trim(fmt),IOSTAT=io) reac1,reac2,prod1,prod2,prod3,label

    ! if photolysis reaction involving species
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
      call addpnt(x2,qy,kdata+4,kdata,x2(1)*(1.d0-deltax),zero)
      call addpnt(x2,qy,kdata+4,kdata,               zero,zero)
      call addpnt(x2,qy,kdata+4,kdata,x2(kdata)*(1.d0+deltax),zero)
      call addpnt(x2,qy,kdata+4,kdata,            biggest,zero)
      call inter2(nw+1,wavl,yq1,kdata,x2,qy,ierr)
      if (ierr .NE. 0) THEN
        print*, "failed while reading in quantum yields for ", species
      endif

      do iw = 1, nw
        do i = 1, nz
          sq(jn,i,iw) = yg1(iw)*yq1(iw)
        enddo
      enddo

      jn = jn+1

      deallocate(qy)
      deallocate(x2)

    endif
  enddo
  close(12)

  !deallocate
  deallocate(yy)
  deallocate(x1, y1)
  deallocate(dumby)

end subroutine


! program test
!   implicit none
!   integer, parameter :: kj = 61,nz = 200, nw = 100, kw = 1000
!   double precision,dimension(kj,nz,nw) :: sq
!   double precision,dimension(nw+1) :: wavl
!   integer j, i
!
!   sq = 0.d0
!   wavl = 0.d0
!   j = 1
!   call linspace(1210.d0,8000.d0,wavl,nw+1)
!
!   call Xsections_general('CH4     ','../input/templates/Archean+haze/reactions.rx',kj,nz,nw,kw,wavl,j,sq)
!   print*,sq(1,nz,:)
!
!
! end program
!
!
! subroutine linspace(from, to, array,n)
!   real(8), intent(in) :: from, to
!   real(8), intent(out) :: array(n)
!   real(8) :: range
!   integer :: n
!   integer :: i
!   range = to - from
!
!   if (n == 0) return
!
!   if (n == 1) then
!     array(1) = from
!     return
!   end if
!
!   do i=1, n
!     array(i) = from + range * (i - 1) / (n - 1)
!   end do
! end subroutine
