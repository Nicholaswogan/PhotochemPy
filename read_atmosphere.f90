
      subroutine read_atmosphere(atmosphere_txt)
        implicit none

        ! module variables
        ! character(len=8), allocatable, dimension(:) :: ISPEC
        ! integer :: nsp2
        ! integer :: nz ! number of vertical grid points
        ! real*8, allocatable, dimension(:,:) :: usol_init

        ! local variables
        character(len=*), intent(in) :: atmosphere_txt
        character(len=10000) :: line
        character(len=8), dimension(1000) :: arr1
        character(len=8),allocatable, dimension(:) :: arr2
        real*8,allocatable, dimension(:) :: temp
        integer :: i, n, io, j, k, ii, iii

        open(4, file=trim(atmosphere_txt),status='OLD',iostat=io)
        read(4,'(A)') line
        n = 0
        do i=1,1000
          read(line,*,iostat=io) arr1(1:i)
          if (io==-1) exit
          n = n+1
        enddo
        allocate(arr2(n))
        allocate(temp(n))
        read(line,*) (arr2(i),i=1,n)

        ! reads in mixing ratios
        iii = 0
        do i=1,nsp2
          do j=1,n
            if (arr2(j).eq.ispec(i)) then
              iii= iii+1
              do k=1,nz
                read(4,*) (temp(ii),ii=1,n)
                usol_init(i,k) = temp(j)
              enddo
              rewind(4) ! rewind!
              read(4,*) ! skip first line
              exit
            endif
          enddo
        enddo

        if (iii.ne.nq) then
          print*,'Error when reading in initial atmosphere.'
        endif

        rewind(4)
        read(4,*)
        ! reads in temperature
        do j=1,n
          if (trim(arr2(j)).eq.'temp') then
            do k=1,nz
              read(4,*) (temp(ii),ii=1,n)
              T(k) = temp(j)
            enddo
          endif
        enddo

        rewind(4)
        read(4,*)
        ! reads in density
        do j=1,n
          if (trim(arr2(j)).eq.'density') then
            do k=1,nz
              read(4,*) (temp(ii),ii=1,n)
              DEN(k) = temp(j)
            enddo
          endif
        enddo

        rewind(4)
        read(4,*)
        ! reads in eddy diffusion?
        do j=1,n
          if (trim(arr2(j)).eq.'eddy') then
            do k=1,nz
              read(4,*) (temp(ii),ii=1,n)
              EDD(k) = temp(j)
            enddo
          endif
        enddo

        close(4)

      end subroutine
