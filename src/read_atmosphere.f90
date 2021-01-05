
      subroutine read_atmosphere(atmosphere_txt)
        implicit none

        ! module variables
        ! character(len=8), allocatable, dimension(:) :: ISPEC
        ! integer :: nsp2
        ! integer :: nz ! number of vertical grid points
        ! real*8, allocatable, dimension(:,:) :: usol_init
        ! integer :: np ! number of particles

        ! local variables
        character(len=*), intent(in) :: atmosphere_txt
        character(len=10000) :: line
        character(len=8), dimension(1000) :: arr1
        character(len=24),allocatable, dimension(:) :: arr2
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
        do i=1,nq
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
          print*, iii, nq
        endif

        do i=1,nq
          if (trim(ispec(i)).eq.'CO2') then
            FCO2 = usol_init(i,1)
          endif
        enddo
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

        rewind(4)
        read(4,*)
        ! reads in aersol parameters
        if (np.gt.0) then
        do j=1,n
          do i=1,np
            if (trim(arr2(j)).eq.trim(ISPEC(size(ISPEC)-(nsp2-nq)-np+i))//'_AERSOL') then
              do k=1,nz
                read(4,*) (temp(ii),ii=1,n)
                aersol(k,i) = temp(j)
                aersol_init(k,i) = temp(j)
              enddo
              rewind(4)
              read(4,*)
              exit
            else if (trim(arr2(j)).eq.trim(ISPEC(size(ISPEC)-(nsp2-nq)-np+i))//'_WFALL') then
              do k=1,nz
                read(4,*) (temp(ii),ii=1,n)
                wfall(k,i) = temp(j)
                wfall_init(k,i) = temp(j)
              enddo
              rewind(4)
              read(4,*)
              exit
            else if (trim(arr2(j)).eq.trim(ISPEC(size(ISPEC)-(nsp2-nq)-np+i))//'_RPAR') then
              do k=1,nz
                read(4,*) (temp(ii),ii=1,n)
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

      end subroutine
