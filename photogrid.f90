
      subroutine photogrid(top_atmosphere)
        implicit none

        ! module variables
        ! integer :: nz ! number of vertical grid points
        ! real*8, allocatable, dimension(:) :: z
        ! real*8, allocatable, dimension(:) :: dz

        ! local variables
        real*8, intent(in) :: top_atmosphere
        real*8 :: dzgrid
        integer :: i

        dzgrid = top_atmosphere/nz

        do I=1,NZ
         Z(I) = (I - 0.5)*dzgrid
        enddo


        DZ(1)=Z(1)+0.5*dzgrid
        do I=2,NZ
          DZ(I)=Z(I)-Z(I-1)
        enddo

      end subroutine
