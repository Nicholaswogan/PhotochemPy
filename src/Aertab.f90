
      SUBROUTINE AERTAB
        use photochem_data, only: nz, nf, nt, VH2O, VH2SO4, ftab
        use photochem_vars, only: rootdir, T
      IMPLICIT NONE

      ! module variables
      ! real*8, allocatable, dimension(:,:) :: VH2O
      ! real*8, allocatable, dimension(:,:) :: VH2SO4
      ! real*8, dimension(nf) :: ftab

      ! local variables
      real*8 :: ttab(nt) , ph2o(nt,nf) , ph2so4(nt,nf)
      integer i, j, is, is1, k

      real*8 fr, h2ol, h2ol1, h2so4l, h2sol1

!
!   THIS SUBROUTINE READS A TABLE OF SULFURIC ACID AND H2O VAPOR
!   PRESSURES AS FUNCTIONS OF TEMPERATURE AND CONCENTRATION OF
!   H2SO4 IN THE PARTICLES.  THEN IT PRODUCES A NEW TABLE IN WHICH
!   THE LOG OF THE VAPOR PRESSURES IS STORED AT EACH VERTICAL GRID
!   POINT OF THE MODEL.
!
!  to replace this i want a simple function
!
!   READ DATAFILE (VAPOR PRESSURES IN MM HG)
      open(2, file=trim(rootdir)//'DATA/aerosol.table',status='OLD')
      DO i = 1 , nt
         DO j = 1 , nf
            READ (2,99001) ph2o(i,j) , ph2so4(i,j) , ttab(i) , ftab(j)
99001       FORMAT (E13.5,2x,E13.5,2x,F4.0,1x,F6.2)
         ENDDO
      ENDDO
      close(2)

!
!   CONVERT VAPOR PRESSURES TO BARS
      DO k = 1 , nf
         DO j = 1 , nt
            ph2o(j,k) = ph2o(j,k)*1.013/760.
                                            !ACK - check if OK w.r.t pressure change - Kevin keeps this the same for Mars
            ph2so4(j,k) = ph2so4(j,k)*1.013/760.
         ENDDO
      ENDDO
!
!   INTERPOLATE TABLE TO TEMPERATURE AT EACH VERTICAL GRID POINT
      DO j = 1 , nz
         DO i = 1 , nt
            is = i
            IF ( ttab(i).GT.T(j) ) GOTO 50
         ENDDO
 50      is1 = MAX0(is-1,1)
!   T(J) LIES BETWEEN TTAB(IS) AND TTAB(IS1)
         fr = 1.
         IF ( is.GT.is1 ) fr = (T(j)-ttab(is1))/(ttab(is)-ttab(is1))
!
!   INTERPOLATE PH2O AND PH2SO4 LOGARITHMICALLY
         DO k = 1 , nf
            h2ol = LOG(ph2o(is,k))
            h2ol1 = LOG(ph2o(is1,k))
            h2so4l = LOG(ph2so4(is,k))
            h2sol1 = LOG(ph2so4(is1,k))
            VH2O(k,j) = fr*h2ol + (1.-fr)*h2ol1
            VH2SO4(k,j) = fr*h2so4l + (1.-fr)*h2sol1
         ENDDO
      ENDDO

      END
