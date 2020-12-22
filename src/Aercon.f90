
      SUBROUTINE AERCON(usol, nq, nz)

      IMPLICIT NONE


      ! module variables
      ! real*8, allocatable, dimension(:) :: FSULF
      ! real*8, allocatable, dimension(:) :: H2SO4S
      ! real*8, allocatable, dimension(:) :: S8S
      ! integer nf


      ! local variables
      integer, intent(in) :: nq, nz
      real*8, dimension(nq,nz), intent(in) :: usol

      real*8 zero, one, fr, h2so4l, psats8, plog10
      real*8 ph2ol(nz)
      integer j, k, kt, kt1

!
!   THIS SUBROUTINE FINDS THE WEIGHT PERCENT OF H2SO4 IN THE
!   PARTICLES AND THE H2SO4 VAPOR PRESSURE, GIVEN THE TEMPERATURE
!   AND H2O CONCENTRATION AT EACH ALTITUDE.
!          IT ALSO CALCULATES THE SATURATION VAPOR PRESSURE OVER
!     SOLID (ORTHORHOMBIC) S8.
!
      zero = 0.0
      one = 1.0

      DO j = 1 , nz
         ph2ol(j) = LOG(USOL(lh2o,j)*P(j)) ! P is pressure in dynes
      ENDDO

      DO j = 1 , nz
         DO k = 1 , nf
            kt = k
            IF ( VH2O(k,j).LT.ph2ol(j) ) GOTO 50
         ENDDO
 50      kt1 = MAX(kt-1,1)
!
!   PH2OL(J) LIES BETWEEN VH2O(KT,J) AND VH2O(KT1,J)
         fr = 1.
         IF ( kt.GT.kt1 ) fr = (ph2ol(j)-VH2O(kt1,j))                   &
                             & /(VH2O(kt,j)-VH2O(kt1,j))
         fr = MAX(fr,zero)
         fr = MIN(fr,one)
         FSULF(j) = (fr*FTAB(kt)+(1.-fr)*FTAB(kt1))*0.01
         h2so4l = fr*VH2SO4(kt,j) + (1.-fr)*VH2SO4(kt1,j)
         H2SO4S(j) = EXP(h2so4l)/P(j)
      ENDDO
!
!   FIND THE SATURATION MIXING RATIO OF S8
      DO j = 1 , nz
         plog10 = (11.664-5166./T(j))*0.981
         psats8 = 10.**plog10*1.013/760.
         S8S(j) = psats8/P(j)
      ENDDO
!
      END
