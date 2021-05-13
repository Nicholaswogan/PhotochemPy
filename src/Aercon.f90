
subroutine aercon(nq, nz, usol, P, T, H2SO4S, S8S, fsulf)
  use photochem_data, only: nf, lh2o, vH2O, ftab, VH2SO4
  implicit none

  ! local variables
  integer, intent(in) :: nq, nz
  real*8, dimension(nq,nz), intent(in) :: usol
  real(8), intent(in) :: T(nz), P(nz)
  real(8), intent(out) :: H2SO4S(nz), S8S(nz), fsulf(nz)

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
  zero = 0.0d0
  one = 1.0d0

  do j = 1 , nz
    ph2ol(j) = log(USOL(lh2o,j)*P(j)) ! P is pressure in dynes
  enddo

  do j = 1 , nz
    do k = 1 , nf
      kt = k
      if ( VH2O(k,j) < ph2ol(j) ) then
        kt1 = max(kt-1,1)
        exit
      endif
    enddo

  !   PH2OL(J) LIES BETWEEN VH2O(KT,J) AND VH2O(KT1,J)
    fr = 1.
    if ( kt.GT.kt1 ) fr = (ph2ol(j)-VH2O(kt1,j)) &
                          /(VH2O(kt,j)-VH2O(kt1,j))
    fr = max(fr,zero)
    fr = min(fr,one)
    FSULF(j) = (fr*FTAB(kt)+(1.-fr)*FTAB(kt1))*0.01
    h2so4l = fr*VH2SO4(kt,j) + (1.-fr)*VH2SO4(kt1,j)
    H2SO4S(j) = EXP(h2so4l)/P(j)
  enddo

  !   FIND THE SATURATION MIXING RATIO OF S8
  do j = 1 , nz
    plog10 = (11.664-5166./T(j))*0.981
    psats8 = 10.**plog10*1.013/760.
    S8S(j) = psats8/P(j)
  enddo

end subroutine
