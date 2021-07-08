
subroutine densty(nz, mubar_z, T, den, P, press)
  use photochem_data, only: grav_z, p0, &
                            dz
  implicit none
  
  integer, intent(in) :: nz
  real(8), intent(in) :: mubar_z(nz), T(nz)
  real(8), intent(out) :: den(nz), P(nz), press(nz)
  
  ! local varaibles
  real*8 rgas, bk
  real*8 roverm
  real*8 t0, p1, ha, tav
  integer i
!   THIS SUBROUTINE CALCULATES ATMOSPHERIC NUMBER DENSITIES, ASSUM-
!   ING HYDROSTATIC EQUILIBRIUM
  ! g0 = g
  rgas = 8.3143E7
  bk = 1.38054E-16

  
  roverm = rgas/mubar_z(1)

  t0 = T(1) + (T(1)-T(2))/2.
  ha = roverm*0.5*(t0+T(1))/grav_z(1)
  p1 = p0*1E6*EXP(-0.5*DZ(1)/ha)
  DEN(1) = p1/(bk*T(1))
!
! ***** FIND DENSITY FROM HYDROSTATIC EQUILIBRIUM *****
  DO i = 2 , nz
    roverm = rgas/mubar_z(i)
    ! r = r0 + Z(i)
    ! gz = g0*(r0/r)*(r0/r)
    tav = 0.5*(T(i)+T(i-1))
    ha = roverm*tav/grav_z(i)
    DEN(i) = DEN(i-1)*EXP(-DZ(i)/ha)*T(i-1)/T(i)
  ENDDO
! ***** FIND PRESSURE FROM THIS DENSITY *********
  DO i = 1 , nz
     PRESS(i) = DEN(i)*bk*T(i) ! dynes
  ENDDO
  ! pressure in bars
  do i = 1,nz
    P(i) =  DEN(i)*1.38E-16*T(i) * 1.0e-6 !bars
  enddo

end subroutine

subroutine mean_molecular_weight(nq, usol_layer, mass, background_mu, mubar)
  implicit none
  
  integer, intent(in) :: nq
  real(8), intent(in) :: usol_layer(nq)
  real(8), intent(in) :: mass(nq)
  real(8), intent(in) :: background_mu
  
  real(8), intent(out) :: mubar
  integer :: j
  real(8) :: f_background

  mubar = 0.d0
  do j = 1, nq
    mubar = mubar + usol_layer(j) * mass(j)
  enddo
  f_background = 1 - sum(usol_layer)
  mubar = mubar + f_background * background_mu
  
end subroutine




