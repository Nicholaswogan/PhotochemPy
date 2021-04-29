
subroutine densty
  use photochem_data, only: nz, far, g, p0, r0, &
                            lo2, dz, z
  use photochem_vars, only: den, T, P, Press, fco2, usol_init
  implicit none

  ! local varaibles
  real*8 g0, rgas, bk
  real*8 wt, ft, roverm
  real*8 t0, p1, ha, r, tav, gz
  integer i
!   THIS SUBROUTINE CALCULATES ATMOSPHERIC NUMBER DENSITIES, ASSUM-
!   ING HYDROSTATIC EQUILIBRIUM
  g0 = g
  rgas = 8.3143E7
  bk = 1.38054E-16

  ft = usol_init(lo2,1) + fco2 + far
  wt = usol_init(lo2,1)*32. + fco2*44. + (1.-ft)*28. + far*40.  !assuming O2,CO2,N2 and Ar are the main players and pN2=1-pO2+pCO2+pAr
  roverm = rgas/wt

  t0 = T(1) + (T(1)-T(2))/2.
  ha = roverm*0.5*(t0+T(1))/g0
  p1 = p0*1E6*EXP(-0.5*DZ(1)/ha)
  DEN(1) = p1/(bk*T(1))
!
! ***** FIND DENSITY FROM HYDROSTATIC EQUILIBRIUM *****
  DO i = 2 , nz
    r = r0 + Z(i)
    gz = g0*(r0/r)*(r0/r)
    tav = 0.5*(T(i)+T(i-1))
    ha = roverm*tav/gz
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

subroutine mean_molecular_weight(nq, nz, usol, mass, mu_background, mubar)
  implicit none
  
  integer, intent(in) :: nq, nz 
  real(8), intent(in) :: usol(nq,nz)
  real(8), intent(in) :: mass(nq)
  real(8), intent(in) :: mu_background
  
  real(8), intent(out) :: mubar(nz)
  integer :: i, j
  real(8) :: f_background

  mubar = 0.d0
  do i = 1, nz
    do j = 1, nq
      mubar(i) = mubar(i) + usol(j,i) * mass(j)
    enddo
    f_background = 1 - sum(usol(:,i))
    mubar(i) = mubar(i) + f_background * mu_background
  enddo
  
end subroutine




