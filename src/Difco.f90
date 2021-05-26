
subroutine DIFCO(nq, nz, mubar_z, T, den, edd, &
                 hscale, tauedd, DK, H_atm, bx1x2, scale_H)
  use photochem_data, only: grav_z, mass, &
                            background_mu
  implicit none
  
  integer, intent(in) :: nq, nz
  real(8), intent(in) :: mubar_z(nz), T(nz), den(nz), edd(nz)
  
  real(8), intent(out) :: hscale(nz), tauedd(nz), DK(nz)
  real(8), intent(out) :: H_atm(nz)
  real(8), intent(out) :: bx1x2(nq,nz)
  real(8), intent(out) :: scale_H(nq,nz)
  
  
  ! local variables
  real*8 bkmg, eddav, denav
  real*8 h, tav
  integer i,j

  ! ***** DK(I) = K*N AT GRID STEP I+1/2 *****
  DO i = 1 , nz - 1
    eddav = SQRT(EDD(i)*EDD(i+1))
                                 !average eddy diffusion in grid center
    denav = SQRT(DEN(i)*DEN(i+1))
                               !average density at grid center
    DK(i) = eddav*denav
  ENDDO
  !
  !   COMPUTE DIFFUSION LIFETIME AT EACH HEIGHT (H*H/K)
  DO i = 1 , nz
    bkmg = 1.38E-16/(1.67E-24*mubar_z(i)*grav_z(i))  
    h = bkmg*T(i)
    HSCALE(i) = h
    TAUEDD(i) = h*h/EDD(i)
  ENDDO

  DO i = 1 , nz - 1 
    bkmg = 1.38E-16/(1.67E-24*mubar_z(i)*grav_z(i))  
    tav = SQRT(T(i)*T(i+1)) !average temperature at grid center
    H_ATM(i) = bkmg*tav
    !  compute scale heights of all the species in the atmosphere at all
    !   heights
    DO j = 1 , nq
      SCALE_H(j,i) = bkmg*tav*mubar_z(i)/MASS(j)
    ENDDO
    do j = 1,nq
      bx1x2(j,i) = 1.52d18*(1.d0/mass(j)+1.d0/mubar_z(i))**0.5d0*(tav**0.5)
    enddo
    
  ENDDO
  
  ! top layer of atmosphere
  bkmg = 1.38E-16/(1.67E-24*mubar_z(nz)*grav_z(nz))  
  H_ATM(nz) = bkmg*tav
  DO j = 1 , nq
     SCALE_H(j,nz) = bkmg*T(nz)*mubar_z(nz)/MASS(j)
  ENDDO
  do j = 1,nq
    bx1x2(j,NZ) = 1.52d18*(1.d0/mass(j)+1.d0/mubar_z(nz))**0.5d0*(T(nz)**0.5d0)
  enddo

end subroutine
