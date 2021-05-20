module photochem_clima
  implicit none
  private
  
  public :: eta_clima, pahlevan_H2_clima, zahnle_eddy
  
  integer, parameter :: real_kind = kind(1.d0)
  ! globals
  real(real_kind) :: mubar_g, grav_g, z_g, P_surf_g
  real(real_kind), parameter :: k_boltz = 1.3807d-16 ! boltzmann constant (cgs units)
  real(real_kind), parameter :: n_avo = 6.022d23 ! avo's number
  ! 3rd order polynomial fit to T(log10P) in the troposphere
  ! from Pahlevan et al. 2019
  real(real_kind), parameter :: poly(4) = [295.74510905, 76.02965729, &
                                           -3.60269118, 9.88900334] 
                                   
contains

  subroutine eta_clima(T_trop, P_surf, P_trop, eta, mubar, grav, P_top, nz, &
                       zout, Tout, ztrop)
    implicit none
    real(real_kind), intent(in) :: T_trop, P_surf, P_trop, eta
    real(real_kind), intent(in) :: mubar, grav, P_top
    integer, intent(in) :: nz
    
    real(real_kind), intent(out) :: zout(nz), Tout(nz), ztrop
    
    real(real_kind) :: T_surf, m_slope, z_top, dz
    
    integer :: i
    
    T_surf = T_trop *(P_surf/P_trop)**eta

    ! tropopause altitude (cm)
    ztrop = -k_boltz/(mubar/N_avo*grav)* &
            (T_trop-T_surf)*(1.d0/dlog(T_trop/T_surf)) &
            *dlog(P_trop/P_surf)
            
    m_slope = ((T_trop-T_surf)/ztrop)
    
    z_top = ztrop + (dlog(P_top/P_surf) + mubar* grav/(N_avo*k_boltz*m_slope)*dlog(T_trop/T_surf)) &
          *(-N_avo*k_boltz*T_trop)/(mubar*grav)
          
    ! compute z
    dz = z_top/nz
    zout(1) = dz/2.d0
    do i = 2,nz
      zout(i) = zout(i-1) + dz
    enddo
    
    do i = 1,nz
      if (zout(i) <= ztrop) then
        Tout(i) = m_slope*zout(i)+T_surf
      elseif (zout(i) > ztrop) then
        Tout(i) = T_trop
      endif
    enddo    

  end subroutine
  
  subroutine pahlevan_H2_clima(P_surf, mubar, grav, P_top, nz, &
                               zout, Pout, Tout, ztrop, err)
    
    real(real_kind), intent(in) :: P_surf
    real(real_kind), intent(in) :: mubar, grav, P_top
    integer, intent(in) :: nz
    
    real(real_kind), intent(out) :: zout(nz), Pout(nz), Tout(nz), ztrop
    
    character(len=1000), intent(out) :: err

    real(real_kind), parameter :: T_trop = 235.d0
    real(real_kind), parameter :: P_trop = 0.21d0
    
    real(real_kind) :: log10P, z_top, dz
    integer :: i
    
    !!! minpack
    integer, parameter :: neqs = 1
    real(real_kind) :: x(neqs), fvec(neqs)
    real(real_kind), parameter :: tol = 2.d-8
    integer :: info
    integer, parameter :: lwa = (neqs*(3*neqs+13))/2+2
    real(real_kind) :: wa(lwa)
    !!! minpack

    err = ""
    
    if (P_surf < 0.75d0) then
      err = "Surface pressure must be > 0.75 bar in H2 climate model"
      return
    elseif (P_surf > 200.0d0) then
      err = "Surface pressure must be < 200.00 bar in H2 climate model"
      return
    endif
             
    call alt_trop(P_trop, P_surf, mubar, grav, ztrop)
    call alt_above_trop(P_top, P_trop, T_trop, ztrop, mubar, grav, z_top)
    
    ! compute z
    dz = z_top/nz
    zout(1) = dz/2.d0
    do i = 2,nz
      zout(i) = zout(i-1) + dz
    enddo
    
    ! compute T
    mubar_g = mubar
    grav_g = grav
    P_surf_g = P_surf
    x(1) = P_surf*0.999d0
    do i = 1,nz
      if (zout(i) <= ztrop) then
        ! compute P
        z_g = zout(i)
        call hybrd1(fcn_pahlevan,neqs,x,fvec,tol,info,wa,lwa)
        if (info /= 1) then
          err = "hybrd failed to find a solution in pahlevan_H2_clima"
          return 
        endif
        ! now compute T
        log10P = dlog10(x(1))
        Tout(i) = poly(1) + poly(2)*log10P +  &
                  poly(3)*log10P**2.d0 + poly(4)*log10P**3.d0
        ! pressure
        Pout(i) = x(1)
      elseif (zout(i) > ztrop) then
        Tout(i) = T_trop
        call P_above_trop(zout(i), P_trop, T_trop, ztrop, mubar, grav, Pout(i))
      endif
    enddo
    
  end subroutine
  
  ! given z, compute eddy (zahnle et al. 2006, pg. 272)
  subroutine zahnle_eddy(nz, z, P_trop, T_trop, ztrop, mubar, grav, edd)
    integer, intent(in) :: nz
    real(real_kind), intent(in) :: z(nz), P_trop, T_trop, ztrop, mubar, grav
    real(real_kind), intent(out) :: edd(nz)
     
    real(real_kind) :: P
    integer :: i
     
    do i = 1,nz
      if (z(i) <= ztrop) then
        edd(i) = 1.d5
      elseif (z(i) > ztrop) then
        call P_above_trop(z(i), P_trop, T_trop, ztrop, mubar, grav, P)
        edd(i) = 1.d3*(P_trop/P)**0.5d0
        if (edd(i) > 1.d6) then
          edd(i) = 1.d6
        endif
      endif
    enddo
    
  end subroutine
  
  ! compute z for a given P in troposphere
  subroutine alt_trop(P, P_surf, mubar, grav, z)
    real(real_kind), intent(in) :: P, P_surf, mubar, grav
    real(real_kind), intent(out) :: z
    
    real(real_kind) :: term1, term2
    
    term1 = -(poly(1)*dlog(P) + poly(2)*dlog(P)**2/dlog(100.d0) + &
              poly(3)*dlog(P)**3/(3.d0*dlog(10.d0)**2) + &
              poly(4)*dlog(P)**4/(4.d0*dlog(10.d0)**3))
    term2 = poly(1)*dlog(P_surf) + poly(2)*dlog(P_surf)**2/dlog(100.d0) + &
            poly(3)*dlog(P_surf)**3/(3.d0*dlog(10.d0)**2) + &
            poly(4)*dlog(P_surf)**4/(4.d0*dlog(10.d0)**3)
    z = (term1+term2)*k_boltz*n_avo/(mubar*grav)
  end subroutine
  
  ! compute z for a given P above tropopause (isothermal)
  subroutine alt_above_trop(P, P_trop, T_trop, ztrop, mubar, grav, z)
    real(real_kind), intent(in) :: P, P_trop, T_trop, ztrop, mubar, grav
    real(real_kind), intent(out) :: z
    z = dlog(P/P_trop)*(-k_boltz*n_avo*T_trop)/(mubar*grav) + ztrop
  end subroutine
  
  ! compute P for a given z above tropopause (isothermal)
  subroutine P_above_trop(z, P_trop, T_trop, ztrop, mubar, grav, P)
    real(real_kind), intent(in) :: z, P_trop, T_trop, ztrop, mubar, grav
    real(real_kind), intent(out) :: P
    P = P_trop*dexp((-mubar*grav)/(k_boltz*n_avo*T_trop)*(z-ztrop))
  end subroutine
  
  ! for solving for P at a given z in troposphere
  subroutine fcn_pahlevan(n,x,fvec,iflag)
    integer, intent(in) :: n
    real(real_kind), intent(in) :: x(n)
    real(real_kind), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    real(real_kind) :: P, term1, term2, term3
    P = x(1)
    term1 = -(poly(1)*dlog(P) + poly(2)*dlog(P)**2/dlog(100.d0) + &
              poly(3)*dlog(P)**3/(3.d0*dlog(10.d0)**2) + &
              poly(4)*dlog(P)**4/(4.d0*dlog(10.d0)**3))
    term2 = poly(1)*dlog(P_surf_g) + poly(2)*dlog(P_surf_g)**2/dlog(100.d0) + &
            poly(3)*dlog(P_surf_g)**3/(3.d0*dlog(10.d0)**2) + &
            poly(4)*dlog(P_surf_g)**4/(4.d0*dlog(10.d0)**3)
    term3 = -(mubar_g*grav_g*z_g)/(k_boltz*n_avo)
    fvec(1) = term1+term2+term3
  end subroutine
  
end module
