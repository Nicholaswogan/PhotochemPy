
module photochem_lightning
  implicit none
  integer, parameter :: real_kind = kind(1.d0)
  ! integer, parameter :: err_len = 1024
  
  ! Equilibrium constants for 3500 K, pg. 240 (Kasting 1979)
  real(real_kind), parameter :: K1 = 0.103d0
  real(real_kind), parameter :: K2 = 0.619d0
  real(real_kind), parameter :: K3 = 5.3d0
  real(real_kind), parameter :: K4 = 0.22d0

  ! for Kasting_equilibrate
  real(real_kind) :: CT, H2T, O2T
  real(real_kind) :: K1p, K2p, K3p, K4p, KT
  
  ! for wogan_equilibrate
  real(real_kind) :: Ptot_global
  real(real_kind) :: nN2T, nCT, nH2T, nO2T

contains
  
  subroutine lightning_scaled2NO(Px, Ptot, prod_NO, prod_rates, err)
    real(real_kind), intent(in) ::  Px(8), Ptot
    real(real_kind), intent(in) :: prod_NO 
    real(real_kind), intent(out) ::  prod_rates(8) ! molecules/cm2/s
    character(len=1000), intent(out) :: err
    real(real_kind) ::  Px_out(8)
    integer :: i
    
    err = ''
    if (prod_NO < 0.d0) then
      err = "NO production can not be negative (lightning_scaled2NO)"
      return
    endif
    call wogan_equilibrate(Px, Ptot, Px_out, err) 
    if (len_trim(err) /= 0) return
    
    do i = 1,8
      prod_rates(i) = prod_NO*(Px_out(i)/Px_out(8)) - &
                      prod_NO*(Px(i)/Px_out(8))
    enddo

  end subroutine

  subroutine lightning_generalized(Px, Ptot, light_disp_rate, prod_rates, err)
    real(real_kind), intent(in) ::  Px(8), Ptot
    real(real_kind), intent(in) :: light_disp_rate ! energy disipated by lightning (Watts/cm2)
    real(real_kind), intent(out) ::  prod_rates(8) ! molecules/cm2/s
    character(len=1000), intent(out) :: err
    
    real(real_kind) ::  Px_out(8)
    ! Average of range of constants determined
    ! in Chameides et al. (1981), Equation 4
    real(real_kind), parameter :: Chameides_const = 6.095d21
    real(real_kind), parameter :: T_freeze = 3500.d0 ! Kelvin. Freeze-out temperature of NO.
    
    integer :: i
    
    err = ''
    if (light_disp_rate < 0.d0) then
      err = "lighting dissipation rate can not be negative (lightning_generalized)"
      return
    endif
    
    call wogan_equilibrate(Px, Ptot, Px_out, err) 
    if (len_trim(err) /= 0) return
    do i = 1,8
      prod_rates(i) = ((Px_out(i)/Ptot)/T_freeze)*Chameides_const*light_disp_rate - &
                      ((Px(i)/Ptot)/T_freeze)*Chameides_const*light_disp_rate
    enddo
    
  end subroutine

  subroutine kasting_equilibrate(Px, Px_out, err)
    real(real_kind), intent(in) ::  Px(8)
    real(real_kind), intent(out) ::  Px_out(8)
    character(len=1000), intent(out) :: err
    
    real(real_kind) :: PN2, PCO2, PCO, PH2, PH2O, PO2, PO, PNO
    real(real_kind) :: PN2_out, PCO2_out, PCO_out, PH2_out, &
                       PH2O_out, PO2_out, PO_out, PNO_out
    real(real_kind) :: x1, x2, x3, x4, x5
    
    integer, parameter :: neqs = 1
    real(real_kind) :: x(neqs), fvec(neqs)
    real(real_kind), parameter :: tol = 2.d-8
    integer :: info
    integer, parameter :: lwa = (neqs*(3*neqs+13))/2+2
    real(real_kind) :: wa(lwa)
    ! This subroutine equilibrates gas mixtures of
    ! N2, CO2, CO, H2, H2O, O2, O, and NO. Following
    ! James Kasting Thesis, Appendix C, pg 233.
    ! "The Evolution of Oxygen and Ozone in Earth's atmosphere" (1979)
    ! However, this approach does not conserve mass!
    
    err = ''
    
    PN2 = Px(1)
    PCO2 = Px(2)
    PCO = Px(3)
    PH2 = Px(4)
    PH2O = Px(5)
    PO2 = Px(6)
    PO = Px(7)
    PNO = Px(8)
    
    if (PNO*1.d2 > PN2) then
      err = "NO is compariable in concentrations to N2. "// &
            "The lightning routine assumes this is not true. So turn off lightning!"
      return
    endif
    if (any(Px <= 0.d0)) then
      err = "You have input a negative partial pressure into 'kasting_equilibrate'."// &
            " This is not OK!"
      return
    endif

    CT = PCO2 + PCO
    H2T = PH2 + PH2O
    O2T = PO2 + PCO2 + 0.5d0*(PCO + PH2O + PO + PNO)
    
    K1p = K1*PN2**0.5d0*O2T**(-0.5d0)
    K2p = K2*O2T**0.5d0
    K3p = K3*O2T**0.5d0
    K4p = K4*O2T**(-0.5d0)
    KT = 0.5d0*(K1p + K4p)
    
    x(1) = log(0.5d0)
    ! hybrd from Minpack
    call hybrd1(fcn_kasting,neqs,x,fvec,tol,info,wa,lwa)
    if (info /= 1 .and. sum(abs(fvec)) > 1.d-4) then
      err = "hybrd failed to converged in Kasting lightning routine."
      return
    endif
    x3 = exp(x(1))
    if (x3 > 1.d0 .or. x3 < 0.d0) then
      err = "hybrd failed to converged to a "// &
            "physical solution in Kasting lighting routine."
      return 
    endif
    ! pg. 236
    x5 = K1p*x3**0.5d0
    x1 = 1.d0/(1.d0 + K2p*x3**0.5d0)
    x2 = 1.d0/(1.d0 + K3p*x3**0.5d0)
    x4 = K4p*x3**0.5d0
    
    PCO_out = x1*CT
    PH2_out = x2*H2T
    PO2_out = x3*O2T
    PO_out = x4*O2T
    PNO_out = x5*O2T
    PCO2_out = CT - PCO_out
    PH2O_out = H2T - PH2_out
    PN2_out = PN2
    
    Px_out = [PN2_out, PCO2_out, PCO_out, PH2_out, &
              PH2O_out, PO2_out, PO_out, PNO_out]
  
  end subroutine
  
  subroutine wogan_equilibrate(Px, Ptot, Px_out, err)
    real(real_kind), intent(in) ::  Px(8)
    real(real_kind), intent(in) :: Ptot
    real(real_kind), intent(out) ::  Px_out(8)
    character(len=1000), intent(out) :: err
    
    real(real_kind) :: nN2, nCO2, nCO, nH2, nH2O, nO2, nO, nNO, nT

    integer, parameter :: neqs = 9
    real(real_kind) :: x(neqs), fvec(neqs)
    real(real_kind), parameter :: tol = 2.d-8
    integer :: info
    integer, parameter :: lwa = (neqs*(3*neqs+13))/2+2
    real(real_kind) :: wa(lwa)
    ! This subroutine equilibrates gas mixtures of
    ! N2, CO2, CO, H2, H2O, O2, O, and NO to 3500 K
    ! to represent lightning.
    ! I adjusted James Kasting's approach because his doesn't 
    ! quite conserve atoms.

    err = ''
    if (Ptot <= 0.d0) then
      err = "The lighting routine recieved a negative total pressure."// &
            " The total pressure must be positive."
      return 
    endif
    if (Ptot < sum(Px)) then
      err = "The partial pressures sum to greater than the total pressure."
      return 
    endif
    Ptot_global = Ptot
    
    nN2 = Px(1)/Ptot
    nCO2 = Px(2)/Ptot
    nCO = Px(3)/Ptot
    nH2 = Px(4)/Ptot
    nH2O = Px(5)/Ptot
    nO2 = Px(6)/Ptot
    nO = Px(7)/Ptot
    nNO = Px(8)/Ptot
    nT = sum(Px)/Ptot
    
    nN2T = nN2 + 0.5d0*nNO
    nCT = nCO2 + nCO
    nH2T = nH2O + nH2
    nO2T = nO2 + nCO2 + 0.5d0*(nCO + nH2O + nO + nNO)
    
    call kasting_equilibrate(Px, Px_out, err)
    if (len_trim(err) /= 0) return
    
    x(1) = log(Px_out(1)/Ptot)
    x(2) = log(Px_out(2)/Ptot)
    x(3) = log(Px_out(3)/Ptot)
    x(4) = log(Px_out(4)/Ptot)
    x(5) = log(Px_out(5)/Ptot)
    x(6) = log(Px_out(6)/Ptot)
    x(7) = log(Px_out(7)/Ptot)
    x(8) = log(Px_out(8)/Ptot)
    x(9) = log(sum(Px_out)/Ptot)

    ! hybrd from Minpack
    call hybrd1(fcn_wogan,neqs,x,fvec,tol,info,wa,lwa)
    if (info /= 1 .and. sum(abs(fvec)) > 1.d-4) then
      err = "hybrd failed to converged in Wogan lightning routine."
      return 
    endif
    
    Px_out = exp(x(1:8))*Ptot

  end subroutine
  
  subroutine fcn_kasting(n,x,fvec,iflag)
    integer, intent(in) :: n
    real(real_kind), intent(in) :: x(n)
    real(real_kind), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    
    ! Equation C18, pg. 238.
    fvec(1) = exp(x(1))+KT*exp(x(1))**0.5d0 + &
              (1.d0 - 1.d0/(2.d0*(1.d0+ K2p*exp(x(1))**0.5d0)))*CT/O2T + &
              0.5d0*(1.d0-(1.d0/(1.d0+K3P*exp(x(1))**0.5d0)))*(H2T/O2T) - 1.d0
    
  end subroutine
  
  subroutine fcn_wogan(n,x,fvec,iflag)
    integer, intent(in) :: n
    real(real_kind), intent(in) :: x(n)
    real(real_kind), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    
    real(real_kind) :: nN2, nCO2, nCO, nH2, nH2O, nO2, nO, nNO, nT
    
    nN2 = exp(x(1))
    nCO2 = exp(x(2))
    nCO = exp(x(3))
    nH2 = exp(x(4))
    nH2O = exp(x(5))
    nO2 = exp(x(6))
    nO = exp(x(7))
    nNO = exp(x(8))
    nT = exp(x(9))
    
    fvec(1) = nNO - K1*nN2**0.5d0*nO2**0.5d0
    fvec(2) = nT**0.5d0*nCO2 - K2*nCO*nO2**0.5d0*Ptot_global**0.5d0
    fvec(3) = nT**0.5d0*nH2O - K3*nH2*nO2**0.5d0*Ptot_global**0.5d0
    fvec(4) = Ptot_global**0.5d0*nO - K4*nO2*nT**0.5d0
    fvec(5) = (nN2 + 0.5d0*nNO) - nN2T
    fvec(6) = (nCO2 + nCO) - nCT
    fvec(7) = (nH2 + nH2O) - nH2T
    fvec(8) = (nO2 + nCO2 + 0.5d0*(nCO + nH2O + nO + nNO)) - nO2T
    fvec(9) = sum(exp(x(1:8))) - nT
    
  end subroutine
  
end module


