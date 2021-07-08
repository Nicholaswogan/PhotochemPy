
subroutine ltning(nq, nz, usol, Psurf, light_disp_rate, prod_rates, err)
  use photochem_data, only: background_spec, &
                            LN2, LCO2, LCO, LH2, LH2O, LO2, LO, LNO
  use photochem_lightning, only: lightning_generalized
  implicit none
  integer, intent(in) :: nq,nz
  real(8), intent(in) :: usol(nq,nz)
  real(8), intent(in) :: Psurf, light_disp_rate
  real(8), intent(out) :: prod_rates(8)
  character(len=err_len), intent(out) :: err
  
  real(8) :: Px(8)
  
  err = ''

  ! PN2, PCO2, PCO, PH2, PH2O, PO2, PO, PNO
  if (background_spec == "N2") then
    Px(1) = (1.d0 - sum(usol(:,1)))*Psurf
  else
    Px(1) = usol(LN2,1)*Psurf
  endif
  if (background_spec == "CO2") then
    Px(2) = (1.d0 - sum(usol(:,1)))*Psurf
  else
    Px(2) = usol(LCO2,1)*Psurf
  endif
  Px(3) = usol(LCO,1)*Psurf
  if (background_spec == "H2") then
    Px(4) = (1.d0 - sum(usol(:,1)))*Psurf
  else
    Px(4) = usol(LH2,1)*Psurf
  endif
  Px(5) = usol(LH2O,1)*Psurf
  Px(6) = usol(LO2,1)*Psurf
  Px(7) = usol(LO,1)*Psurf
  Px(8) = usol(LNO,1)*Psurf

  call lightning_generalized(Px, Psurf, light_disp_rate, prod_rates, err)
  if (len_trim(err)/=0) return
  
end subroutine