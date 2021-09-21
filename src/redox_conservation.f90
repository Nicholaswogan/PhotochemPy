subroutine redox_conservation(flow,fup,sr)
  use photochem_data, only: nq, redoxstate, atomsS
  use photochem_vars, only: redox_factor, sulfur_factor, distflux, lbound
  implicit none
  real(8), dimension(nq), intent(in) :: FLOW ! flux at lower boundary
  real(8), dimension(nq), intent(in) :: FUP ! flux at upper boundary
  real(8), dimension(nq), intent(in) :: SR ! flux at

  real(8) :: oxid_in_new, oxid_out_new, red_in_new
  real(8) :: red_out_new, red_rain_new, oxy_rain_new
  real(8) :: redox_new, redox_arr(6)
  
  real(8) :: sulfur_in, sulfur_out, sulfur_con
  
  integer i
  
  !!! redox conservation !!!
  
  oxid_in_new=0.0d0
  oxid_out_new=0.0d0
  red_in_new=0.0d0
  red_out_new=0.0d0
  red_rain_new=0.0d0
  oxy_rain_new=0.0d0

  ! lower boundary, upper boundary, and rainout.
  do i=1,nq
    if (redoxstate(I) .GT. 0.) then
      oxid_in_new=oxid_in_new + FLOW(I)*redoxstate(I)
      oxid_out_new=oxid_out_new + FUP(I)*redoxstate(I)
      oxy_rain_new=oxy_rain_new + SR(I)*redoxstate(I)
    else if (redoxstate(I) .LT. 0) then
      red_in_new=red_in_new+ FLOW(I)*redoxstate(I)*(-1.0)
      red_out_new=red_out_new + FUP(I)*redoxstate(I)*(-1.0)
      red_rain_new=red_rain_new + SR(I)*redoxstate(I)*(-1.0)
    endif
  enddo

  ! distributed fluxes
  do i = 1,nq
    if (lbound(i) == 3) then
      if (redoxstate(i) > 0) then
        oxid_in_new = oxid_in_new + distflux(i)*redoxstate(i)
      elseif (redoxstate(i) < 0) then
        red_in_new = red_in_new + distflux(i)*redoxstate(i)*(-1.d0)
      endif
    endif
  enddo
  
  redox_new = red_in_new - red_out_new &
              - oxid_in_new + oxid_out_new & 
              - red_rain_new + oxy_rain_new
  
  redox_arr = [red_in_new, red_out_new, oxid_in_new, &
               oxid_out_new, red_rain_new, oxy_rain_new]
  
  ! redox_new should be small compared to largest redox fluxes
  redox_factor = redox_new/maxval(abs(redox_arr))
  
  !!! Sulfur conservation !!!
  sulfur_in = 0.d0
  sulfur_out = 0.d0
  
  do i = 1,nq
    ! lower boundary
    if (flow(i) > 0) then
      sulfur_in = sulfur_in + flow(i)*atomsS(i)
    elseif (flow(i) < 0) then
      sulfur_out = sulfur_out + flow(i)*atomsS(i)  
    endif
    ! upper boundary
    if (fup(i) < 0) then
      sulfur_in = sulfur_in - fup(i)*atomsS(i)
    elseif (fup(i) > 0) then
      sulfur_out = sulfur_out - fup(i)*atomsS(i)  
    endif
    ! rainout
    sulfur_out = sulfur_out - sr(i)*atomsS(i)
    ! distributed flux
    if (lbound(i) == 3) then
      if (distflux(i) > 0) then
        sulfur_in = sulfur_in + distflux(i)*atomsS(i)
      elseif (distflux(i) < 0) then
        sulfur_out = sulfur_out + distflux(i)*atomsS(i)  
      endif
    endif
  enddo
  
  sulfur_con = sulfur_in + sulfur_out
  sulfur_factor = sulfur_con/max(abs(sulfur_in),abs(sulfur_out))
  
end subroutine
