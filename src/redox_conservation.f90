subroutine redox_conservation(flow,fup,sr)

  double precision, dimension(nq), intent(in) :: FLOW
  double precision, dimension(nq1), intent(in) :: FUP
  double precision, dimension(nq1), intent(in) :: SR


  double precision :: oxid_in_new, oxid_out_new, red_in_new
  double precision :: red_out_new, red_rain_new, oxy_rain_new
  integer i
  double precision :: redox_new, oxid_rain_new




  oxid_in_new=0.0d0
  oxid_out_new=0.0d0
  red_in_new=0.0d0
  red_out_new=0.0d0
  red_rain_new=0.0d0
  oxy_rain_new=0.0d0

  do i=1,NQ1
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

  !distributed fluxes
  !ACK - hardcoding - the below needs to be wrapped in an IF loop on the LBOUND...

  oxid_in_new=oxid_in_new + 2.0*distflux(LO2)
  red_in_new=red_in_new + distflux(LCO) + distflux(LH2) + 3.*distflux(LH2S)! +1.5*distflux(LHCL)    !ACK

  oxid_rain_new = oxid_out_new
  red_rain_new = red_rain_new

  !ok this needs to be finished up.  I also need to make sure that the boundary conditions are actually working as intended.
  !check in particular the distributed fluxes.  In general, this should be ready to go.
  !I have updated the sulfur balance, so will probably be able to follow the same scheme.
  redox_new = red_in_new - red_out_new -oxid_in_new + oxid_out_new &
              -red_Rain_new + oxy_rain_new

  ! print*,red_in_new , red_out_new ,oxid_in_new , oxid_out_new, red_Rain_new , oxy_rain_new
  ! print*,redox_new, redox_new/oxid_in_new   !mc - for ease in debugging lightning print to screen

  redox_factor = redox_new/oxid_in_new

end subroutine
