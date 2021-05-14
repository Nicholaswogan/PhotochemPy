
subroutine read_planet(planet_dat,err)
  use photochem_data, only: G, Fscale, Alb, ztrop,far,R0,P0, planet
  implicit none

  ! local variables
  character(len=*), intent(in) :: planet_dat
  character(len=err_len), intent(out) :: err
  real*8 :: timeGa ! not going to actually read in
  integer :: msun ! same here
  integer :: ihzscale ! same
  integer :: ireset ! same
  integer :: io
  err = ''

  ! open reactions.rx
  open(9, file=trim(planet_dat),status='old')

502  FORMAT(F7.1/,F7.2/,F7.3/,E7.1/,F7.3/,E8.3/,F8.3/,A8/,F4.2/,I1/ &
  ,I2/,I1/,I3/,E10.4)

  READ(9,502,iostat = io) G,FSCALE,ALB,ZTROP,FAR,R0,P0,PLANET,TIMEGA,IRESET, &
                          msun, ihzscale
  if (io /= 0) then
    err = 'Problem reading '//trim(planet_dat)
    return
  endif
  close(9)


end subroutine
