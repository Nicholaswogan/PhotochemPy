

program main
  use photochem, only: setup, integrate, cvode_equilibrium
  use photochem_data, only: ispec
  use photochem_vars, only: usol_out, rootdir
  implicit none
  logical success
  character(len=1000) :: err
  character(len=:), allocatable :: template
  integer :: ind(1)
  
  err = ''

  rootdir = '../../PhotochemPy/'
  template = "Archean2Proterozoic"

  call setup('../../input/templates/'//template//'/species.dat', &
             '../../input/templates/'//template//'/reactions.rx', &
             '../../input/templates/'//template//'/settings.yaml', &
             '../../input/templates/'//template//'/atmosphere.txt', &
             '../../input/templates/'//template//'/Sun_2.7Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  ! find photochemical equilibrium with backward euler
  call integrate(10000,success,err) 
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  ! or find photochemical equilibrium with CVODE BDF
  call cvode_equilibrium(1.d-3, 1.d-30, .true., success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  ! The output atmosphere is in usol_vars
  ind = findloc(ispec,'CH4')
  print*,'Surface CH4 mixing ratio = ', usol_out(ind(1),1)
  
end program
