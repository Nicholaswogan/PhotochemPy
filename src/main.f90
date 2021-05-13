

program main
  use photochem, only: setup, integrate, allocate_memory, cvode_equilibrium
  use photochem_vars, only: rootdir,  max_cvode_steps, redox_factor!, usol_init
  implicit none
  logical converged
  logical success
  character(len=1000) :: err
  character(len=:), allocatable :: template
  err = ''

  rootdir = '../PhotochemPy/'
  template = "Hadean+HCN"

  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/planet.dat', &
             '../input/templates/'//template//'/input_photchem.dat', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_4.0Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  
  ! call integrate(10000,converged,err)
  ! max_cvode_steps = 100000
  ! print*,max_cvode_steps
  call cvode_equilibrium(1.d-3, 1.d-27, .true., success, err)
  
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  
  print*,redox_factor
  
end program
