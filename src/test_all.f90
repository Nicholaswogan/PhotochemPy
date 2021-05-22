

program test_all
  use photochem, only: setup, integrate, cvode_equilibrium, steam2photochem
  use photochem_data, only: nq, np
  use photochem_vars, only: rootdir, max_cvode_steps
  implicit none
  logical success
  character(len=1000) :: err
  character(len=:), allocatable :: template
  logical :: memtest
  integer :: nsteps
  
  ! for steam2photochem
  integer :: nz_in
  real(8) :: P_surf, P_top
  real(8), allocatable :: usol_layer(:), rpar_layer(:)
  
  err = ''
  rootdir = '../PhotochemPy/'
  
  print*,'nsteps:'
  read*, nsteps
  max_cvode_steps = nsteps
  
  template = "ModernEarth"
  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/settings.yaml', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_now.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  call integrate(nsteps, success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  call cvode_equilibrium(1.d-3, 1.d-27, .true., success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  template = "Archean2Proterozoic"
  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/settings.yaml', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_2.7Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  call integrate(nsteps, success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  call cvode_equilibrium(1.d-3, 1.d-27, .true., success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  template = "Archean+haze"
  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/settings.yaml', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_2.7Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  call integrate(nsteps, success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  call cvode_equilibrium(1.d-3, 1.d-27, .true., success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  template = "Hadean+HCN"
  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/settings.yaml', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_4.0Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  call integrate(nsteps, success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  call cvode_equilibrium(1.d-3, 1.d-27, .true., success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  
  
  
  template = "Hadean+HCN_test"
  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/settings.yaml', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_4.0Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    stop
  endif
  nz_in = 200
  P_surf = 7.d0
  P_top = 2.1427d-08
  allocate(usol_layer(nq))
  allocate(rpar_layer(np))
  usol_layer = 1.d-30
  usol_layer(65) = 0.05d0
  usol_layer(66) = 0.05d0
  rpar_layer = [1.d-5, 1.d-5, 1.d-7, 1.d-7]
  call steam2photochem(nq, np, nz_in, P_surf, P_top, usol_layer, rpar_layer, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  
  call cvode_equilibrium(1.d-3, 1.d-30, .true., success, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  
  
  
end program
