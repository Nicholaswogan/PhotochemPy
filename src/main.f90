

program main
  use photochem, only: setup, integrate, allocate_memory, cvode_equilibrium, steam2photochem!, cvode_save!, read_settings
  use photochem_data, only: nq, np, z, nz, P0, LH2, LH, nz, dz
  use photochem_vars, only: rootdir,  max_cvode_steps, redox_factor, usol_init, veff, usol_out,den,flow
  use photochem_wrk, only: surf_radiance
  implicit none
  logical success
  character(len=1000) :: err
  character(len=:), allocatable :: template
  integer :: i
  ! real(8) :: t_eval(100)
  ! integer :: num_t_eval = 100
  
  ! for steam2photochem
  integer :: nz_in
  real(8) :: P_surf, P_top
  real(8), allocatable :: usol_layer(:), rpar_layer(:)
  
  err = ''

  rootdir = '../PhotochemPy/'
  template = "ModernEarth"

  call setup('../input/templates/'//template//'/species.dat', &
             '../input/templates/'//template//'/reactions.rx', &
             '../input/templates/'//template//'/settings.yaml', &
             '../input/templates/'//template//'/atmosphere.txt', &
             '../input/templates/'//template//'/Sun_now.txt', err)
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
