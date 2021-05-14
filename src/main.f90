

program main
  use photochem, only: setup, integrate, allocate_memory, cvode_equilibrium!, read_settings
  ! use photochem_data, only: wavl, nw, flux
  use photochem_vars, only: rootdir,  max_cvode_steps, redox_factor!, usol_init
  use photochem_wrk, only: surf_radiance
  implicit none
  logical converged
  logical success
  character(len=1000) :: err
  character(len=:), allocatable :: template
  ! integer :: i
  err = ''
  
  
  ! call read_settings('../input/templates/Hadean+HCN/settings.yaml',err)
  ! if (len_trim(err) /= 0) then
  !   print*,trim(err)
  !   print*,'error yaml'
  !   stop
  ! endif
  ! stop
  

  rootdir = '../PhotochemPy/'
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
  
  call integrate(10000,converged,err)
  ! max_cvode_steps = 100000
  ! print*,max_cvode_steps
  ! call cvode_equilibrium(1.d-3, 1.d-27, .true., success, err)
  
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  
  ! do i = 1,nw
    ! print*,wavl(i)/10.d0,surf_radiance(i)
  ! enddo
  
end program
