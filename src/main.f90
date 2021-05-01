

program main
  use photochem, only: setup, integrate, allocate_memory, cvode_equilibrium
  use photochem_vars, only: rootdir
  implicit none
  integer converged
  logical success
  character(len=1000) :: err

  rootdir = '../PhotochemPy/'

  call setup('../input/templates/Archean2Proterozoic/species.dat', &
             '../input/templates/Archean2Proterozoic/reactions.rx', &
             '../input/templates/Archean2Proterozoic/planet.dat', &
             '../input/templates/Archean2Proterozoic/input_photchem.dat', &
             '../input/templates/Archean2Proterozoic/atmosphere.txt', &
             '../input/templates/Archean2Proterozoic/Sun_2.7Ga.txt', err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error worked properly'
    stop
  endif
  
  ! call integrate(1,converged)
  call cvode_equilibrium(1.d-3, 1.d-27, .true., success)
  ! call right_hand_side(usol_flat,rhs,neq)
  
end program
