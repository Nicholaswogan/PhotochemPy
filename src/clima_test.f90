
program main
  use photochem_clima
  
  implicit none
  real(8) :: T_trop, P_surf, P_trop, mubar, grav, eta, Ptop
  integer, parameter :: nz = 200
  real(8) :: zout(nz), Tout(nz), ztrop, edd(nz), Pout(nz)
  
  integer :: i
  character(len=1000) :: err
  err = ""
  
  T_trop = 180.d0
  P_surf = 1.00d0
  P_trop = 0.1d0
  mubar = 7.d0
  grav = 981.d0
  eta = 0.14d0
  Ptop = 2.1427d-08

  call eta_clima(T_trop, P_surf, P_trop, eta, mubar, grav, Ptop, nz, &
                zout, Tout, ztrop)
                
  call pahlevan_H2_clima(P_surf, mubar, grav, Ptop, nz, &
                         zout, Pout, Tout, ztrop, err)
  if (len_trim(err) /= 0) return
  
  call zahnle_eddy(nz, zout, 0.21d0, 235.d0, ztrop, mubar, grav, edd)
                
  do i = 1,nz            
    print*,zout(i)/1.d5, Tout(i), Pout(i)
  enddo
  print*,ztrop/1.d5

end program