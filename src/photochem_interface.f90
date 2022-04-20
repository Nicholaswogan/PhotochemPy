module photochem_interface
  implicit none
  
  abstract interface
    function lbound_fcn_interface(tn) result(flux)
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: tn
      real(c_double) :: flux
    end function
    
    subroutine time_dependent_flux_fcn(tn, nw, photon_flux)
      use iso_c_binding, only: c_double, c_int
      real(c_double), value, intent(in) :: tn
      integer(c_int), value, intent(in) :: nw
      real(c_double), intent(out) :: photon_flux(nw)
    end subroutine
    
  end interface
  
end module