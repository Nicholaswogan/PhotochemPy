module photochem_interface
  implicit none
  
  abstract interface
    function lbound_fcn_interface(tn) result(flux)
      real(8), value, intent(in) :: tn
      real(8) :: flux
    end function
  end interface
  
end module