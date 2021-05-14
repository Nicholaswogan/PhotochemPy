

subroutine read_settings(set_file, err)
  use yaml, only : parse, error_length
  use yaml_types, only : type_node, type_dictionary, type_error
  implicit none
  
  ! will be global
  real(8) :: g, fscale, alb, ztrop, r0, p0
  character(len=20) :: planet
  
  character(len=*), intent(in) :: set_file
  character(len=1000), intent(out) :: err 
  
  character(error_length) :: error
  class (type_node), pointer :: root
  class (type_dictionary), pointer :: planet, settings
  type (type_error), pointer :: io_err
  err = ''
  
  ! parse yaml file
  root => parse(set_file,unit=100,error=error)
  if (error /= '') then
    err = trim(error)
    return
  end if
  select type (root)
    class is (type_dictionary)
      planet => root%get_dictionary('planet',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      
      planet%get_real('surface-gravity',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      

      call root%finalize()
      deallocate(root)
    class default
      err = "settings file must have dictionaries at root level"
      return
  end select

end subroutine

! subroutine planet_settings(planet)
!   class(type_dictionary), intent(in) :: planet
! 
!   implicit none
! 
! 
! end subroutine




program main
  implicit none
  
  character(len=1000) :: err 
  err = ''
  
  call read_settings('../input/templates/Hadean+HCN/settings.yaml', err)
  
  if (len_trim(err) /= 0) then
    print*,trim(err)
    print*,'error yaml'
    stop
  endif
  
end program