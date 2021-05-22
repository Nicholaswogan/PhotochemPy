

subroutine read_settings(set_file, err)
  use yaml, only : parse, error_length
  use yaml_types, only : type_node, type_dictionary, type_error
  use photochem_data, only: grav_surf, fscale, alb, ztrop, r0, p0, planet, &
                            AGL, EPSJ, prono, hcdens, zy, Lgrid, &
                            IO2, ino, frak, ihztype, lightning, &
                            np, top_atmos, bottom_atmos, nz, rainout_on
  implicit none
  
  character(len=*), intent(in) :: set_file
  character(len=err_len), intent(out) :: err 
  
  character(error_length) :: error
  class (type_node), pointer :: root
  class (type_dictionary), pointer :: planet_dict, set_dict
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
      planet_dict => root%get_dictionary('planet',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      
      grav_surf = planet_dict%get_real('surface-gravity',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      fscale = planet_dict%get_real('fscale',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      alb = planet_dict%get_real('surface-albedo',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      ztrop = planet_dict%get_real('tropopause-altitude',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      r0 = planet_dict%get_real('planet-radius',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      p0 = planet_dict%get_real('surface-pressure',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      planet = planet_dict%get_string('planet',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      
      set_dict => root%get_dictionary('settings',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      agl = set_dict%get_real('diurnal-averaging-factor',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      lgrid = set_dict%get_integer('lgrid',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      io2 = set_dict%get_integer('io2',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      ino = set_dict%get_integer('ino',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      epsj = set_dict%get_real('epsj',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      prono = set_dict%get_real('prono',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      frak = set_dict%get_integer('frak',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      hcdens = set_dict%get_real('hcdens',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      ihztype = set_dict%get_integer('ihztype',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      zy = set_dict%get_real('solar-zenith-angle',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      lightning = set_dict%get_logical('lightning',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      rainout_on = set_dict%get_logical('rainout-on',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      nz = set_dict%get_integer('number-layers',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      bottom_atmos = set_dict%get_real('bottom-atmosphere',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      top_atmos = set_dict%get_real('top-atmosphere',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
  
      call root%finalize()
      deallocate(root)
    class default
      err = "settings file must have dictionaries at root level"
      return
  end select
  
  if ((frak == 1) .and. (np < 3)) then
    err = 'frak must be 0 if no hydrocarbons are present. See '//trim(set_file)
    return
  endif
  if (IO2 /= 2) then
    err = "IO2 = 2 is the only option. See "//trim(set_file)
    return
  endif
  if ((iNO /= 0) .and. (iNO /= 1)) then
    err = "INO = 0 or 1 are the only options. See "//trim(set_file)
    return
  endif
  
  if ((planet /= "EARTH") .and. (planet /= "MARS")) then
    err = "Planet must be set to EARTH or MARS. See "//trim(set_file)
    return
  endif
  
  if (nz < 50) then
    err = "nz can not be less than 50. See "//trim(set_file)
    return
  endif
  if (nz > 800) then
    err = "nz can not be >800 See "//trim(set_file)
    return
  endif
  
  if (bottom_atmos < 0.d0) then
    err = "bottom of the atmosphere can not be less than 0. See "//trim(set_file)
    return
  endif
  
  if (bottom_atmos > top_atmos) then
    err = "bottom of the atmosphere must be lower than the top. See "//trim(set_file)
    return
  endif
  
  if ((.not.rainout_on) .and. (lightning)) then
    err = "rainout must be on if lightning is on. See "//trim(set_file)
    return
  endif
  

end subroutine

