

subroutine read_settings(set_file, err)
  use yaml, only : parse, error_length
  use yaml_types, only : type_node, type_dictionary, type_error
  use photochem_data, only: grav_surf, fscale, alb, ztrop, r0, p0, &
                            AGL, EPSJ, light_disp_rate, hcdens, zy, Lgrid, &
                            IO2, ino, frak, ihztype, lightning, &
                            np, top_atmos, bottom_atmos, nz, rainout_on, &
                            H2O_strat_condensation, fix_water_in_troposphere, &
                            relative_humidity, use_manabe, confac, rhcold, &
                            estimate_CO2_photo_above_grid
  implicit none
  
  character(len=*), intent(in) :: set_file
  character(len=err_len), intent(out) :: err 
  
  character(error_length) :: error
  class (type_node), pointer :: root
  class (type_dictionary), pointer :: planet_dict, set_dict
  type (type_error), pointer :: io_err
  character(len=30) :: temp_char
  integer :: io
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
      ! planet = planet_dict%get_string('planet',error = io_err)
      ! if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      
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
      light_disp_rate = set_dict%get_real('lightning-dissipation-rate',error = io_err)
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
      
      H2O_strat_condensation = set_dict%get_logical('H2O-stratosphere-condensation', .true.,error = io_err)
      if (H2O_strat_condensation) then ! then we need rh of cold trap
        rhcold = set_dict%get_real('relative-humidity-cold-trap',error = io_err)
        if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      endif
      
      confac = set_dict%get_real('condensation-factor',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      
      fix_water_in_troposphere = set_dict%get_logical('fix-water-in-troposphere', .true.,error = io_err)
      relative_humidity = -1.d0
      if (fix_water_in_troposphere) then
        temp_char = set_dict%get_string('relative-humidity',error = io_err)
        if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
        
        read(temp_char,*,iostat = io) relative_humidity
        if (io /= 0) then
          ! it isn't a float
          if (trim(temp_char) == "manabe") then
            use_manabe = .true.
          else
            err = '"relative-humidity" can only be a number between 0 and 1, or "manabe". See '//trim(set_file)
            return 
          endif
        else
          use_manabe = .false.
        endif  
      endif

      nz = set_dict%get_integer('number-layers',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      bottom_atmos = set_dict%get_real('bottom-atmosphere',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      top_atmos = set_dict%get_real('top-atmosphere',error = io_err)
      if (associated(io_err)) then; err = trim(set_file)//trim(io_err%message); return; endif
      
      ! no error possible  
      estimate_CO2_photo_above_grid = set_dict%get_logical('estimate-CO2-photo-above-grid', &
                                      .true.,error = io_err)
    class default
      err = "settings file must have dictionaries at root level"
      return
  end select
  
  call root%finalize()
  deallocate(root)
  
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
  
  ! if ((planet /= "EARTH") .and. (planet /= "MARS")) then
  !   err = "Planet must be set to EARTH or MARS. See "//trim(set_file)
  !   return
  ! endif
  
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
  
  if (H2O_strat_condensation .and. .not. fix_water_in_troposphere) then
    err = "fix-water-in-troposhere must be on if stratospheric H2O condensation is on. See "//trim(set_file)
    return
  endif
  
  if (rainout_on .and. .not. fix_water_in_troposphere) then
    err = "fix-water-in-troposhere must be on if rainout is on. See "//trim(set_file)
    return
  endif
  

end subroutine

