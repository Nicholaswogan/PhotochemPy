

integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                              result(ierr) bind(c, name='RhsFn')
  use, intrinsic :: iso_c_binding
  use fcvode_mod
  use fsundials_nvector_mod
  use photochem_data, only: neq
  use photochem_vars, only: max_cvode_steps, verbose
  use photochem_wrk, only: nsteps_previous, global_err, cvode_mem
  implicit none
  ! calling variables
  real(c_double), value :: tn        ! current time
  type(N_Vector)        :: sunvec_y  ! solution N_Vector
  type(N_Vector)        :: sunvec_f  ! rhs N_Vector
  type(c_ptr), value    :: user_data ! user-defined dat
  
  integer(c_long) :: nsteps(1)
  integer(c_int) :: loc_ierr
  real(c_double) :: hcur(1)
  
  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: yvec(:)
  real(c_double), pointer :: fvec(:)
  
  ierr = 0
  
  ! get data arrays from SUNDIALS vectors
  yvec(1:neq) => FN_VGetArrayPointer(sunvec_y)
  fvec(1:neq) => FN_VGetArrayPointer(sunvec_f)
  
  ! fill RHS vector
  call right_hand_side(yvec, fvec, neq, global_err)
  if (len_trim(global_err) /= 0) then
    if (trim(global_err) == "Mixing ratios sum to > 1.0. "//&
                          "Atmosphere is probably in a run-away state.") then
      ierr = -1
    else
      print*,trim(global_err)
      ierr = 1
    endif
  endif
  loc_ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
  
  if ((nsteps(1) /= nsteps_previous) .and. (verbose)) then
    loc_ierr = FCVodeGetCurrentStep(cvode_mem, hcur)
    print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
         nsteps, tn, hcur(1),maxval(abs(fvec))
  endif

  nsteps_previous = nsteps(1)
  
  return
end function

integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                              tmp1, tmp2, tmp3) &
                              result(ierr) bind(C,name='JacFn')
  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fcvode_mod
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use fsunmatrix_band_mod
  use fsundials_matrix_mod

  use photochem_data, only: neq, lda
  use photochem_wrk, only: global_err
  implicit none
  
  ! calling variables
  real(c_double), value :: tn        ! current time
  type(N_Vector)        :: sunvec_y  ! solution N_Vector
  type(N_Vector)        :: sunvec_f
  type(SUNMatrix)        :: sunmat_J  ! rhs N_Vector
  type(c_ptr), value    :: user_data ! user-defined data
  type(N_Vector)        :: tmp1, tmp2, tmp3

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: yvec(:)
  real(c_double), pointer :: djac(:,:)

  ierr = 0

  yvec(1:neq) => FN_VGetArrayPointer(sunvec_y)
  djac(1:lda,1:neq) => FSUNBandMatrix_Data(sunmat_J)

  call jacobian(yvec, djac, neq, lda, global_err)
  if (len_trim(global_err) /= 0) then
    if (trim(global_err) == "Mixing ratios sum to > 1.0. "//&
                          "Atmosphere is probably in a run-away state.") then
      ierr = -1
    else
      print*,trim(global_err)
      ierr = 1
    endif
  endif
  return
end function
