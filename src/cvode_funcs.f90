

integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                              result(ierr) bind(c, name='RhsFn')
  use, intrinsic :: iso_c_binding
  use fcvode_mod
  use fsundials_nvector_mod
  use photochem_data, only: neq
  use photochem_vars, only: max_cvode_steps, verbose
  use photochem_wrk, only: cvode_stepper, time_prev, global_err
  implicit none
  ! calling variables
  real(c_double), value :: tn        ! current time
  type(N_Vector)        :: sunvec_y  ! solution N_Vector
  type(N_Vector)        :: sunvec_f  ! rhs N_Vector
  type(c_ptr), value    :: user_data ! user-defined dat
  
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
    print*,trim(global_err)
    ierr = -1
  endif
  
  if ((tn .ne. time_prev) .and. (verbose)) then
    print"(1x,'N =',i6,3x,'Time = ',es20.14,3x,'max(df/dt) = ',es10.3)",cvode_stepper,tn,maxval(abs(fvec))
    cvode_stepper = cvode_stepper + 1
  endif
  
  if (cvode_stepper >= max_cvode_steps) then
    global_err = 'max steps'
    ierr = -1
  endif
  time_prev = tn
  
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
    ierr = -1
  endif
  return
end function