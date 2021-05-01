SUBROUTINE FCVFUN(TT, U, UDOT, IPAR, RRPAR, IER)
  use photochem, only: right_hand_side
  use photochem_data, only: neq
  use photochem_vars, only: max_cvode_steps, verbose
  use photochem_wrk, only: cvode_stepper, time_prev, global_err
  implicit none
  INTEGER*8 IPAR(*)
  INTEGER*4 IER
  DOUBLE PRECISION TT, U(*), UDOT(neq), RRPAR(*)
  ier = 0
  
  call right_hand_side(U,UDOT,neq,global_err)
  if (len_trim(global_err) /= 0) then
    ier = -1
  endif
  
  if ((TT .ne. time_prev) .and. (verbose)) then
    print"(1x,'N =',i6,3x,'Time = ',es20.14,3x,'max(df/dt) = ',es10.3)",cvode_stepper,TT,maxval(UDOT)
    cvode_stepper = cvode_stepper + 1
  endif
  
  if (cvode_stepper >= max_cvode_steps) then
    global_err = 'max steps'
    ier = -1
  endif
  time_prev = TT
end subroutine

SUBROUTINE FCVBJAC(N, MU, ML, MDIM, TT, U, FU, &
                  BJAC, HH, IPAR, RRPAR, V1, V2, V3, IER)
  use photochem, only: jacobian
  use photochem_data, only: neq, nq
  use photochem_wrk, only: global_err
  implicit none
  ! The following declaration specification should match C type long int.
  INTEGER*8 N, MU, ML, MDIM, IPAR(*)
  INTEGER*4 IER
  DOUBLE PRECISION TT, U(*), FU(*), BJAC(MDIM,neq), HH, RRPAR(*)
  DOUBLE PRECISION V1(*), V2(*), V3(*)
  ier = 0
  call jacobian(U,BJAC(1:nq+nq+1,:),neq,nq+nq+1,global_err)
  if (len_trim(global_err) /= 0) then
    ier = -1
  endif
end subroutine
