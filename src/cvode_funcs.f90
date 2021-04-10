SUBROUTINE FCVFUN(TT, U, UDOT, IPAR, RRPAR, IER)
  use Photochem
  implicit none
  INTEGER*8 IPAR(*)
  INTEGER*4 IER
  DOUBLE PRECISION TT, U(*), UDOT(neq), RRPAR(*)
  integer i
  ier = 0
  call right_hand_side(U,UDOT,neq)
  if ((TT .ne. time_prev) .and. (verbose)) then
    print"(1x,'N =',i6,3x,'Time = ',es20.14,3x,'max(df/dt) = ',es10.3)",cvode_stepper,TT,maxval(UDOT)
    if (cvode_stepper >= max_cvode_steps) ier = -1
    cvode_stepper = cvode_stepper + 1
  endif
  time_prev = TT
end subroutine

SUBROUTINE FCVBJAC(N, MU, ML, MDIM, TT, U, FU, &
                  BJAC, HH, IPAR, RRPAR, V1, V2, V3, IER)
  use Photochem
  implicit none
  ! The following declaration specification should match C type long int.
  INTEGER*8 N, MU, ML, MDIM, IPAR(*)
  INTEGER*4 IER
  DOUBLE PRECISION TT, U(*), FU(*), BJAC(MDIM,neq), HH, RRPAR(*)
  DOUBLE PRECISION V1(*), V2(*), V3(*)
  integer i
  BJAC = 0.d0
  call jacobian(U,BJAC(1:nq+nq+1,:),neq,nq+nq+1)
  ier = 0
end subroutine
