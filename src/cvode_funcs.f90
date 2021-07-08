SUBROUTINE FCVFUN(TT, U, UDOT, IPAR, RRPAR, IER)
  use photochem, only: right_hand_side
  use photochem_data, only: neq, z, ispec, nq
  use photochem_vars, only: max_cvode_steps, verbose
  use photochem_wrk, only: cvode_stepper, time_prev, global_err
  implicit none
  INTEGER*8 IPAR(*)
  INTEGER*4 IER
  DOUBLE PRECISION TT, U(*), UDOT(neq), RRPAR(*)
  
  real(8) :: tmp, mx
  integer :: i, k, j
  ier = 0
  
  call right_hand_side(U,UDOT,neq,global_err)
  if (len_trim(global_err) /= 0) then
    ier = -1
  endif
  
  if ((TT .ne. time_prev) .and. (verbose)) then
    
    tmp = 0.d0
    mx = tmp
    k = 1
    do i = 1,neq
      tmp = abs(udot(i)/u(i))
      if (tmp > mx .and. U(i) > 1.d-30) then
        mx = tmp
        k = i
      endif
    enddo
    j = k/nq
    i = k-j*nq
    
    ! print"(1x,'N =',i6,3x,'Time = ',es25.10,3x,a10,3x,es12.5,3x,es12.5,3x,es12.5)", &
            ! cvode_stepper, tt, trim(ispec(i)),Udot(k),U(k),z(j+1)/1.d5
    
    print"(1x,'N =',i6,3x,'Time = ',es20.14,3x,'max(df/dt) = ',es10.3)",cvode_stepper,TT,maxval(abs(UDOT))
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
