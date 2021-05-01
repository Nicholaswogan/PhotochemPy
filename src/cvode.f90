
subroutine cvode(T0, usol_start, nnq, nnz, t_eval, num_t_eval, rtol, atol, &
                 use_fast_jacobian, solution, success)
  use photochem_data, only: neq, nq, nz, jtrop
  use photochem_vars, only: verbose, lbound, fixedmr, T, den, P, Press
  use photochem_wrk, only: cvode_stepper, rain, raingc
  implicit none

  ! inputs
  double precision, intent(in) :: T0 ! intial time
  double precision, dimension(nnq,nnz) :: usol_start ! initial atmosphere
  integer, intent(in) :: nnq
  integer, intent(in) :: nnz
  double precision, dimension(num_t_eval), intent(in) :: t_eval ! times to evaluate solution
  integer, intent(in) :: num_t_eval ! dimension of t_eval
  double precision, intent(in) :: rtol ! relative tolerance
  double precision, intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian (keep false!)

  ! output
  double precision, dimension(num_t_eval,nnq,nnz),intent(out) :: solution
  logical, intent(out) :: success

  ! other
  integer*4 LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE ! for printing stats
  data LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, &
      LNJE/17/
  integer*4 :: IER, METH, IATOL, ITASK ! some options
  integer i, j, k, ii

  ! The following declaration specification should match C type long int.
  integer*8 nneq, IOUT(25), IPAR(2), mu, ml
  double precision TT
  double precision U(neq), ROUT(10), RRPAR(1)

  cvode_stepper = 0

  nneq = neq ! number of equations
  METH = 2 ! CVODE BDF method
  IATOL = 1 ! scalar absolute tolerance
  mu = nq ! upper diagonal
  ml = nq ! lower diagonal
  ITASK = 1 ! for output

  ! begin stuff that needs to be inizialized
  call densty(nq, nz, usol_start, T, den, P, press) 
  call rainout(.true.,Jtrop,usol_start,nq,nz, T,den, rain, raingc) 
  ! end stuff that needs to be inizialized

  ! initial conditions
  DO I=1,NQ
    DO J=1,NZ
      K = I + (J-1)*NQ
      U(K) = usol_start(i,j)
    enddo
  enddo

  CALL FNVINITS(1, nneq, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FNVINITS'
    stop
  endif

  ! initialize banded matrix module
  call FSUNBANDMATINIT(1, nneq, mu, ml, mu+ml, ier)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FSUNBANDMATINIT'
    stop
  endif

  ! initialize banded linear solver module
  call FSUNBANDLINSOLINIT(1, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FSUNBANDLINSOLINIT'
    stop
  endif

  call FCVMALLOC(T0, U, METH, IATOL, RTOL, ATOL, &
                 IOUT, ROUT, IPAR, RRPAR, IER)
  if (ier .ne. 0) then
   print*,'SUNDIALS_ERROR: FCVMALLOC'
   stop
  endif

  ! attach matrix and linear solver modules to CVLs interface
  CALL FCVLSINIT(IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVLSINIT'
    CALL FCVFREE
    stop
  endif

  CALL FCVSETIIN('MAX_NSTEPS', 1000000, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    stop
  endif

  CALL FCVSETIIN('MAX_ERRFAIL', 15, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    stop
  endif

  ! start with tiny step
  call FCVSETRIN('INIT_STEP',1.d-6,ier)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVSETRIN'
    CALL FCVFREE
    stop
  endif

  if (use_fast_jacobian) then
    CALL FCVBANDSETJAC(1, IER)
    if (ier .ne. 0) then
      print*,'SUNDIALS_ERROR: FCVBANDSETJAC'
      CALL FCVFREE
      stop
    endif
  endif

  do ii=1,num_t_eval

    CALL FCVODE(t_eval(ii), TT, U, ITASK, IER)
    if (ier .ne. 0) then
      print*,'SUNDIALS_ERROR: FCVODE ', ier
      ! CALL FCVFREE
      ! stop
    endif

    ! save the solution
    do I=1,NQ
      do J=1,NZ
        K = I + (J-1)*NQ
        solution(ii,i,j) = U(k)
      enddo
    enddo

    ! fix lbound = 1
    do i=1,nq
      if (lbound(i) .eq.1) then
        solution(ii,i,1) = fixedmr(i)
      endif
    enddo

  enddo

  if (verbose) then
    WRITE(6,80) IOUT(LNFE), IOUT(LNJE), IOUT(LNSETUP), &
               IOUT(LNNI), IOUT(LNCF), IOUT(LNETF)
80  FORMAT(//'Final statistics:'// &
          ' No. f-s = ', I4, &
          '  No. J-s = ', I4, '   No. LU-s = ', I4/ &
          ' No. nonlinear iterations = ', I4/ &
          ' No. nonlinear convergence failures = ', I4/ &
          ' No. error test failures = ', I4)
  endif

  if (ier == 0) then
    success = .true.
  else
    success = .false.
  endif
  call fcvfree

end subroutine

subroutine cvode_save(T0, usol_start, nnq, nnz, t_eval, num_t_eval, rtol, atol, &
                      use_fast_jacobian, outfilename, success)
  use photochem_data, only: neq, nq, nz, jtrop, ispec
  use photochem_vars, only: verbose, lbound, fixedmr, T, den, P, Press
  use photochem_wrk, only: cvode_stepper, rain, raingc
  implicit none

  ! inputs
  double precision, intent(in) :: T0 ! intial time
  double precision, dimension(nnq,nnz) :: usol_start ! initial atmosphere
  integer, intent(in) :: nnq
  integer, intent(in) :: nnz
  double precision, dimension(num_t_eval), intent(in) :: t_eval ! times to evaluate solution
  integer, intent(in) :: num_t_eval ! dimension of t_eval
  double precision, intent(in) :: rtol ! relative tolerance
  double precision, intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian (keep false!)
  character(len=*), intent(in) :: outfilename ! were to save the solution

  ! output
  logical, intent(out) :: success

  ! other
  double precision, dimension(nnq,nnz) :: solution_temp
  character(len=25) :: tmp0, tmp1 ! for io
  integer*4 LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE ! for printing stats
  data LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, &
      LNJE/17/
  integer*4 :: IER, METH, IATOL, ITASK ! some options
  integer i, j, k, ii

  ! The following declaration specification should match C type long int.
  integer*8 nneq, IOUT(25), IPAR(2), mu, ml
  double precision TT
  double precision U(neq), ROUT(10), RRPAR(1)

  cvode_stepper = 0

  nneq = neq ! number of equations
  METH = 2 ! CVODE BDF method
  IATOL = 1 ! scalar absolute tolerance
  mu = nq ! upper diagonal
  ml = nq ! lower diagonal
  ITASK = 1 ! for output

  ! file prep
  open(1,file=outfilename,status='replace')
  write(1,'(a11)',advance='no') ispec(1)
  do i=2,nq
    write(1,'(a25)',advance='no') ispec(i)
  enddo
  close(1)

  ! begin stuff that needs to be inizialized
  call densty(nq, nz, usol_start, T, den, P, press) 
  call rainout(.true.,Jtrop,usol_start,nq,nz, T,den, rain, raingc) 
  ! end stuff that needs to be inizialized

  ! initial conditions
  DO I=1,NQ
    DO J=1,NZ
      K = I + (J-1)*NQ
      U(K) = usol_start(i,j)
    enddo
  enddo

  CALL FNVINITS(1, nneq, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FNVINITS'
    stop
  endif

  ! initialize banded matrix module
  call FSUNBANDMATINIT(1, nneq, mu, ml, mu+ml, ier)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FSUNBANDMATINIT'
    stop
  endif

  ! initialize banded linear solver module
  call FSUNBANDLINSOLINIT(1, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FSUNBANDLINSOLINIT'
    stop
  endif

  call FCVMALLOC(T0, U, METH, IATOL, RTOL, ATOL, &
                 IOUT, ROUT, IPAR, RRPAR, IER)
  if (ier .ne. 0) then
   print*,'SUNDIALS_ERROR: FCVMALLOC'
   stop
  endif

  ! attach matrix and linear solver modules to CVLs interface
  CALL FCVLSINIT(IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVLSINIT'
    CALL FCVFREE
    stop
  endif

  CALL FCVSETIIN('MAX_NSTEPS', 1000000, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    stop
  endif

  CALL FCVSETIIN('MAX_ERRFAIL', 15, IER)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    stop
  endif

  ! start with tiny step
  call FCVSETRIN('INIT_STEP',1.d-6,ier)
  if (ier .ne. 0) then
    print*,'SUNDIALS_ERROR: FCVSETRIN'
    CALL FCVFREE
    stop
  endif

  if (use_fast_jacobian) then
    CALL FCVBANDSETJAC(1, IER)
    if (ier .ne. 0) then
      print*,'SUNDIALS_ERROR: FCVBANDSETJAC'
      CALL FCVFREE
      stop
    endif
  endif

  do ii=1,num_t_eval

    CALL FCVODE(t_eval(ii), TT, U, ITASK, IER)
    if (ier .ne. 0) then
      print*,'SUNDIALS_ERROR: FCVODE ', ier
      ! CALL FCVFREE
      ! stop
    endif

    ! save the solution
    do I=1,NQ
      do J=1,NZ
        K = I + (J-1)*NQ
        solution_temp(i,j) = U(k)
      enddo
    enddo

    ! fix lbound = 1
    do i=1,nq
      if (lbound(i) .eq.1) then
        solution_temp(i,1) = fixedmr(i)
      endif
    enddo

    ! save to file
    open(1,file=outfilename,status='old',position="append")
    write(1,*)
    write(tmp0,"(es25.15)") t_eval(ii)
    write(tmp1,"(i25)") ii
    write(1,"(a, 25a, a, 25a)") 't = ',adjustl(tmp0), 'n = ', adjustl(tmp1)
    do i = 1,nz
      do j = 1,nq ! write each line
        write(1,'(es25.15e3)',advance='no') solution_temp(j,i)
      enddo
      write(1,*) ! new line
    enddo
    close(1)

  enddo

  if (verbose) then
    WRITE(6,80) IOUT(LNFE), IOUT(LNJE), IOUT(LNSETUP), &
               IOUT(LNNI), IOUT(LNCF), IOUT(LNETF)
80  FORMAT(//'Final statistics:'// &
          ' No. f-s = ', I4, &
          '  No. J-s = ', I4, '   No. LU-s = ', I4/ &
          ' No. nonlinear iterations = ', I4/ &
          ' No. nonlinear convergence failures = ', I4/ &
          ' No. error test failures = ', I4)
  endif

  if (ier == 0) then
    success = .true.
  else
    success = .false.
  endif
  call fcvfree

end subroutine


subroutine cvode_equilibrium(rtol, atol, use_fast_jacobian, success)
  use photochem_data, only: nr, nq, np, nz1, nq1, isl, nz, jtrop, &
                            lh, lh2, lh2o, dz, nsp2, background_spec, jtrop
  use photochem_vars, only: verbose, usol_init, usol_out, &
                            den, fluxo, flow, T, P, den, press
  use photochem_wrk, only: cvode_stepper, wfall, dk, scale_h, h_atm, yp, yl, &
                           bh2n2, bhn2, raingc, A, rain
  implicit none
  ! input
  double precision, intent(in) :: rtol ! relative tolerance
  double precision, intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian

  ! output
  logical, intent(out) :: success

  ! other
  double precision :: T0 = 0.d0 ! intial time
  integer :: num_t_eval = 1
  double precision, dimension(1) :: t_eval
  double precision, dimension(1,nq,nz) :: solution
  double precision, dimension(nq,nz) :: usol
  integer cr, cm, c1, c2
  integer i,j,k,ll
  real*8,dimension(nq1) :: SR, FUP
  double precision :: wnf
  real*8, dimension(nq,nz) :: Fval
  real*8, dimension(nsp2,nz) :: D

  cvode_stepper = 0
  ! begin stuff that needs to be inizialized
  call densty(nq, nz, usol_init, T, den, P, press) 
  call rainout(.true.,Jtrop,usol_init,nq,nz, T,den, rain, raingc) 
  ! end stuff that needs to be inizialized

  t_eval(1) = 1.d17
  call system_clock(count = c1, count_rate = cr, count_max = cm)
  call cvode(T0, usol_init, nq, nz, t_eval, num_t_eval, rtol, atol, &
             use_fast_jacobian, solution, success)
  call system_clock(count = c2)
  usol = solution(1,:,:)

  if (success) then
    if (verbose) then
      print"('Time to find equilibrium =',f10.3,' seconds')", &
            (c2-c1)/real(cr)
    endif
  ! the output
    do i=1,nq
      do j=1,nz
        usol_out(i,j) = usol(i,j)
      enddo
    enddo

    ! the flux lots of stuff below
    call dochem(0, nr, nsp2, nq, nz, usol, A, isl, jtrop, D, fval)
    DO K=1,NQ
      DO I=1,NZ
      D(K,I) = USOL(K,I)*DEN(I)
      enddo
      DO I=1,NZ-1
        FLUXO(K,I) = - DK(I)*(USOL(K,I+1) - USOL(K,I))/DZ(I)
      enddo
    enddo

    if(NP.GT.0) THEN
      do k=1,np
        do I=1,NZ1
          LL=NQ-NP+K
          WNF = 0.5*(WFALL(I,K)*DEN(I)*USOL(LL,I) + WFALL(I+1,K)*DEN(I+1) &
            *USOL(LL,I+1))
          FLUXO(LL,I) = FLUXO(LL,I) - WNF
        enddo
      enddo
    endif

    if (background_spec /= 'H2') then
      do i=1,NZ-1
        fluxo(LH,i) = fluxo(LH,i) &
          - bHN2(i)*(usol(LH,i+1) - usol(LH,i))/dz(i) &
          + bHN2(i)*0.5*(usol(LH,i) + usol(LH,i+1)) &
          *(1./H_atm(i) - 1./scale_H(LH,i))
        fluxo(LH2,i) =fluxo(LH2,i) &
          - bH2N2(i)*(usol(LH2,i+1) - usol(LH2,i))/dz(i) &
          + bH2N2(i)*0.5*(usol(LH2,i) + usol(LH2,i+1)) &
          *(1./H_atm(i) - 1./scale_H(LH2,i))
      enddo
    endif

    DO K=1,NQ
      FLOW(K) = FLUXO(K,1) - (YP(K,1) - YL(K,1)*D(K,1))*DZ(1)
      FUP(K) = FLUXO(K,NZ1) + (YP(K,NZ) - YL(K,NZ)*D(K,NZ))*DZ(NZ)
    enddo
    FLOW(LH2O) = FLUXO(LH2O,jtrop)

    DO I=1,NQ
      SR(I) = 0.
      DO J=1,JTROP
        SR(I) = SR(I) + RAINGC(I,J)*USOL(I,J)*DEN(J)*DZ(J)
      enddo
    enddo

    ! redox state
    call redox_conservation(FLOW,FUP,SR)
  endif ! end if converged

end subroutine
