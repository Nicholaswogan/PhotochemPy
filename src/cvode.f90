
subroutine cvode(T0,usol_start,nnq,nnz,t_eval,num_t_eval,rtol,atol,use_fast_jacobian,outfilename,solution)
  implicit none

  ! inputs
  double precision, intent(in) :: T0 ! intial time
  double precision, dimension(nnq,nnz) :: usol_start ! initial atmosphere
  integer, intent(in) :: num_t_eval
  integer, intent(in) :: nnq
  integer, intent(in) :: nnz
  double precision, dimension(num_t_eval), intent(in) :: t_eval ! times to evaluate solution
  double precision, intent(in) :: rtol ! relative tolerance
  double precision, intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian (keep false!)
  character(len=*), intent(in) :: outfilename ! were to save the solution

  double precision, dimension(num_t_eval,nnq,nnz),intent(out) :: solution

  character(len=25) :: tmp0, tmp1 ! for io

  integer*4 LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE ! for printing stats
  data LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, &
      LNJE/17/
  integer*4 :: IER, METH, IATOL, ITASK ! some options
  integer i, j, k, ii

  ! The following declaration specification should match C type long int.
  integer*8 nneq, IOUT(25), IPAR(2), mu, ml
  double precision TT, TOUT
  double precision U(neq), ROUT(10), RRPAR(1)

  nneq = neq
  METH = 2 ! CVODE BDF method
  IATOL = 1 ! scalar absolute tolerance
  mu = nq ! upper diagonal
  ml = nq ! lower diagonal
  ITASK = 1 ! for output

  if (outfilename .ne. "None") then
    ! file prep
    open(1,file=outfilename,status='replace')
    write(1,'(a11)',advance='no') ispec(1)
    do i=2,nq
      write(1,'(a25)',advance='no') ispec(i)
    enddo
    close(1)
  endif

  ! so rainout works
  call rainout(Jtrop,0,usol_start,nq,nz)

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

  ! no bound on number of steps
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

  ! keep this false! Let CVODE calculate the jacobian
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
      CALL FCVFREE
      stop
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

    if (outfilename .ne. "None") then
      ! save to file
      open(1,file=outfilename,status='old',position="append")
      write(1,*)
      write(tmp0,"(es25.15)") t_eval(ii)
      write(tmp1,"(i25)") ii
      write(1,"(a, 25a, a, 25a)") 't = ',adjustl(tmp0), 'n = ', adjustl(tmp1)
      do i = 1,nz
        do j = 1,nq ! write each line
          write(1,'(es25.15e3)',advance='no') solution(ii,j,i)
        enddo
        write(1,*) ! new line
      enddo
      close(1)
    endif

  enddo

  if (verbose) then
    WRITE(6,80) IOUT(LNST), IOUT(LNFE), IOUT(LNJE), IOUT(LNSETUP), &
               IOUT(LNNI), IOUT(LNCF), IOUT(LNETF)
80  FORMAT(//'Final statistics:'// &
          ' No. steps = ', I4, '  No. f-s = ', I4, &
          '  No. J-s = ', I4, '   No. LU-s = ', I4/ &
          ' No. nonlinear iterations = ', I4/ &
          ' No. nonlinear convergence failures = ', I4/ &
          ' No. error test failures = ', I4)
  endif

end subroutine
