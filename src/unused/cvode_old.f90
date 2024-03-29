
subroutine cvode99(T0, usol_start, nnq, nnz, t_eval, num_t_eval, rtol, atol, &
                 use_fast_jacobian, solution, success, err)
  use photochem_data, only: neq, nq, nz, jtrop, np, mass, background_mu, rainout_on
  use photochem_vars, only: verbose, lbound, fixedmr, T, den, P, Press, max_cvode_steps, &
                            rpar_init
  use photochem_wrk, only: cvode_stepper, rain, raingc, global_err, rpar
  
  
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
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian
  character(len=1000), intent(out) :: err
  character(len=10) :: message

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
  real(8) :: mubar_z(nz)

  cvode_stepper = 0
  err = ''

  nneq = neq ! number of equations
  METH = 2 ! CVODE BDF method
  IATOL = 1 ! scalar absolute tolerance
  mu = nq ! upper diagonal
  ml = nq ! lower diagonal
  ITASK = 1 ! for output

  ! begin stuff that needs to be inizialized
  if (np.gt.0) then
    rpar = rpar_init
  endif
  do i = 1,nz
    call mean_molecular_weight(nq, usol_start(:,i), mass, background_mu, mubar_z(i))
  enddo
  call densty(nz, mubar_z, T, den, P, press) 
  if (rainout_on) then
    call rainout(.true.,Jtrop,usol_start,nq,nz, T,den, rain, raingc, err)
    if (len_trim(err) /= 0) return
  else
    rain = 0.d0
    raingc = 0.d0
  endif
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
    err = 'SUNDIALS_ERROR: FNVINITS'
    return
  endif

  ! initialize banded matrix module
  call FSUNBANDMATINIT(1, nneq, mu, ml, mu+ml, ier)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FSUNBANDMATINIT'
    return
  endif

  ! initialize banded linear solver module
  call FSUNBANDLINSOLINIT(1, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FSUNBANDLINSOLINIT'
    return
  endif

  call FCVMALLOC(T0, U, METH, IATOL, RTOL, ATOL, &
                 IOUT, ROUT, IPAR, RRPAR, IER)
  if (ier .ne. 0) then
   err = 'SUNDIALS_ERROR: FCVMALLOC'
   return
  endif

  ! attach matrix and linear solver modules to CVLs interface
  CALL FCVLSINIT(IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVLSINIT'
    CALL FCVFREE
    return
  endif

  CALL FCVSETIIN('MAX_NSTEPS', 1000000, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    return
  endif

  CALL FCVSETIIN('MAX_ERRFAIL', 15, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    return
  endif

  ! start with tiny step
  call FCVSETRIN('INIT_STEP',1.d-6,ier)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVSETRIN'
    CALL FCVFREE
    return
  endif

  if (use_fast_jacobian) then
    CALL FCVBANDSETJAC(1, IER)
    if (ier .ne. 0) then
      err = 'SUNDIALS_ERROR: FCVBANDSETJAC'
      CALL FCVFREE
      return
    endif
  endif

  do ii=1,num_t_eval

    CALL FCVODE(t_eval(ii), TT, U, ITASK, IER)
    if (ier < 0) then ! then there is a problem
      if ((global_err == 'max steps') .and. (ier == -8)) then ! reached max steps
        write(message,'(i10)')  max_cvode_steps 
        print*,'CVODE stopped because it reached the maximum number of of specified steps: '//trim(message)
        CALL FCVFREE
        return
      else if ((global_err /= 'max steps') .and. (ier == -8)) then ! err in my rhs
        err = global_err
        CALL FCVFREE
        return
      else
        write(message,'(i4)') ier
        err = 'SUNDIALS_ERROR: FCVODE '//trim(message)
        CALL FCVFREE
        return
      endif
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

  if (ier >= 0) then
    success = .true.
  else
    success = .false.
  endif
  call fcvfree

end subroutine

subroutine cvode_save123(T0, usol_start, nnq, nnz, t_eval, num_t_eval, rtol, atol, &
                      use_fast_jacobian, outfilename, amount2save, success, err)
  use photochem_data, only: neq, nq, nz, jtrop, ispec, np, nr, background_spec, &
                            nw, wavl, z, jchem, kj, ks, photoreac, photonums, photospec, &
                            jchem, nmax, iprod, iloss,nump,numl, nsp, background_mu, mass, rainout_on, &
                            Flux
  use photochem_vars, only: verbose, lbound, fixedmr, T, den, P, Press, max_cvode_steps, &
                            rpar_init, edd
  use photochem_wrk, only: cvode_stepper, rain, raingc, global_err, rpar, surf_radiance, A, yp, yl, D
  
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
  integer, intent(in) :: amount2save
  

  ! output
  logical, intent(out) :: success
  character(len=1000), intent(out) :: err
  character(len=10) :: message

  ! other
  double precision, dimension(nnq,nnz) :: solution_temp
  integer*4 LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE ! for printing stats
  data LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, &
      LNJE/17/
  integer*4 :: IER, METH, IATOL, ITASK ! some options
  integer i, j, k, ii

  ! The following declaration specification should match C type long int.
  integer*8 nneq, IOUT(25), IPAR(2), mu, ml
  double precision TT
  double precision U(neq), ROUT(10), RRPAR(1)
  real(8), allocatable :: UDOT(:), loss(:,:)
  real(8), allocatable :: photorates(:,:)
  integer :: ind(1)
  real(8) :: mubar_z(nz)
  
  allocate(UDOT(neq))
  allocate(loss(nq,nz))
  allocate(photorates(ks,nz))
  
  err = ''

  cvode_stepper = 0

  nneq = neq ! number of equations
  METH = 2 ! CVODE BDF method
  IATOL = 1 ! scalar absolute tolerance
  mu = nq ! upper diagonal
  ml = nq ! lower diagonal
  ITASK = 1 ! for output
  
  ! begin stuff that needs to be inizialized
  if (np.gt.0) then
    rpar = rpar_init
  endif
  do i = 1,nz
    call mean_molecular_weight(nq, usol_start(:,i), mass, background_mu, mubar_z(i))
  enddo
  call densty(nz, mubar_z, T, den, P, press) 
  if (rainout_on) then
    call rainout(.true.,Jtrop,usol_start,nq,nz, T,den, rain, raingc, err)
    if (len_trim(err) /= 0) return
  else
    rain = 0.d0
    raingc = 0.d0
  endif
  ! end stuff that needs to be inizialized

  ! initial conditions
  DO I=1,NQ
    DO J=1,NZ
      K = I + (J-1)*NQ
      U(K) = usol_start(i,j)
    enddo
  enddo
  call right_hand_side(U,UDOT,neq,err)
  
  
  ! file prep
  open(2,file=outfilename,status='replace',form="unformatted")
  write(2) nq 
  write(2) np
  write(2) nw
  write(2) kj
  write(2) ks
  write(2) nz
  write(2) background_spec ! len=8
  write(2) ispec
  write(2) photospec
  
  ! new
  write(2) nsp
  write(2) nr
  write(2) nmax
  write(2) jchem
  write(2) iprod
  write(2) iloss
  write(2) nump
  write(2) numl
  write(2) A
  ! new
  
  write(2) z
  write(2) T
  write(2) edd
  write(2) rpar_init
  write(2) wavl
  write(2) flux
  
  write(2) num_t_eval
  write(2) amount2save
  close(2)

  CALL FNVINITS(1, nneq, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FNVINITS'
    return
  endif

  ! initialize banded matrix module
  call FSUNBANDMATINIT(1, nneq, mu, ml, mu+ml, ier)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FSUNBANDMATINIT'
    return
  endif

  ! initialize banded linear solver module
  call FSUNBANDLINSOLINIT(1, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FSUNBANDLINSOLINIT'
    return
  endif

  call FCVMALLOC(T0, U, METH, IATOL, RTOL, ATOL, &
                 IOUT, ROUT, IPAR, RRPAR, IER)
  if (ier .ne. 0) then
   err = 'SUNDIALS_ERROR: FCVMALLOC'
   return
  endif

  ! attach matrix and linear solver modules to CVLs interface
  CALL FCVLSINIT(IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVLSINIT'
    CALL FCVFREE
    return
  endif

  CALL FCVSETIIN('MAX_NSTEPS', 1000000, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    return
  endif

  CALL FCVSETIIN('MAX_ERRFAIL', 15, IER)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVSETIIN'
    CALL FCVFREE
    return
  endif

  ! start with tiny step
  call FCVSETRIN('INIT_STEP',1.d-6,ier)
  if (ier .ne. 0) then
    err = 'SUNDIALS_ERROR: FCVSETRIN'
    CALL FCVFREE
    return
  endif

  if (use_fast_jacobian) then
    CALL FCVBANDSETJAC(1, IER)
    if (ier .ne. 0) then
      err = 'SUNDIALS_ERROR: FCVBANDSETJAC'
      CALL FCVFREE
      return
    endif
  endif

  do ii=1,num_t_eval

    CALL FCVODE(t_eval(ii), TT, U, ITASK, IER)
    if (ier < 0) then ! then there is a problem
      if ((global_err == 'max steps') .and. (ier == -8)) then ! reached max steps
        write(message,'(i10)')  max_cvode_steps 
        print*,'CVODE stopped because it reached the maximum number of of specified steps: '//trim(message)
        CALL FCVFREE
        return
      else if ((global_err /= 'max steps') .and. (ier == -8)) then ! err in my rhs
        err = global_err
        CALL FCVFREE
        return
      else
        write(message,'(i4)') ier
        err = 'SUNDIALS_ERROR: FCVODE '//trim(message)
        CALL FCVFREE
        return
      endif
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
    
    call right_hand_side(U,UDOT,neq,err)
    photorates = 0.d0
    do j=1,kj
      i = photoreac(j)   
      ind = findloc_integer(size(photospec),photospec,i)
      photorates(ind(1),:) = photorates(ind(1),:) + A(photonums(j),:)*solution_temp(i,:)*den
    enddo
    do j = 1,nq
      loss(j,:) = yl(j,:)*solution_temp(j,:)*den
    enddo
    
    
    open(2,file=outfilename,status='old',form="unformatted", position="append")
    write(2) 999
    write(2) t_eval(ii)
    ! write(2) solution_temp
    write(2) D
    write(2) den
    if (amount2save == 1) then
      write(2) P
      write(2) surf_radiance
      write(2) photorates
      write(2) yp
      write(2) loss
    endif
    close(2)

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

  if (ier >= 0) then
    success = .true.
  else
    success = .false.
  endif
  call fcvfree
  deallocate(UDOT)
  deallocate(loss)
  deallocate(photorates)

end subroutine


subroutine cvode_equilibrium(rtol, atol, use_fast_jacobian, success, err)
  use photochem_data, only: nr, nq, np, nz1, nq1, isl, nz, jtrop, &
                            lh2o, dz, nsp2, jtrop, background_mu, mass, rainout_on, &
                            fix_water_in_troposphere
  use photochem_vars, only: verbose, usol_init, usol_out, &
                            den, fluxo, flow, T, P, den, press, rpar_init, &
                            equilibrium_time
  use photochem_wrk, only: cvode_stepper, wfall, dk, scale_h, h_atm, yp, yl, &
                           bx1x2, raingc, A, rain, rpar
  implicit none
  ! input
  double precision, intent(in) :: rtol ! relative tolerance
  double precision, intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian

  ! output
  logical, intent(out) :: success
  character(len=1000), intent(out) :: err

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
  real(8) :: mubar_z(nz)
  err = ''
  cvode_stepper = 0
  ! begin stuff that needs to be inizialized
  if (np.gt.0) then
    rpar = rpar_init
  endif
  do i = 1,nz
    call mean_molecular_weight(nq, usol_init(:,i), mass, background_mu, mubar_z(i))
  enddo
  call densty(nz, mubar_z, T, den, P, press) 
  if (rainout_on) then
    call rainout(.true.,Jtrop,usol_init,nq,nz, T,den, rain, raingc, err)
    if (len_trim(err) /= 0) return
  else
    rain = 0.d0
    raingc = 0.d0
  endif
  ! end stuff that needs to be inizialized

  t_eval(1) = equilibrium_time
  call system_clock(count = c1, count_rate = cr, count_max = cm)
  call cvode(T0, usol_init, nq, nz, t_eval, num_t_eval, rtol, atol, &
             use_fast_jacobian, solution, success, err)
  if (len_trim(err) /= 0) return
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
    
    do i=1,NZ-1
      do j=1,(NQ-NP)
        fluxo(j,i) = fluxo(j,i) &
                    - bX1X2(j,i)*(usol(j,i+1) - usol(j,i))/dz(i) &
                    + bX1X2(j,i)*0.5*(usol(j,i) + usol(j,i+1)) &
                    *(1./H_atm(i) - 1./scale_H(j,i))
      enddo
    enddo

    DO K=1,NQ
      FLOW(K) = FLUXO(K,1) - (YP(K,1) - YL(K,1)*D(K,1))*DZ(1)
      FUP(K) = FLUXO(K,NZ1) + (YP(K,NZ) - YL(K,NZ)*D(K,NZ))*DZ(NZ)
    enddo
    if (fix_water_in_troposphere) then
      FLOW(LH2O) = FLUXO(LH2O,jtrop)
    endif

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
