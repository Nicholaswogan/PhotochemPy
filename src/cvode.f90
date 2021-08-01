
subroutine cvode(t0, usol_start, nnq, nnz, t_eval, num_t_eval, rtol, atol, &
                     use_fast_jacobian, solution, success, err)
  use photochem_data, only: neq, nq, nz, jtrop, np, mass, background_mu, rainout_on
  use photochem_vars, only: lbound, fixedmr, T, den, P, Press, &
                            rpar_init, &
                            max_cvode_steps, initial_dt, max_err_test_failures, max_order
  use photochem_wrk, only: cvode_stepper, rain, raingc, global_err, rpar
  
  use, intrinsic :: iso_c_binding
  use fcvode_mod, only: CV_BDF, CV_NORMAL, FCVodeInit, FCVodeSStolerances, &
                        FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                        FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                        FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd
  use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
  use fnvector_serial_mod, only: FN_VMake_Serial   
  use fsunmatrix_band_mod, only: FSUNBandMatrix
  use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
  use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
  use fsunlinsol_band_mod, only: FSUNLinSol_Band
  
  implicit none
  
  real(8), intent(in) :: t0 ! initial time
  real(8), dimension(nnq,nnz) :: usol_start ! initial atmosphere
  integer, intent(in) :: nnq
  integer, intent(in) :: nnz
  real(8), dimension(num_t_eval), intent(in) :: t_eval ! times to evaluate solution
  integer, intent(in) :: num_t_eval ! dimension of t_eval
  real(8), intent(in) :: rtol ! relative tolerance
  real(8), intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian
  
  real(8), dimension(num_t_eval,nnq,nnz),intent(out) :: solution
  logical, intent(out) :: success
  character(len=1000), intent(out) :: err
  
  ! local
  real(c_double) :: tcur(1)    ! current time
  integer(c_int) :: ierr       ! error flag from C functions
  type(c_ptr)    :: cvode_mem  ! CVODE memory
  type(N_Vector), pointer :: sunvec_y ! sundials vector
  real(c_double) :: yvec(neq)
  integer(c_long) :: neq_long
  integer(c_long) :: mu, ml
  integer(c_long) :: mxsteps_ ! unused
  type(SUNMatrix), pointer :: sunmat
  type(SUNLinearSolver), pointer :: sunlin
  
  real(8) :: mubar_z(nz) ! needed for initialization
  character(len=10) :: message
  integer :: i, j, k, ii
  
  cvode_stepper = 0
  err = ''
  
  mxsteps_ = max_cvode_steps
  neq_long = neq
  tcur   = t0
  mu = nq
  ml = nq
  
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
  do i = 1,nq
    do j = 1,nz
      k = i + (j-1)*nq
      yvec(k) = usol_start(i,j)
    enddo
  enddo
  
  ! create SUNDIALS N_Vector
  sunvec_y => FN_VMake_Serial(neq_long, yvec)
  if (.not. associated(sunvec_y)) then
    err = "CVODE setup error."
    return
  end if
  
  ! create CVode memory
  cvode_mem = FCVodeCreate(CV_BDF)
  if (.not. c_associated(cvode_mem)) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), t0, sunvec_y)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  sunmat => FSUNBandMatrix(neq_long, mu, ml)
  sunlin => FSUNLinSol_Band(sunvec_y,sunmat)
  
  ierr = FCVodeSetLinearSolver(cvode_mem, sunlin, sunmat)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  if (use_fast_jacobian) then
    ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
  endif
  
  ierr = FCVodeSetMaxNumSteps(cvode_mem, mxsteps_)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSetInitStep(cvode_mem, initial_dt)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSetMaxErrTestFails(cvode_mem, max_err_test_failures)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSetMaxOrd(cvode_mem, max_order)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  do ii = 1, num_t_eval
    ierr = FCVode(cvode_mem, t_eval(ii), sunvec_y, tcur, CV_NORMAL)
    if (ierr /= 0) then
      success = .false.
      if ((global_err == 'max steps') .and. (ierr == -8)) then ! reached max steps
        write(message,'(i10)')  max_cvode_steps 
        print*,'CVODE stopped because it reached the maximum number of of specified steps: '//trim(message)
        call FN_VDestroy(sunvec_y)
        call FCVodeFree(cvode_mem)
        ierr = FSUNLinSolFree(sunlin)
        call FSUNMatDestroy(sunmat)
        return
      else if ((global_err /= 'max steps') .and. (ierr == -8)) then ! err in my rhs
        err = global_err
        call FN_VDestroy(sunvec_y)
        call FCVodeFree(cvode_mem)
        ierr = FSUNLinSolFree(sunlin)
        call FSUNMatDestroy(sunmat)
        return
      else
        write(message,'(i4)') ierr
        err = 'SUNDIALS_ERROR: FCVODE '//trim(message)
        call FN_VDestroy(sunvec_y)
        call FCVodeFree(cvode_mem)
        ierr = FSUNLinSolFree(sunlin)
        call FSUNMatDestroy(sunmat)
        return
      endif
    else
      success = .true.

      ! save the solution
      do I=1,NQ
        do J=1,NZ
          K = I + (J-1)*NQ
          solution(ii,i,j) = yvec(k)
        enddo
      enddo

      ! fix lbound = 1
      do i=1,nq
        if (lbound(i) .eq.1) then
          solution(ii,i,1) = fixedmr(i)
        endif
      enddo
  
    endif
    
  enddo
    
  ! free memory
  call FN_VDestroy(sunvec_y)
  call FCVodeFree(cvode_mem)
  ierr = FSUNLinSolFree(sunlin)
  if (ierr /= 0) then
    err = "CVODE deallocation error"
    return
  end if
  call FSUNMatDestroy(sunmat)

end subroutine


subroutine cvode_save(t0, usol_start, nnq, nnz, t_eval, num_t_eval, rtol, atol, &
                      use_fast_jacobian, outfilename, amount2save, success, err)
  use photochem_data, only: neq, nq, nz, jtrop, ispec, np, nr, background_spec, &
                            nw, wavl, z, jchem, kj, ks, photoreac, photonums, photospec, &
                            jchem, nmax, iprod, iloss,nump,numl, nsp, background_mu, mass, rainout_on, &
                            Flux                         
  use photochem_vars, only: lbound, fixedmr, T, den, P, Press, &
                            rpar_init, edd, &
                            max_cvode_steps, initial_dt, max_err_test_failures, max_order
  use photochem_wrk, only: cvode_stepper, rain, raingc, global_err, rpar, surf_radiance, A, yp, yl, D

  use, intrinsic :: iso_c_binding
  use fcvode_mod, only: CV_BDF, CV_NORMAL, FCVodeInit, FCVodeSStolerances, &
                        FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                        FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                        FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd
  use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
  use fnvector_serial_mod, only: FN_VMake_Serial   
  use fsunmatrix_band_mod, only: FSUNBandMatrix
  use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
  use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
  use fsunlinsol_band_mod, only: FSUNLinSol_Band

  implicit none

  real(8), intent(in) :: t0 ! intial time
  real(8), dimension(nnq,nnz) :: usol_start ! initial atmosphere
  integer, intent(in) :: nnq
  integer, intent(in) :: nnz
  real(8), dimension(num_t_eval), intent(in) :: t_eval ! times to evaluate solution
  integer, intent(in) :: num_t_eval ! dimension of t_eval
  real(8), intent(in) :: rtol ! relative tolerance
  real(8), intent(in) :: atol ! absolute tolerance
  logical, intent(in) :: use_fast_jacobian ! if True then use my jacobian (keep false!)
  character(len=*), intent(in) :: outfilename ! were to save the solution
  integer, intent(in) :: amount2save
  
  logical, intent(out) :: success
  character(len=1000), intent(out) :: err
  
  ! local
  real(c_double) :: tcur(1)    ! current time
  integer(c_int) :: ierr       ! error flag from C functions
  type(c_ptr)    :: cvode_mem  ! CVODE memory
  type(N_Vector), pointer :: sunvec_y ! sundials vector
  real(c_double) :: yvec(neq)
  integer(c_long) :: neq_long
  integer(c_long) :: mu, ml
  integer(c_long) :: mxsteps_ ! unused
  type(SUNMatrix), pointer :: sunmat
  type(SUNLinearSolver), pointer :: sunlin
  
  real(8), allocatable :: UDOT(:), loss(:,:), solution_temp(:,:)
  real(8), allocatable :: photorates(:,:)
  integer :: ind(1)
  real(8) :: mubar_z(nz) ! needed for initialization
  character(len=10) :: message
  integer :: i, j, k, ii
  
  allocate(UDOT(neq))
  allocate(loss(nq,nz))
  allocate(photorates(ks,nz))
  allocate(solution_temp(nnq, nnz))
  
  cvode_stepper = 0
  err = ''
  
  mxsteps_ = max_cvode_steps
  neq_long = neq
  tcur   = t0
  mu = nq
  ml = nq
  
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
  do i = 1,nq
    do j = 1,nz
      k = i + (j-1)*nq
      yvec(k) = usol_start(i,j)
    enddo
  enddo
  call right_hand_side(yvec, udot, neq, err)
  if (len_trim(err) /= 0) return
  
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
  
  ! create SUNDIALS N_Vector
  sunvec_y => FN_VMake_Serial(neq_long, yvec)
  if (.not. associated(sunvec_y)) then
    err = "CVODE setup error."
    return
  end if
  
  ! create CVode memory
  cvode_mem = FCVodeCreate(CV_BDF)
  if (.not. c_associated(cvode_mem)) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), t0, sunvec_y)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  sunmat => FSUNBandMatrix(neq_long, mu, ml)
  sunlin => FSUNLinSol_Band(sunvec_y,sunmat)
  
  ierr = FCVodeSetLinearSolver(cvode_mem, sunlin, sunmat)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  if (use_fast_jacobian) then
    ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
  endif
  
  ierr = FCVodeSetMaxNumSteps(cvode_mem, mxsteps_)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSetInitStep(cvode_mem, initial_dt)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSetMaxErrTestFails(cvode_mem, max_err_test_failures)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FCVodeSetMaxOrd(cvode_mem, max_order)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  
  do ii = 1, num_t_eval
    ierr = FCVode(cvode_mem, t_eval(ii), sunvec_y, tcur, CV_NORMAL)
    if (ierr /= 0) then
      success = .false.
      if ((global_err == 'max steps') .and. (ierr == -8)) then ! reached max steps
        write(message,'(i10)')  max_cvode_steps 
        print*,'CVODE stopped because it reached the maximum number of of specified steps: '//trim(message)
        call FN_VDestroy(sunvec_y)
        call FCVodeFree(cvode_mem)
        ierr = FSUNLinSolFree(sunlin)
        call FSUNMatDestroy(sunmat)
        return
      else if ((global_err /= 'max steps') .and. (ierr == -8)) then ! err in my rhs
        err = global_err
        call FN_VDestroy(sunvec_y)
        call FCVodeFree(cvode_mem)
        ierr = FSUNLinSolFree(sunlin)
        call FSUNMatDestroy(sunmat)
        return
      else
        write(message,'(i4)') ierr
        err = 'SUNDIALS_ERROR: FCVODE '//trim(message)
        call FN_VDestroy(sunvec_y)
        call FCVodeFree(cvode_mem)
        ierr = FSUNLinSolFree(sunlin)
        call FSUNMatDestroy(sunmat)
        return
      endif
    else
      success = .true.

      ! save the solution
      do I=1,NQ
        do J=1,NZ
          K = I + (J-1)*NQ
          solution_temp(i,j) = yvec(k)
        enddo
      enddo

      ! fix lbound = 1
      do i=1,nq
        if (lbound(i) .eq.1) then
          solution_temp(i,1) = fixedmr(i)
        endif
      enddo
      
      call right_hand_side(yvec, udot, neq, err)
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
  
    endif
    
  enddo
  
  ! free memory
  call FN_VDestroy(sunvec_y)
  call FCVodeFree(cvode_mem)
  ierr = FSUNLinSolFree(sunlin)
  if (ierr /= 0) then
    err = "CVODE deallocation error"
    return
  end if
  call FSUNMatDestroy(sunmat)

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