  subroutine integrate(nsteps,converged,err)
    use photochem_data, only: nw, nr, nz, nz1, nq, nq1, kj, np, lda, neq, isl, &
                              agl, epsj, frak, hcdens, ino, io2, zy, &
                              lco, lh, lh2, lhcaer, lhcaer2, lh2o, lo, &
                              ls4, ls8aer, lso4aer, &
                              z, dz, jtrop, ispec, photoreac, photonums, nsp2, &
                              background_spec, lightning, mass, background_mu, rainout_on, &
                              fix_water_in_troposphere, P0, light_disp_rate, &
                              estimate_CO2_photo_above_grid
    use photochem_vars, only: verbose, usol_init, usol_out, rpar_init, wfall_init, aersol_init, &
                              lbound, fixedmr, vdep, vdep0, veff, veff0, smflux, sgflux, &
                              distheight, distflux, mbound, T, den, edd, fluxo, flow, H2Osat, P, &
                              press, equilibrium_time
    use photochem_wrk, only: rpar, wfall, aersol, hscale, scale_h, h_atm, bx1x2, &
                             A, yl, yp, rain, raingc, &
                             adl, add, adu, dl, dd, du, dk, &
                             prod_rates, tauedd, &
                             H2SO4S, S8S, fsulf, surf_radiance, D
    implicit none
    ! module variables
    ! all of them?

    ! local variables
    integer,intent(in) :: nsteps
    logical,intent(out) :: converged
    character(len=1000), intent(out) :: err

    real*8, dimension(nq,nz) :: usol
    real*8, dimension(nq,nz) :: usaveold
    real*8, dimension(nq,nz) :: usave
    ! real*8, dimension(nq,nz) :: usolprev
    real*8 Dt,dtinv, time, tstop
    integer nn, i, j, k, jj, nparti
    integer JO2_O1D,JCO2,JO2, JCO2_O1D
    real*8 absorbers(kj,nz)
    real*8,dimension(nz) :: H2O, O2, CO2
    real*8,dimension(kj,nz) :: prates
    integer NPSO4,NPS8, NPHC, NPHC2, MM
    real*8 VCO2
    real*8, dimension(lda,neq) :: djac
    real*8, dimension(neq) :: RHS
    integer, dimension(neq) :: ipvt
    real*8, dimension(nz) :: R
    real*8, dimension(nz,np) :: conver
    real*8, dimension(nq) :: U
    integer LB, MB, L, is, m, n, LL
    real*8 disth, dtsave, emax, erel, smallest, RMAX, UMAX
    integer jdisth, js, kd, ku, kl, info
    integer j25, j70, j50
    integer lcountso4, lcounts8, lcountHC, lcountHC2, lpolyscount
    real*8 ZTOP, ZTOP1, WNF
    real*8 DPU(NZ,NP),DPL(NZ,NP)
    real*8, dimension(nq,nz) :: Fval, fv
    integer cr, cm, c1, c2
    real*8,dimension(nq1) :: SR, FUP
    ! real*8, dimension(nsp2,nz) :: D
    real(8) :: mubar_z(nz)

    converged = .true.
    err = ''
    info = 0

    ! set atmosphere to inital atmosphere
    usol = usol_init
    ! aersol parameters to initial parameters
    if (np.gt.0) then
      rpar = rpar_init
      wfall = wfall_init
      aersol = aersol_init
    endif
    ! begin stuff that needs to be inizialized
    do i = 1,nz
      call mean_molecular_weight(nq, usol_init(:,i), mass, background_mu, mubar_z(i), err)
      if (len_trim(err) > 0) return
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

    do i=1,nq
      if (LBOUND(i).EQ.1) USOL(i,1)=fixedmr(i)
    enddo
    
    DT = 1.D-6
    DTINV = 1./DT
    TIME = 0.
    TSTOP = equilibrium_time
    ! NSTEPS = 50000
    nn = 0
    KD = 2*NQ + 1
    KU = KD - NQ
    KL = KD + NQ

    NPSO4 = 0
    NPS8 = 0
    NPHC = 0
    NPHC2 = 0
    is = 1
    emax = 0.d0

    ! start the time-stepping loop
    call system_clock(count = c1, count_rate = cr, count_max = cm)
    do n = 1,nsteps
      if (verbose) then
        print"(1x,'N =',i6,3x,'Time = ',es12.5,3x,'DT = ',es12.5"// &
        ",3x,'emax = ',es12.5,3x,'for ',a8)", &
        n,time,dt,emax,ispec(is)
      endif
      TIME = TIME + DT
      nn = nn+1
      
      do i = 1,nz
        call mean_molecular_weight(nq, usol(:,i), mass, background_mu, mubar_z(i), err)
        if (len_trim(err) > 0) return
      enddo
      call densty(nz, mubar_z, T, den, P, press) 
      call difco(nq, nz, mubar_z, T, den, edd, &
                hscale, tauedd, DK, H_atm, bx1x2, scale_H)
      if (fix_water_in_troposphere) then
        call photsatrat(nz, T, P, den, Jtrop, H2Osat, H2O) ! H2o mixing ratio
        DO J=1,JTROP
          USOL(LH2O,J) = H2O(J) 
        enddo
      endif
      if (lightning) then 
        call ltning(nq, nz, usol, P0, light_disp_rate, prod_rates, err)
        
      endif
      call diffusion_coeffs(nq, nz, den, dz, DK, bx1x2, scale_H, H_atm, &
                           DU, DL, DD, ADU, ADL, ADD)
                                    
      ! below is H escape
      if (background_spec /= 'H2') then ! then we can consider its upper and lower boundary
       if (mbound(LH2) == 0) then
         Veff(LH2) = 1.0*bx1x2(lh2,nz)/DEN(NZ)*(1./Hscale(nz) &
         - 1./scale_H(LH2,nz))
       endif
      endif  
      if (mbound(lH) == 0) then
       Veff(LH) = 1.0*bx1x2(LH,nz)/DEN(NZ)*(1./Hscale(nz) &
        - 1./scale_H(LH,nz))
      endif
      
      call rates(nz, nr, T, den, A, err)
      if (len_trim(err) > 0) return
      call dochem(-1, nr, nsp2, nq, nz, usol, A, isl, jtrop, D, fval)

      lpolyscount = 0
      do k=1,kj
        do i=1,nz
          if (ISPEC(INT(photoreac(k))).eq.'S8      ') then
            absorbers(k,i)=ABS(usol(INT(photoreac(k)),i))
            if (lpolyscount .eq. 2) absorbers(k,i)=0.0
            if (i.eq.nz) then
              lpolyscount=lpolyscount+1
            endif
          else
            absorbers(k,i)=ABS(USOL(INT(photoreac(k)),i))
          endif
        enddo
      enddo

      JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2    ')
      JCO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CO2    ')
      JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')
      do I=1,NZ
        O2(I) =  absorbers(JO2_O1D,I)
        CO2(I) = absorbers(JCO2,I)
      enddo

      if (rainout_on) then
        call rainout(.false.,Jtrop,Usol,nq,nz, T,den, rain, raingc, err)
        if (len_trim(err) /= 0) return
      else
        rain = 0.d0
        raingc = 0.d0
      endif

      call aercon(nq, nz, usol, P, T, H2SO4S, S8S, fsulf)

      if (NP.GT.0) then   !particles in main loop
        CALL SEDMNT(frak,HCDENS,nz,np,conver, .true.)

        do J=1,NZ
          do JJ=1, NP
            if (JJ.eq.1)  nparti = LSO4AER
            if (JJ.eq.2)  nparti = LS8AER
            if (JJ.eq.3)  nparti = LHCAER
            if (JJ.eq.4)  nparti = LHCAER2
            AERSOL(J,JJ) = USOL(nparti,J)*DEN(J)/(CONVER(J,JJ))
          enddo
        enddo
        do J=1,NP
          DPU(1,J) = WFALL(2,J)*DEN(2)/DEN(1)/(2.*DZ(1))
          DPL(NZ,J) = WFALL(NZ1,J)*DEN(NZ1)/DEN(NZ)/(2.*DZ(NZ))
          DO I=2,NZ1
            DPU(I,J) = WFALL(I+1,J)*DEN(I+1)/DEN(I)/(2.*DZ(I))
            DPL(I,J) = WFALL(I-1,J)*DEN(I-1)/DEN(I)/(2.*DZ(I))
          enddo
        enddo
      ! else
        ! print*,'Note: Since NP = 0, did not call SDMNT...'
      endif
      
      ! TIME-DEPENDENT BOUNDARY CONDITIONS
      if(NP.EQ.0) then
!-mab Time-dependent boundary conditions for particles set to 0
        NPSO4 = 0
        NPS8 = 0
        NPHC = 0
        NPHC2 = 0
      elseif(NP.GT.0) then
        NPSO4 = LSO4AER - NQ + NP
        NPS8 = LS8AER - NQ + NP
        NPHC = LHCAER - NQ + NP
        NPHC2 = LHCAER2 - NQ + NP

        VDEP(LSO4AER) = VDEP0(LSO4AER) + WFALL(1,NPSO4)
        VEFF(LSO4AER) = VEFF0(LSO4AER) + WFALL(NZ,NPSO4)

        VDEP(LS8AER) = VDEP0(LS8AER) + WFALL(1,NPS8)
        VEFF(LS8AER) = VEFF0(LS8AER) + WFALL(NZ,NPS8)

        if(NP.ge.3) then
          VDEP(LHCAER) = VDEP0(LHCAER) + WFALL(1,NPHC)
          VEFF(LHCAER) = VEFF0(LHCAER) + WFALL(NZ,NPHC)

          if (NP.eq.4) then
            VDEP(LHCAER2) = VDEP0(LHCAER2) + WFALL(1,NPHC2)
            VEFF(LHCAER2) = VEFF0(LHCAER2) + WFALL(NZ,NPHC2)
          endif
        endif
      endif
      
      call photo(zy, agl, io2, ino, usol, mubar_z, D, nsp2, nw, &
                 nq, nz, kj, prates, surf_radiance)
      do j=1,kj
        do i=1,nz
          A(photonums(j),i)=prates(j,i)
        enddo
      enddo
      
! estimate CO2 photolysis above the top of the grid
! and return CO + O to the upper grid point
! this behavior is turned on and off by setting MBOUND=2 for CO2 in species.dat
      JCO2_O1D = JCO2+1
      VCO2 = (prates(JCO2,NZ) + prates(JCO2_O1D,NZ) ) * HSCALE(NZ)
      if (estimate_CO2_photo_above_grid) then
        SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ)
        SMFLUX(LCO) = SMFLUX(LO)
      endif

! ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      djac = 0.d0
      rhs = 0.d0

      call dochem(0, nr, nsp2, nq, nz, usol, A, isl, jtrop, D, fval)

      DO I=1,NQ
        DO J=1,NZ
          K = I + (J-1)*NQ
          RHS(K) = FVAL(I,J)
          USAVEOLD(I,J) = USOL(I,J)
!-mc testing 4/29/06  used to revert if timestep is too big.
          USAVE(I,J) = USOL(I,J)
        enddo
      enddo

      !$OMP PARALLEL PRIVATE(i,j,R,usol,mm,m,k,fv,D)
      DO I=1,NQ
        DO J=1,NZ
          USOL(I,J) = USAVE(I,J)
        enddo
      enddo


      ! print*,usol(1,1)


      !$OMP DO
      DO I=1,NQ
      ! Loop through all vertical atmospheric layers
        DO J=1,NZ
          R(J) = EPSJ * USOL(I,J)
      !   as it was - USOL should be positive here anyway,
          IF(R(J).LT.1.D-100) R(J) = 1.D-100
      !   Add perturbing quantity to mixing ratio
          USOL(I,J) = USAVE(I,J) + R(J)
        enddo
      !   Call the photochemistry routine
      !      FV has dimension (NQ,NZ) and holds gas densities
        call dochem(0, nr, nsp2, nq, nz, usol, A, isl, jtrop, D, fv)



        DO M=1,NQ
          MM = M - I + KD
          DO J=1,NZ
            K = I + (J-1)*NQ
      !  -J since its orig - perturbed
            DJAC(MM,K) = (FVAL(M,J) - FV(M,J))/R(J)
          enddo
        enddo



        DO J=1,NZ
          USOL(I,J) = USAVE(I,J)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      ! stop


!     COMPUTE TRANSPORT TERMS AT INTERIOR GRID POINTS
      DO I = 1,NQ
        DO J=2,NZ1
          K = I + (J-1)*NQ
          RHS(K) = RHS(K) - DD(i,J)*USOL(I,J)-ADD(i,j)*USOL(I,J) &
            + DU(i,J)*USOL(I,J+1) + ADU(i,j)*USOL(i,j+1) &
            + DL(i,J)*USOL(I,J-1) + ADL(i,J)*USOL(I,J-1)
          DJAC(KD,K) = DJAC(KD,K) + DTINV + DD(i,J) + ADD(i,j)
          DJAC(KU,K+NQ) = - DU(i,J) - ADU(i,j)
          DJAC(KL,K-NQ) = - DL(i,J) - ADL(i,j)
        enddo
      enddo



      if(NP.gt.0) then  !particles in main loop
! C ack - these need to be abstracted... WARNING
! C   ADD ADVECTION TERMS FOR PARTICLES
        do L=1,NP
          if(L.EQ.1) I = LSO4AER
          if(L.EQ.2) I = LS8AER
          if(L.EQ.3) I = LHCAER
          if(L.EQ.4) I = LHCAER2
          do J=2,NZ1
            K = I + (J-1)*NQ
            RHS(K) = RHS(K) + DPU(J,L)*USOL(I,J+1) - DPL(J,L)*USOL(I,J-1)
            DJAC(KU,K+NQ) = DJAC(KU,K+NQ) - DPU(J,L)
            DJAC(KL,K-NQ) = DJAC(KL,K-NQ) + DPL(J,L)
          enddo
        enddo
      endif

! ***** LOWER BOUNDARY CONDITIONS *****
      DO K=1,NQ
        U(K) = USOL(K,1)
!   OK as long as we don't model atmospheric Boron (WARNING?)
        LB = LBOUND(K)
        if (LB.eq.0 .OR. LB.eq.3) then
!       CONSTANT DEPOSITION VELOCITY/SPECIES WITH DISTRIBUTED FLUXES
          RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K)) &
          + ADU(k,1)*(USOL(K,2) + U(K)) - VDEP(K)*U(K)/DZ(1)
          DJAC(KD,K) = DJAC(KD,K) +DTINV +DU(k,1) -ADU(k,1) &
          + VDEP(K)/DZ(1)
          DJAC(KU,K+NQ) = - DU(k,1) - ADU(k,1)
!  is this right for particles??
        else if (LB .eq. 1) then
!       CONSTANT MIXING RATIO
          RHS(K) = 0.d0
          do M=1,NQ
            MM = KD + K - M
            DJAC(MM,M) = 0.d0
          enddo
          DJAC(KU,K+NQ) = 0.d0
          DJAC(KD,K) = DTINV + DU(k,1) - ADU(k,1)
        else if (LB == 2) then
!       CONSTANT UPWARD FLUX
          RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K)) &
          + ADU(k,1)*(USOL(K,2) + U(K)) + SGFLUX(K)/DEN(1)/DZ(1)
          DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(k,1) - ADU(k,1)
          DJAC(KU,K+NQ) = - DU(k,1) - ADU(k,1)
        else if (LB == -1) then 
          ! Julie Moses 2000 boundary condition
          RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K)) &
          + ADU(k,1)*(USOL(K,2) + U(K)) - (edd(1)/h_atm(1))*U(K)/DZ(1)
          DJAC(KD,K) = DJAC(KD,K) +DTINV +DU(k,1) -ADU(k,1) &
          + (edd(1)/h_atm(1))/DZ(1)
          DJAC(KU,K+NQ) = - DU(k,1) - ADU(k,1)
          
        endif
      enddo

! ***** UPPER BOUNDARY CONDITIONS *****
      DO I=1,NQ
        U(I) = USOL(I,NZ)
        K = I + NZ1*NQ
        MB = MBOUND(I)
        if (MB.eq.0) then
!       CONSTANT EFFUSION VELOCITY
          RHS(K) = RHS(K) + DL(i,NZ)*(USOL(I,NZ1) - U(I)) &
          + ADL(i,NZ)*(USOL(I,NZ1) + U(I)) - VEFF(I)*U(I)/DZ(NZ)
          DJAC(KD,K) = DJAC(KD,K) +DTINV +DL(i,NZ) -ADL(i,NZ) &
          + VEFF(I)/DZ(NZ)
          DJAC(KL,K-NQ) = - DL(i,NZ) -ADL(i,NZ)
        else if (MB.eq. 1) then
!       Constant mixing ratio at the top. not debugged WARNING
          RHS(K) = 0.d0
          do M=1,NQ
            MM = KD + K - M
            DJAC(MM,M) = 0.d0
          enddo
          DJAC(KU,K+NQ) = 0.d0
          DJAC(KD,K) = DTINV + DL(i,NZ) -ADL(i,NZ)
        else
!   CONSTANT UPWARD FLUX
          RHS(K) = RHS(K) + DL(i,NZ)*(USOL(I,NZ1) - U(I)) &
          + ADL(i,NZ)*(USOL(I,NZ1) + U(I)) - SMFLUX(I)/DEN(NZ)/DZ(NZ)
          DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(i,NZ) - ADL(i,NZ)
          DJAC(KL,K-NQ) = - DL(i,NZ) - ADL(i,NZ)
        endif
      enddo

!   HOLD H2O AND S8 CONSTANT BELOW ZTROP
      if (fix_water_in_troposphere) then
        L = LH2O
        DO J=1,JTROP
          K = L + (J-1)*NQ
          RHS(K) = 0.
          DO M=1,NQ
            MM = M - L + KD
            DJAC(MM,K) = 0.
          enddo
          DJAC(KD,K) = DTINV
          DJAC(KU,K+NQ) = 0.
          IF(J.NE.1) then
            DJAC(KL,K-NQ) = 0.
          endif
        enddo
      endif



!  distributed (volcanic) sources

      do i=1,nq
        if (LBOUND(i).eq.3) then
!     convert to cm
          disth=distheight(i)*1.e5
!     height index (-1 given the staggered grid)
!     the 1 in the second postion tells minloc to return a scalar
          jdisth=minloc(Z,1, Z .ge. disth)-1
          if (jdisth < 2) then
            jdisth = 2
          endif

          ZTOP=Z(jdisth)-Z(1)
          ZTOP1=Z(jdisth)+0.5*DZ(jdistH)

!     distribute from second level to distheight
          do j=2,jdisth
            K = i + (J-1)*NQ
            rhs(k) = rhs(k) + 2.*distflux(i)*(ZTOP1-Z(j))/(Den(j)*ZTOP**2)
          enddo
        endif
      enddo


      CALL SGBFA(DJAC,LDA,NEQ,NQ,NQ,IPVT,INFO)
      IF(INFO.NE.0) then
        ! print *, N,INFO
        ! print *, 'ERROR in SGBFA'
        err = 'Linear solve failed in subroutine integrate'
        return
      endif

!-mc  we have set up RHS so that DJAC*X=RHS
!-mc DJAC is now upper triangular...

      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)

      J25=minloc(Z,1, Z/1.e5 .ge. 25)-1
      J70=minloc(Z,1, Z/1.e5 .ge. 70)-1
      J50=minloc(Z,1, Z/1.e5 .ge. 50)-1

!   COMPUTE NEW CONCENTRATIONS (IGNORE ERRORS IN SEVERAL SPECIES
!    THAT VIRTUALLY DISAPPEAR UP HIGH)
      EMAX = 0.d0
      DO I=1,NQ
        DO J=1,NZ
          K = I + (J-1)*NQ
          IF (I.EQ.LSO4AER.AND. J.GT.J25) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          ELSEIF(I.EQ.LS8AER.AND. J.GT.J25) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          ELSEIF(I.EQ.LHCAER.AND. J.GT.J25) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          ELSEIF(I.EQ.LHCAER2.AND. J.GT.J25) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          ELSEIF(I.EQ.LS4) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          ELSEIF(USOL(I,J).LT. 1.E-20) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          else
            EREL = ABS(RHS(K)/USOL(I,J))
            EMAX = max(EMAX,EREL)
            IF(EREL.LT.EMAX) THEN
              USOL(I,J) = USOL(I,J) + RHS(K)
!       store info on species with largest error
            ELSE
              IS = I
!       mc -this label is OK, because S will never have a photolysis reaction
              JS = J
              UMAX = USOL(I,J)
              RMAX = RHS(K)
              USOL(I,J) = USOL(I,J) + RHS(K)
            endif
          endif
        enddo
      enddo

      ! make sure H2O doesn't change below trop
      if (fix_water_in_troposphere) then
        DO J=1,JTROP
          USOL(LH2O,J) = H2O(J)
        enddo
      endif

      if (NP.GT.0) then
! switch around main loop particles
! 1.e-38 is the smallest number for single precision. We should upgrade
! this to double precision at some point. -Shawn D-G  WARNING
        smallest = 1.e-38
        lcountso4=0
        lcounts8=0
        lcountHC=0
        lcountHC2=0


        DO J=1,NZ
          if(USOL(LSO4AER,J).LT.0) lcountso4=lcountso4+1

          if(lcountso4.gt.0) then
            if (J.GT.1) then
              USOL(LSO4AER,J)= &
              USOL(LSO4AER,J-1)*EXP(-WFALL(J,NPSO4)*DZ(J)/EDD(J))
            else
              USOL(LSO4AER,J)=-1
            endif
          endif

          if(USOL(LS8AER,J).LT.0) lcounts8=lcounts8+1
          if(lcounts8.gt.0) then
            if(J.GT.1) then
              USOL(LS8AER,J) = &
              USOL(LS8AER,J-1) * EXP(-WFALL(J,NPS8)*DZ(J)/EDD(J))
            else
              USOL(LS8AER,J)=-1
            endif
          endif

          if (NP.GE.3) then
            if(USOL(LHCAER,J).LT.0) lcountHC=lcountHC+1

            if(lcountHC.gt.0) then
              if (J.GT.1) then
                USOL(LHCAER,J) = &
                USOL(LHCAER,J-1)*EXP(-WFALL(J,NPHC)*DZ(J)/EDD(J))
              else
                USOL(LHCAER,J)=-1
              endif
            endif

            if (NP.EQ.4) then
              if(USOL(LHCAER2,J).LT.0) lcountHC2=lcountHC2+1

              if(lcountHC2.gt.0) then
                if (J.GT.1) then
                  USOL(LHCAER2,J)= &
                  USOL(LHCAER2,J-1)*EXP(-WFALL(J,NPHC2)*DZ(J)/EDD(J))
                else
                  USOL(LHCAER2,J)=-1
                endif
              endif
            endif
          endif

          USOL(LSO4AER,J) = max (USOL(LSO4AER,J),smallest)
          USOL(LS8AER,J) = max(USOL(LS8AER,J),smallest)
          if (NP.GE.3) then
            USOL(LHCAER,J) = max (USOL(LHCAER,J),smallest)
            USOL(LHCAER2,J) = max (USOL(LHCAER2,J),smallest)
          endif
        enddo
! end particle switching
      endif


      !test
      do i=1,nq
        do j=1,nz
          USOL(i,j)=abs(USOL(i,j))
        enddo
      enddo

!   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
      IF(EMAX.LT.0.15)  DT = 1.1*DTSAVE
      IF(EMAX.LT.0.07)  DT = 1.2*DTSAVE
      IF(EMAX.LT.0.01)  DT = 1.4*DTSAVE
      IF(EMAX.LT.0.008)  DT = 2.0*DTSAVE
      IF(EMAX.LT.0.004)  DT = 3.0*DTSAVE
      IF(EMAX.LT.0.001) DT = 4.0*DTSAVE
      IF(EMAX.LT.0.0005) DT = 5.*DTSAVE
      DTINV = 1./DT

!   RETRY TIME STEP IF EMAX EXCEEDS 25 PERCENT
      IF (EMAX.GT.0.25) THEN
        DT = 0.7*DTSAVE
        TIME = TIME - DTSAVE
        do I=1,NQ
          do J=1,NZ
            USOL(I,J) = USAVEOLD(I,J)
          enddo
        enddo
!   valid timestep, so update USOLPREV vector
      ! ELSE
      !   do i=1,nq
      !     do j=1,nz
      !       USOLPREV(I,J)=USAVEOLD(I,J)
      !     enddo
      !   enddo
      endif

      IF (INFO.NE.0) STOP
      IF (NN.EQ.NSTEPS) then
        converged = .false.
        if (verbose) then
          print"('Photochemical model did not converge in ',i6,' steps')",&
                  nsteps
        endif
        exit
      endif

      IF (TIME.GT.TSTOP) then
        converged = .true.
        exit
      endif
    ! should never finish the loop if converged
    converged = .false.
    enddo
! ***** END THE TIME-STEPPING LOOP *****
    call system_clock(count = c2)

    if (converged) then
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

    ! the flux
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
