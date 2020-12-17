  subroutine integrate
    implicit none
    ! module variables
    ! all of them?

    ! local variables
    real*8, dimension(nq,nz) :: usol
    real*8, dimension(nq,nz) :: usaveold
    real*8, dimension(nq,nz) :: usave
    real*8, dimension(nq,nz) :: usolprev
    real*8, dimension(nq,nz) :: REL
    real*8 Ha
    real*8 Dt,dtinv, time, tstop
    integer nn, i, j, k, jj, nparti, nsteps
    integer jh2o,JO2_O1D,JO3_O1D,JCO2,JO1D,JO2, JCO2_O1D
    real*8 absorbers(kj,nz)
    real*8,dimension(nz) :: H2O, O2, O3, CO2
    real*8,dimension(kj,nz) :: prates
    integer NPSO4,NPS8, NPHC, NPHC2, MM
    real*8 VCO2
    real*8, dimension(lda,neq) :: djac
    real*8, dimension(neq) :: RHS
    integer, dimension(neq) :: ipvt
    real*8, dimension(nz) :: R
    real*8, dimension(nz,np) :: conver
    real*8, dimension(nq) :: U
    integer LB, MB, L, is, m, n
    real*8 disth, dtsave, emax, erel, smallest, RMAX, UMAX
    integer jdisth, js, kd, ku, kl, info
    integer j15, j25, j70, j50
    integer lcountso4, lcounts8, lcountHC, lcountHC2, lpolyscount
    real*8 ZTOP, ZTOP1
    real*8 DPU(NZ,NP),DPL(NZ,NP)
    real*8, dimension(nq,nz) :: Fval, fv


    do i=1,nq
      do j=1,nz
        usol(i,j) = usol_init(i,j)
      enddo
    enddo
    do i=1,nq
      if (LBOUND(i).EQ.1) USOL(i,1)=fixedmr(i)
    enddo
    call photsatrat(jtrop,nz,h2o)
    DO J=1,JTROP
      USOL(LH2O,J) = H2O(J) ! problme!!
    enddo

    HA = HSCALE(NZ)
    DT = 1.D-6
    DTINV = 1./DT
    TIME = 0.
    TSTOP = 1.D17
    NSTEPS = 50000
    nn = 0
    KD = 2*NQ + 1
    KU = KD - NQ
    KL = KD + NQ




    ! start the time-stepping loop
    do n = 1,nsteps
      print"(2x,'N =',i6,3x,'Time = ',es10.2,3x,'DT = ',es10.2)", &
      n,time,dt
      TIME = TIME + DT
      nn = nn+1

      call rates

      lpolyscount = 0
      do k=1,kj
        do i=1,nz
      !     this gets any SL/IN species that have photorates
          if (photoreac(k).gt.nq) then
            absorbers(k,i)=ABS(SL(INT(photoreac(k)),i))/DEN(I)
      !     quasi-hardcoded S8 behavior WARNING
          else if (ISPEC(INT(photoreac(k))).eq.'S8      ') then
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

      JH2O=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'H2O    ')
      JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2    ')
      JO3_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O3    ')
      JCO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CO2    ')
      JO1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O1D    ')
      JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')
      do I=1,NZ
        H2O(I) = absorbers(JH2O,I)
        O2(I) =  absorbers(JO2_O1D,I)
        O3(I) = absorbers(JO3_O1D,I)
        CO2(I) = absorbers(JCO2,I)
      enddo

      call photo(zy, agl, io2, ino, usol, nq, nz, kj, prates)


      if (n.eq.1) then
        call rainout(Jtrop,0,Usol,nq,nz)
      else
        call rainout(Jtrop,1,Usol,nq,nz)
      endif


      call aercon(usol,nq,nz)

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

! estimate CO2 photolysis above the top of the grid
! and return CO + O to the upper grid point
! this behavior is turned on and off by setting MBOUND=2 for CO2 in species.dat
      JCO2_O1D = JCO2+1
      VCO2 = (prates(JCO2,NZ) + prates(JCO2_O1D,NZ) ) * HA
      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ)
      SMFLUX(LCO) = SMFLUX(LO)

      if (NP.GT.0) then   !particles in main loop
        CALL SEDMNT(frak,HCDENS,ihztype,nz,np,conver)

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
      else
        print*,'Note: Since NP = 0, did not call SDMNT...'
      endif

! ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****

      DO J=1,LDA
        DO K=1,NEQ
          DJAC(J,K) = 0.d0
        enddo
      enddo
      DO K=1,NEQ
        RHS(K) = 0.d0
      enddo



      call dochem(Fval,0,jtrop,isl,usol,nq,nz)



      ! print*,fval


      DO I=1,NQ
        DO J=1,NZ
          K = I + (J-1)*NQ
          RHS(K) = FVAL(I,J)
          USAVEOLD(I,J) = USOL(I,J)
!-mc testing 4/29/06  used to revert if timestep is too big.
          USAVE(I,J) = USOL(I,J)
        enddo
      enddo

      !$OMP PARALLEL PRIVATE(i,j,R,usol,mm,m,k,fv)
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
        call dochem(Fv,0,jtrop,isl,usol,nq,nz)



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
        else
!       CONSTANT UPWARD FLUX
          RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K)) &
          + ADU(k,1)*(USOL(K,2) + U(K)) + SGFLUX(K)/DEN(1)/DZ(1)
          DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(k,1) - ADU(k,1)
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



!  distributed (volcanic) sources

      do i=1,nq
        if (LBOUND(i).eq.3) then
!     convert to cm
          disth=distheight(i)*1.e5
!     height index (-1 given the staggered grid)
!     the 1 in the second postion tells minloc to return a scalar
          jdisth=minloc(Z,1, Z .ge. disth)-1

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
        print *, N,INFO
        print *, 'ERROR in SGBFA'
        stop
      endif

!-mc  we have set up RHS so that DJAC*X=RHS
!-mc DJAC is now upper triangular...

      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)

      J15=minloc(Z,1, Z/1.e5 .ge. 15)-1
      J25=minloc(Z,1, Z/1.e5 .ge. 25)-1
      J70=minloc(Z,1, Z/1.e5 .ge. 70)-1
      J50=minloc(Z,1, Z/1.e5 .ge. 50)-1

!   COMPUTE NEW CONCENTRATIONS (IGNORE ERRORS IN SEVERAL SPECIES
!    THAT VIRTUALLY DISAPPEAR UP HIGH)
      EMAX = 0.
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
            REL(I,J) = RHS(K)/USOL(I,J)
            EREL = ABS(REL(I,J))
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
      DO J=1,JTROP
        USOL(LH2O,J) = H2O(J)
      enddo

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
      ELSE
        do i=1,nq
          do j=1,nz
            USOLPREV(I,J)=USAVEOLD(I,J)
          enddo
        enddo
      endif

      IF (INFO.NE.0) STOP
      IF (NN.EQ.NSTEPS) then
        exit
      endif

      IF (TIME.GT.TSTOP) then
        exit
      endif

    enddo
! ***** END THE TIME-STEPPING LOOP *****

  ! the output
  do i=1,nq
    do j=1,nz
      usol_out(i,j) = usol(i,j)
    enddo
  enddo


  end subroutine
