
subroutine dochem(N, nr, nsp2, nq, nz, usol, A, nshort, jtrop, D, fval)
  use photochem_data, only: nsp, lco, &
                            lh2, lh2o, &
                            lh2so4, lno, lo, lo2, lso4aer, ln2, lco2, &
                            planet, lightning, H2O_strat_condensation, &
                            background_spec, z
                            
  use photochem_vars, only: den, H2OSAT
  use photochem_wrk, only: prod_rates, raingc, &
                           yp, yl, H2SO4S

  implicit none

  ! local variables
  integer, intent(in) :: n
  integer, intent(in) :: jtrop
  integer, intent(in) :: nshort
  integer, intent(in) :: nq, nsp2
  integer, intent(in) :: nz, nr
  real*8, dimension(nq,nz),intent(in) :: usol
  real(8), intent(in) :: A(nr,nz)
  real*8, dimension(nq,nz),intent(out) :: Fval
  real*8, dimension(nsp2,nz), intent(out) :: D
  ! real*8, dimension(nq,nz) :: YP, YL

  integer i, j, ll, lla
  real*8 xp(nz), xl(nz), conso4(nz)
  real*8 xlj, confac
  integer jt1
  real*8 rhcold, h2ocrt, zap, CONDEN1

!   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.  PHOTO-
!   CHEMICAL EQUILIBRIUM SPECIES ARE DONE FIRST.  THESE MUST CON-
!   TAIN NO NONLINEARITIES (SUCH AS S4 REACTING WITH ITSELF TO FORM
!   S8) AND MUST BE DONE IN THE PROPER ORDER (I.E. IF SPECIES A
!   REACTS TO FORM B, THEN A MUST BE FOUND FIRST).  LONG-LIVED
!   SPECIES CAN BE DONE IN ANY ORDER.
      
  D = 0.d0

!compute number densities for long-lived and particle species
  do i = 1 , nq
    do j = 1 , nz
      d(i,j) = usol(i,j)*den(j)
    enddo
  enddo
!now do the last 4 inert species that are the same in both codes
  do j = 1 , nz
    D(nsp,j) = Den(j) * (1.d0 - sum(usol(:,j))) ! density of the background atmosphere
    d(nsp2-1,j) = 1.d0         ! HV has density of 1 for photorate calculations
    d(nsp2,j) = DEN(j)     ! M - background density for three body reactions
  enddo
  IF ( N.GE.0 ) THEN
! short lived densities
  do i = nq + 1 , nq + Nshort
    CALL CHEMPL(A,d,xp,xl,i)
    do j = 1 , nz
      d(i,j) = xp(j)/xl(j)
    enddo
  enddo

! ***** LONG-LIVED SPECIES CHEMISTRY *****
  do i = 1 , nq
    CALL CHEMPL(A, d,xp,xl,i)
    DO j = 1 , nz
      xlj = xl(j) + RAINGC(i,j)
      Fval(i,j) = xp(j)/DEN(j) - xlj*USOL(i,j)
      YP(i,j) = xp(j)
      YL(i,j) = xlj
    ENDDO
  ENDDO
  
  IF ( planet.EQ.'EARTH' ) THEN
    confac = 1.6D-5     !condensation factor
  ELSEIF ( planet.EQ.'MARS' ) THEN
    confac = 1.6D-5*10. ! reduce supersaturation of stratosphere
  ELSEIF ( planet.EQ.'DRY' ) THEN
    confac = 1.6D-5*10.  ! reduce supersaturation of stratosphere
  ENDIF
!   ZERO OUT H2O TERMS IN THE TROPOSPHERE AND INCLUDE LIGHTNING
!   PRODUCTION OF NO AND O2
  ! changel = 1     !mc temp var for testing lightning changes versus OLD JFK method
                 !changeL=1 uses new code, changeL=0 uses old code
  jt1 = Jtrop + 1     ! same as NH1
  do j = 1 , Jtrop
    Fval(lh2o,j) = 0.
    YP(lh2o,j) = 0.
    YL(lh2o,j) = 0.
  enddo
  if (lightning) then
    do j = 1,jtrop
      ! PN2, PCO2, PCO, PH2, PH2O, PO2, PO, PNO
      if (background_spec == "N2") then
        ! pass
      else
        zap = (prod_rates(1)/z(jtrop))
        Fval(lN2,j) = Fval(lN2,j) + zap/DEN(j)
        YP(lN2,j) = YP(lN2,j) + zap
      endif
      if (background_spec == "CO2") then
        ! pass
      else
        zap = (prod_rates(2)/z(jtrop))
        Fval(lCO2,j) = Fval(lCO2,j) + zap/DEN(j)
        YP(lCO2,j) = YP(lCO2,j) + zap
      endif
      zap = (prod_rates(3)/z(jtrop))
      Fval(lCO,j) = Fval(lCO,j) + zap/DEN(j)
      YP(lCO,j) = YP(lCO,j) + zap
      if (background_spec == "H2") then
        ! pass
      else
        zap = (prod_rates(4)/z(jtrop))
        Fval(lH2,j) = Fval(lH2,j) + zap/DEN(j)
        YP(lH2,j) = YP(lH2,j) + zap
      endif
      zap = (prod_rates(5)/z(jtrop))
      Fval(lH2O,j) = Fval(lH2O,j) + zap/DEN(j)
      YP(lH2O,j) = YP(lH2O,j) + zap
      
      zap = (prod_rates(6)/z(jtrop))
      Fval(lO2,j) = Fval(lO2,j) + zap/DEN(j)
      YP(lO2,j) = YP(lO2,j) + zap
      
      zap = (prod_rates(7)/z(jtrop))
      Fval(lO,j) = Fval(lO,j) + zap/DEN(j)
      YP(lO,j) = YP(lO,j) + zap
      
      zap = (prod_rates(8)/z(jtrop))
      Fval(lNO,j) = Fval(lNO,j) + zap/DEN(j)
      YP(lNO,j) = YP(lNO,j) + zap
      
      ! scale = RAIN(j)/RAIN(1)
      ! zap = zapno*scale
      ! Fval(lno,j) = Fval(lno,j) + zap/DEN(j)
      ! YP(lno,j) = YP(lno,j) + zap
      ! IF ( changel.EQ.1 ) THEN
    !-mc 4/28/06      making NO requires subtracting 1/2 O2   (1/2 N2 + 1/2 O2 <-> NO)
      ! Fval(lo2,j) = Fval(lo2,j) - 0.5*zap/DEN(j)
      ! YP(lo2,j) = YP(lo2,j) - 0.5*zap
      ! ELSE
    !-mc   4/28/06 the following is no longer needed as the code computes how much of each reductant is produced
    !-mc   i use these numbers to compute the amount of O2 produced. This is fine because CO2 and H2O reservoirs are effectivly infinte
    !-mc  "un-commenting" these for test versus JFK's original code.  in the else statement
      ! zap = zapo2*scale
      ! Fval(lo2,j) = Fval(lo2,j) + zap/DEN(j)
      ! YP(lo2,j) = YP(lo2,j) + zap
      ! ENDIF
      
      ! IF ( changel.EQ.1 ) THEN
    ! - mc adding lightning CO/H2 into chemistry for better redox conservation
    !- kevin's addition of 3-20-06
      ! zap = zaph2*scale
      ! Fval(lh2,j) = Fval(lh2,j) + zap/DEN(j)
      ! YP(lh2,j) = YP(lh2,j) + zap
  !-mc 4/28/06 adding H2 also requires adding 1/2 O2    (H2O <-> H2 + 1/2 O2)
      ! Fval(lo2,j) = Fval(lo2,j) + 0.5*zap/DEN(j)
      ! YP(lo2,j) = YP(lo2,j) + 0.5*zap
  !-mc
      ! zap = zapco*scale
      ! Fval(lco,j) = Fval(lco,j) + zap/DEN(j)
      ! YP(lco,j) = YP(lco,j) + zap
  !-mc 4/28/06  adding CO also requires adding 1/2 O2
      ! Fval(lo2,j) = Fval(lo2,j) + 0.5*zap/DEN(j)
      ! YP(lo2,j) = YP(lo2,j) + 0.5*zap
  !-mc now figure out how much O is produced: add to O, subtract 1/2 O2  (1/2 O2 <-> O)
      ! zap = zapo*scale
      ! Fval(lo,j) = Fval(lo,j) + zap/DEN(j)
      ! YP(lo,j) = YP(lo,j) + zap
      ! Fval(lo2,j) = Fval(lo2,j) - 0.5*zap/DEN(j)
      ! YP(lo2,j) = YP(lo2,j) - 0.5*zap
      ! ENDIF

    !        stop
    !-end 3-20-06 addition
    ENDDO
  endif
! ACK - this may be part of the reason the time-dependent code is having trouble
!   H2O CONDENSATION IN THE STRATOSPHERE
!   (RHCOLD IS THE ASSUMED RELATIVE HUMIDITY AT THE COLD TRAP)
! dunno what to do here, I'll take it to be small
  rhcold = 0.1d0
  IF ( planet.EQ.'EARTH' ) THEN
    rhcold = 0.1
! Jim had 0.1 ; what needs to be here is something that will give the right stratospheric H2O 3ppm
  ELSEIF ( planet.EQ.'MARS' ) THEN
!      RHCOLD = 0.4  ! Jim had 0.1 ?  my standard is 0.4  <-Kevin words (mc - this seems wrong)
    rhcold = 0.17
 ! from Kevin's Mars paper
!       RHCOLD = 0.10  ! Jim had 0.1 ?  my standard is 0.4  <-Kevin words (mc - this seems wrong)
  ENDIF
  if (H2O_strat_condensation) then
    DO j = jt1 , nz
      h2ocrt = rhcold*H2OSAT(j)
      IF ( USOL(lh2o,j).GE.h2ocrt ) THEN
        ! CONDEN(j) = confac*(USOL(lh2o,j)-h2ocrt)
        CONDEN1 = confac*(USOL(lh2o,j)-h2ocrt)
                                               !this is saved in SATBLK to be printed out in output file
        Fval(lh2o,j) = Fval(lh2o,j) - CONDEN1
      ENDIF
    ENDDO
  endif
  
!   H2SO4 CONDENSATION
  ll = lh2so4
  lla = lso4aer
  DO j = 1 , nz
    CONSO4(j) = confac*(USOL(ll,j)-H2SO4S(j))
    IF ( CONSO4(j).GT.0 ) THEN
      ! This exponent is the reduce stiffness, and increase stability. 
      ! An "if" statement discontinuity is very stiff!
      CONSO4(j) = CONSO4(j)*(-exp(-CONSO4(j)/1.d-20) + 1)
                        !dont allow artifical evaporation
      Fval(ll,j) = Fval(ll,j) - CONSO4(j)
      Fval(lla,j) = Fval(lla,j) + CONSO4(j)
      YL(ll,j) = YL(ll,j) + confac
      YP(ll,j) = YP(ll,j) + confac*H2SO4S(j)*DEN(j)
      YP(lla,j) = YP(lla,j) + CONSO4(j)*DEN(j)
    ENDIF
  ENDDO
!    !what follows is not in Jim's code
!   S8 CONDENSATION (this is needed if we every want to deal with 'hot air' - s8 stays in the vapor phase
         ! skips8 = 1
!          IF ( skips8.EQ.0 ) THEN
!
!             ll = ls8
!             lla = ls8aer
!
!             confc2 = 1.E-2
!                        !whats this?  CONFAC = 1.6E-5   whatever that is
!
!             DO j = 1 , nz
!                   !(fixes an error in our earlier codes where S8 didn't condense in the troposphere)
! !      DO 15 J=JT1,NZ  !(fixes an error in our earlier codes where S8 didn't condense in the troposphere) (temp return..)
! !      EVAPS8(J) = 0.!what is this? it appears to do nothing,and is printed out as 0 in the output file...
!
!                CONS8(j) = confc2*(USOL(ll,j)-S8S(j))
!
!                IF ( CONS8(j).GT.0. ) THEN
!                                !dont allow artifical evaporation
!
!                   Fval(ll,j) = Fval(ll,j) - CONS8(j)
!                   Fval(lla,j) = Fval(lla,j) + CONS8(j)
!                                                            !else do RHS in tri-diag
!                   YL(ll,j) = YL(ll,j) + confc2
!                   YP(ll,j) = YP(ll,j) + confc2*S8S(j)*DEN(j)
!                   YP(lla,j) = YP(lla,j) + CONS8(j)*DEN(j)
!
! !       print *,j,CONS8(J),CONFC2*S8S(J)*DEN(J),CONS8(J)*DEN(J)
! !       print *, j, USOL(LL,J),S8S(J),(USOL(LL,J) - S8S(J))
!
!                ENDIF
!             ENDDO
!
!          ENDIF
             !end S8 skip loop
  ENDIF  !end normal operation loop (i.e. if IDO = 0 or 1)
  
end subroutine

