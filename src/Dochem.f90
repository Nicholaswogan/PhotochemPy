
      SUBROUTINE DOCHEM(Fval,N,Jtrop,Nshort,Usol,nq,nz)

      IMPLICIT NONE

      ! module variables
      ! SL

      ! local variables
      integer, intent(in) :: n
      integer, intent(in) :: jtrop
      integer, intent(in) :: nshort
      integer, intent(in) :: nq
      integer, intent(in) :: nz
      real*8, dimension(nq,nz),intent(in) :: usol
      real*8, dimension(nq,nz),intent(out) :: Fval
      ! real*8, dimension(nq,nz) :: YP, YL

      integer i, j, iss4sl, ll, lla
      real*8, dimension(nsp2,nz) :: D
      real*8 xp(nz), xl(nz), conso4(nz)
      real*8 aq, bq, cq, dls4, xlj, confac, scale
      integer jt1, changel
      real*8 rhcold, h2ocrt, zap, CONDEN1




      ! TAUo2 , TAUch4 , TAUso2

!
!   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.  PHOTO-
!   CHEMICAL EQUILIBRIUM SPECIES ARE DONE FIRST.  THESE MUST CON-
!   TAIN NO NONLINEARITIES (SUCH AS S4 REACTING WITH ITSELF TO FORM
!   S8) AND MUST BE DONE IN THE PROPER ORDER (I.E. IF SPECIES A
!   REACTS TO FORM B, THEN A MUST BE FOUND FIRST).  LONG-LIVED
!   SPECIES CAN BE DONE IN ANY ORDER.

      do i=1,nz
        do j=1,nsp2
          D(j,i) = 0.d0
        enddo
      enddo
      iss4sl = 0




!compute number densities for long-lived and particle species
!these are re-computed below...
      DO i = 1 , nq1
         DO j = 1 , nz

!         if(j.eq.1)print *, i, ISPEC(I),' long-lived densities'
            d(i,j) = USOL(i,j)*DEN(j)
!          print  *, usol(1,j)


         ENDDO
      ENDDO

!now do the last 4 inert species that are the same in both codes
      DO j = 1 , nz



        if (CO2_inert.eq.1) D(LCO2,j) = FCO2*Den(j) !If CO2 is inert, then put it in D
        if (N2_inert.eq.1) D(NSP,j) = (1. - USOL(LO2,J) - FCO2 - USOL(LCO,J))* DEN(J) !if N2 is inert, then it is rest of atmosphere.



         d(nsp2-1,j) = 1.         ! HV has density of 1 for photorate calculations
         d(nsp2,j) = DEN(j)     ! M - background density for three body reactions
      ENDDO



      IF ( N.GE.0 ) THEN
                        !normal operation mode (-1 just fills up D and SL for first timestep)
!
! ***** SOLVE FOR THE PHOTOCHEMICAL EQUILIBRIUM SPECIES *****
!

         ! sls4 = 0
              !checking if S4 is the the short-lived loop

         DO i = nq1 + 1 , nq1 + Nshort
                               !loop through the equilibrium species...
!         print *, I,ISPEC(I),' short-lived'
            IF ( ISPEC(i).EQ.'S4' ) iss4sl = i

            CALL CHEMPL(d,xp,xl,i)

            DO j = 1 , nz
               d(i,j) = xp(j)/xl(j)
            ENDDO
         ENDDO




!-mab Uncomment below to help with short-lived species XP/XL/D debugging
!      print*,'NQ1+1 = ',NQ1+1
!      DO J=1,NZ
!        IF(J.EQ.1) print*,"Densities from XP/XL for species no.",NQ1+1
!        IF(J.EQ.1) print*,"(printing within chempl....)"
!        print*,"Layer no.,D(I,J),XP,XL",J,D(NQ1+1,J),XP(J),XL(J)
!        print*,'---'
!	  ENDDO
!!!!!!!!!!!!!!!!!!!
!      pause 5

!   SOLVE QUADRATIC FOR S4, if S4 is in the SL lived loop
! equation is production=loss, which turns into a quadratic equation in S4 density
! production (C) is from S+S3 and S2+S2
! losses are S4+S4 ->S8AER (A) and S4 + Hv -> S2+S2 (B)
! so we have C=A*DenS4^2 + B*DENS4*DENHv, where DENHv=1 by definition.

         IF ( iss4sl.GT.0 ) THEN


! add a loop over all reactions here, just dipping into to do loop if S4 is invloved
!an use the chempl stuff, but need to figure out K, the species number of S4
            ! nps4 = NUMP(iss4sl)
            ! nls4 = NUML(iss4sl)
            DO j = 1 , nz

               aq = 2.*A(148,j)
                            !S4+S4 -> S8Aer loss term
               bq = A(149,j)
                            ! S4+ Hv -> S2 + S2 loss term
               cq = A(146,j)*d(ls2,j)*d(ls2,j) + A(147,j)*d(ls,j)*d(ls3,j)                                         ! production terms

!        CQ=0.0

!ack need to abstract the above - below is a start, altough something is wrong.  check carefully against code in CHEMPL and shawn's code
!the problem is the N - it is overwriting one of the inputs... anyway.

!ACK - also down at the bottom of the file is an s8col printout that is suppressed.  Also a diagnostic printout in out.so2 in output

        !do the production terms generically
!        do I=1,nps4
!          L=IPROD(ISS4SL,I)   !reaction number
!          M = JCHEM(1,L)   !reactant 1 for reaction number J
!          N = JCHEM(2,L)   !reactant 2 for reaction number J
!          CQ=CQ + A(L,J)*D(M,J)*D(N,J) !rate*density1*density2
!        enddo

               dls4 = (SQRT(bq*bq+4.*aq*cq)-bq)/(2.*aq)
               d(ls4,j) = MAX(dls4,1.D-99)
            ENDDO

         ENDIF



!
! ***** LONG-LIVED SPECIES CHEMISTRY *****
         DO i = 1 , nq
!         print *, I,ISPEC(I),' long-lived'
            CALL CHEMPL(d,xp,xl,i)
            !
            ! IF ( ISPEC(i).EQ.'O2' ) TAUo2 = 1/xl(1)
            ! IF ( ISPEC(i).EQ.'CH4' ) TAUch4 = 1/xl(1)
            ! IF ( ISPEC(i).EQ.'SO2' ) TAUso2 = 1/xl(1)


            DO j = 1 , nz
               xlj = xl(j) + RAINGC(i,j)
               Fval(i,j) = xp(j)/DEN(j) - xlj*USOL(i,j)

!      print *, usol(i,j)
               YP(i,j) = xp(j)
               YL(i,j) = xlj
!      IF (ISPEC(I).EQ.'H2SO4') print *, J,XL(J),RAINGC(I,J),XP(j),XLJ
            ENDDO

         ENDDO


!          IF ( Usetd.EQ.1 ) THEN
! ! ***** TRIDIAGONAL SPECIES  *****
!             DO i = nq + 1 , nq1
! !         print *, I,ISPEC(I),' tri-diag'
!                CALL CHEMPL(d,xp,xl,i)
!                DO j = 1 , nz
! !        YL(I,J) = XL(J) + RAINGC(LH2SO4,J)   !this seems wrong for S8 but is status quo
! !        YL(I,J) = XL(J) + RAINGC(LSO4AER,J)   !this seems wrong for S8 but is status quo
!
!
! ! - these two below are the behavior that I am abstracting, but may be wrong
! !         if (ISPEC(I).EQ.'SO4AER') YL(I,J) = XL(J) + RAINGC(LSO4AER,J)
! !         if (ISPEC(I).EQ.'S8AER') YL(I,J) = XL(J) + RAINGC(LS8AER,J)
! ! - Jim/Kevin's code originally had the behavoir where H2SO4 was used for all particles.  Is this the source of the sulfur redox balance issues?
!
! !-mc - OK, if I come back to this later I should remember two things here:
! !- there is a mix of two different codes here. In one, I took RAINGC , HEFF, and some of the others to NQ+NP, and computed them directly.  If we are going to do the tridiagonal, they should be removed completely and I should go back to the original behavior, which is jsut uing H2SO4 everywhere.  At the same time, Kevin had S8 not being rained out at all because it isn't soluble. Either way, what we have in this section is both wrong and confusing...
!
!
! !         YL(I,J) = XL(J) + RAINGC(I,J)
!                   YL(i,j) = xl(j) + RAINGC(lh2so4,j)
!                   YP(i,j) = xp(j)
! !         if (ISPEC(I).EQ.'SO4AER') YL(I,J) = YL(I,J) + RAINGC(LH2SO4,J)
!                ENDDO
!             ENDDO
!
!          ENDIF


         IF ( planet.EQ.'EARTH' ) THEN
            confac = 1.6D-5     !condensation factor
         ELSEIF ( planet.EQ.'MARS' ) THEN
            confac = 1.6D-5*10. ! reduce supersaturation of stratosphere
         ELSEIF ( planet.EQ.'DRY' ) THEN
                                       !added by EWS 9/14/2015
            confac = 1.6D-5*10.  ! reduce supersaturation of stratosphere
         ENDIF



!   ZERO OUT H2O TERMS IN THE TROPOSPHERE AND INCLUDE LIGHTNING
!   PRODUCTION OF NO AND O2

         changel = 1     !mc temp var for testing lightning changes versus OLD JFK method
                         !changeL=1 uses new code, changeL=0 uses old code

!-mab: Let's NOT zero-out the H2O terms or include lightning for giant templates...
!-mab: (Basing this one on FH2-based distinction...) !mc - could we use the PLANET variable here?
         ! IF ( fh2.LT.0.50 ) THEN
          jt1 = Jtrop + 1     ! same as NH1

          DO j = 1 , Jtrop
             Fval(lh2o,j) = 0.
             YP(lh2o,j) = 0.
             YL(lh2o,j) = 0.


             scale = RAIN(j)/RAIN(1)

             zap = zapno*scale
             Fval(lno,j) = Fval(lno,j) + zap/DEN(j)
             YP(lno,j) = YP(lno,j) + zap


             IF ( changel.EQ.1 ) THEN
  !-mc 4/28/06      making NO requires subtracting 1/2 O2   (1/2 N2 + 1/2 O2 <-> NO)
                Fval(lo2,j) = Fval(lo2,j) - 0.5*zap/DEN(j)
                YP(lo2,j) = YP(lo2,j) - 0.5*zap
             ELSE

  !-mc   4/28/06 the following is no longer needed as the code computes how much of each reductant is produced
  !-mc   i use these numbers to compute the amount of O2 produced. This is fine because CO2 and H2O reservoirs are effectivly infinte
  !-mc  "un-commenting" these for test versus JFK's original code.  in the else statement
                zap = zapo2*scale
                Fval(lo2,j) = Fval(lo2,j) + zap/DEN(j)
                YP(lo2,j) = YP(lo2,j) + zap
             ENDIF




             IF ( changel.EQ.1 ) THEN

  ! - mc adding lightning CO/H2 into chemistry for better redox conservation
  !- kevin's addition of 3-20-06
                zap = zaph2*scale
                Fval(lh2,j) = Fval(lh2,j) + zap/DEN(j)
                YP(lh2,j) = YP(lh2,j) + zap

  !-mc 4/28/06 adding H2 also requires adding 1/2 O2    (H2O <-> H2 + 1/2 O2)
                Fval(lo2,j) = Fval(lo2,j) + 0.5*zap/DEN(j)
                YP(lo2,j) = YP(lo2,j) + 0.5*zap
  !-mc

                zap = zapco*scale
                Fval(lco,j) = Fval(lco,j) + zap/DEN(j)
                YP(lco,j) = YP(lco,j) + zap

  !-mc 4/28/06  adding CO also requires adding 1/2 O2
                Fval(lo2,j) = Fval(lo2,j) + 0.5*zap/DEN(j)
                YP(lo2,j) = YP(lo2,j) + 0.5*zap
  !-mc


  !-mc now figure out how much O is produced: add to O, subtract 1/2 O2  (1/2 O2 <-> O)
                zap = zapo*scale
                Fval(lo,j) = Fval(lo,j) + zap/DEN(j)
                YP(lo,j) = YP(lo,j) + zap
                Fval(lo2,j) = Fval(lo2,j) - 0.5*zap/DEN(j)
                YP(lo2,j) = YP(lo2,j) - 0.5*zap
  !-mc

             ENDIF

  !        stop
  !-end 3-20-06 addition

          ENDDO
         ! ENDIF
              !end hot jupiter skip loop



! ACK - this may be part of the reason the time-dependent code is having trouble
!
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
         DO j = jt1 , nz
            h2ocrt = rhcold*H2OSAT(j)
            IF ( USOL(lh2o,j).GE.h2ocrt ) THEN
               ! CONDEN(j) = confac*(USOL(lh2o,j)-h2ocrt)
               CONDEN1 = confac*(USOL(lh2o,j)-h2ocrt)
                                                       !this is saved in SATBLK to be printed out in output file
               Fval(lh2o,j) = Fval(lh2o,j) - CONDEN1
            ENDIF
         ENDDO


!
!   H2SO4 CONDENSATION

         ll = lh2so4
         lla = lso4aer

         DO j = 1 , nz
            CONSO4(j) = confac*(USOL(ll,j)-H2SO4S(j))

            IF ( CONSO4(j).GT.0 ) THEN
                                !dont allow artifical evaporation
               Fval(ll,j) = Fval(ll,j) - CONSO4(j)
               Fval(lla,j) = Fval(lla,j) + CONSO4(j)
                                                           !else handled in RHS of tri-diag in main code
               YL(ll,j) = YL(ll,j) + confac
               YP(ll,j) = YP(ll,j) + confac*H2SO4S(j)*DEN(j)
               YP(lla,j) = YP(lla,j) + CONSO4(j)*DEN(j)
!      print *,j,CONSO4(J),CONFAC*H2SO4S(J)*DEN(J),CONSO4(J)*DEN(J)
!      print *,j,USOL(LL,J),H2SO4S(J),USOL(LL,J) - H2SO4S(J)
            ENDIF
!      endif

         ENDDO

!      stop

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


!
! what is this?  its hardwired for case where O3 is trace
!     DO 7 J=1,NZ
!  7  O3(J) = D(LO3,J)/DEN(J)


! ***** SAVE THESE DENSITIES FOR PRINTOUT *****
!-mc and for allowing short-lived species to photolyze...?
!-mc modifying to contain all species.
! orig      DO 9 I=NQ1,NSP
      DO i = 1 , nsp
         DO j = 1 , nz
            SL(i,j) = d(i,j)
         ENDDO
      ENDDO



!       IF ( N.LT.1 ) RETURN           !on final time step compute TP and TL
!
! !
! ! ***** CALCULATE COLUMN-INTEGRATED PRODUCTION AND LOSS *****
!       o3col = 0.
!       h2scol = 0.
!       so2col = 0.
!       s2col = 0.
!       s4col = 0.
!       s8col = 0.
! !
!       DO l = 1 , nr
!          DO j = 1 , nz
!             REACRAT(l,j) = 0
!          ENDDO
!          RAT(l) = 0.
!       ENDDO
!
!
!
!       DO k = 1 , nq1   !num species
!          TP(k) = 0.
!          TL(k) = 0.
!       ENDDO
! !
!
!       IF ( fh2.LT.0.50 ) THEN
! !-mab: Disabling below for the giant planet template(s).
!          DO j = 1 , nz
!             RELH(j) = d(lh2o,j)/DEN(j)/H2OSAT(j)
!                                             !this gets H2O mixing ratio in a more general way
!             o3col = o3col + d(lo3,j)*DZ(j)
!
!             h2scol = h2scol + d(lh2s,j)*DZ(j)
!             so2col = so2col + d(lso2,j)*DZ(j)
!             s2col = s2col + d(ls2,j)*DZ(j)
!             s4col = s4col + d(ls4,j)*DZ(j)
! !       S8COL = S8COL + D(LS8,J)*DZ(J)   !ACK need to deal if S8 in gas phase
!
!
!          ENDDO
!       ENDIF   !end hot jupiter skip loop
! !
!       DO l = 1 , nr
!          m = JCHEM(1,l)      !identifies first reactant of equation L
!          k = JCHEM(2,l)      !second reactant
!          DO j = 1 , nz
!             REACRAT(l,j) = A(l,j)*d(m,j)*d(k,j)
!                                                !reaction rate*densities
! ! cm^3/mol/s * (mol/cm^3)^2 ->  mol/cm^3/s (i.e. rate units)
!
! !-mab:      IF(L.EQ.179) THEN For reaction rate debugging...
! !-mab: Uncomment below to track intermediate rates for chosen L
! !      	IF(J.EQ.1)print*, 'L =',L
! !      	print*, 'J,A,D(M),D(K)',J,A(L,J),D(M,J),D(K,J)
! !      ENDIF
!             RAT(l) = RAT(l) + REACRAT(l,j)*DZ(j)
!          ENDDO
!       ENDDO
! !      mol/cm^3/s * cm ->  mol/cm^2/s (i.e. RAT is in height integrated flux units)
!
! !-mab: Below is bebugging help...
! !-mab: Uncomment below to get the integrated rates at the end of each time step N
! !      print*,'Integrated Rates	J:'
! !      DO L=1,NR
! !      	print*, RAT(L),L
! !      ENDDO
!       DO i = 1 , nq1
!          XLG(i) = YL(i,1)
!          DO j = 1 , nz
!             TP(i) = TP(i) + YP(i,j)*DZ(j)
!             TL(i) = TL(i) + YL(i,j)*d(i,j)*DZ(j)
!          ENDDO
!       ENDDO
! !
      END
