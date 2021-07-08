
      SUBROUTINE RAINOUT(initialize,Jtrop,Usol,nq,nz, T,den, rain, raingc, err)
        use photochem_data, only: naq, lh2co, lh2s, lso2, lso4aer, &
                                  lh2so4, ispec, z, lco2, background_spec
        use photochem_wrk, only: H, xsave

      implicit none

      logical, intent(in) :: initialize
      integer, intent(in) :: jtrop
      integer, intent(in) :: nq
      integer, intent(in) :: nz
      real*8, dimension(nq,nz), intent(in) :: usol
      real(8), intent(in) :: T(nz), den(nz)
      
      ! out
      real(8), intent(out) :: rain(nz), raingc(nq,nz)
      character(len=err_len), intent(out) :: err

      ! local
      real*8 HEFF(NQ),IPVT(NAQ),DJAC(NAQ,NAQ)
      real*8 F(NAQ),FP(NAQ),TAQ(NZ),X(NAQ)
      integer LSO2g,LH2COg,LSO2aq,LH2COaq,LHCO3_,LCO3_2,LHSO3_
      integer LSO3_2,LH2COSO3,LOH_
      real*8 eps, GAM15, GAM8, AV, wl, R, zkm, xs, gamma
      integer NH, nh1, i, in, inewt, info, j, k, l1, l2, ltest
      real*8 ch2oh2, dx, fac, fhso3_, fz, h2cog
      real*8 hso3_, qj, rkj, so2aq, so2g, so3_2
      real*8 t_triple, temp, test, tfac, wh2o, y, fco2

      real*8, dimension(nz) :: HCO2
      real*8, dimension(nz) :: R4, R5, R6, R7, R8, R9
      real*8, dimension(6,nz) :: RRRR ! for passing Rs to aqueous
      real*8, dimension(nz) :: PH

      DATA lso2g , lh2cog , lso2aq , lh2coaq , lhco3_ , lco3_2 ,        &
         & lhso3_ , lso3_2 , lh2coso3 , loh_/1 , 2 , 3 , 4 , 5 , 6 , 7 ,&
         & 8 , 9 , 10/
      real(8) :: hplus
      real*8 alpharain, co2aq, h2cog0, hh2co
      real*8 hso2, so2g0, so4_2
      real(8) :: rain_parameters(7)
      err = ''



!  I am wanting to adjust rainout downwards for another planet
!  or for lower temperatures generally
!  Giorgi and Chameides use WH2O for the current rainout rate.
!   they quote 3.3e-6 g/cm2/s  as the global average when WH2O is
!  integrated over altitude
!
!      THIS SUBROUTINE CALCULATES RAINOUT RATES USING GIORGI AND
!   CHAMEIDES (1985) MODEL.  FIRST, IT CALCULATES THE NORMAL HENRY'S
!   LAW COEFFICIENTS.  THEN IT SOLVES A SYSTEM OF NAQ AQUEOUS PHASE
!   REACTIONS TO FIND EFFECTIVE HENRY'S LAW COEFFICIENTS.  THE
!   REACTIONS ARE:
!
!     1)  (SO2)g + ALPHARAIN*[(SO2)aq  +  HSO3-  +  SO3=  +  CH2OHSO3-]
!                  =  (SO2)go
!     2)  (SO2)g   =  (SO2)aq
!     3)  (H2CO)g  =  CH2(OH)2
!     4)  (CO2)aq  =  HCO3-  +  H+
!     5)  (SO2)aq  =  HSO3-  +  H+
!     6)   HCO3-   =  CO3=   +  H+
!     7)   HSO3-   =  SO3=   +  H+
!     8)  CH2(OH)2  +  HSO3-  =  H2O  +  CH2OHSO3-
!     9)   H2O     =  H+  +  OH-
!    10)  (H2CO)g + ALPHARAIN*[CH2(OH)2  +  CH2OHSO3-]  =  (H2CO)go

!   ALONG WITH
!         (CO2)g   =  (CO2)aq
!         (H2SO4)g =  2H+  +  SO4=
!         H+  =  OH-  +  HCO3-  +  HSO3-  +  CH2OHSO3-  +
!                2*[CO3=  +  SO3=]
!
!     THE VARIABLES IN THE NEWTON STEP ARE:
!     1)  X(1)  =  (SO2)g
!     2)  X(2)  =  (H2CO)g
!     3)  X(3)  =  (SO2)aq
!     4)  X(4)  =  CH2(OH)2
!     5)  X(5)  =  HCO3-
!     6)  X(6)  =  CO3=
!     7)  X(7)  =  HSO3-
!     8)  X(8)  =  SO3=
!     9)  X(9)  =  CH2OHSO3-
!    10)  X(10) =  OH-
!
!   FIRST DEFINE RELEVANT CONSTANTS\



      eps = 1.E-7         !Shawn has this as 1e-4
      inewt = 30
      gam15 = 8.64E+05/2.0
      gam8 = 7.0E+06/2.0
      av = 6.02E+23
      wl = 1.0
      r = 1.36E-22
      nh = Jtrop  !height index of troposphere is called NH in this subroutine...

      nh1 = nh + 1



!   MODIFY TEMPERATURE PROFILE TO DO AQUEOUS CHEMISTRY
      t_triple = 273.15
      DO i = 1 , nh
         taq(i) = MAX(T(i),t_triple)
      ENDDO
!   CALCULATE EQUILIBRIUM CONSTANTS FOR AQUEOUS PHASE REACTIONS
      DO i = 1 , nh
!   4)  (CO2)aq  =  HCO3-  +  H+
        R4(i) = 4.3E-7*EXP(-913.*(1./taq(i)-1./298.))
        RRRR(1,i) = R4(i)
!   5)  (SO2)aq  =  HSO3-  +  H+
        R5(i) = 1.7E-2*EXP(2090.*(1./taq(i)-1./298.))
        RRRR(2,i) = R5(i)
!   6)   HCO3-   =  CO3=   +  H+
        R6(i) = 5.6E-11
        RRRR(3,i) = R6(i)
!   7)   HSO3-   =  SO3=   +  H+
        R7(i) = 6.E-8*EXP(1120.*(1./taq(i)-1./298.))
        RRRR(4,i) = R7(i)
!   8)  CH2(OH)2  +  HSO3-  =  H2O  +  CH2OHSO3-
        R8(i) = 1.E5
        RRRR(5,i) = R8(i)
!   9)   H2O     =  H+  +  OH-
        R9(i) = 1.E-14*EXP(-6716.*(1./taq(i)-1./298.))
        RRRR(6,i) = R9(i)
      ENDDO
      DO i = 1 , nh
         tfac = (1./taq(i)-1./298.)
         HCO2(i) = 3.5E-2*EXP(2400.*tfac)
      enddo


!
!   CALCULATE NORMAL HENRY'S LAW COEFFICIENTS (PHYSICAL DISSOLUTION ONLY)
!   lets play guess the units.  looks like mol/liter/atm
      IF (initialize) THEN

         DO i = 1 , nh
            tfac = (1./taq(i)-1./298.)
            HCO2(i) = 3.5E-2*EXP(2400.*tfac)
                                          ! updated


            DO j = 1 , nq
               IF ( ISPEC(j).EQ.'O' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'H2O' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'H' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'HCO' ) H(j,i) = 0.
                                         ! I think jim made this up
               IF ( ISPEC(j).EQ.'CH3' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'S' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'S2' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'S4' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'S8' ) H(j,i) = 0.
      !infinite (CHECK THIS against PAVLOV/KEVIN)
               IF ( ISPEC(j).EQ.'S8AER' ) H(j,i) = 7.E11
               IF ( ISPEC(j).EQ.'HNO' ) H(j,i) = 7.E11
       ! I think jim made the below up
               IF ( ISPEC(j).EQ.'HS' ) H(j,i) = 1.E5
       ! I think jim made the below up
               IF ( ISPEC(j).EQ.'SO' ) H(j,i) = 1.9E-3
       ! updated, infinite
               IF ( ISPEC(j).EQ.'SO3' ) H(j,i) = 7.E11
               IF ( ISPEC(j).EQ.'H2SO4' ) H(j,i) = 7.E11
       ! I think jim made the below up
               IF ( ISPEC(j).EQ.'HSO' ) H(j,i) = 9.E3
       !infinite for particle species
               IF ( ISPEC(j).EQ.'SO4AER' ) H(j,i) = 7.E11
               IF ( ISPEC(j).EQ.'HCAER' ) H(j,i) = 7.E11
                                                !infinite for particle species
               IF ( ISPEC(j).EQ.'HCAER2' ) H(j,i) = 7.E11
                                                 !infinite for particle species
               IF ( ISPEC(j).EQ.'O2' ) H(j,i) = 1.3E-3*EXP(1500.*tfac)
                                                             ! updated
               IF ( ISPEC(j).EQ.'OH' ) H(j,i) = 30.*EXP(4500.*tfac)
                                                            ! updated
               IF ( ISPEC(j).EQ.'HO2' ) H(j,i) = 4.E3*EXP(5900.*tfac)
                                                              ! updated
               IF ( ISPEC(j).EQ.'H2O2' ) H(j,i) = 8.3E4*EXP(7400.*tfac)
                                                               ! updated
               IF ( ISPEC(j).EQ.'H2' ) H(j,i) = 7.8E-4*EXP(500.*tfac)
                                                              ! updated
               IF ( ISPEC(j).EQ.'CO' ) H(j,i) = 1.0E-3*EXP(1300.*tfac)
                                                               ! updated
               IF ( ISPEC(j).EQ.'H2CO' ) H(j,i) = 3.2E3*EXP(6800.*tfac)
                                                                ! updated
               IF ( ISPEC(j).EQ.'CH4' ) H(j,i) = 1.4E-3*EXP(1600.*tfac)
                                                                ! updated
               IF ( ISPEC(j).EQ.'C2H6' ) H(j,i) = 1.9E-3*EXP(2300.*tfac)
                                                                 ! updated
               IF ( ISPEC(j).EQ.'NO' ) H(j,i) = 1.9E-3*EXP(1500.*tfac)
                                                               ! updated
               IF ( ISPEC(j).EQ.'NO2' ) H(j,i) = 1.2E-2*EXP(2500.*tfac)
                                                                ! updated
               IF ( ISPEC(j).EQ.'HNO2' ) H(j,i) = 50.*EXP(4900.*tfac)
                                                              ! added, updated
               IF ( ISPEC(j).EQ.'HNO3' ) H(j,i) = 2.1E5*EXP(8700.*tfac)
                                                                ! added, updated
               IF ( ISPEC(j).EQ.'HO2NO2' ) H(j,i)                       &
                  & = 1.2E4*EXP(6900.*tfac)                      ! added, updated
               IF ( ISPEC(j).EQ.'H2S' ) H(j,i) = 0.1*EXP(2000.*tfac)
                               ! updated
               IF ( ISPEC(j).EQ.'SO2' ) H(j,i) = 1.4*EXP(2900.*tfac)
                                ! updated
               IF ( ISPEC(j).EQ.'OCS' ) H(j,i) = 0.022*EXP(2100.*tfac)
                                  ! updated
               IF ( ISPEC(j).EQ.'CH3SH' ) H(j,i) = 0.2*EXP(2800.*tfac)
               IF ( ISPEC(j).EQ.'C2H6S' ) H(j,i) = 0.48*EXP(3100.*tfac)
               IF ( ISPEC(j).EQ.'C2H6S2' ) H(j,i) = 0.96*EXP(4000.*tfac)
               IF ( ISPEC(j).EQ.'CS2' ) H(j,i) = 0.055*EXP(2800.*tfac)
! no info found for CS, CH3S, HCS. CO solubility is small, as is HCO, so 0's are probably OK. no info on CH3O either.

!none of the carbon species have any rainout terms in Shawn's model - how valid is this?
               IF ( ISPEC(j).EQ.'CH' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C2H2' ) H(j,i) = 0.

               IF ( ISPEC(j).EQ.'CH3O2' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'CH3O' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'CH2CO' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'CH3CO' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'CH3CHO' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C2H3' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C2H4' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C2H2OH' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C2H4OH' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C3H8' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C3H7' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C3H6' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C2H5HCO' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'C3H5' ) H(j,i) = 0.
               IF ( ISPEC(j).EQ.'CH2CCH2' ) H(j,i) = 0.

! added by Nick
!       if(ISPEC(J).EQ.'HNCO')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'NH3')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'NH2')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'HCN')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'NCO')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'CN')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'N2H3')   H(J,I) = 0.
!       if(ISPEC(J).EQ.'N2H4')   H(J,I) = 0.



!-mc should go through the henry table and look for chlorine as well as hazy species

            ENDDO
         ENDDO


!   NOW ESTIMATE INITIAL CONCENTRATIONS AT GRID STEP 1 ON THE
!     FIRST CALL
         ! get CO2 mixing ratio
         if (background_spec == "CO2") then
           fco2 = 1 - sum(usol(:,1))
         else
           fco2 = usol(lco2,1)
         endif

         alpharain = wl*1.E-9*6.02E23/DEN(1)
         co2aq = fco2*DEN(1)*HCO2(1)*r*taq(1)
         hplus = SQRT(R4(1)*co2aq)
         x(lhco3_) = R4(1)*co2aq/hplus     !5
         x(lco3_2) = R6(1)*x(lhco3_)/hplus !6
         x(loh_) = R9(1)/hplus             !10
!      print *, 'rain 4'

         so2g = USOL(lso2,1)
         so2aq = H(lso2,1)*so2g

         hso3_ = R5(1)*so2aq/hplus
         so3_2 = R7(1)*hso3_/hplus
         fac = so2g/(so2g+alpharain*(so2aq+hso3_+so3_2))

         x(lso2g) = so2g*fac        !1
         x(lso2aq) = so2aq*fac      !3
         x(lhso3_) = hso3_*fac      !7
         x(lso3_2) = so3_2*fac      !8
!

         h2cog = USOL(lh2co,1)
         ch2oh2 = H(lh2co,1)*h2cog

         fhso3_ = R8(1)*ch2oh2*x(lhso3_)
         fac = h2cog/(h2cog+alpharain*(ch2oh2+fhso3_))
         x(lh2cog) = h2cog*fac   !2
         x(lh2coaq) = ch2oh2*fac !4
         x(lh2coso3) = fhso3_*fac
                                 !9

      ENDIF


!
! ***** LOOP OVER ALTITUDE *****
      DO i = 1 , nh           !this is a big loop (NH is the tropopause height index JTROP elsewhere)

         IF (.not. initialize) THEN


            DO k = 1 , naq
               x(k) = XSAVE(k,i)
            ENDDO
         ENDIF
         
         if (background_spec == "CO2") then
           fco2 = 1 - sum(usol(:,i))
         else
           fco2 = usol(lco2,i)
         endif

         alpharain = wl*1.E-9*6.02E23/DEN(i)

         ! IF ( Usetd.EQ.1 ) THEN
         ! so4_2 = (USOL(lh2so4,i)+PARTICLES(i,lso4aer-nq))/alpharain
         ! ELSE
          so4_2 = (USOL(lh2so4,i)+USOL(lso4aer,i))/alpharain ! this will not work
         ! ENDIF

         ! SO4SAV(i) = so4_2
         co2aq = fco2*DEN(i)*HCO2(i)*r*taq(i)

         so2g0 = USOL(lso2,i)
         h2cog0 = USOL(lh2co,i)
         hso2 = H(lso2,i)
         hh2co = H(lh2co,i)
         rain_parameters = [alpharain, co2aq, h2cog0, hh2co, &
                            hso2, so2g0, so4_2]




!
!   START NEWTON ITERATION

         DO in = 1 , inewt
            ! CALL AQUEOUS(x,f,i)
            CALL AQUEOUS(x,naq,f,i,RRRR,nz, rain_parameters, hplus) ! changed input a little
!
            DO j = 1 , naq
               xs = x(j)
               dx = eps*x(j)
               x(j) = x(j) + dx
               CALL AQUEOUS(x,naq,fp,i,RRRR,nz, rain_parameters, hplus)
!
               DO k = 1 , naq
                  djac(k,j) = (fp(k)-f(k))/dx
               ENDDO
               x(j) = xs
            ENDDO
!

            CALL SGEFA(djac,naq,naq,ipvt,info)
            IF ( info.NE.0 ) THEN
!                PRINT 99001 , info , i
! 99001          FORMAT (//1X,'NEWTON SOLVER FAILED IN AQUEOUS'/,5X,      &
!                       &'INFO =',I3,'  GRID STEP =',I3)
               err = 'Newton solve failed in rainout calculation.'
               return
            ELSE
               CALL SGESL(djac,naq,naq,ipvt,f,0)
            ENDIF
!

            ltest = 0
            DO j = 1 , naq
               x(j) = x(j) - f(j)
               test = ABS(f(j)/x(j))
               IF ( test.GT.1.E-5 ) ltest = 1
            ENDDO
            IF ( ltest.EQ.0 ) GOTO 50
                                  ! loop has converged; jump out of loop
         ENDDO
!
         ! PRINT 99002 , i 
! 99002    FORMAT (//1X,'NEWTON SOLVER FAILED TO CONVERGE IN AQUEOUS'/,5X,&
                ! &'GRID STEP =',I3)
         ! STOP
         err = 'Newton solve failed to converge in rainout calculation.'
         return

!
!   CALCULATE EFFECTIVE HENRY'S LAW COEFFICIENTS (INCLUDING AQUEOUS
!      PHASE REACTIONS)
 50      DO j = 1 , nq
            heff(j) = H(j,i) + 1.E-30
         ENDDO
!
         heff(lh2co) = (x(lh2coaq)+x(lh2coso3))/x(lh2cog)

         l1 = lso2
         l2 = lh2s

         heff(l1) = (x(lso2aq)+x(lhso3_)+x(lso3_2)+x(lh2coso3))/x(lso2g)
         heff(l2) = H(l2,i)*(1.+1.1E-7/hplus)
         PH(i) = -LOG10(hplus)

!
!   SAVE DENSITIES AND CALCULATE ENHANCEMENTS

         DO j = 1 , naq
            XSAVE(j,i) = x(j)
         ENDDO

!   NOW BEGIN GIORGI AND CHAMEIDES FORMULATION FOR RAINOUT RATES
         zkm = Z(i)/1.E5
                        !convert vertical grid into kilometers

! this mod for high CO2
!      eleven = 11.0
!      ZKM = min(ZKM,eleven)
         temp = T(i)
!
!  Find appropriate GAMMA
         IF ( zkm.LE.1.51 ) THEN
            gamma = gam15
         ELSEIF ( zkm.LT.8. ) THEN
            gamma = gam15 + (gam8-gam15)*((zkm-1.5)/6.5)
         ELSE
            gamma = gam8
         ENDIF
!
!  Find WH2O
         IF ( zkm.LE.1. ) THEN
            y = 11.35 + 0.1*zkm
         ELSE
            y = 11.5444 - 0.085333*zkm - 9.1111E-03*zkm*zkm
         ENDIF
         wh2o = 10.**y

         ! IF ( planet.EQ.'EARTH' ) THEN
            ! wh2o = 1.0*wh2o
                          ! nominal earthly rain
!       wh2o=1e-9*wh2o  !ATACAMA
         ! ELSEIF ( planet.EQ.'MARS' ) THEN
            ! wh2o = 1E-9*wh2o
                          !turn off the rain
         ! ELSEIF ( planet.EQ.'DRY' ) THEN
            ! wh2o = 1E-9*wh2o
                          !turn off rain for dry planet - EWS 9/14/2015
         ! ENDIF

!  Find F(Z)
         IF ( zkm.LE.1.51 ) THEN
            fz = 0.1
         ELSE
            fz = 0.16615 - 0.04916*zkm + 3.37451E-3*zkm*zkm
         ENDIF
!
!  Loop over species
         DO j = 1 , nq
            rkj = wh2o/55./(av*wl*1.E-9+1./(heff(j)*r*temp))
            qj = 1. - fz + fz/(gamma*rkj)*(1.0-EXP(-rkj*gamma))
            RAINGC(j,i) = (1.-EXP(-rkj*gamma))/(gamma*qj)
         ENDDO
      ENDDO
! ***** END ALTITUDE LOOP *****

!
!-mc set rainout rates to zero above the tropopause (nh=jtrop)

      DO i = nh1 , nz
         DO j = 1 , nq
            RAINGC(j,i) = 0.d0
         ENDDO
      ENDDO


!
! ***** OLD (FISHMAN AND CRUTZEN) RAINOUT RATE *****
!   (USED FOR SCALING THE VERTICAL DISTRIBUTION OF LIGHTNING)
!   what is the 2.4e-6 constant?
      DO i = 1 , nh
         zkm = Z(i)/1.E5
         RAIN(i) = 2.4E-6*EXP((6.-zkm)/2.42)
         IF ( zkm.LT.6. ) RAIN(i) = 2.4E-6
      ENDDO

      DO i = nh1 , nz
         RAIN(i) = 0.
      ENDDO

      end subroutine


