
      SUBROUTINE SEDMNT(Frak,Hcdens,nz,np, conver,change_radius)
        use photochem_data, only: nq, g, lhcaer, lhcaer2, ls8aer, lso4aer
        use photochem_vars, only: T, den
        use photochem_wrk, only: rpar, wfall, aersol, raingc, hscale, tauedd, &
                                 fsulf

      IMPLICIT NONE


      ! module variables
      ! Fsulf
      ! aersol
      ! wfall
      ! rpar


      ! local variables
      integer, intent(in) :: nz, np, frak
      real*8, intent(in) :: hcdens
      real*8, dimension(nz,np), intent(out) :: conver
      logical, intent(in) :: change_radius

      real*8 a, b, c, bk, pi, adensity, e_minus_5
      real*8 alph, e_hc, f1, factor, r, rf, gpermolec
      real*8, dimension(nz,np) :: tauc, RFRAC, TAUSED
      integer i,j,k,l,ll, nz1


      real*8 tautrn(nz) , rhop(nz)
      real*8 tauran(nz,np) , alam(nz) , taucpk(nz,np) , eta(nz)


      real*8 cuning(nz,np) , amass(nz,np) , thermsp(nz,np) ,         &
              & taurelaxc(nz,np) , taurelax(nz,np) , afpl(nz,np) ,      &
              & delta(nz,np) , betaf(nz,np)

      real*8 nmon , rmon , df
                         !frachack
      factor = 0.d0

!
!
!   THIS SUBROUTINE CALCULATES FALL VELOCITIES AND ESTIMATES PARTICLE SIZE
!   BASED ON THEIR COAGULATION LIFETIMES
!
!   CONSTANTS FOR STOKES-CUNNINGHAM EQUATION (KASTEN, 1968)
!-AP      A = 1.249
!-AP      B = 0.42
!-AP      C = 0.87
!   CONSTANTS from Pruppacher and Klett page 450
      a = 1.257
      b = 0.4
      c = 1.1
!-AP
!-AP       A = 0.866
!-AP       B = 0.29
!-AP       C = 1.25
!-AP
      bk = 1.38*1.D-16
      pi = 3.14159
      nz1 = nz - 1
!


      IF ( Frak.EQ.1 ) THEN
! **********************************************
!-EW  IMPLEMENTATION OF FRACTAL MICROPHYSICS
!-EW  ONLY FOR HYDROCARBONS, K=3 and 4
!
!-EW  IMPLEMENTATION USES A SIZE BIN DEPENDENT FRACTAL
!-EW  DIMENSION TO PARAMETERIZE AGGREGATE RESTRUCTURING
!-EW  PERMEABILITY EFFECTS ARE NOT TREATED HERE
!
!-EW RPAR = equal mass spherical radii
!-EW RFRAC = fractal aggregate radii
!
!-EW  THIS IS THE MONOMER RADIUS [cm]
!       if (ihztype.eq.0.) then RMON = 50.E-7
!       if (ihztype.eq.1.) then RMON = 10.E-7
!       if (ihztype.eq.2.) then RMON = 20.E-7
!       if (ihztype.eq.3.) then RMON = 70.E-7
!       if (ihztype.eq.4.) then RMON = 10.E-6
!       if (ihztype.eq.5.) then RMON = 50.E-7
!       if (ihztype.eq.6.) then RMON = 50.E-7

         rmon = 50.D-7
         DO k = 3 , 4
            DO j = 1 , nz
               nmon = (RPAR(j,k)/rmon)**3.
               IF ( nmon.LE.1. ) THEN
                  df = 3.
               ELSE
                  df = 2.4 - 0.9*EXP(-nmon/500.)
               ENDIF
               RFRAC(j,k) = RPAR(j,k)**(3./df)*rmon**(1.-3./df)
            ENDDO
         ENDDO


!-EW *******************************************
!
      ENDIF

      DO j = 1 , nz
         alam(j) = 1.63D14/DEN(j)
         eta(j) = ABS((1.718+0.0049*(T(j)-273.)-1.2*(1.D-5)*(T(j)-273.)*&
                & (T(j)-273.))*1.D-4)
      ENDDO


      DO k = 1 , np
!   (1 = SULFATE, 2 = S8, 3 = HYDROCARBON, 4=HCAER2)

         l = lso4aer
                    !so using the sulfate aersol for all particles?

                            !Jim's code with particles in tri-diag uses H2SO4 rather than the aerosol
         ! IF ( Usetd.EQ.1 ) l = lh2so4
                      !using h2so4 for all particles? (sure why not - it's infinite rainout...)

!
!-AP ESTIMATION OF THE AEROSOL FREE PATH LENGTH
         DO j = 1 , nz
            IF ( k.GE.3 .AND. Frak.EQ.1 ) THEN
                                       !ACK hardcoded particle number (getting HCAER and HCAER2)
!-EW  HYDROCARBONS: USE FRACTAL MICROPHYSICS
               alph = a + b*EXP(-c*RFRAC(j,k)/alam(j))
                                               !giada - related to particle diffusion I think
               cuning(j,k) = 1 + alph*alam(j)/RFRAC(j,k)
!-AP Here we assume that the density of aerosol is 1 g/cm3
!-Giada: actually, density of hc aerosols should be 0.63 g/cm3
!        see Trainer et al (2006)
!-AP Notation is similar Fusch 1964
! THERMSP = thermal velocity of molecule
!       adensity = 0.64
               amass(j,k) = (4./3.)*pi*RPAR(j,k)**3*Hcdens
               thermsp(j,k) = SQRT((8*bk*T(j))/(pi*amass(j,k)))
                                                        !giada - thermal velocity
               taurelaxc(j,k) = 2*RPAR(j,k)**3/(9*eta(j)*RFRAC(j,k))
               taurelax(j,k) = taurelaxc(j,k)*cuning(j,k)
               afpl(j,k) = thermsp(j,k)*taurelax(j,k)
                                              !giada - this is particle mean free path?
                                                        !giada - delta is mean distance from ctr of sphere
                                                        !reached by particle leaving the surface and traveling distance
               delta(j,k) = (((2*RFRAC(j,k)+afpl(j,k))**3-(4*RFRAC(j,k)*&
                          & RFRAC(j,k)+afpl(j,k)*afpl(j,k))**1.5)       &
                          & /(6*RFRAC(j,k)*afpl(j,k))-2*RFRAC(j,k))     &
                          & *SQRT(2.)                   !equal to the mean free path (AFPL)
            ELSE
!-EW  S8,SO4,HC if frak=0: USE SPHERICAL MICROPHYSICS
               alph = a + b*EXP(-c*RPAR(j,k)/alam(j))
               cuning(j,k) = 1 + alph*alam(j)/RPAR(j,k)
!-AP Here we assume that the density of aerosol is 1 g/cm3
!-AP Notation is similar Fusch 1964
               adensity = Hcdens
               amass(j,k) = (4./3.)*pi*RPAR(j,k)**3*adensity
               thermsp(j,k) = SQRT((8*bk*T(j))/(pi*amass(j,k)))
               taurelaxc(j,k) = 2*RPAR(j,k)*RPAR(j,k)/(9*eta(j))
               taurelax(j,k) = taurelaxc(j,k)*cuning(j,k)
               afpl(j,k) = thermsp(j,k)*taurelax(j,k)
               delta(j,k) = (((2*RPAR(j,k)+afpl(j,k))**3-(4*RPAR(j,k)*  &
                          & RPAR(j,k)+afpl(j,k)*afpl(j,k))**1.5)        &
                          & /(6*RPAR(j,k)*afpl(j,k))-2*RPAR(j,k))       &
                          & *SQRT(2.)
            ENDIF
            !end frak loop
         ENDDO

!-AP Calculation of the correction to the coagulation kernel
         DO j = 1 , nz
            IF ( k.GE.3 .AND. Frak.EQ.1 ) THEN
                                         !get HCAER and HCAER2
!-EW  HYDROCARBONS: USE FRACTAL MICROPHYSICS - giada: betaf is effective radii?
               betaf(j,k) = 1/(RFRAC(j,k)/(RFRAC(j,k)+delta(j,k)/2)     &
                          & +pi*afpl(j,k)/(2*SQRT(2.)*RFRAC(j,k)))
            ELSE
!-EW  S8,SO4,HC if frak=0: USE SPHERICAL MICROPHYSICS
               betaf(j,k) = 1/(RPAR(j,k)/(RPAR(j,k)+delta(j,k)/2)+pi*   &
                          & afpl(j,k)/(2*SQRT(2.)*RPAR(j,k)))
            ENDIF
         ENDDO
!-AP ******************************************************



!   ESTIMATE COAGULATION AND SEDIMENTATION LIFETIMES (TOON AND FARLOW, 1981
         DO i = 1 , nz
! this commented one was from Toon and Farlow 1981; giada-new methodology from Pavlov 2001
!     TAUC(I,K) = 1.E6/(AERSOL(I,K)*SQRT(RPAR(I,K)))      !T&F (81) p.41 - e-folding lifetime against coagulation
                                                          ! 1/tauc = 1/N DN/DT

! recompute tauc using new methodology (giada- from Pavlov 2001)
            taucpk(i,k) = 3*eta(i)/(4*AERSOL(i,k)*bk*T(i)*cuning(i,k))
            TAUC(i,k) = taucpk(i,k)/betaf(i,k)
                                        !ok for fractals as BETAF has new effective radii if frak=1

            tauran(i,k) = 1./(RAINGC(l,i)+1.D-20)
                                                 ! this is really wrong for S8 !!  !where does this come from?

!it looks like I never completly fixed the rainout thing.  TAURAN is using LSO4AER for both.  but isn't RAINGC really small for S8?

                                                 ! what is this? rainout lifetime, i presume. should check magnitude

            TAUSED(i,k) = HSCALE(i)/WFALL(i,k)
         ENDDO
!
!   FIND MINIMUM OF DIFFUSION AND SEDIMENTATION LIFETIMES, THEN SCALE PRTICLE SIZES
         e_minus_5 = 1.D-5
         e_hc = 1.3D-7
         if (change_radius) then
           DO i = 1 , nz
              tautrn(i) = MIN(TAUSED(i,k),TAUEDD(i))            !find the minimum of the three destruction timescales
              tautrn(i) = MIN(tautrn(i),tauran(i,k))            !lifetime against eddy diffusion is H*H/K
                                                                !where K is eddy diffusion coefficient, H is scale height
              RPAR(i,k) = RPAR(i,k)*(tautrn(i)/TAUC(i,k))**0.25 !particle growth depends on
              IF ( k.GE.3 ) THEN
                 RPAR(i,k) = MAX(RPAR(i,k),e_hc)
                                          !largest HC particles are smaller?
              ELSE
                 RPAR(i,k) = MAX(RPAR(i,k),e_minus_5)            !largest particles are 1 micron
              ENDIF
           ENDDO
         
!
!   DON'T ALLOW PARTICLES TO DECREASE IN SIZE AT LOW ALTITUDES
           DO i = 1 , nz1
              j = nz - i
              RPAR(j,k) = MAX(RPAR(j,k),RPAR(j+1,k))
           ENDDO
           
         endif
!
!   COMPUTE PARTICLE-TO-GAS CONVERSION FACTORS AND DENSITIES
! - it would nice to document where these come from - the 1989 sulfur/UV paper, I imagine.
!- no hints in Shawn's code
         DO i = 1 , nz
            r = RPAR(i,k)
            ll = nq - np + k
                    !ACK - assuming particles are last LL elements
            ! IF ( Usetd.EQ.1 ) ll = k + nq
                                 !particles in tri-diag

            IF ( ll.EQ.ls8aer ) THEN
               rhop(i) = 2.07
               factor = 2.03D7
            ENDIF

            IF ( ll.EQ.lso4aer ) THEN
               rhop(i) = 1. + 0.8*Fsulf(i)
               factor = 4.6D7*Fsulf(i)
               if (.not. change_radius) then ! 
                 rhop(i) = 1.6d0
                 factor = 4.6d7
               endif
            ENDIF

            IF ( ll.EQ.lhcaer .OR. ll.EQ.lhcaer2 ) THEN
!         RHOP(I) = 1.4
!         RHOP(I) = 1.0  !from feng's code
               rhop(i) = Hcdens
!         factor=7.06E7 !giada - this was computed using 1.4 g/cm3
               IF ( ll.EQ.lhcaer ) gpermolec = 2*1.66D-24 +             &
                  & 12*4*1.66D-24                 !grams per molecule for HCAER

               IF ( ll.EQ.lhcaer2 ) gpermolec = 4*1.66D-24 +            &
                  & 12*5*1.66D-24                 !grams per molecule for HCAER2
               factor = (4./3.)*pi*(1.E-5)**3*Hcdens/gpermolec

            ENDIF
            CONVER(i,k) = factor*(r/1.D-5)**3


! - giada - factor is the NUMBER OF MOLECULES PER 0.1UM SPHERE (don't know why .1um was chosen, but that is how it is...)

!rho_p is particle density (g/cm3)
!factor is related to density somehow (jim's words) (giada - yeah, this was all I had to go off of to figue
!out the mysterious 'factor'...)

!conver is the number of molecules/particle, so that main calcuations are done in molecule space
!conver is used in the main code to calculate aersol (= number density of aerosols)!


         ENDDO

!   NOW COMPUTE FALL VELOCITIES
         DO j = 1 , nz

            IF ( k.GE.3 .AND. Frak.EQ.1 ) THEN
!-EW  HYDROCARBONS: USE FRACTAL MICROPHYSICS
               r = RPAR(j,k)
               rf = RFRAC(j,k)
               f1 = 2./9.*rhop(j)*r*r*r*g/eta(j)/rf
                                              !from stokes law F1 - settling velocity
       !giada -  it's computing terminal velocity  when frictional and buoyant forces
       !are equal to gravitational force
               alph = a + b*EXP(-c*rf/alam(j))
                                       ! I think this is related to particle resistance to motion
               WFALL(j,k) = f1*(1.+alam(j)*alph/rf)
                                              !wfall = fall velocity
                             !this term (alam*alph/rf) is particle diffusion?  maybe?
            ELSE
!-EW  S8, SO4,HC if frak=0: USE SPHERICAL MICROPHYSICS
               r = RPAR(j,k)
!-AP      ETA = 1.77E-4 * SQRT(T(J)/288.)
!-AP  From Prupacher & Klett
               f1 = 2./9.*rhop(j)*r*r*g/eta(j)
               alph = a + b*EXP(-c*r/alam(j))
               WFALL(j,k) = f1*(1.+alam(j)*alph/r)
            ENDIF
         ENDDO


      ENDDO

!


!      print *, 'stopping in Sedmnt'
!      stop
      END
