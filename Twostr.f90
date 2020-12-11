
      SUBROUTINE TWOSTR(Sigr,U0,Sq,Wav,Ll,nzp1,nz2,absorbers,S)

!note i haven't yet done the work to fully abstract particles in this subroutine.
!hydrocarbons are using optical properties files, while sulfate and sulfur are hardcoded.

! in call, from Photo, send WAV(L) in the WAV position, and L in the LL position
!IKN is a printing option, NN is just N as sent to Photo
!S is the returned source function

!   This is my version of the Toon et al. 2-stream code.  (Ref.: JGR
!   94, 16287, 1989).  It vectorizes over height, rather than wavelength,
!   and is designed to work with PRIMS3 and its companion photochemical
!   models.
!
!   For now, at least, it is hardwired as the quadrature approximation.
!   NP is the number of different types of particles

      IMPLICIT NONE

      ! module variables
      ! ...

      ! local varaibles
      real*8, intent(in) :: sigr(nz)
      real*8, intent(in) :: U0
      real*8, intent(in) :: sq(kj,nz,kw)
      real*8, intent(in) :: wav
      integer, intent(in) :: Ll
      integer, intent(in) :: nzp1
      integer, intent(in) :: nz2
      real*8, intent(in) :: absorbers(kj,nz)
      real*8, intent(out) :: S(nz)

      real*8 tau(NZ) , tauctstr(NZP1) , gt(NZ) , gam1(NZ) , gam2(NZ) &
              & , gam3(NZ) , gam4(NZ) , alam(NZ) , cgam(NZ) , e1(NZ) ,  &
              & e2(NZ) , e3(NZ) , e4(NZ) , cp0(NZ) , cpb(NZ) , cm0(NZ) ,&
              & cmb(NZ) , y1(NZ) , y2(NZ) , w0(NZ) , tausg(NZ) ,        &
              & tausp(NZ) , direct(NZP1) , amean(NZP1) , taug(NZ) ,     &
              & fmt(NZ)
      real*8 a(NZ2) , b(NZ2) , d(NZ2) , e(NZ2) , fup(NZP1) ,         &
              & fdn(NZP1)
      real*8 w0p(np) , qext(np) , tauscat_part(np,NZ) , taup(NZ)
      real*8 tau_part(np,NZ)
      real*8 denom, emlt, et0, etb, facm, facp, gp, gpnew
      integer i, j, k, l, l1, mz2, n,nflag, nzm1
      real*8 pi, sq3, ssfc, taua, taup1, u0m, u0m2, u1, u1m
      integer O

      O = 0

!this FUP is different from the lower boundary flux terms in Output.f.  Should be fine here.


!     U1 = 0.5  (Eddington value)
      sq3 = SQRT(3.)
      pi = 3.14159
      u1 = 1./sq3
      u0m = 1./U0
      u0m2 = u0m*u0m
      u1m = 1./u1
      nzm1 = NZ - 1
      mz2 = NZ2
      gp = 0.8  !old

!
!   Calculate the optical depths of the different layers.  TAUA is absorption,
!   TAUSG is scattering by gases, TAUSP is scattering by particles, TAUG is
!   extinction due to gases, TAU is total extinction due to gases and
!   particles.
!   Note that the grid for this subroutine is numbered from top to bottom,
!   whereas the main program is numbered from bottom to top.
!   First do gases

      DO i = 1 , NZ
         n = NZP1 - i
         tausp(n) = 0.
         taup(n) = 0.
         tausg(n) = Sigr(i)*DEN(i)*DZ(i)

! Can below be removed?
!orig      TAUA = ( SO3(I)*O3(I) + SO2(I)*O2(I) + SCO2*CO2(I) + SH2O*H2O(I)
!orig     2   + SSO2*FSO2(I) + SS2*S2(I) + SH2S*H2S(I) + SNH3*FNH3(I) )
!orig     3   * DEN(I)*DZ
!orig   TAUG(N) = TAUA + TAUSG(N)

!Sspecies(I) is the species cross section at the given height
!Fspecies are absolute values of species mixing ratios
!so, in my scheme where everything absorbs, this would be:
!a loop over kj (to get all species that absorb (at each height)) -
!of sq*absorbers

         taua = 0.0
         DO j = 1 , kj
            ! IF ( planet.EQ.'WASP12B' ) THEN
            !    IF ( j.LT.9 ) taua = taua + Sq(j,i,Ll)*ABSORBERS(j,i)
            ! ELSE
               taua = taua + Sq(j,i,Ll)*ABSORBERS(j,i)
            ! ENDIF
         ENDDO

         taug(n) = taua*DEN(i)*DZ(i) + tausg(n)

      ENDDO

!   Now do particles.  Must combine their scattering optical depths into a
!   single array in order to calculate W0 and G.  (TAUSP(N))
!   Don't need any arrays for pure absorption.
!   Scale optical depth, W0, and G for the particles using the Delta-
!   Eddington approximation of Joseph et al. (Ref: J. Atmos. Sci. 33,
!   2452, 1976)

!gna - this is the correction for forward scattering: delta eddington approximation

      IF ( np.GT.o ) THEN
                        !REMOVE THE PARTICLE STUFF, YET AGAIN....
!   Particle 1 is sulfate, 2 is S8, 3 is HCAER
         DO k = 1 , np
            w0p(k) = 0.0
            IF ( k.LT.3 ) w0p(k) = 1.
            IF ( k.EQ.2 .AND. Wav.LE.3500. ) w0p(k) = 0.5
                                                  !ACK hardcoded grid
 !      print*, 'K,W0P=', K,W0P
         ENDDO
!here, Kevin uses a formula to compute WOPS and QEXT (fom G0=0.8 and WOP)

         DO j = 1 , np
            DO i = 1 , NZ
               n = NZP1 - i
!-MC Here, we are including correct value of Qext and W0P for hydrocarbons
               IF ( j.GE.3 ) THEN
       !hardcoded for HCAER and HCAER2
                  qext(j) = QEXTT(Ll,i,j)
                  w0p(j) = W0T(Ll,i,j)
               ELSE
       !for sulfate and S8 : valid for large particles only, but this is status quo for now...
                  qext(j) = 2.
       !W0P set above
       !for sulfate and S8, set assymetry factor to GP=0.8 (alt and wl independent)
                  GFT(Ll,i,j) = 0.8
               ENDIF

               taup1 = qext(j)*pi*RPAR(i,j)*RPAR(i,j)*AERSOL(i,j)*DZ(i)
      !particle extinction for each particle at each height
      !not used, but might as well keep
               tau_part(j,i) = taup1
      !particle scattering for each particle at each height
               tauscat_part(j,i) = w0p(j)*taup1
      !TAUSP contains total particle scattering at each height
               tausp(n) = tausp(n) + tauscat_part(j,i)
      !TAUP contains total particle extinction at each height
               taup(n) = taup(n) + taup1
            ENDDO
         ENDDO
      ENDIF

!   Calculate W0 and G by averaging over Rayleigh and Mie scatterers.
!   (scattering due to gases vs. particles)
!   Avoid letting W0 equal exactly 1.
      DO n = 1 , NZ
      !TAU is total extinction (Gas + Rayleigh + Particle)
         tau(n) = taug(n) + taup(n)
      !w0 = total scattering/total extinction
         w0(n) = (tausg(n)+tausp(n))/tau(n)
         w0(n) = AMIN1(w0(n),0.99999)
!-mab: Was "0.999" before - extra 99s needed for WASP12B convergence

      !GFT is still bottom to top, so needs a switch
         i = NZP1 - n

         gpnew = 0.0
         IF ( np.GT.o ) THEN
                        !REMOVE THE PARTICLE STUFF, YET AGAIN....
            DO k = 1 , np
      !asymmetry factors weighted by particle scattering/total scattering
               gpnew = gpnew + GFT(Ll,i,k)*tauscat_part(k,i)            &
                     & /(tausp(n)+tausg(n))
            ENDDO

      !dont let assymetry factor get larger than 1
            gt(n) = AMIN1(gpnew,0.99999)
         ELSE
            gt(n) = 0.0
         ENDIF

      ENDDO


!   Delta-Eddington scaling
!-AP I used approximation from Joseph et al. 1976
      DO n = 1 , NZ
         fmt(n) = gt(n)*gt(n)
         tau(n) = tau(n)*(1.-w0(n)*fmt(n))
         w0(n) = w0(n)*(1.-fmt(n))/(1.-w0(n)*fmt(n))
         gt(n) = gt(n)/(1.+gt(n))
      ENDDO
!-AP**************************************************
!
!   Calculate the gamma's, lambda's, and e's
      DO n = 1 , NZ
!     GAM1(N) = (7. - W0(N)*(4.+3.*GT(N)))/4.
!     GAM2(N) = - (1. - W0(N)*(4.-3.*GT(N)))/4.
!     GAM3(N) = (2. - 3.*GT(N)*U0)/4.
!   (Eddington values above; quadrature values below)
!
         gam1(n) = sq3*(2.-w0(n)*(1.+gt(n)))/2.
         gam2(n) = sq3*w0(n)*(1.-gt(n))/2.
         gam3(n) = (1.-sq3*gt(n)*U0)/2.
         IF ( np.EQ.0 ) gam3(n) = (1.-sq3*U0)/2.
         gam4(n) = 1. - gam3(n)

         alam(n) = SQRT(gam1(n)*gam1(n)-gam2(n)*gam2(n))
         cgam(n) = (gam1(n)-alam(n))/gam2(n)
         emlt = EXP(-alam(n)*tau(n))

         e1(n) = 1. + cgam(n)*emlt
         e2(n) = 1. - cgam(n)*emlt
         e3(n) = cgam(n) + emlt
         e4(n) = cgam(n) - emlt
      ENDDO

!   Calculate A, B, and D, i.e. the coefficients of the tridiagonal matrix
!   Top of atmosphere
      a(1) = 0.
      b(1) = e1(1)
      d(1) = -e2(1)

!   Odd coefficients
      DO n = 1 , nzm1
         l = 2*n + 1
         a(l) = e2(n)*e3(n) - e4(n)*e1(n)
         b(l) = e1(n)*e1(n+1) - e3(n)*e3(n+1)
         d(l) = e3(n)*e4(n+1) - e1(n)*e2(n+1)
      ENDDO

!   Even coefficients
      DO n = 1 , nzm1
         l = 2*n
         a(l) = e2(n+1)*e1(n) - e3(n)*e4(n+1)
         b(l) = e2(n)*e2(n+1) - e4(n)*e4(n+1)
         d(l) = e1(n+1)*e4(n+1) - e2(n+1)*e3(n+1)
      ENDDO

!   Bottom of atmosphere
      a(NZ2) = e1(NZ) - alb*e3(NZ)
      b(NZ2) = e2(NZ) - alb*e4(NZ)
      d(NZ2) = 0.

!   Now, set up the RHS of the equation:
!   TAUCTSTR(N) is the optical depth above layer N
      tauctstr(1) = 0.
      DO n = 2 , NZP1
         tauctstr(n) = tauctstr(n-1) + tau(n-1)
      ENDDO
! On last call,
! Print out TAUCTSTR(N), TAU(N), W0(N), GT(N), TAUSG(N), TAUSP(N).
! Also print TAUC at the ground for all wavelengths.

!       IF ( Nn.EQ.1 .AND. Ikn.EQ.1 ) THEN
!       !ACK - check LUN
!          WRITE (22,99001) Wav , tauctstr(NZP1)
! 99001    FORMAT (1X,F6.1,2X,1PE10.3)
!       !ACK - what are these hardcoded wavelengths???
!          IF ( Wav.EQ.1860.5 .OR. Wav.EQ.2010. .OR. Wav.EQ.2116.5 .OR.   &
!             & Wav.EQ.2211. .OR. Wav.EQ.2312.5 .OR. Wav.EQ.2516. .OR.    &
!             & Wav.EQ.3007.5 .OR. Wav.EQ.3900. .OR. Wav.EQ.4500. ) THEN
!       !ACK - check LUN
!             WRITE (20,99002) Wav , U0 , alb
! 99002       FORMAT ('# WAV = ',F6.1,2X,'U0 = ',F6.4,2X,'Rsfc = ',F5.3,/ &
!              &'# TWOSTR: TAUCTSTR(N) is the optical depth above layer N'&
!             & ,/'#  Z',4X,'TAUCTSTR(N)',4X,'TAU(N)',5X,'W0(N)',6X,      &
!             & 'G(N)',6X,'TAUSG(N)',3X,'TAUSP(N)',3X,'N')
!             DO n = 1 , NZ
!                i = NZP1 - n
!       !ACK - check LUN
!                WRITE (20,99003) i , tauctstr(n) , tau(n) , w0(n) ,      &
!                               & gt(n) , tausg(n) , tausp(n) , n
! 99003          FORMAT (1X,I3,1P6E11.3,1X,I3)
!             ENDDO
!       !ACK - check LUN
!             WRITE (20,99004) tauctstr(NZP1)
! 99004       FORMAT ('#  0',1PE11.3)
!          ENDIF
!       ENDIF

!   DIRECT(N) is the direct solar flux at the top of layer N.  Values
!   are normalized to unity.  DIRECT(NZP1) is the direct flux at the ground.
      direct(1) = U0
      DO n = 1 , NZ
         et0 = EXP(-tauctstr(n)/U0)
         etb = et0*EXP(-tau(n)/U0)
         direct(n+1) = etb*U0
         denom = alam(n)*alam(n) - u0m2
         facp = w0(n)*((gam1(n)-u0m)*gam3(n)+gam4(n)*gam2(n))
         facm = w0(n)*((gam1(n)+u0m)*gam4(n)+gam2(n)*gam3(n))

         cp0(n) = et0*facp/denom
         cpb(n) = etb*facp/denom
         cm0(n) = et0*facm/denom
         cmb(n) = etb*facm/denom
      ENDDO
      ssfc = alb*direct(NZP1)

!   Odd coefficients
      e(1) = -cm0(1)
      DO n = 1 , nzm1
         l = 2*n + 1
         e(l) = (cp0(n+1)-cpb(n))*e3(n) + (cmb(n)-cm0(n+1))*e1(n)
      ENDDO

!   Even coefficients
      DO n = 1 , nzm1
         l = 2*n
         e(l) = (cp0(n+1)-cpb(n))*e2(n+1) - (cm0(n+1)-cmb(n))*e4(n+1)
      ENDDO
      e(NZ2) = ssfc - cpb(NZ) + alb*cmb(NZ)

!   Call the tridiagonal solver (from LINPACK).  E is the RHS of the matrix
!   equation on input and is the solution vector Y on output
      CALL SGTSL(mz2,a,b,d,e,nflag)
      IF ( nflag.NE.0 ) PRINT 99005 , nflag
99005 FORMAT (/1X,'Tridiagonal solver failed in TWOSTR, NFLAG =',I4)

      DO n = 1 , NZ
         l = 2*n
         l1 = l - 1
         y1(n) = e(l1)
         y2(n) = e(l)
      ENDDO

!   Calculate the mean intensities, AMEAN(N), at the boundaries between
!   the layers.  AMEAN(N) is the intensity at the top of layer N.
      amean(1) = u1m*(y1(1)*e3(1)-y2(1)*e4(1)+cp0(1)) + 1.
      DO n = 1 , NZ
         amean(n+1) = u1m*(y1(n)*(e1(n)+e3(n))+y2(n)*(e2(n)+e4(n))      &
                    & +cpb(n)+cmb(n)) + direct(n+1)/U0
      ENDDO

!   Reset any AMEAN values that may go negative.  Check error file
!   to be sure this only happens near the ground where AMEAN ~ 0.
      DO n = 1 , NZP1
         IF ( amean(n).LT.0.0 ) THEN
         !ACK - check LUN
            WRITE (13,99006) Wav , n , amean(n)
99006       FORMAT ('WAVE =',F6.1,' AMEAN(',I3,')=',1PE11.3)
            amean(n) = ABS(amean(n))
         ENDIF
      ENDDO

!  Calculate upward and downward fluxes.

      fup(1) = ((y1(1)*e3(1)-y2(1)*e4(1))+cp0(1))
      fdn(1) = direct(1)
      DO n = 1 , NZ
         fup(n+1) = (y1(n)*e1(n)+y2(n)*e2(n)+cpb(n))
         fdn(n+1) = (y1(n)*e3(n)+y2(n)*e4(n)+cmb(n)) + direct(n+1)
      ENDDO

!   Convert back to main program grid.  S(I) is the mean intensity at the
!   midpoint of layer I.
      DO i = 1 , NZ
         n = NZP1 - i
         S(i) = SQRT(amean(n)*amean(n+1))
      ENDDO

!  Print out the results at a few wavelengths

!       IF ( Nn.EQ.1 .AND. Ikn.EQ.1 ) THEN
!          IF ( Wav.EQ.1860.5 .OR. Wav.EQ.2010. .OR. Wav.EQ.2116.5 .OR.   &
!             & Wav.EQ.2211. .OR. Wav.EQ.2312.5 .OR. Wav.EQ.2516. .OR.    &
!             & Wav.EQ.3007.5 .OR. Wav.EQ.3900. .OR. Wav.EQ.4500. ) THEN
!       !ACK - check LUN
!             WRITE (21,99007) Wav , U0 , alb
! 99007       FORMAT (/'# WAV = ',F6.1,2X,'U0 = ',F6.4,2X,'Rsfc = ',F5.3, &
!                    &/'#  Z',4X,'S(Z)',7X,'Y1',9X,'Y2',9X,'DIRECT',5X,   &
!                    &'AMEAN',6X,'FUP',8X,'FDN',7X,'N')
!             DO n = 1 , NZ
!                i = NZP1 - n
!       !ACK - check LUN
!                WRITE (21,99008) i , S(i) , y1(n) , y2(n) , direct(n) ,  &
!                               & amean(n) , fup(n) , fdn(n) , n
! 99008          FORMAT (1X,I3,1P7E11.3,1X,I3)
!             ENDDO
!       !ACK - check LUN
!             WRITE (21,99009) direct(NZP1) , amean(NZP1) , fup(NZP1) ,   &
!                            & fdn(NZP1)
! 99009       FORMAT ('#',36X,1P4E11.3)
!          ENDIF
!       ENDIF

      END
!-PK ***************************************
