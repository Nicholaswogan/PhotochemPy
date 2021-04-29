
      SUBROUTINE DIFCO
        use photochem_data, only: nq, nz, nz1, g, planet, mass, &
                                  background_mu
        use photochem_vars, only: usol_init, edd, den, T
        use photochem_wrk, only: hscale, tauedd, dk, H_atm, bhn2, bh2n2, scale_h


      IMPLICIT NONE

      ! module variables
      ! HSCALE
      ! tauedd
      ! H_ATM, SCALE_H, BHN2, BH2N2
      ! real*8, dimension(nz) :: DK


      ! local variables
      real*8 wt
      real*8 bkmg, eddav, denav
      real*8 h, tav
      integer i,j



      call mean_molecular_weight(nq, usol_init(:,1), mass, background_mu, wt)
      bkmg = 1.38E-16/(1.67E-24*wt*g)    !a good pressure
      tav = 0.d0
!
! ***** DK(I) = K*N AT GRID STEP I+1/2 *****
!
      DO i = 1 , nz1
         eddav = SQRT(EDD(i)*EDD(i+1))
                                     !average eddy diffusion in grid center
         denav = SQRT(DEN(i)*DEN(i+1))
                                     !average density at grid center
         DK(i) = eddav*denav
      ENDDO
!
!   COMPUTE DIFFUSION LIFETIME AT EACH HEIGHT (H*H/K)
      DO i = 1 , nz
         h = bkmg*T(i)
         HSCALE(i) = h
         TAUEDD(i) = h*h/EDD(i)
      ENDDO

      DO i = 1 , nz1
         tav = SQRT(T(i)*T(i+1)) !average temperature at grid center
         H_ATM(i) = bkmg*tav
!  compute scale heights of all the species in the atmosphere at all
!   heights
         DO j = 1 , nq
            SCALE_H(j,i) = bkmg*tav*wt/MASS(j)
         ENDDO

!          IF ( fh2.GE.0.50 ) THEN
! !-mab: Proceed with this loop for giant planets instead...
! !-mab: transferring the "b" portion of the "DI" expression from Ravi's version since DI = b/n
! !-mab: where b is the binary molecular diffusion parameter and DI is the diffusion coefficient
! !-mab: got relation from Jim's new book Chapter 5.
! ! (b has a constant*T^(some power) form for terrestrial gases in Chapter 5. Not established for giant planets?_?)
!             DO j = 1 , nq
!                BX1X2(j,i) = 1.52E18*(1./MASS(j)+1./wt)**0.5*(tav**0.5)
!             ENDDO
!-mab: One of the below options to be used for terrestrial planets.
         IF ( planet.EQ.'EARTH' ) THEN
            BHN2(i) = 2.7E19*(tav/200.)**0.75
                                             ! correct for N2
            BH2N2(i) = 1.4E19*(tav/200.)**0.75
                                             ! correct for N2
!        bXN2(i) = 4.0D18*(TAV/200.)**0.75
         ELSEIF ( planet.EQ.'MARS' ) THEN
            BHN2(i) = 0.8*1.8*1.4E19*(tav/200.)**0.75
                                                    ! correct for CO2
            BH2N2(i) = 0.8*1.4E19*(tav/200.)**0.75
                                                 ! correct for CO2
!        bXN2(i) = 4.0D18*(TAV/200.)**0.75
         ELSEIF ( planet.EQ.'DRY' ) THEN ! assume O2-dominated - this is the same as N2 for now
            BHN2(i) = 2.7E19*(tav/200.)**0.75
                                             ! correct for N2
            BH2N2(i) = 1.4E19*(tav/200.)**0.75
                                             ! correct for N2
! Note: there will be a separate expression in output.f for this scenario as well
         ENDIF
!        bXN2(i) = 4.0D18*(TAV/200.)**0.75
!-mab: Should we put an error check for cases that don't fall above or below?
!-mab: What are we doing for Titan?
      ENDDO

      H_ATM(nz) = bkmg*tav
      DO j = 1 , nq
         SCALE_H(j,nz) = bkmg*T(nz)*wt/MASS(j)
      ENDDO

      ! IF ( fh2.GE.0.50 ) THEN
      !    DO j = 1 , nq
      !       BX1X2(j,nz) = 1.52E18*(1./MASS(j)+1./wt)**0.5*(T(nz)**0.5)
      !    ENDDO
!-mab: One of the below options to be used for terrestrial planets.
      IF ( planet.EQ.'EARTH' ) THEN
         BHN2(nz) = 2.7E19*(T(nz)/200.)**0.75   ! correct for N2
         BH2N2(nz) = 1.4E19*(T(nz)/200.)**0.75  ! correct for N2
!        bXN2(nz) = 4.0D18*(T(nz)/200.)**0.75
      ELSEIF ( planet.EQ.'MARS' ) THEN
         BHN2(nz) = 0.8*1.8*1.4E19*(T(nz)/200.)**0.75  ! correct for CO2
         BH2N2(nz) = 0.8*1.4E19*(T(nz)/200.)**0.75     ! correct for CO2
!        bXN2(nz) = 4.0D18*(T(nz)/200.)**0.75
      ELSEIF ( planet.EQ.'DRY' ) THEN    ! assume O2-dominated - this is the same as N2 for now
         BHN2(i) = 2.7E19*(tav/200.)**0.75   ! correct for N2
         BH2N2(i) = 1.4E19*(tav/200.)**0.75  ! correct for N2
      ENDIF     !output.f impacted



      END
