
      SUBROUTINE PHOTSATRAT(Jtrop, nz, H2O)

      IMPLICIT NONE

      ! module variables
      ! h2osat


      ! local variables
      integer, intent(in) :: jtrop
      integer, intent(in) :: nz
      real*8, dimension(nz) ,intent(out) :: H2O ! H2O mixing ratio!

      real*8 t0, pzero, amv, vapl, subl, r, a0, bk, ps
      integer j
      real*8 hl(nz) , a(nz)
      real*8 p1, p2, pv, rel

!
!   THIS SUBROUTINE CALCULATES THE SATURATION VAPOR PRESSURE OF WATER
!   FROM MAGNUS' EQUATION (SEE RUNAWAY GREENHOUSE PAPER IN ICARUS -
!   APPENDIX A).  IT THEN FIXES TROPOSPHERIC H2O USING A MANABE/
!   WETHERALD RELATIVE HUMIDITY DISTRIBUTION.
!
      t0 = 273.15
      pzero = 6.103E-3          !vapor pressure at T0 I assume
      amv = 18.
      vapl = 597.3
      subl = 677.
      r = 1.9872
      a0 = 0.553
      bk = 1.38E-16
      ps = 1.E-6*DEN(1)*bk*T(1)      !this is actually a good measure of pressure
!
      DO j = 1 , nz
         hl(j) = subl
         a(j) = 0.
         IF ( T(j).GE.t0 ) THEN
            hl(j) = vapl + a0*t0
            a(j) = a0
         ENDIF
      ENDDO
!
!   FIND SATURATION VAPOR PRESSURE
      DO j = 1 , nz
         p1 = pzero*(t0/T(j))**(amv*a(j)/r)
         p2 = EXP(amv*hl(j)/r*(1./t0-1./T(j)))
         pv = p1*p2
         ! P(j) = 1.E-6*DEN(j)*bk*T(j)   ! bars !already defined
         H2OSAT(j) = pv/P(j)
      ENDDO
!
!   CALCULATE TROPOSPHERIC H2O CONCENTRATIONS
      DO j = 1 , Jtrop
         IF ( planet.EQ.'EARTH' ) THEN
            rel = 0.77*(P(j)/ps-0.02)/0.98  !manabe formula


!       REL = 1.0 ! Saturated troposphere - eventual should abstract this and ATACAMA

!       REL = 0.17   !test ATACAMA
         ELSEIF ( planet.EQ.'MARS' ) THEN
!      REL = 0.12                         ! 7 microns
            rel = 0.17                    ! this is good for 9.5 micrometers
! this makes the model very dry  - warning nonstandard!!!

         ELSEIF ( planet.EQ.'DRY' ) THEN
                                       !Dry planets. Implemented originally for post-runaway, O2-rich atmospheres
            rel = 1.E-8
                  ! For H2O-free atmospheres. Set above zero to prevent floating point errors. - EWS - 9/14/2015
         ENDIF

         ! RELH(j) = rel
         H2o(j) = rel*H2OSAT(j)


      ENDDO

!


      END
