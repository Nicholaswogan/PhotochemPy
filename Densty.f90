!*==DENSTY.spg  processed by SPAG 6.72Dc at 16:52 on 10 Dec 2020
      SUBROUTINE DENSTY

      IMPLICIT NONE

      ! module variables
      ! some?

      ! local varaibles
      real*8 g0, rgas, bk
      real*8 wt, ft, roverm
      real*8 t0, p1, ha, r, tav, gz
      integer i


!   THIS SUBROUTINE CALCULATES ATMOSPHERIC NUMBER DENSITIES, ASSUM-
!   ING HYDROSTATIC EQUILIBRIUM
      g0 = g
      rgas = 8.3143E7
      bk = 1.38054E-16


      ft = usol_init(lo2,1) + fco2 + far
      wt = usol_init(lo2,1)*32. + fco2*44. + (1.-ft)*28. + far*40.  !assuming O2,CO2,N2 and Ar are the main players and pN2=1-pO2+pCO2+pAr


      roverm = rgas/wt

!      PG = 1.013E6    !primary specification of pressure in the model
!      P0 = PG  (P0 now specified in PLANET.dat)
!     P0 = PG GIVES YOU A ONE BAR ATMOSPHERE

!-mc      DZ = Z(2) - Z(1)   !ACK

      t0 = T(1) + (T(1)-T(2))/2.
      ha = roverm*0.5*(t0+T(1))/g0
      p1 = p0*1E6*EXP(-0.5*DZ(1)/ha)
      DEN(1) = p1/(bk*T(1))
!
! ***** FIND DENSITY FROM HYDROSTATIC EQUILIBRIUM *****
      DO i = 2 , nz
!-mc      DZ = Z(I) - Z(I-1)
         r = r0 + Z(i)
         gz = g0*(r0/r)*(r0/r)
         tav = 0.5*(T(i)+T(i-1))
         ha = roverm*tav/gz
         DEN(i) = DEN(i-1)*EXP(-DZ(i)/ha)*T(i-1)/T(i)
      ENDDO
!


! ***** FIND PRESSURE FROM THIS DENSITY *********
      DO i = 1 , nz
         PRESS(i) = DEN(i)*bk*T(i) ! dynes
      ENDDO

      ! pressure in bars
      do i = 1,nz
        P(i) =  DEN(i)*1.38E-16*T(i) * 1.0e-6 !bars
      enddo


      ! PRINT * , "Molecular weight of atmosphere = " , wt

!
      END
