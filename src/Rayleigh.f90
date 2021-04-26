!*==RAYLEIGH.spg  processed by SPAG 6.72Dc at 18:36 on 11 Dec 2020
      SUBROUTINE RAYLEIGH(Wleff,Ncomp,Icomp,Volmix,nz,Sigr2)

!ccccccccccccccccccccccccccccc  r a y l e i g h cccccccccccccccccccccccc
!c                                                                    cc
!c    p u r p o s e :                                                 cc
!c                                                                    cc
!c    this subroutine computes the rayleigh scattering cross section  cc
!c    per molecule (meters**-2)  for atmospheres composed of (1) air, cc
!c    (2) co2, (3) n2, (4) o2, (5), h2o, (6) h2, (7) he, or any       cc
!c    combination of these gases.                                     cc
!c    to find the rayleigh scattering optical depth at any level of   cc
!c    the atmosphere, these cross sections must be multiplied by      cc
!c    the pathlength-integrated number density, n(z)*dz.  in a        cc
!c    hydrostatic atmosphere, n(z)*dz = a0*dp/(u0*grav), where u0     cc
!c    is the mean molecular mass (kg/mole), and a0 is avagadro's      cc
!c    number(kg/kmole)                                                cc
!c                                                                    cc
!c    r e f e r e n c e s :                                           cc
!c                                                                    cc
!c    e.j. mccartney, optics of the atmosphere, wiley, p. 187-215,    cc
!c          1976.                                                     cc
!c    a.t. young, revised depolarization corrections for atmospheric  cc
!c          extinction, appl. opt. 19, 3427-3428, 1980.               cc
!c    c.w. allen, astrophysical quantities, athlone press, p. 87,     cc
!c         1964.                                                      cc
!c                                                                    cc
!c    i n p u t :                                                     cc
!c                                                                    cc
!c    ncomp = number of major atmosphere constituents                 cc
!c    icomp = atmosphere constituent index                            cc
!c            (1) air  (2) co2  (3) n2  (4) o2  (5) h2o (6) h2 (7) he cc
!c   volmix = volume mixing ratio of each constituent                 cc
!c    wleff = effective wavelength (microns)                          cc
!c                                                                    cc
!c    o u t p u t :                                                   cc
!c                                                                    cc
!c    sigr2  - raleigh scattering cross section per molecule (cm**2)  cc
!c                                                                    cc
!ccccccccccccccccccccccccccccc  r a y l e i g h  ccccccccccccccccccccccc

      IMPLICIT NONE
      ! global variables
      ! ...

      ! local variables
! RAYLEIGH(Wleff,Ncomp,Icomp,Volmix,Sigr2)
      integer, intent(in) :: nz
      real*8, intent(in) :: wleff
      integer, intent(in) :: ncomp(nz)
      integer, intent(in) :: icomp(10,nz)
      real*8, intent(in) :: volmix(10,nz)

      real*8, intent(out) :: sigr2(nz)
      INTEGER i , j
      ! REAL*8 Sigr2 , Volmix , Wleff
      REAL*8 delta(7) , a(7) , b(7) , wl2i , r , aniso , r2 ,    &
           & sum , cnst
      ! real*8 Volmix(10,nz) , Ncomp(nz) , Icomp(10,nz) , Sigr2(nz)

!       depolarization factors for air, co2, n2, o2 (young, 1980)
!-mab: h2o depolarization from Marshall & Smith (1990).
      DATA delta/0.0279E0 , 0.078E0 , 0.021E0 , 0.058E0 , 0.17E0 ,      &
         & 0.0000 , 0.0000/

!***    wavelength dependence coefficients for the refractive index
!       (allen, 1964) (note: wavelengths must be in microns)
!-mab: Tabulating h2o A, B as 0.85*(air value) as per suggested in
!-mab: Kopparapu et al. 2013 (Habitable zones) and P. von Paris 2013.

      DATA a/2.871E-04 , 4.39E-04 , 2.906E-04 , 2.663E-04 , 2.44E-04 ,  &
         & 1.358E-4 , 3.48E-5/
      DATA b/5.67E-3 , 6.4E-3 , 7.7E-3 , 5.07E-3 , 4.82E-3 , 7.52E-3 ,  &
         & 2.3E-3/

!****   Define the constant multiplier, cnst.  this
!	constant is equal to 24*pi**3/(1.e-24*L**2),
!	where L = loschmidt's number (mks units) ,
!	(L = 2.687e25 molecules/m3) and the factor 1.e-24 is
!	needed to convert wl**4 from microns to meters.

      DATA cnst/1.031E-24/

      wl2i = 1.0E0/(Wleff*Wleff)

      DO j = 1 , nz
         sum = 0.
!      print *, ncomp(j)
         DO i = 1 , Ncomp(j)
!          print *, icomp(i,j)
            r = 1.0E0 + a(Icomp(i,j))*(1.0E0+b(Icomp(i,j))*wl2i)
            aniso = (6.0E0+3.0E0*delta(Icomp(i,j)))                     &
                  & /(6.0E0-7.0E0*delta(Icomp(i,j)))
            r2 = r*r + 2.0E0
            sum = sum + Volmix(i,j)*aniso*((r*r-1.)/r2)**2

         ENDDO
       !1d4 converts from m^2 to cm^2 for photohcemical model
         Sigr2(j) = cnst*wl2i*wl2i*sum*1E4
      ENDDO

      END
