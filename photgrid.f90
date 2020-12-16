
      subroutine photgrid(top_atmosphere)
        implicit none

        ! module variables
        ! integer :: nz ! number of vertical grid points
        ! real*8, allocatable, dimension(:) :: z
        ! real*8, allocatable, dimension(:) :: dz

        ! local variables
        real*8, intent(in) :: top_atmosphere !
        real*8 :: dzgrid
        integer :: i

        dzgrid = top_atmosphere/nz

        do I=1,NZ
         Z(I) = (I - 0.5)*dzgrid
        enddo


        DZ(1)=Z(1)+0.5*dzgrid
        do I=2,NZ
          DZ(I)=Z(I)-Z(I-1)
        enddo

        JTROP=minloc(Z,1, Z .ge. ztrop)-1

      end subroutine


      subroutine gridw(nw,wl,wc,wu,LGRID)
!*-----------------------------------------------------------------------------*
!*=  PURPOSE:                                                                 =*
!*=  Create the wavelength grid for all interpolations and radiative transfer =*
!*=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
!*=  No gaps are allowed within the wavelength grid.                          =*
!*-----------------------------------------------------------------------------*
!*=  PARAMETERS:                                                              =*
!*=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
!*=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
!*=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
!*=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
!*=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
!*-----------------------------------------------------------------------------*
!*=  EDIT HISTORY:                                                            =*
!*=  Original                                                                 =*
!C R. F. Esswein 020214 Change all REAL declarations to REAL*8
!C R. Esswein 020221 Change file references to full path names.
!C R. Esswein 030407 Choose extension to lower wavelengths.
!C M. Claire  060802  Integrating into Kevin's code
!*-----------------------------------------------------------------------------*
!*= This program is free software;  you can redistribute it and/or modify     =*
!*= it under the terms of the GNU General Public License as published by the  =*
!*= Free Software Foundation;  either version 2 of the license, or (at your   =*
!*= option) any later version.                                                =*
!*= The TUV package is distributed in the hope that it will be useful, but    =*
!*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!*= License for more details.                                                 =*
!*= To obtain a copy of the GNU General Public License, write to:             =*
!*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!*-----------------------------------------------------------------------------*
      implicit none

      ! module variables
      ! none!

      ! local variables
      integer, parameter :: kw = 2900
      integer, intent(in) :: LGRID ! input:
      real*8, dimension(kw), intent(out) :: wl, wc, wu ! output
      integer, intent(out) :: nw
      character(len=60) :: string
      real*8 wincr
      integer iw, I, l, kin
      logical ok
      integer idum
      real*8 dum
      integer mopt

      DO iw = 1, kw ! zero out
        wl(iw) = 0.
        wu(iw) = 0.
        wc(iw) = 0.
      ENDDO

      !**** chose wavelengths

      !* some pre-set options
      !*     mopt = 1    equal spacing
      !*     mopt = 2    Isaksen's grid
      !*     mopt = 3    combined Kockarts/Isaksen grid + Lyman-Alpha
      !*     mopt = 4    user-defined
      !*     mopt = 5    (Notes below)
      ! grid from Kevin's code,used in Zahnle.flx/.grid
      ! (This is also Allen Grid for S-R + old JPL grid)
      !C-mab The present stellar flux files for Hot Jupiters also use this grid.
      !*     mopt = 6    grid from Jim's climate code  (entirly a hack right now for interpolative purposes)
      !*     mopt = 7    Jim's old grid, but high resolution from 175-220

      !note - before using, make sure that nw is the number of wavelengths to loop over and
      !that wl(nw+1)=wu(nw)
      ! this is needed for the interpolations to work correctly
      if (LGRID.eq.0) mopt = 5
      ! if (LGRID.eq.1) mopt = 7

      ! IF (mopt .EQ. 1) GO TO 1
      ! IF (mopt .EQ. 2) GO TO 2
      ! IF (mopt .EQ. 3) GO TO 3
      ! IF (mopt .EQ. 4) GO TO 4
      IF (mopt .EQ. 5) GO TO 5 ! main option!
      ! IF (mopt .EQ. 6) GO TO 6
      ! IF (mopt .EQ. 7) GO TO 7
 5    CONTINUE
      ! read in wogan's grid here. THis is the same as kevins but extended down
      ! to 121.0 nm instead of 121.6 nm
      nw = 118
      OPEN(kin, file='DATA/GRIDS/wogan.grid',status='old')

      DO i = 1,2
        READ(kin,*)  !skip header
      ENDDO

      do L=1,nw
        READ(kin,*) WL(L),WU(L)
        wc(L) = ( wl(L) + wu(L) )/2.
      enddo
      wl(nw+1) = wu(nw)  !final point for interpolative array
      CLOSE (kin)
      GO TO 9
 9    CONTINUE

      !c-mc should probably print these out to output file rather than screen
      ! print *, 'NW = ',nw,'   WAVELENGTH GRID:',wl(1),' -',wu(nw), &
      !        ' Angstroms'
      !* check grid for assorted improprieties:
      CALL gridck(kw,nw,wl,ok)

      IF (.NOT. ok) THEN
        print *,'STOP in GRIDW:  The w-grid does not make sense'
        STOP
      ENDIF

      end subroutine gridw


      SUBROUTINE readflux(flux_txt,nw,wl,f)

!*-----------------------------------------------------------------------------*
!*=  PURPOSE:                                                                 =*
!*=  Read and re-grid extra-terrestrial flux data.                            =*
!*-----------------------------------------------------------------------------*
!*=  PARAMETERS:                                                              =*
!*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!*=           wavelength grid                                                 =*
!*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!*=           working wavelength grid                                         =*
!*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
!*=           each specified wavelength                                       =*
!*=  timega - solar system age in Ga    nope!                                 =*
!*-----------------------------------------------------------------------------*
!*=  EDIT HISTORY:                                                            =*
!*=  02/97  Changed offset for grid-end interpolation to relative number      =*
!*=         (x * (1 +- deltax))                                               =*
!*=  05/96  Put in different preset options                                   =*
!C R. F. Esswein 020214 Change all REAL declarations to REAL*8
!C    Use parameters "zero" and "largest" as
!C      parameters to "addpnt" calls.
!C R. Esswein 020221 Change file references to full path names.
!c       M. Claire       091306 integrating into Kevin's code
!*-----------------------------------------------------------------------------*
!*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!*-----------------------------------------------------------------------------*

      implicit none

      ! module variables
      ! none!

      ! local variables
      real*8, PARAMETER :: deltax = 1.E-4, biggest=1.E+36, zero=0.0
      integer, parameter :: kdata = 26500, kw=2900
      ! * input: (wavelength grid)
      INTEGER, intent(in) :: nw
      REAL*8, intent(in) :: wl(kw)
      character(len=*), intent(in) :: flux_txt
      ! * output: (extra terrestrial solar flux)
      REAL*8, intent(out) :: f(kw)
      ! * work arrays for input data files:
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata), y3(kdata)
      INTEGER nhead, n, i, ierr, kin, io, n3
      INTEGER iw
      ! * data gridded onto wl(kw) grid:
      REAL*8 yg1(kw)
      REAL*8 yg2(kw)
      REAL*8 yg3(kw)
      REAL*8, parameter :: hc = 6.62E-34 * 2.998E8

      nhead = 0
      ierr = 0
      OPEN(UNIT=kin, &
      file=flux_txt, &
      STATUS='old')

      n = 0
      DO i=1,kdata
        READ(kin,*,iostat=io) x1(i), y1(i)  !this flux in mw/m2/nm, but is sampled at subangstom resolution
        if (io/=0) exit
        x1(i)=x1(i)*10e0   !convert wavelength from nm to Angstoms
        x2(i)=x1(i)      ! x2 also angstroms
        x3(i)=x1(i)      ! x3 also angstroms
        y3(i)=y1(i)/10e0   ! for y3, convert thuillier flux to mw/m2/A
        n = n + 1
      ENDDO
      CLOSE (kin)




      n3=n
      ierr=0

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,          zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,        biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,ierr)  !inter2 is points to bins

      !so yg3 is flux on the model wavelength grid

        !error check for call to inter2
       IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr,'  Something wrong in Grid.f/readflux'
        STOP
       ENDIF

! NOTE: This explicitly assumes the correct Lyman-Alpha flux *at the planet* has been included
!      in the input spectrum. It would benefit you to ensure that this is really the case!
      DO iw = 1, nw
        f(iw) = yg3(iw)*(wl(iw+1)-wl(iw))*5.039e8*wl(iw)/10. !convert to photons/cm2/s
        ! print *, wl(iw), yg3(iw)
      ENDDO

      END subroutine


      SUBROUTINE gridck(k,n,x,ok)

!*-----------------------------------------------------------------------------*
!*=  PURPOSE:                                                                 =*
!*=  Check a grid X for various improperties.  The values in X have to comply =*
!*=  with the following rules:                                                =*
!*=  1) Number of actual points cannot exceed declared length of X            =*
!*=  2) Number of actual points has to be greater than or equal to 2          =*
!*=  3) X-values must be non-negative                                         =*
!*=  4) X-values must be unique                                               =*
!*=  5) X-values must be in ascending order                                   =*
!*-----------------------------------------------------------------------------*
!*=  PARAMETERS:                                                              =*
!*=  K  - INTEGER, length of X as declared in the calling program          (I)=*
!*=  N  - INTEGER, number of actual points in X                            (I)=*
!*=  X  - REAL, vector (grid) to be checked                                (I)=*
!*=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
!*=                .FALSE.-> X violates at least one of 1)-5)                 =*
!*-----------------------------------------------------------------------------*
!*=  EDIT HISTORY:                                                            =*
!*=  Original                                                                 =*
!C R. F. Esswein 020214 Change all REAL declarations to REAL*8
!C       M.C.            060802  Integrated into Kevin's code
!*-----------------------------------------------------------------------------*
!*= This program is free software;  you can redistribute it and/or modify     =*
!*= it under the terms of the GNU General Public License as published by the  =*
!*= Free Software Foundation;  either version 2 of the license, or (at your   =*
!*= option) any later version.                                                =*
!*= The TUV package is distributed in the hope that it will be useful, but    =*
!*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!*= License for more details.                                                 =*
!*= To obtain a copy of the GNU General Public License, write to:             =*
!*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!*-----------------------------------------------------------------------------*
!*= To contact the authors, please mail to:                                   =*
!*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!*= send email to:  sasha@ucar.edu                                            =*
!*-----------------------------------------------------------------------------*
!*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!*-----------------------------------------------------------------------------*
      IMPLICIT NONE

!* input:
      INTEGER k, n
      REAL*8 x(k)

!* output:
      LOGICAL ok

!* local:
      INTEGER i
!*_______________________________________________________________________

      ok = .TRUE.

!* check if dimension meaningful and within bounds

      IF (n .GT. k) THEN
        ok = .false.
        print *,'Number of data exceeds dimension'
        print *, k,n
        RETURN
       ENDIF

      IF (n .LT. 2) THEN
        ok = .FALSE.
        print *, 'Too few data, number of data points must be >= 2'
        RETURN
      ENDIF

!* disallow negative grid values

      IF(x(1) .LT. 0.) THEN
        ok = .FALSE.
        print *,'Grid cannot start below zero'
        RETURN
      ENDIF

!* check sorting
      DO i = 2, n
        IF( x(i) .LE. x(i-1)) THEN
          ok = .FALSE.
          print *,'Grid is not sorted or contains multiple values'
          print *, i, x(i),x(i-1)
          RETURN
        ENDIF
      enddo

      END subroutine



      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )

!*-----------------------------------------------------------------------------*
!*=  PURPOSE:                                                                 =*
!*=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!*=  ascending order                                                          =*
!*-----------------------------------------------------------------------------*
!*=  PARAMETERS:                                                              =*
!*=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!*=  Y    - REAL vector of length LD, y-values                            (IO)=*
!*=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
!*=         program                                                           =*
!*=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
!*=         N < LD.  On exit, N is incremented by 1.                          =*
!*=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!*=  YNEW - REAL, y-value of point to be added                             (I)=*
!*-----------------------------------------------------------------------------*

      IMPLICIT NONE

! calling parameters

      INTEGER ld, n
      REAL*8 x(ld), y(ld)
      REAL*8 xnew, ynew
      INTEGER ierr

!C local variables
      INTEGER insert
      INTEGER i



!* initialize error flag

      ierr = 0

!* check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
        print *, '>>> ERROR (ADDPNT) <<<  Cannot expand array '
        print *, '                        All elements used.'
        STOP
      ENDIF

      insert = 1
      i = 2

!* check, whether x is already sorted.
!* also, use this loop to find the point at which xnew needs to be inserted
!* into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
          print *, '>>> ERROR (ADDPNT) <<<  x-data must be '// &
                     'in ascending order!', i,x(i),x(i-1)
          STOP
      ELSE
        IF (xnew .GT. x(i)) insert = i + 1
      ENDIF
      i = i+1
      GOTO 10
      ENDIF

!* if <xnew,ynew> needs to be appended at the end, just do so,
!* otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN

        x(n+1) = xnew
        y(n+1) = ynew

      ELSE

!* shift all existing points one index up

        DO i = n, insert, -1
          x(i+1) = x(i)
          y(i+1) = y(i)
        ENDDO

!* insert new point

        x(insert) = xnew
        y(insert) = ynew

      ENDIF

!* increase total number of elements in x, y

      n = n+1

      END subroutine

      SUBROUTINE inter2(ng,xg,yg,n,x,y,ierr)

!*-----------------------------------------------------------------------------*
!*=  PURPOSE:                                                                 =*
!*=  Map input data given on single, discrete points onto a set of target     =*
!*=  bins.                                                                    =*
!*=  The original input data are given on single, discrete points of an       =*
!*=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!*=  of target bins.  In general, this is the case for most of the weighting  =*
!*=  functions (action spectra, molecular cross section, and quantum yield    =*
!*=  data), which have to be matched onto the specified wavelength intervals. =*
!*=  The average value in each target bin is found by averaging the trapezoi- =*
!*=  dal area underneath the input data curve (constructed by linearly connec-=*
!*=  ting the discrete input values).                                         =*
!*=  Some caution should be used near the endpoints of the grids.  If the     =*
!*=  input data set does not span the range of the target grid, an error      =*
!*=  message is printed and the execution is stopped, as extrapolation of the =*
!*=  data is not permitted.                                                   =*
!*=  If the input data does not encompass the target grid, use ADDPNT to      =*
!*=  expand the input array.                                                  =*
!*-----------------------------------------------------------------------------*
!*=  PARAMETERS:                                                              =*
!*=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!*=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!*=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!*=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!*=        bin i (i = 1..NG-1)                                                =*
!*=  N   - INTEGER, number of points in input grid                         (I)=*
!*=  X   - REAL, grid on which input data are defined                      (I)=*
!*=  Y   - REAL, input y-data                                              (I)=*
!*-----------------------------------------------------------------------------*

      IMPLICIT NONE
!* input:
      INTEGER ng, n
      REAL*8 x(n), y(n), xg(ng)
      INTEGER ierr
!* output:
      REAL*8 yg(ng)

!* local:
      REAL*8 area, xgl, xgu
      REAL*8 darea, slope
      REAL*8 a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart


!*_______________________________________________________________________

!*  test for correct ordering of data, by increasing value of x

      DO 10, i = 2, n
         IF (x(i) .LE. x(i-1)) THEN

            ierr = 1
            WRITE(*,*)'data not sorted', i,x(i),x(i-1)
            RETURN
         ENDIF
   10 CONTINUE


      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
           ierr = 2
          WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
          RETURN
        ENDIF
      ENDDO




!* check for xg-values outside the x-range

      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
          WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span '// &
                    'grid.  '
          WRITE(0,*) '                        Use ADDPNT to '// &
                    'expand data and re-run.'
          STOP
      ENDIF

!*  find the integral of each grid interval and use this to
!*  calculate the average y value for the interval
!*  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1
      ngintv = ng - 1
      DO 50, i = 1,ngintv

!* initalize:

            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)

!*  discard data before the first grid interval and after the
!*  last grid interval
!*  for internal grid intervals, start calculating area by interpolating
!*  between the last point which lies in the previous interval and the
!*  first point inside the current interval

            k = jstart

            IF (k .LE. n-1) THEN

!*  if both points are before the first grid, go to the next point
   30         CONTINUE
                IF (x(k+1) .LE. xgl) THEN
                   jstart = k - 1
                   k = k+1
                   IF (k .LE. n-1) GO TO 30
                ENDIF

!*  if the last point is beyond the end of the grid, complete and go to the next
!*  grid
   40         CONTINUE
                 IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN

                    jstart = k-1

!* compute x-coordinates of increment

                    a1 = MAX(x(k),xgl)
                    a2 = MIN(x(k+1),xgu)

!*  if points coincide, contribution is zero

                    IF (x(k+1).EQ.x(k)) THEN
                       darea = 0.e0
                    ELSE
                       slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                       b1 = y(k) + slope*(a1 - x(k))
                       b2 = y(k) + slope*(a2 - x(k))
                       darea = (a2 - a1)*(b2 + b1)/2.
!c                       print *,a2,a1,k,y(k),slope,b2,b1,darea
                    ENDIF


!*  find the area under the trapezoid from a1 to a2

                    area = area + darea

!* go to next point

                    k = k+1
                    GO TO 40

                ENDIF

            ENDIF

!*  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)


   50 CONTINUE
!*_______________________________________________________________________


      END subroutine


      SUBROUTINE inter3(ng,xg,yg, n,x,y, FoldIn)

!*-----------------------------------------------------------------------------*
!*=  PURPOSE:                                                                 =*
!*=  Map input data given on a set of bins onto a different set of target     =*
!*=  bins.                                                                    =*
!*=  The input data are given on a set of bins (representing the integral     =*
!*=  of the input quantity over the range of each bin) and are being matched  =*
!*=  onto another set of bins (target grid).  A typical example would be an   =*
!*=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
!*=  vals, that has to be matched onto the working wavelength grid.           =*
!*=  The resulting area in a given bin of the target grid is calculated by    =*
!*=  simply adding all fractional areas of the input data that cover that     =*
!*=  particular target bin.                                                   =*
!*=  Some caution should be used near the endpoints of the grids.  If the     =*
!*=  input data do not span the full range of the target grid, the area in    =*
!*=  the "missing" bins will be assumed to be zero.  If the input data extend =*
!*=  beyond the upper limit of the target grid, the user has the option to    =*
!*=  integrate the "overhang" data and fold the remaining area back into the  =*
!*=  last target bin.  Using this option is recommended when re-gridding      =*
!*=  vertical profiles that directly affect the total optical depth of the    =*
!*=  model atmosphere.                                                        =*
!*-----------------------------------------------------------------------------*
!*=  PARAMETERS:                                                              =*
!*=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
!*=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
!*=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
!*=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
!*=           y-value for bin i (i = 1..NG-1)                                 =*
!*=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
!*=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
!*=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
!*=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
!*=           y-value for bin i (i = 1..N-1)                                  =*
!*=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
!*=           FoldIn = 0 -> No folding of "overhang" data                     =*
!*=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
!*=                         last target bin                                   =*
!*-----------------------------------------------------------------------------*
!*=  EDIT HISTORY:                                                            =*
!*=  06/96  Added FoldIn switch                                               =*
!C R. F. Esswein 020214 Change all REAL declarations to REAL*8
!C                               Use generic names of intrinsic functions.
!C       M. Claire       060806  integrated into Kevin's code
!*-----------------------------------------------------------------------------*
!*= This program is free software;  you can redistribute it and/or modify     =*
!*= it under the terms of the GNU General Public License as published by the  =*
!*= Free Software Foundation;  either version 2 of the license, or (at your   =*
!*= option) any later version.                                                =*
!*= The TUV package is distributed in the hope that it will be useful, but    =*
!*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!*= License for more details.                                                 =*
!*= To obtain a copy of the GNU General Public License, write to:             =*
!*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!*-----------------------------------------------------------------------------*
!*= To contact the authors, please mail to:                                   =*
!*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!*= send email to:  sasha@ucar.edu                                            =*
!*-----------------------------------------------------------------------------*
!*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!*-----------------------------------------------------------------------------*

      IMPLICIT NONE

!* input:
      INTEGER n, ng
      REAL*8 xg(ng)
      REAL*8 x(n), y(n)

      INTEGER FoldIn

!* output:
      REAL*8 yg(ng)

!* local:
      REAL*8 a1, a2, sum
      REAL*8 tail
      INTEGER jstart, i, jl, k
!*_______________________________________________________________________

!* check whether flag given is legal
      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
         WRITE(0,*) '                        Must be 0 or 1'
         STOP
      ENDIF

!* do interpolation

      jstart = 1

      DO 30, i = 1, ng - 1

         yg(i) = 0.
         sum = 0.
         jl = jstart

         IF (jl .LE. n-1) THEN

   20      CONTINUE

             IF (x(jl+1) .LT. xg(i)) THEN
                jstart = jl
                jl = jl+1
                IF (jl .LE. n-1) GO TO 20
             ENDIF

   25      CONTINUE

             IF ((x(jl) .LE. xg(i+1)) .AND. (jl .LE. n-1)) THEN

                a1 = MAX(x(jl),xg(i))
                a2 = MIN(x(jl+1),xg(i+1))

                sum = sum + y(jl) * (a2-a1)/(x(jl+1)-x(jl))
                jl = jl+1
                GO TO 25

             ENDIF

           yg(i) = sum

         ENDIF

   30 CONTINUE


!* if wanted, integrate data "overhang" and fold back into last bin

      IF (FoldIn .EQ. 1) THEN

         jl = jl-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(jl+1)     ! upper limit of last input bin considered

!*        do folding only if grids don't match up and there is more input
         IF ((a2 .GT. a1) .OR. (jl+1 .LT. n)) THEN
           tail = y(jl) * (a2-a1)/(x(jl+1)-x(jl))
           DO k = jl+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF
!*_______________________________________________________________________
      !
      ! RETURN
      END
