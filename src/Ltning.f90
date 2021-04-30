
subroutine ltning(nq, nz, usol, &
                  zapNO, zapO2, proNOP, zapCO, zapH2, zapO)
  use photochem_data, only: lCH4, lco, lh2, lo2, lh2o, P0, prono, &
                            ln2, lco2, background_spec
  !  3-21-06   This subroutine does not work
  !  it does not conserve redc  it needs to be replaced by something that isn't wrong
  !-mc 4-29-06  it conserves redox now and is only as wrong as chamedies, which is probably wrong...
  implicit none
  
  integer, intent(in) :: nq,nz
  real*8, dimension(nq,nz), intent(in) :: usol
  
  real(8), intent(out) :: zapNO, zapO2, proNOP, zapCO, zapH2, zapO

  real*8 ak1,ak2,ak3,ak4
  real*8 pbar, xs, x, err, alpha, b, c, ct, fco
  real*8 fh2o, h2t, o2t, pch4, pco, pco2, ph2, ph2o
  real*8 pn2, pno, pnonow, po, prnox, x2, x_orig, zapco2
  real*8 zaph2o, po2, a, beta

  real*8 fx, fpx
  real(8) :: fn2, fco2
  integer n, ns

  !     THIS SUBROUTINE CALCULATES LIGHTNING PRODUCTION RATES FOR O2
  !     AND NO IN AN N2-O2-CO2-H2O ATMOSPHERE USING THE EQUATIONS IN
  !     THESIS APPENDIX C.
  !     this is incomprehensible without more hints
  !
  !     EQUILIBRIUM CONSTANTS AT 3500 K
  ak1 = .103                ! = pNo/sqrt(pN2*pO2)  dimensionless
  ak2 = .619                ! = pCO2/PCO/sqrt(pO2) bars
  ak3 = 5.3                 ! = pH20/pH2/sqrt(pO2) bars
  ak4 = 0.22                ! = pO/sqrt(pO2)

  !      Pbar = P0/1.013E6   !converting from pascals to bars
  pbar = p0/1.013   !converting from pascals to bars

  fco = USOL(lco,1)
  ph2o = USOL(lh2o,1)*pbar
  ph2 = USOL(lh2,1)*pbar     ! I am having problems here if PH2 gets too big
  pch4 = USOL(lch4,1)*pbar

  fh2o = ph2o/pbar
  if (background_spec == 'N2') then
    fn2 = 1 - sum(usol(:,1))
  else
    fn2 = usol(ln2,1)
  endif 
  if (background_spec == 'CO2') then
    fco2 = 1 - sum(usol(:,1))
  else
    fco2 = usol(lco2,1)
  endif

  pn2 = fn2*pbar
  po2 = USOL(lo2,1)*pbar
  pco2 = fco2*pbar
  pco = fco*pbar
  o2t = po2 + pco2 + 0.5*(ph2o+pco)
  h2t = ph2 + ph2o
  ct = pco2 + pco

  alpha = ak2*SQRT(o2t)
  beta = ak3*SQRT(o2t)
  a = (ak1*SQRT(pn2)+ak4)/(2.*SQRT(o2t))      !K_T in jim's thesis
  b = 0.5*ct/o2t
  c = 0.5*h2t/o2t

  !     INITIAL GUESS FOR XO2 AT 3500 K
  x = 0.1 + 0.9*po2/o2t + 0.2*pco2/o2t ! this initial guess seems to be causing trouble sometimes
  x_orig = x

  !
  !     NEWTON STEP
  do n = 1 , 20
  50      ns = n
     xs = x
     x2 = SQRT(x)


     fx = x + a*x2 - b/(1.+alpha*x2) - c/(1.+beta*x2) &
     + 2.*b + c - 1.
   ! ^^^ C18) in JFK thesis

     fpx = 1. + (a+alpha*b/(1.+alpha*x2)**2+beta*c/ &
     (1.+beta*x2)**2)/(2.*x2)
     x = x - fx/fpx
                  !GNA - bad initial guess for X (if X too big) is making this occasionally go negative for me, which make PO2 negative, which wrecks havoc in the sqrt(PN2*PO2) below...
     if ( x.lt.0 ) then
        x = x_orig*0.95
                       !decrease initial guess slightly until small enough
        x_orig = x
        PRINT * , 'fixing x in ltning.f. new x = ' , x
        PRINT * , 'old x = ' , x/0.95
        GOTO 50
     endif

     err = ABS((x-xs)/x)

     IF ( err.LT.1.E-5 ) GOTO 100
  enddo
100  po2 = x*o2t
  pno = ak1*SQRT(pn2*po2)
  pco = ct/(1.+ak2*SQRT(po2))     !added 3-20-06
  ph2 = h2t/(1.+ak3*SQRT(po2))     !added 3-20-06

  !-mc
  po = ak4*SQRT(po2)   !added 4-28-06
  !-mc

  !      print *, PO2, 0.5*(PCO+PH2-PO-PNO)

  !      stop
  !
  !     SCALE AGAINST ESTIMATED COLUMN PRODUCTION OF NO IN THE PRESENT-
  !     DAY TROPOSPHERE.  DISTRIBUTE PRODUCTION OVER LOWEST 6 KM.
  !
  !   COLUMN-INTEGRATED NO PRODUCTION RATE IN PRESENT ATMOSPHERE IS
  !   EQUAL TO PRONO
  prnox = prono/7.942E5
  pnonow = 3.574E-2
  pronop = prono*pno/pnonow

  zapno = prnox*pno/pnonow
  zapo2 = zapno*po2/pno

  zapco = zapno*pco/pno   ! added 3-20-06
  zaph2 = zapno*ph2/pno   ! added 3-20-06

  zapo = zapno*po/pno     ! added 4-28-06
  zapco2 = zapno*pco2/pno
  zaph2o = zapno*ph2o/pno



  !      print *, ZAPNO, ZAPO2,ZAPCO,ZAPH2,ZAPO,ZAPCO,ZAPH2O
  !      stop

  !  - add below
  !       WRITE (14,99001) po2 , pno , zapo2 , zapno , prono , pronop ,     &
  !                      & zapco , zaph2
  ! 99001 FORMAT (/1X,'PO2=',1PE10.3,2X,'PNO=',1PE10.3,2X,'ZAPO2=',E10.3,2X,&
  !              &'ZAPNO =',E10.3,2X,'PRONO =',E10.3,2X,'PRONOP =',E10.3,2X,&
  !              &'ZAPCO =',E10.3,2X,'ZAPH2 =',E10.3)
  !       WRITE (14,99002) ns , x , err
  ! 99002 FORMAT (1X,'N=',I2,5X,'X=',1PE12.5,2X,'ERR=',1PE12.5)
end subroutine
