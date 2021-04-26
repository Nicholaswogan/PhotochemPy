
      SUBROUTINE AQUEOUS(X, naq, F, I, RRRR, nz)
        use photochem_wrk, only: alpharain, co2aq, h2cog0, hh2co, &
                                 hplus, hso2, so2g0, so4_2
      implicit none

      ! module variables

      ! local variables
      integer, intent(in) :: naq
      real*8, dimension(naq), intent(in) :: X
      integer, intent(in) :: I
      real*8, dimension(6,nz), intent(in) :: RRRR
      integer,intent(in) :: nz

      real*8, dimension(naq),intent(out) :: F
      integer LSO2g,LH2COg,LSO2aq,LH2COaq,LHCO3_,LCO3_2,LHSO3_
      integer LSO3_2,LH2COSO3,LOH_, k
      real*8, dimension(nz) :: R4, R5, R6, R7, R8, R9


      r4(i) = RRRR(1,i)
      r5(i) = RRRR(2,i)
      r6(i) = RRRR(3,i)
      r7(i) = RRRR(4,i)
      r8(i) = RRRR(5,i)
      r9(i) = RRRR(6,i)



      DATA lso2g , lh2cog , lso2aq , lh2coaq , lhco3_ , lco3_2 ,        &
         & lhso3_ , lso3_2 , lh2coso3 , loh_/1 , 2 , 3 , 4 , 5 , 6 , 7 ,&
         & 8 , 9 , 10/
!
!     THIS SUBROUTINE DOES THE AQUEOUS PHASE CHEMISTRY DESCRIBED IN
!     SUBROUTINE RAINOUT
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

      IF ( naq.GT.0 ) THEN
         hplus = X(lhco3_) + X(lhso3_) + X(lh2coso3) + X(loh_)          &
               & + 2.*(X(lco3_2)+X(lso3_2)+so4_2)
         DO k = 1 , naq
            F(k) = 0.0
            IF ( k.EQ.1 ) F(k) = X(lso2g) - so2g0 +                     &
                               & alpharain*(X(lso2aq)+X(lhso3_)         &
                               & +X(lso3_2)+X(lh2coso3))
            IF ( k.EQ.2 ) F(k) = X(lso2aq) - hso2*X(lso2g)
            IF ( k.EQ.3 ) F(k) = X(lh2coaq) - hh2co*X(lh2cog)
            IF ( k.EQ.4 ) F(k) = X(lhco3_)*hplus - R4(I)*co2aq
            IF ( k.EQ.5 ) F(k) = X(lhso3_)*hplus - R5(I)*X(lso2aq)
            IF ( k.EQ.6 ) F(k) = X(lco3_2)*hplus - R6(I)*X(lhco3_)
            IF ( k.EQ.7 ) F(k) = X(lso3_2)*hplus - R7(I)*X(lhso3_)
            IF ( k.EQ.8 ) F(k) = X(lh2coso3) - R8(I)*X(lh2coaq)         &
                               & *X(lhso3_)
            IF ( k.EQ.9 ) F(k) = X(loh_)*hplus - R9(I)
            IF ( k.EQ.10 ) F(k) = X(lh2cog) - h2cog0 +                  &
                                & alpharain*(X(lh2coaq)+X(lh2coso3))
         ENDDO
!
      ENDIF
      END
