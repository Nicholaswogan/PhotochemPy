
      SUBROUTINE CHEMPL(D,Xp,Xl,K)

      IMPLICIT NONE

      ! module variables


      ! local variables
      real*8, dimension(nz), intent(out) :: XP
      real*8, dimension(nz), intent(out) :: XL
      real*8, dimension(nsp2,nz), intent(in) :: D
      integer, intent(in) :: k

      integer i, l, m, n, j, nl, nprx





!
!   THIS SUBROUTINE CALCULATES CHEMICAL PRODUCTION AND LOSS RATES
!   USING THE INFORMATION IN THE MATRICES JCHEM, ILOSS, AND IPROD.
!   CALLED BY SUBROUTINE DOCHEM.
!
      DO i = 1 , nz
         Xp(i) = 0.d0
         Xl(i) = 0.d0
      ENDDO
!
!   LOSS FREQUENCY XL
      nl = NUML(K)        !chempl is called with given species (K)- NUML is how many reactions involve K

      DO l = 1 , nl
         j = ILOSS(1,K,l) !reaction number for loss process
         m = ILOSS(2,K,l) !reactant number for loss process
         DO i = 1 , nz
            Xl(i) = Xl(i) + A(j,i)*D(m,i)
                                    !rate*density of reactant 1
         ENDDO

      ENDDO


!XL has units of s^-1
!
!   PRODUCTION RATE XP
      nprx = NUMP(K)

      DO l = 1 , nprx
         j = IPROD(K,l)
                       !reaction number
         m = JCHEM(1,j)
                       !reactant 1 for reaction number J
         n = JCHEM(2,j)
                       !reactant 2 for reaction number J
         DO i = 1 , nz !loop over height (also for each reaction)
            Xp(i) = Xp(i) + A(j,i)*D(m,i)*D(n,i)


                                           !rate*density1*density2
         ENDDO
      ENDDO
!HV has a density of 1, so rates make up for this
!XP in units of mol/cm^3/s



      END
