
      SUBROUTINE INITPHOTO(flux_txt)
      implicit none !!!!!

      ! Module variables
      ! real*8, allocatable, dimension(:) :: Flux ! Solar flux photons/(cm2 s)
      ! real*8, allocatable, dimension(:) :: wavl, wav, wavu ! wavelength bins
      ! real*8, dimension(17,4) :: alphap ! this
      ! real*8, dimension(17,4) :: beta  ! this
      ! integer, dimension(17) :: nk !this
      ! real*8, dimension(2900) :: SO2HZ ! this

      ! local variables
      character(len=*),intent(in) :: flux_txt
      Integer j,i,k,l



      REAL*8 dum, EFAC, ha
      real*8 columndepth(KJ,NZ)
      real*8 absorbers(kj,nz)
      integer lpolyscount, m
      real*8 rgas, wt, rmg


      ! first thing is to allocate

      ! LGRID = 0


      ! this function defines the wavelength grid and returns
      CALL gridw(nw,wavl,wav,wavu,LGRID)
      ! nw, the number of points on the grid
      ! wavl, a vector of lower grid points  !note wavl(nw+1)=wavu(nw) i.e. final grid point for interpolation
      ! wavu, a vector of upper grid points (i.w. wavl+delta)
      ! wav, a vector of centered values (i.e. (wavl + wavu)/2 )


      ! this subroutine returns the flux data interpolated to the wavelength grid along with
      CALL readflux(flux_txt,nw,wavl,flux)


!C *****  ***** ***** READ THE PHOTOLYSIS DATAFILE ***** ***** *****
!c eventually, this will die...


      OPEN (3,FILE='DATA/photo.dat',STATUS='OLD')

      READ (3,99001)
99001 FORMAT (/)
      DO l = 1 , 108
        READ (3,*)
      ENDDO

      READ (3,99006)
      DO l = 1 , 17
          READ (3,*)
      ENDDO

      READ (3,99006)
      DO l = 1 , 17
        READ (3,*)
      ENDDO

      !c Below SO2HZ is still used. what is the best way to deal?
      READ (3,99006)
      DO l = 1 , 35
        READ (3,99002) so2hz(l) , dum , dum , dum , dum , dum , &
              dum,dum,dum,dum,dum
      ENDDO
99002 FORMAT (5X,11(E8.1,1X))


      READ (3,99007)
      DO l = 1 , 68
        READ (3,*)
      ENDDO

      READ (3,99007)
      DO l = 1 , 68
        READ (3,*)
      ENDDO


      !these lines might need to be kept...
      !what I don't know is where the coefficients come from or how they would scale with a varying wavelength grid...


      do i=1,17
        do j=1,4
          alphap(i,j) = 0.
          beta(i,j) = 0.
        enddo
      enddo

      READ (3,99003)
99003 FORMAT (//////)
      DO l = 1 , 17
        READ (3,99004) nk(l) , (alphap(l,k),k=1,4)
  99004 FORMAT (6X,I1,1X,4F13.5)
        READ (3,99005) (beta(l,k),k=1,4)
  99005 FORMAT (8X,4E13.4/)
      ENDDO


      READ (3,99007)
      DO l = 1 , 68
        READ (3,*)
      ENDDO
99006 FORMAT (////)
99007 FORMAT (//)
      CLOSE (3)

! C *****  ***** ***** END READING PHOTOLYSIS DATAFILE ***** ***** *****
! c at some point do away with all of the above, just taking what's relevant into single datafiles

      ! Here calculate column depth
      do i=1,kj
        do j = 1,nz
          columndepth(i,j) = 0.0
        enddo
      enddo

      RGAS = 8.3143D7             !erg/mol/K
      WT = usol_init(LO2,1)*32. + FCO2*44. + (1.-usol_init(LO2,1)-FCO2)*28. + 0.4 !(g) mean molar mass
      RMG = RGAS/(WT*G)           !gm cm^2/s^2/mol/K  / g *s^2/cm ->  cm/mol/K

      do k=1,kj
        do i=1,nz

!     this gets any SL/IN species that have photorates
          if (photoreac(k).gt.nq) then

            absorbers(k,i)=ABS(SL(INT(photoreac(k)),i))/DEN(I)

!     quasi-hardcoded S8 behavior WARNING
          else if (ISPEC(INT(photoreac(k))).eq.'S8      ') then

            absorbers(k,i)=ABS(USOL_init(INT(photoreac(k)),i))
            if (lpolyscount .eq. 2) absorbers(k,i)=0.0
!           S8R doesn't really exist
!           S8L doesn't really either, but this is how we are rolling for now...
!           we are filling up absorbers for S8R and S8, but S8 doesn't
!           have a cross section
!           so only S8R is actually absorbing photons here in the RT scheme.
            if (i.eq.nz) then
!            print *, k,photoreac(k),ISPEC(INT(photoreac(k))),lpolyscount
              lpolyscount=lpolyscount+1
            endif
          else
            absorbers(k,i)=ABS(USOL_init(INT(photoreac(k)),i))
          endif
        enddo
      enddo

      do k=1,kj
        DO  M=1,NZ1
          I = NZ - M      !run through heights from the top down.
          HA = RMG*0.5*(T(I) + T(I+1))  !scale height RT/MG
! c-mc        DZ = Z(I+1) - Z(I)  !ACK - this is good, but should already exist as a vector
! c-mc in our new scheme DZ(I)=Z(I)-Z(I-1) so DZ(I+1)=Z(I+1)-Z(I)
! C ACK - may have to return when I take this to a variable grid
          EFAC = (1. - EXP(-DZ(I+1)/HA))*DEN(I)*HA     !column depth of each layer
          ! TTOT(I) = TTOT(I+1) + EFAC     !total column depth above height level I
          columndepth(k,I)= columndepth(k,I+1) &
                         + EFAC*SQRT(absorbers(k,i)*absorbers(k,i+1)) !species column depth above level I
        enddo
      enddo


! C read in relevant cross sections and set the Jnumbers,photospec,and photolabel
      do i=1,kj
         photolabel(i)=''
      enddo

!
      j=1
      do i=1,ks   !number of photolysis species
        CALL XS(ISPEC(INT(photospec(i))),nw,wavl,wav,T,DEN,j, &
              columndepth,zy,IO2)
      enddo

      ! fill up stuff....
! c       stop
      ! RETURN
      END
