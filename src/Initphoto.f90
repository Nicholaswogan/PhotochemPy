
      SUBROUTINE INITPHOTO(err)
        use photochem_data, only: kj, ks, nw, &
                                  photospec, ispec, &
                                  wavl, wav
        use photochem_vars, only: T
      implicit none !!!!!

      ! local variables
      ! character(len=*),intent(in) :: flux_txt
      character(len=err_len) :: err
      Integer j,i



      ! REAL*8 dum !, EFAC, ha
      ! real*8 columndepth(KJ,NZ)
      ! real*8 absorbers(kj,nz)
      ! integer lpolyscount
      err = ''


      ! lpolyscount = 0


      ! first thing is to allocate

      ! LGRID = 0


      ! this function defines the wavelength grid and returns
      ! CALL gridw(nw,wavl,wav,wavu,LGRID)
      ! nw, the number of points on the grid
      ! wavl, a vector of lower grid points  !note wavl(nw+1)=wavu(nw) i.e. final grid point for interpolation
      ! wavu, a vector of upper grid points (i.w. wavl+delta)
      ! wav, a vector of centered values (i.e. (wavl + wavu)/2 )


      ! this subroutine returns the flux data interpolated to the wavelength grid along with
      ! CALL readflux(flux_txt,nw,wavl,flux)


!C *****  ***** ***** READ THE PHOTOLYSIS DATAFILE ***** ***** *****
!c eventually, this will die...


      ! OPEN (3,FILE=trim(rootdir)//'DATA/photo.dat',STATUS='OLD')
! 
!       READ (3,99001)
! 99001 FORMAT (/)
!       DO l = 1 , 108
!         READ (3,*)
!       ENDDO
! 
!       READ (3,99006)
!       DO l = 1 , 17
!           READ (3,*)
!       ENDDO
! 
!       READ (3,99006)
!       DO l = 1 , 17
!         READ (3,*)
!       ENDDO

      !c Below SO2HZ is still used. what is the best way to deal?
!       READ (3,99006)
!       DO l = 1 , 35
!         READ (3,99002) so2hz(l) , dum , dum , dum , dum , dum , &
!               dum,dum,dum,dum,dum
!       ENDDO
! 99002 FORMAT (5X,11(E8.1,1X))
! 
! 
!       READ (3,99007)
!       DO l = 1 , 68
!         READ (3,*)
!       ENDDO
! 
!       READ (3,99007)
!       DO l = 1 , 68
!         READ (3,*)
!       ENDDO


      !these lines might need to be kept...
      !what I don't know is where the coefficients come from or how they would scale with a varying wavelength grid...


!       do i=1,17
!         do j=1,4
!           alphap(i,j) = 0.
!           beta(i,j) = 0.
!         enddo
!       enddo
! 
!       READ (3,99003)
! 99003 FORMAT (//////)
!       DO l = 1 , 17
!         READ (3,99004) nk(l) , (alphap(l,k),k=1,4)
!   99004 FORMAT (6X,I1,1X,4F13.5)
!         READ (3,99005) (beta(l,k),k=1,4)
!   99005 FORMAT (8X,4E13.4/)
!       ENDDO
! 
! 
!       READ (3,99007)
!       DO l = 1 , 68
!         READ (3,*)
!       ENDDO
! 99006 FORMAT (////)
! 99007 FORMAT (//)
!       CLOSE (3)

! C *****  ***** ***** END READING PHOTOLYSIS DATAFILE ***** ***** *****
! c at some point do away with all of the above, just taking what's relevant into single datafiles

! C read in relevant cross sections and set the Jnumbers,photospec,and photolabel

!
      j=1
      do i=1,ks   !number of photolysis species
        CALL XS(ISPEC(INT(photospec(i))),nw,wavl,wav,T,j,err)
        if (len_trim(err) /= 0) return
      enddo
      if (j /= kj + 1) then
        err = 'Problem reading in xsections. It is hard to tell what is '// &
              'wrong because cross sections are handled in a terrible way (not my fault).'
        return
      endif
      

      ! fill up stuff....
! c       stop
      ! RETURN
      END
