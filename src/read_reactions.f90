
subroutine read_reactions(reactions_rx, err)
  use photochem_data, only: nq, kj, ks, nr, nsp, nsp2, nmax, &
                            reactype, rateparams, chemj, jchem, &
                            nump, numl, iprod, iloss, ispec, &
                            photospec, photoreac, photonums, &
                            atomsO, atomsH, atomsC, atomsS, atomsN, atomsCl
  implicit none

  ! local variables
  character(len=*), intent(in) :: reactions_rx
  character(len=err_len), intent(out) :: err
  character(len=50) :: fmat1, fmat11, temp1, temp2, temp3, temp4
  character(len=50) :: fmat2, fmat22, temp
  integer :: i, m, j, n, k, io
  integer :: rhOcount, rhHcount, rhCcount, rhScount
  integer :: rhNcount, rhCLcount, numprod, bad
  integer, allocatable, dimension(:) :: testvec
  integer :: jcount, jrcount, juniq
  err = ''

  ! allocate memeory
  allocate(testvec(nr))

  ! open reactions.rx
  open(9, file=trim(reactions_rx),status='OLD')
  ! formats
  fmat1 = '(A8,2X,A8,2X,A8,2X,A8,2X,A8)'
  fmat2 = '(48X,A5)'
  do i=1,nr
    read (9,fmat1,iostat=io) CHEMJ(:,i)
    if (io /= 0) then
      write(temp,'(i5)') i
      err = 'Problem reading in reaction at line '//trim(adjustl(temp))// &
            ' in the file '//trim(reactions_rx)
      return
    endif
  enddo
  rewind(9)
  do i=1,nr
    read(9,fmat2,iostat=io) REACTYPE(i)
    if (io /= 0) then
      write(temp,'(i5)') i
      err = 'Problem reading in reaction type at line '//trim(adjustl(temp))// &
            ' in the file '//trim(reactions_rx)
      return
    endif
  enddo
  rewind(9)

  ! zero out! Always.
  rateparams = 0.d0
  ! now the rate parameters.
  fmat1 = '(58X,E9.2,3X,F8.2)'
  fmat11 = '(58X,A9,3X,A8)'
  fmat2 = '(58X,E9.2,3X,E9.2,2X,2F5.2)'
  fmat22 = '(58X,A9,3X,A9,2X,2A5)'
  do i=1,nr
    if (REACTYPE(i) .EQ. '2BODY') then
      read(9,fmat1,iostat = io) rateparams(1,i), rateparams(2,i)
      backspace(9)
      read(9,fmat11,iostat = io) temp1, temp2
      if ((len_trim(temp1) == 0) .or. (len_trim(temp2) == 0) .or. &
          (io /= 0)) then
        write(temp,'(i5)') i
        err = 'Problem reading in reaction rate parameters at line '// &
              trim(adjustl(temp))//' in the file '//trim(reactions_rx)
        return
      endif
    else if (REACTYPE(i) .EQ. '3BODY') then
      read(9,fmat2,iostat = io) rateparams(1,i), rateparams(2,i), &
                                rateparams(3,i), rateparams(4,i)
      backspace(9)
      read(9,fmat22,iostat = io) temp1, temp2, temp3, temp4
      
      if ((len_trim(temp1) == 0) .or. (len_trim(temp2) == 0) .or. &
          (len_trim(temp3) == 0) .or. (len_trim(temp4) == 0) .or. &
          (io /= 0)) then
        write(temp,'(i5)') i
        err = 'Problem reading in reaction rate parameters at line '// &
              trim(adjustl(temp))//' in the file '//trim(reactions_rx)
        return
      endif
    else if (REACTYPE(i) .EQ. 'ELEMT') then ! elementary reaction from Cantera
      read(9,'(58X,E9.2,3X,F8.2,3X,F8.2)',iostat = io) rateparams(1,i), rateparams(2,i), rateparams(3,i)
      backspace(9)
      read(9,'(58X,A9,3X,A8,3X,A8)',iostat = io) temp1, temp2, temp3
      if ((len_trim(temp1) == 0) .or. (len_trim(temp2) == 0) .or. &
          (len_trim(temp3) == 0) .or. &
          (io /= 0)) then
        write(temp,'(i5)') i
        err = 'Problem reading in reaction rate parameters at line '// &
              trim(adjustl(temp))//' in the file '//trim(reactions_rx)
        return
      endif
    else if (REACTYPE(i) .EQ. 'WEIRD') then
      ! If it is weird, then there is nothing to read. Just skip line
      read(9,*)
    else if (REACTYPE(i) .EQ. 'PHOTO') then
      ! if photo then skip this reaction as well
      read(9,*)
    else
      write(temp,'(i5)') i
      err = 'Reaction type '//trim(reactype(i))//' on line '//trim(adjustl(temp))// &
            ' in the file '//trim(reactions_rx)//' is not a valid reaction type'
      return
    endif
  enddo
  close(9)

  ! ****  JCHEM has species numbers; CHEMJ is corresponding characters
  ! ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
  DO J=1,NR
    DO M=1,5
      JCHEM(M,J) = 0 ! this is important!
      IF (CHEMJ(M,J).EQ.' ') then
        cycle
      else
        DO I=1,NSP2
          if (CHEMJ(M,J).eq.ISPEC(I)) then
            JCHEM(M,J) = I
            exit
          elseif (i.eq.nsp2) then
            write(temp,'(i5)') j
            err = 'Species '//trim(CHEMJ(M,J))//' in reaction '//trim(adjustl(temp))// &
                  ' is not in the list of species.'
            return
          endif
        enddo
      endif
    enddo
  enddo

  ! very important to zero out!
  nump = 0
  numl = 0
  iprod = 0
  iloss = 0

  ! ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
  DO M=1,2
  ! so N=2, then 1
    N = 3-M
    DO J=1,NR
      !so I = JCHEM(1,NR) then JCEHM(2,NR)
      I = JCHEM(M,J)
      !skips 0 (i.e. nothing) and NSP1 (HV)
      IF(I.LT.1.OR.I.GT.NSP) then
        cycle
      else
        !counter for how many reactions species I is involved with
        NUML(I) = NUML(I) + 1
        K = NUML(I)
        IF(NUML(I).GT.NMAX) then
          err = trim(ispec(I))//' is involved in too many reactions.'
          return
        else
          !ILOSS(1,species in pos 1, species reac#)
          ! -then ILOSS(1,spec in pos 2, reac#)= global reaction #
          ILOSS(1,I,K) = J
          ! ILOSS(1,species in pos 1, species reac#)
          ! then ILOSS(1,spec in pos 2, reac#)= other species
          ILOSS(2,I,K) = JCHEM(N,J)
        endif
      endif
    enddo
  enddo

  DO M=3,5
    DO J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) then
        cycle
      else
        NUMP(I) = NUMP(I) + 1
        K = NUMP(I)
        IF(NUML(I).GT.NMAX) then
          err = trim(ispec(I))//' is involved in too many reactions.'
          return
        else
          IPROD(I,K) = J
        endif
      endif
    enddo
  enddo

  ! check mass balance of chemical reactions
  do i=1,nr

    rhOcount=0
    rhHcount=0
    rhCcount=0
    rhScount=0
    rhNcount=0
    rhCLcount=0

    !assume 3 products unless..
    numprod=3
    if (JCHEM(5,i).EQ.0) numprod=2
    if (JCHEM(4,i).EQ.0) numprod=1

    ! This loop counts up the mass on the right hand side of the .rx
    do j=0,numprod-1
     rhOcount=rhOcount+atomsO(JCHEM(3+j,i))
     rhHcount=rhHcount+atomsH(JCHEM(3+j,i))
     rhCcount=rhCcount+atomsC(JCHEM(3+j,i))
     rhScount=rhScount+atomsS(JCHEM(3+j,i))
     rhNcount=rhNcount+atomsN(JCHEM(3+j,i))
     rhCLcount=rhCLcount+atomsCL(JCHEM(3+j,i))
    enddo


    bad=0
    if (rhOcount.ne.atomsO(JCHEM(1,i))+atomsO(JCHEM(2,i))) bad=1
    if (rhHcount.ne.atomsH(JCHEM(1,i))+atomsH(JCHEM(2,i))) bad=1
    if (rhCcount.ne.atomsC(JCHEM(1,i))+atomsC(JCHEM(2,i))) bad=1
    if (rhScount.ne.atomsS(JCHEM(1,i))+atomsS(JCHEM(2,i))) bad=1
    if (rhNcount.ne.atomsN(JCHEM(1,i))+atomsN(JCHEM(2,i))) bad=1
    if (rhCLcount.ne.atomsCL(JCHEM(1,i))+atomsCL(JCHEM(2,i))) then
      bad=1
    endif

    if (bad .eq. 1) then
      write(temp,'(i5)') i
      err = 'Bad mass balance in reaction '//trim(adjustl(temp))
      return
   endif
  enddo   !end mass balance check


  ! PROCESS the photolysis reactions
  !
  ! this next little bit creates:
  ! photoreac(kj) - an array of species numbers for each photolysis reaction.
  ! used in absorbers/columndepth computations
  ! photospec(ks) - the unique photolysis reaction numbers
  ! (i.e. unique elements of photoreac)
  ! used in Initphoto.f to fill up sq, the cross section vector
  ! photonums(kj) - the reaction number of each photolysis reaction
  ! used in Photo.f to fill up the A vector of rates


  ! A vector of length nr with 1's in the location where photolysis reactions are
  do i=1,nr
    testvec(i) = 0
  enddo
  do i=1,ks
    photospec(i) = 0
  enddo
  do i=1,kj
    photoreac(i) = 0
    photonums(i) = 0
  enddo
  testvec=INDEX(REACTYPE,'PHOTO')

  jcount=1
  jrcount=1
  juniq=1
  do i=1,nr
    if (testvec(i).eq.1.) then
      if ((JCHEM(1,i) > nq) .and. ((CHEMJ(1,i) /= 'CO2') .and. (CHEMJ(1,i) /= 'N2'))) then
        err = 'Short-lived species '//trim(CHEMJ(1,i))//' can not have photolysis reactions. '// &
              'Only long-lived species can have photolysis reactions.'
        return
      endif
      ! Captures the species number of each photo reaction
      photoreac(jrcount)=JCHEM(1,i)
      ! Captures the reaction number of each photoreaction
      photonums(jrcount)=i
      jrcount=jrcount+1
      ! Captures the unique photolysis species in photospec
      if (juniq.eq.1) then
        photospec(juniq)=JCHEM(1,i)
        juniq=juniq+1
      else
        bad=0
        do m=1,ks
          if (JCHEM(1,i).eq.photospec(m)) bad=1
        enddo
        if (bad.eq.0) then
          photospec(juniq)=JCHEM(1,i)
          juniq=juniq+1
        endif
      endif
    endif
  jcount=jcount+1
  enddo


  if (juniq-1.ne.ks) then
    err = 'There is some problem with the photolysis reactions. The problem is not clear.'
    return
  endif
  !
  if (SUM(INDEX(REACTYPE,'PHOTO')) .NE. kj) then
    err = 'There is some problem with the photolysis reactions. The problem is not clear.'
    return
  endif


end subroutine
