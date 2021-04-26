
      subroutine read_reactions(reactions_rx)
        use photochem_data, only: kj, ks, nr, nsp, nsp2, nmax, &
                                  reactype, rateparams, chemj, jchem, &
                                  nump, numl, iprod, iloss, ispec, &
                                  photospec, photoreac, photonums, &
                                  atomsO, atomsH, atomsC, atomsS, atomsN, atomsCl
        implicit none

        ! module variables
        ! Character(len=8), allocatable, dimension(5,:) :: chemj
        ! integer, allocatable, dimension(:,:) :: jchem
        ! Character(len=8), allocatable, dimension(:) :: reactype
        ! real*8, allocatable, dimension(:,:) :: rateparams ! a new one
        ! integer, allocatable, dimension(:,:,:) :: iloss
        ! integer, allocatable, dimension(:,:) :: iprod
        ! integer, allocatable, dimension(:) :: photoreac
        ! integer, allocatable, dimension(:) :: photonums
        ! integer, allocatable, dimension(:) :: photospec
        ! integer, allocatable, dimension(:) :: NUML, NUMP

        ! local variables
        character(len=*) :: reactions_rx
        character(len=50) :: fmat1
        character(len=50) :: fmat2
        integer :: i, m, j, l, ierr, n, k
        integer :: rhOcount, rhHcount, rhCcount, rhScount
        integer :: rhNcount, rhCLcount, numprod, bad
        integer, allocatable, dimension(:) :: testvec
        integer :: jcount, jrcount, juniq
!f2py   intent(in) :: reactions_rx

        ! allocate memeory
        allocate(testvec(nr))

        ! open reactions.rx
        open(9, file=trim(reactions_rx),status='OLD')
        ! formats
        fmat1 = '(A8,2X,A8,2X,A8,2X,A8,2X,A8)'
        fmat2 = '(48X,A5)'
        read (9,fmat1) CHEMJ
        rewind(9)
        read(9,fmat2) REACTYPE
        rewind(9)

        ! zero out! Always.
        do i=1,nr
          do j=1,4
            rateparams(j,i) = 0.0
          enddo
        enddo
        ! now the rate parameters.
        fmat1 = '(58X,E9.2,3X,F8.2)'
        fmat2 = '(58X,E9.2,3X,E9.2,2X,2F5.2)'
! 667     FORMAT(58X,E9.2,3X,F8.2)            !for two body reaction rates
! 668     FORMAT(58X,E9.2,3X,E9.2,2X,2F5.2)   !for three body reaction rates
        do i=1,nr
          if (REACTYPE(i) .EQ. '2BODY') then
            read(9,fmat1) rateparams(1,i), rateparams(2,i)
          else if (REACTYPE(i) .EQ. '3BODY') then
            read(9,fmat2) rateparams(1,i), rateparams(2,i), &
                          rateparams(3,i), rateparams(4,i)
          else if (REACTYPE(i) .EQ. 'WEIRD') then
            ! If it is weird, then there is nothing to read. Just skip line
            read(9,*)
          else if (REACTYPE(i) .EQ. 'PHOTO') then
            ! if photo then skip this reaction as well
            read(9,*)
          ! else
          !   ! just skip the line!
          !   read(9,*)
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
                  IERR = J
                  print *, ISPEC
                  print *, 'ispec(i)', ISPEC(i)
                  print *, (CHEMJ(L,J),L=1,5)
                  ! quit; error in reactions
                  PRINT '(1X,"error in reaction",I3)',IERR
                  stop
                endif
              enddo
            endif
          enddo
        enddo

        ! very important to zero out!
        do i = 1,nsp
          NUMP(i) = 0
          NUML(i) = 0
          do j = 1,nmax
            iprod(i,j) = 0
            do k = 1,2
              iloss(k,i,j) = 0
            enddo
          enddo
        enddo

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
                print *,'NMAX exceeded ILOSS'
                stop
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
                print *,'NMAX exceeded IPROD'
                stop
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
            print *, 'bad mass balance in reaction',i
            print *, (CHEMJ(j,i),j=1,5)
            print *, numprod
            !the problem is either in the .rx file or the species.dat file
            print *, rhNcount,atomsN(JCHEM(1,i)),atomsN(JCHEM(2,i))
            stop
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
          print *, 'discrepency between unique photolysis reactions/ks'
          print *, juniq-1, ks
          stop
        endif
        !
        if (SUM(INDEX(REACTYPE,'PHOTO')) .NE. kj) then
          print *,'discrepency between number of photo reactions and kj'
          print *, SUM(INDEX(REACTYPE,'PHOTO')), kj
          stop
        endif


        ! deallocate
        ! deallocate(testvec) ! the program doesn't like this at all


      end subroutine
