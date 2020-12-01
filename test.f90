
      program test
        implicit none

        character(len=1000) :: line

667     FORMAT(58X,E9.2,3X,F8.2)            !for two body reaction rates
        open(4, file='input/reactions.rx',status='OLD')

        read(4,*) line



        
        close(4)




      end program
