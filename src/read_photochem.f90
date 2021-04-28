
subroutine read_photochem(input_photochem,err)
  use photochem_data, only: AGL, EPSJ, prono, hcdens, zy, &
                            Lgrid, IO2, ino, frak, ihztype, np
  implicit none

  ! module variables
  ! real*8 :: AGL, EPSJ, prono, hcdens, zy
  ! integer :: Lgrid, IO2, frak, ihztype

  ! local variables
  character(len=*), intent(in) :: input_photochem
  character(len=err_len), intent(out) :: err
  character(len=11) AA ! dummy
  integer :: iseason, IZYO2, icouple

  ! open reactions.rx
  err = ''
  open(231, file=trim(input_photochem),status='OLD')

  READ(231,*)
  READ(231,*)
  READ(231,*)
  READ(231,*)

  READ(231,*)AA, AGL
  READ(231,*)AA, ISEASON
  READ(231,*)AA, IZYO2
  READ(231,*)AA, LGRID
  READ(231,*)AA, IO2
  READ(231,*)AA, INO
  READ(231,*)AA, EPSJ
  READ(231,*)AA, PRONO
  READ(231,*)AA, frak
  READ(231,*)AA, HCDENS
  READ(231,*)AA, ICOUPLE
  READ(231,*)AA, ihztype
  READ(231,*)AA, ZY

  close(231)
  
  if ((frak == 1) .and. (np < 3)) then
    err = 'frak must be 0 if no hydrocarbons are present. See '//trim(input_photochem)
    return
  endif
  
  ! need more errors for io2 and ino but later


end subroutine
