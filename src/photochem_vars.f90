module photochem_vars
  implicit none
  public
  integer, private, parameter :: real_kind = kind(1.0d0)
  ! this module contains stuff that can change, or can potentially change,
  ! between photochemical integrations (boundary conditions)
  
  character(len=500) :: rootdir = 'PhotochemPy/'
  logical :: verbose = .true.
  
  ! Defined in species.dat
  integer, allocatable, dimension(:) :: LBOUND
  real*8, allocatable, dimension(:) :: VDEP0
  real*8, allocatable, dimension(:) :: VDEP
  real*8, allocatable, dimension(:) :: FIXEDMR
  real*8, allocatable, dimension(:) :: distflux
  real*8, allocatable, dimension(:) :: SGFLUX
  real*8, allocatable, dimension(:) :: distheight
  integer, allocatable, dimension(:) :: MBOUND
  real*8, allocatable, dimension(:) :: SMFLUX
  real*8, allocatable, dimension(:) :: VEFF0
  real*8, allocatable, dimension(:) :: VEFF
  
  ! for output
  double precision :: redox_factor
  
  ! needed in read_atmosphere.f90
  real*8, allocatable, dimension(:,:) :: usol_init ! initial atmospheric composition
  real*8, allocatable, dimension(:,:) :: aersol_init ! aersol parameter
  real*8, allocatable, dimension(:,:) :: wfall_init ! aersol parameter
  real*8, allocatable, dimension(:,:) :: rpar_init ! aersol parameter
  real*8, allocatable, dimension(:) :: den ! total number density vs altitude
  real*8, allocatable, dimension(:) :: T ! Temperature vs altitude
  real*8, allocatable, dimension(:) :: EDD ! Eddy diffusion coefficients
  
  ! needed in Densty.f90
  real*8, allocatable, dimension(:) :: Press ! pressure in dynes
  real*8, allocatable, dimension(:) :: P ! pressure in bars
  
  ! needed in PhotSatrat.f90
  real*8, allocatable, dimension(:) :: h2osat
  
  ! needed in integrate.f90
  real*8, allocatable, dimension(:,:) :: usol_out
  real*8, allocatable, dimension(:) :: flow
  real*8, allocatable, dimension(:,:) :: fluxo
  
  ! cvode
  integer :: max_cvode_steps = 100000000
  
  
end module