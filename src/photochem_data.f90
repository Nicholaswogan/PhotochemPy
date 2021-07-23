
module photochem_data
  implicit none
  public
  integer, private, parameter :: real_kind = kind(1.0d0)
  ! this module contains stuff that generally does not change from
  ! one integration to the next. Its the unchanging "data" that is read
  ! in at the very beginning before doing a bunch of simulations.
  
  ! All stuff in here should be considered read-only
  
  integer :: nz ! number of vertical grid points
  integer :: nz1 !nz-1
  integer :: nq ! number of long lived species
  integer :: nq1
  integer :: np ! number of particles
  integer :: nsp ! total number of species
  integer :: nsp2 ! total number of species, including HV and M
  integer :: nr ! number of reactions
  integer, parameter :: nmax = 300 ! max number of reactions a species can be involved in
  integer :: ks ! number of photo species
  integer :: kj ! number of photo reactions
  ! integer, parameter :: kw = 1000 ! max number of wavelength bins
  integer :: nw
  integer, parameter :: naq = 10 !number of aqueous species
  integer, parameter :: nt = 50 !number of temperatures in sulfate/H2O vapor pressure file (DATA/aerosol.table)
  integer, parameter :: nf = 50 !NT=number of pressures per temperature in DATA/aerosol.table
  integer :: lda
  integer :: neq
  
  ! Defined in species.dat
  integer :: iSL ! number of sl species
  real(8) :: background_mu
  character(len=8) :: background_spec
  character(len=8), allocatable, dimension(:) :: ISPEC
  real*8, allocatable, dimension(:) :: redoxstate
  real*8, allocatable, dimension(:) :: mass
  integer LSO2, LH2CO, lh2so4, lso4aer, lh2s ! indexes of a few things
  integer LCO, LH2O, LH2, LCH4, LO2, LH
  integer Ls8aer, Lhcaer, Lhcaer2, ls2, ls3, ls4
  integer lno, lo, LCO2, ln2, ls
  integer, allocatable, dimension(:) :: atomsO
  integer, allocatable, dimension(:) :: atomsH
  integer, allocatable, dimension(:) :: atomsC
  integer, allocatable, dimension(:) :: atomsS
  integer, allocatable, dimension(:) :: atomsN
  integer, allocatable, dimension(:) :: atomsCL
  
  ! Defined in reactions.rx
  Character(len=8), allocatable, dimension(:,:) :: chemj
  integer, allocatable, dimension(:,:) :: jchem
  Character(len=8), allocatable, dimension(:) :: reactype
  real*8, allocatable, dimension(:,:) :: rateparams ! a new one.
  integer, allocatable, dimension(:,:,:) :: iloss
  integer, allocatable, dimension(:,:) :: iprod
  integer, allocatable, dimension(:) :: photoreac
  integer, allocatable, dimension(:) :: photonums
  integer, allocatable, dimension(:) :: photospec
  integer, allocatable, dimension(:) :: NUML, NUMP
  
  ! needed in read_planet.f90
  real*8 :: grav_surf, Fscale, Alb, ztrop,far,R0,P0
  real(8), allocatable :: grav_z(:)
  
  ! needed in read_photochem.f90
  real*8 :: AGL, EPSJ, light_disp_rate, hcdens, zy
  integer :: Lgrid, IO2, ino, frak, ihztype
  logical :: lightning
  logical :: rainout_on
  logical :: H2O_strat_condensation
  real(8) :: confac, rhcold
  logical :: fix_water_in_troposphere
  logical :: use_manabe
  real(8) :: relative_humidity
  logical :: estimate_CO2_photo_above_grid
  
  ! needed in subroutine photgrid (in photgrid.f90)
  real*8, allocatable, dimension(:) :: z ! altitude of middle of grid
  real*8, allocatable, dimension(:) :: dz ! Delta_z of each altitude grid
  real(8) :: top_atmos, bottom_atmos
  integer JTROP
  
  ! needed in initphoto.f90.
  ! real*8, dimension(kw) :: Flux ! Solar flux photons/(cm2 s)
  ! real*8, dimension(kw) :: wavl, wav, wavu ! wavelength bins
  real*8, allocatable, dimension(:) :: Flux ! Solar flux photons/(cm2 s)
  real*8, allocatable, dimension(:) :: wavl, wav, wavu ! wavelength bins
  real*8, allocatable, dimension(:,:,:) :: sq ! cross sections * qy
  
  ! needed in initmie.f90
  real*8, dimension(51) :: Rstand
  ! real*8, dimension(kw,51) :: W0HC
  ! real*8, dimension(kw,51) :: GHC, QEXTHC
  real*8, allocatable, dimension(:,:) :: W0HC
  real*8, allocatable, dimension(:,:) :: GHC, QEXTHC
  
  ! needed in Aertab.f90
  real*8, allocatable, dimension(:,:) :: VH2O
  real*8, allocatable, dimension(:,:) :: VH2SO4
  real*8, dimension(nf) :: ftab
  
  ! atmosphere_txt contents
  real(8), allocatable :: usol_file(:,:), rpar_file(:,:), &
                          wfall_file(:,:), aersol_file(:,:)
  real(8), allocatable :: T_file(:), edd_file(:), z_file(:)
  integer :: nzf

end module