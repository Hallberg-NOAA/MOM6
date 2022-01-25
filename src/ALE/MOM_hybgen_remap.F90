!> This module contains the hybgen remapping routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_remap

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains,         only : pass_var
use MOM_debugging,       only : hchksum, uvchksum
use MOM_EOS,             only : EOS_type, calculate_density, calculate_density_derivs
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_open_boundary,   only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,   only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_string_functions, only : uppercase
use MOM_tracer_registry, only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen remapping
type, public :: hybgen_remap_CS ; private

  !> Hybgen remapper flag (0=PCM,1=PLM,2=PPM,-ve:isoPCM,3=WENO-like), usually 3
  integer :: hybmap

  !> Hybgen velocity remapper flag (0=PCM,1=PLM,3=WENO-like), usually 3
  integer :: hybmap_vel

end type hybgen_remap_CS

!> Documentation for external callers
character(len=256), public :: remappingSchemesDoc = &
                 "PCM         (1st-order accurate)\n"//&
                 "PLM         (2nd-order accurate)\n"//&
                 "PPM_H4      (3rd-order accurate)\n"//&
                 "WENO        (Weighted Essentially Non-Oscillatory)\n"

! The following are private parameter constants
integer, parameter  :: REMAPPING_PCM        = 0 !< O(h^1) remapping scheme
integer, parameter  :: REMAPPING_PLM        = 1 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PPM_H4     = 2 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_WENO       = 15 !< O(h^5) remapping scheme

public hybgen_remap, init_hybgen_remap

contains

!> Initialise a hybgen_remap_CS control structure and store its parameters
subroutine init_hybgen_remap(CS, GV, US, param_file)
  type(hybgen_remap_CS),   pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file

  character(len=40) :: mdl = "MOM_hybgen" ! This module's name.
  character(len=80) :: scheme_name, vel_scheme_name

  if (associated(CS)) call MOM_error(FATAL, "init_hybgen_remap: CS already associated!")
  allocate(CS)

  call get_param(param_file, mdl, "HYBGEN_REMAPPING_SCHEME", scheme_name, &
                 "An integer indicating the remapping scheme to be used by Hygen. "//&
                 "Valid values are: \n"//trim(remappingSchemesDoc), default="WENO")
  CS%hybmap = hybgen_remap_scheme_from_name(scheme_name, "HYBGEN_REMAPPING_SCHEME")

  call get_param(param_file, mdl, "HYBGEN_REMAPPING_VEL_SCHEME", vel_scheme_name, &
                 "An integer indicating the velocity remapping scheme to be used by Hygen. "//&
                 "Valid values are: \n"//trim(remappingSchemesDoc), default=trim(scheme_name))
  CS%hybmap_vel = hybgen_remap_scheme_from_name(vel_scheme_name, "HYBGEN_REMAPPING_VEL_SCHEME")

  if (CS%hybmap_vel == REMAPPING_PPM_H4) then
    call MOM_error(WARNING, "Replacing PPM_H4 scheme for HYBGEN_REMAPPING_VEL_SCHEME with WENO.")
    CS%hybmap_vel = REMAPPING_WENO
  endif

end subroutine init_hybgen_remap

!> This function returns an integer that indicates which remapping scheme to use
integer function hybgen_remap_scheme_from_name(string, param_name)
  character(len=*),   intent(in)    :: string !< String to parse for the remapping scheme
  character(len=*),   intent(in)    :: param_name !< Name of the parameter being interpreted.

  select case ( uppercase(trim(string)) )
    case ("PCM") ; hybgen_remap_scheme_from_name = REMAPPING_PCM
    case ("PLM") ; hybgen_remap_scheme_from_name = REMAPPING_PLM
    case ("PPM_H4") ; hybgen_remap_scheme_from_name = REMAPPING_PPM_H4
    case ("WENO") ; hybgen_remap_scheme_from_name = REMAPPING_WENO
    case default
      call MOM_error(FATAL, "hybgen_remap_scheme: Unrecognized choice for "//&
                            trim(param_name)//" ("//trim(string)//").")
  end select
end function hybgen_remap_scheme_from_name

!> Hybgen_remap remaps the state variables, velocities and any tracers to a new
!! vertical grid using the algorithms imported from the HYCOM ocean model.
subroutine hybgen_remap(G, GV, remap_CS, dp_orig, dp, Reg, OBC, u, v, PCM_cell)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(hybgen_remap_CS),   intent(in)    :: remap_CS !< hybgen remapping control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp_orig !< Thicknesses on the source grid [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp  !< Thicknesses on the destination grid [H ~> m or kg m-2]
  type(tracer_registry_type), pointer    :: Reg !< Tracer registry structure
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(inout) :: u   !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(inout) :: v   !< Meridional velocity field [L T-1 ~> m s-1]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: PCM_cell !< If true, PCM remapping should be used in a cell.

  ! Local variables
  real :: dpu_orig(SZIB_(G),SZJ_(G),SZK_(GV))  ! Original U-point thicknesses [H ~> m or kg m-2]
  real :: dpv_orig(SZI_(G),SZJB_(G),SZK_(GV))  ! Original V-point thicknesses [H ~> m or kg m-2]
  real :: dpu(SZIB_(G),SZJ_(G),SZK_(GV))       ! Updated U-point layer thicknesses [H ~> m or kg m-2]
  real :: dpv(SZI_(G),SZJB_(G),SZK_(GV))       ! Updated v-point layer thicknesses [H ~> m or kg m-2]
  character(len=256) :: mesg  ! A string for output messages
  integer :: j, ntr

  ntr = 0 ; if (associated(Reg)) ntr = Reg%ntr

  call hybgen_remap_tracers(G, GV, remap_CS, Reg, ntr, dp, dp_orig, PCM_cell)

  if (.not.(present(u) .and. present(v))) return

! --- vertical momentum flux across moving interfaces (the s-dot term in the
! --- momentum equation) - required to locally conserve momentum when hybgen
! --- moves vertical coordinates first, store old interface pressures in
! --- -pu-, -pv-

  ! Find and store the previous thicknesses at velocity points.
  call dpudpv(dpu_orig, dpv_orig, dp_orig, G, GV, OBC)

  ! --- update layer thickness at u- and v- points.
  call dpudpv(dpu, dpv, dp, G, GV, OBC)

  !$OMP parallel do default(shared)
  do j=G%jsc,G%jec
    call hybgenbj_u(remap_CS, G, GV, dpu, dpu_orig, u, j)  !update u-velocity
  enddo

  !$OMP parallel do default(shared)
  do J=G%JscB,G%JecB
    call hybgenbj_v(remap_CS, G, GV, dpv, dpv_orig, v, j)  !update v-velocity
  enddo

end subroutine hybgen_remap

!> Apply vertical remapping for all tracers (incluiding temperature and salinity) to the new grid.
subroutine hybgen_remap_tracers(G, GV, CS, Reg, ntracer, dp_new, dp_orig, PCM_cell)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(hybgen_remap_CS),   intent(in)    :: CS  !< hybgen control structure
  type(tracer_registry_type), pointer    :: Reg !< Tracer registry structure. (Tracer values are intent inout.)
  integer,                 intent(in)    :: ntracer !< The number of tracers in the registry, or
                                                !! 0 if the registry is not in use.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp_new  !< New layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp_orig !< Original layer thicknesses [H ~> m or kg m-2]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: PCM_cell !< If true, PCM remapping should be used in cell.

  ! Local variables
  real :: h_src(GV%ke)      ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  real :: h_tgt(GV%ke)      ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: tracer_i_j(GV%ke,max(ntracer,1))  !  Columns of each tracer [Conc]
  logical :: PCM_lay(GV%ke) ! If true for a layer, use PCM remapping for that layer
  real :: dpthin            ! A negligible layer thickness, used to avoid roundoff issues
                            ! or division by 0 [H ~> m or kg m-2]
  integer :: i, j, k, nk, m

  nk = GV%ke
  dpthin = 1.0e-6*GV%m_to_H

  do k=1,nk ; PCM_lay(k) = .false. ; enddo
  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then
    do k=1,nk
      h_src(k) = dp_orig(i,j,k)
      h_tgt(k) = dp_new(i,j,k)
    enddo
    if (present(PCM_cell)) then ; do k=1,nk
      PCM_lay(k) = PCM_cell(i,j,k)
    enddo ; endif
    ! Note that temperature and salinity are among the tracers remapped here.
    do m=1,ntracer ; do k=1,nk
      tracer_i_j(k,m) = Reg%Tr(m)%t(i,j,k)
    enddo ; enddo

    call hybgen_remap_column(CS%hybmap, nk, h_src, h_tgt, ntracer, tracer_i_j, dpthin, PCM_lay)

    ! Note that temperature and salinity are among the tracers remapped here.
    do m=1,ntracer ; do k=1,nk
      Reg%Tr(m)%t(i,j,k) = tracer_i_j(k,m)
    enddo ; enddo
  endif ; enddo ; enddo !i & j

end subroutine hybgen_remap_tracers

!> Vertically remap a column of scalars to the new grid
subroutine hybgen_remap_column(remap_scheme, nk, h_src_in, h_tgt, ntracr, trac_i_j, dpthin, PCM_lay)
  integer,         intent(in)    :: remap_scheme !< A coded integer indicating the remapping scheme to use
  integer,         intent(in)    :: nk           !< Number of layers in this column
  real,            intent(in)    :: h_src_in(nk)    !< Source grid layer thicknesses [H ~> m or kg m-2]
  real,            intent(in)    :: h_tgt(nk)    !< Target grid layer thicknesses [H ~> m or kg m-2]
  integer,         intent(in)    :: ntracr       !< The number of registered tracers (including temperature
                                                 !! and salinity)
  real,            intent(inout) :: trac_i_j(nk, max(ntracr,1)) !< Columns of the tracers [Conc] or [degC] or [ppt]
  real,            intent(in)    :: dpthin       !< A negligibly small thickness for the purpose of cell
                                                 !! reconstructions [H ~> m or kg m-2].
  logical,         intent(in)    :: PCM_lay(nk)  !< If true for a layer, use PCM remapping for that layer

! --- -------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - remap scalars.
! --- -------------------------------------------------------------
  real :: s1d(nk,ntracr)    ! original scalar fields [Conc] or [degC] or [ppt]
  real :: f1d(nk,ntracr)    ! final    scalar fields [Conc] or [degC] or [ppt]
  real :: c1d(nk,ntracr,3)  ! interpolation coefficients
  real :: h_src(nk)         ! original layer thicknesses [H ~> m or kg m-2]
  real :: dpi( nk)          ! original layer thicknesses, >= dpthin [H ~> m or kg m-2]
  real :: pres(nk+1)        ! Source grid interface depths [H ~> m or kg m-2]
  real :: prsf(nk+1)        ! Target grid interface depths [H ~> m or kg m-2]
  integer :: k, ktr, nums1d
  character(len=256) :: mesg  ! A string for output messages

  double precision, parameter :: zp5=0.5 ! for sign function

  nums1d = ntracr

! --- Store the 'old' vertical grid fields and the 'new' vertical grid spacings
  pres(1) = 0.0 ; prsf(1) = 0.0
!### Remove these two lines later
  do k=1,nk ; pres(K+1) = pres(K) + h_src_in(k) ; enddo
  do k=1,nk ; h_src(k) = pres(K+1) - pres(K) ; enddo

  do k=1,nk
    do ktr=1,ntracr
      s1d(k,ktr) = trac_i_j(k,ktr)
    enddo !ktr

    dpi(k) = max(h_src(k), dpthin)
!### Uncomment this line later   pres(K+1) = pres(K) + h_src(k)
    prsf(K+1) = prsf(K) + h_tgt(k)
  enddo !k

!
! --- remap scalar field profiles from the 'old' vertical
! --- grid onto the 'new' vertical grid.
!
  if     (remap_scheme == REMAPPING_PCM) then !PCM
    call hybgen_pcm_remap(s1d, pres, h_src, f1d, prsf, nk, nk, nums1d, dpthin)
  elseif (remap_scheme == REMAPPING_PLM) then !PLM (as in 2.1.08)
    call hybgen_plm_coefs(s1d, h_src, PCM_lay, c1d, nk, nums1d, dpthin)
    call hybgen_plm_remap(s1d, pres, h_src, c1d, f1d, prsf, nk, nk, nums1d, dpthin)
  elseif (remap_scheme == REMAPPING_PPM_H4) then !PPM
    call hybgen_ppm_coefs(s1d, dpi, PCM_lay, c1d, nk, nums1d, dpthin)
    call hybgen_ppm_remap(s1d, pres, h_src, c1d, f1d, prsf, nk, nk, nums1d, dpthin)
  elseif (remap_scheme == REMAPPING_WENO) then !WENO-like
    call hybgen_weno_coefs(s1d, dpi, PCM_lay, c1d, nk, nums1d, dpthin)
    call hybgen_weno_remap(s1d, pres, h_src, c1d, f1d, prsf, nk, nk, nums1d, dpthin)
  endif
  do k=1,nk
    do ktr=1,ntracr
      trac_i_j(k,ktr) = f1d(k,ktr)
    enddo !ktr
  enddo !k

end subroutine hybgen_remap_column

!> Remap a vertical slice of u-velocities
subroutine hybgenbj_u(CS, G, GV, dpu, dpu_orig, u, j)
  type(hybgen_remap_CS),   intent(in)    :: CS  !< hybgen control structure
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dpu !< New layer thicknesses at u-points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)    :: dpu_orig !< Original layer thicknesses at u-points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u   !< Zonal velocities [L T-1 ~> m s-1]
  integer, intent(in) :: j !< The j-slice to work on
!
! --- --------------------------------------------
! --- hybrid grid generator, single j-row at u points (part B).
! --- --------------------------------------------
!
  real :: s1d(GV%ke)      ! Original velocities [L T-1 ~> m s-1]
  real :: f1d(GV%ke)      ! Final velocities [L T-1 ~> m s-1]
  real :: c1d(GV%ke,3)    ! Interpolation coefficients
  real :: dpi(GV%ke)      ! Original layer thicknesses, >= dpthin [H ~> m or kg m-2]
  real :: dprs(GV%ke)     ! Original layer thicknesses [H ~> m or kg m-2]
  real :: pres(GV%ke+1)   ! Original layer interfaces [H ~> m or kg m-2]
  real :: prsf(GV%ke+1)   ! Final    layer interfaces [H ~> m or kg m-2]
  real :: dpthin ! A negligible layer thickness, used to avoid roundoff issues
                 ! or division by 0 [H ~> m or kg m-2]
  real :: onemm  ! one mm in pressure units [H ~> m or kg m-2]
  logical :: lcm(GV%ke)      ! Use PCM for some layers? (always .false.)
  character(len=256) :: mesg  ! A string for output messages
  integer :: i, k, kk

! --- vertical momentum flux across moving interfaces (the s-dot term in the
! --- momentum equation) - required to locally conserve momentum when hybgen
! --- moves vertical coordinates.

  kk = GV%ke
  onemm = 0.001*GV%m_to_H
  dpthin = 1.0e-6*GV%m_to_H
!
! --- always use high order remapping for velocity
  do k=1,kk
    lcm(k) = .false.  !use same remapper for all layers
  enddo !k
!
  do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
!
! ---     store one-dimensional arrays of -u- and -p- for the 'old' vertical grid
    pres(1) = 0.0
    prsf(1) = 0.0
    do k=1,kk
      s1d(k)   = u(I,j,k)
      pres(K+1) = pres(K) + dpu_orig(I,j,k) ! original vertical grid
      prsf(K+1) = prsf(K) + dpu(I,j,k)      ! new vertical grid
      dprs(k)  = dpu_orig(I,j,k)
      dpi( k)  = max(dprs(k), dpthin)
    enddo !k
!
! ---     remap -u- profiles from the 'old' vertical grid onto the
! ---     'new' vertical grid.
!
    if     (CS%hybmap_vel == REMAPPING_PCM) then !PCM
      call hybgen_pcm_remap(s1d, pres, dprs, f1d, prsf, kk, kk, 1, dpthin)
    elseif (CS%hybmap_vel == REMAPPING_PLM) then !PLM (as in 2.1.08)
      call hybgen_plm_coefs(s1d, dprs, lcm, c1d, kk, 1, dpthin)
      call hybgen_plm_remap(s1d, pres, dprs, c1d, f1d, prsf, kk, kk, 1, dpthin)
    else !WENO-like (even if scalar fields are PLM or PPM)
      call hybgen_weno_coefs(s1d, dpi, lcm, c1d, kk, 1, dpthin)
      call hybgen_weno_remap(s1d, pres, dprs, c1d, f1d, prsf, kk, kk, 1, dpthin)
    endif !hybmap_vel
    do k=1,kk
      if ((dpi(k) > dpthin) .or. (prsf(k) <= prsf(kk+1)-onemm)) then
        u(I,j,k) = f1d(k)
      else
      ! --- thin near-bottom layer, zero total current
        u(I,j,k) = 0.0
      endif
!###
!     u(I,j,k) = f1d(k)
    enddo !k

  endif ; enddo !iu
!
end subroutine hybgenbj_u

!> Remap a vertical slice of v-velocities
subroutine hybgenbj_v(CS, G, GV, dpv, dpv_orig, v, j)
  type(hybgen_remap_CS),   intent(in)    :: CS  !< hybgen control structure
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: dpv !< New layer thicknesses at v-points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(in)    :: dpv_orig !< Original layer thicknesses at v-points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v   !< Meridional velocities [L T-1 ~> m s-1]
  integer, intent(in) :: J !< The j-slice to work on
!
! --- --------------------------------------------
! --- hybrid grid generator, single j-row at v points (part B).
! --- --------------------------------------------
!
  real :: s1d(GV%ke)      ! Original velocities [L T-1 ~> m s-1]
  real :: f1d(GV%ke)      ! Final velocities [L T-1 ~> m s-1]
  real :: c1d(GV%ke,3)    ! Interpolation coefficients
  real :: dpi(GV%ke)      ! Original layer thicknesses, >= dpthin [H ~> m or kg m-2]
  real :: dprs(GV%ke)     ! Original layer thicknesses [H ~> m or kg m-2]
  real :: pres(GV%ke+1)   ! Original layer interfaces [H ~> m or kg m-2]
  real :: prsf(GV%ke+1)   ! Final    layer interfaces [H ~> m or kg m-2]
  real :: dpthin ! A negligible layer thickness, used to avoid roundoff issues
                 ! or division by 0 [H ~> m or kg m-2]
  real :: onemm  ! one mm in pressure units [H ~> m or kg m-2]
  logical :: lcm(GV%ke)   ! Use PCM for some layers? (always .false.)
  character(len=256) :: mesg  ! A string for output messages
  integer :: i, k, kk

!
! --- vertical momentum flux across moving interfaces (the s-dot term in the
! --- momentum equation) - required to locally conserve momentum when hybgen
! --- moves vertical coordinates.
!
  kk = GV%ke
  onemm = 0.001*GV%m_to_H
  dpthin = 1.0e-6*GV%m_to_H
!
! --- always use high order remapping for velocity
  do k=1,kk
    lcm(k) = .false.  !use same remapper for all layers
  enddo !k

  do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
!
! ---     store one-dimensional arrays of -v- and -p- for the 'old' vertical grid
    pres(1) = 0.0
    prsf(1) = 0.0
    do k=1,kk
      s1d(k)    = v(i,J,k)
      pres(K+1) = pres(K) + dpv_orig(i,J,k) ! original vertical grid
      prsf(K+1) = prsf(K) + dpv(i,J,k)      ! new vertical grid
      dprs(k)   = dpv_orig(i,J,k)
      dpi( k)   = max(dprs(k), dpthin)
    enddo !k
!
! ---     remap -v- profiles from the 'old' vertical grid onto the
! ---     'new' vertical grid.
!
    if     (CS%hybmap_vel == REMAPPING_PCM) then !PCM
      call hybgen_pcm_remap(s1d, pres, dprs, f1d, prsf, kk, kk, 1, dpthin)
    elseif (CS%hybmap_vel == REMAPPING_PLM) then !PLM (as in 2.1.08)
      call hybgen_plm_coefs(s1d, dprs, lcm, c1d, kk, 1, dpthin)
      call hybgen_plm_remap(s1d, pres, dprs, c1d, f1d, prsf, kk, kk, 1, dpthin)
    else !WENO-like (even if scalar fields are PLM or PPM)
      call hybgen_weno_coefs(s1d, dpi, lcm, c1d, kk, 1, dpthin)
      call hybgen_weno_remap(s1d, pres, dprs, c1d, f1d, prsf, kk, kk, 1, dpthin)
    endif !hybmap_vel
    do k=1,kk
      if ((dpi(k) > dpthin) .or. (prsf(k) <= prsf(kk+1)-onemm)) then
        v(i,J,k) = f1d(k)
      else
        ! --- thin near-bottom layer, zero total current
        v(i,J,k) = 0.0
      endif
!###
!     v(i,J,k) = f1d(k)
    enddo !k

  endif ; enddo !iv
!
end subroutine hybgenbj_v

!> Do piecewise constant remapping for a set of scalars
subroutine hybgen_pcm_remap(si, pi, dpi, so, po, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
!  real,    intent(in) :: dpo(ko+1) !< The output layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: piecewise constant across each input cell
!             the output is the average of the interpolation
!             profile across each output cell.
!
!     PCM (donor cell) is the standard 1st order upwind method.
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       pi    - initial layer interface depths (non-negative)
!                  pi(   1) is the surface
!                  pi(ki+1) is the bathymetry
!                  pi(k+1) >= pi(k)
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       ki    - number of  input layers
!       ko    - number of output layers
!       ks    - number of fields
!       po    - target interface depths (non-negative)
!                  po(   1) is the surface
!                  po(ko+1) is the bathymetry (== pi(ki+1))
!                  po(k+1) >= po(k)
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       so    - scalar fields in po-layer space
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  integer :: i, k, l, lb, lt
  real :: dpb, dpt, xb, xt, zb, zt, zx, o
  real*8 :: sz
  real :: si_min(ks), si_max(ks)
!
! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
!
  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i
!
  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while (pi(lb+1) < zb .and. lb < ki)
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
!       write(mesg,*) 'k,zt,zb = ',k,zt,zb
!       call MOM_mesg(mesg, all_print=.true.)
    lt=lb !top will always correspond to bottom of previous
          !find input layer containing bottom output interface
    do while (pi(lb+1) < zb .and. lb < ki)
      lb = lb+1
    enddo
    if (zb-zt <= thin .or. zt >= zx) then
      if (k /= 1) then
        ! --- thin or bottomed layer, values taken from layer above
        do i=1,ks
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ks
          so(k,i) = si(k,i)
        enddo !i
      endif
    else
!
!     form layer averages.
!     use PPM-like logic (may not have minimum operation count)
!
!     if     (pi(lb) > zt) then
!       write(mesg,*) 'bad lb = ',lb
!       call MOM_error(FATAL, mesg)
!     endif
      if     (lt /= lb) then  !multiple layers
        xt=(zt-pi(lt))/max(dpi(lt),thin)
        xb = (zb-pi(lb))/max(dpi(lb),thin)
        dpt=pi(lt+1)-zt
        dpb = zb-pi(lb)
        do i=1,ks
          o  = si((lt+lb)/2,i)  !offset to reduce round-off
          sz = dpt*(si(lt,i)-o)
          do l=lt+1,lb-1
            sz = sz+dpi(l)*(si(l,i)-o)
          enddo !l
          sz = sz+dpb*(si(lb,i)-o)
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      else  !single layer
        do i=1,ks
          so(k,i) = si(lt,i)
        enddo !i
      endif
    endif !thin:std layer
  enddo !k

end subroutine hybgen_pcm_remap

!> Set up the coefficients for PLM remapping of a set of scalars
subroutine hybgen_plm_coefs(si, dpi, lc, ci, kk, ks, thin)
  integer, intent(in)  :: kk        !< The number of input layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(kk,ks) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(kk)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  logical, intent(in)  :: lc(kk)    !< If true for a layer, use PCM remapping for that layer
  real,    intent(out) :: ci(kk,ks) !< The PLM coefficients (slopes) describing the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: piecewise linear across each input cell with
!             monotonized central-difference limiter.
!
!     van Leer, B., 1977, J. Comp. Phys., 23 276-299.
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       lc    - use PCM for selected layers
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       ci    - coefficients (slopes) for hybgen_plm_remap
!                profile(y) = si+ci*(y-1),  0<=y<=1
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  integer k,i
  real    qcen,zbot,zcen,ztop
!
  do i=1,ks
    ci(1, i) = 0.0
    ci(kk,i) = 0.0
  enddo !i
  do k= 2,kk-1
    if     (lc(k) .or. dpi(k) <= thin) then  !use PCM
      do i=1,ks
        ci(k,i) = 0.0
      enddo !i
    else
! ---     use qcen in place of 0.5 to allow for non-uniform grid
      qcen = dpi(k)/(dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))  !dpi(k)>thin
      do i=1,ks
! ---       PLM (non-zero slope, but no new extrema)
! ---       layer value is si-0.5*ci at top    interface,
! ---                  and si+0.5*ci at bottom interface.
!
! ---       monotonized central-difference limiter (van Leer, 1977,
! ---       JCP 23 pp 276-299).  For a discussion of PLM limiters, see
! ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
        ztop = 2.0*(si(k,  i)-si(k-1,i))
        zbot = 2.0*(si(k+1,i)-si(k,  i))
        zcen = qcen*(si(k+1,i)-si(k-1,i))
        if     (ztop*zbot > 0.0) then !ztop,zbot are the same sign
          ci(k,i) = sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
        else
          ci(k,i) = 0.0  !local extrema, so no slope
        endif
      enddo !i
    endif  !PCM:PLM
  enddo !k

end subroutine hybgen_plm_coefs

!> Do piecewise linear remapping for a set of scalars
subroutine hybgen_plm_remap(si, pi, dpi, ci, so, po, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: ci(ki,ks) !< The PLM coefficients (slopes) describing the scalar fields [A]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
!  real,    intent(in) :: dpo(ko+1) !< The output layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: piecewise linear across each input cell
!             the output is the average of the interpolation
!             profile across each output cell.
!
!     van Leer, B., 1977, J. Comp. Phys., 23 276-299.
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       pi    - initial layer interface depths (non-negative)
!                  pi(   1) is the surface
!                  pi(ki+1) is the bathymetry
!                  pi(k+1) >= pi(k)
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       ci    - coefficients (slopes) from hybgen_plm_coefs
!                profile(y) = si+ci*(y-1),  0<=y<=1
!       ki    - number of  input layers
!       ko    - number of output layers
!       ks    - number of fields
!       po    - target interface depths (non-negative)
!                  po(   1) is the surface
!                  po(ko+1) is the bathymetry (== pi(ki+1))
!                  po(k+1) >= po(k)
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       so    - scalar fields in po-layer space
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  integer :: i, k, l, lb, lt
  real :: c0, qb0, qb1, qt0, qt1, xb, xt, zb, zt, zx, o
  real*8 :: sz
  real :: si_min(ks),si_max(ks)
!
! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
!
  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i
!
  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while ((pi(lb+1) < zb) .and. (lb < ki))
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
!       write(mesg,*) 'k,zt,zb = ',k,zt,zb
!       call MOM_mesg(mesg, all_print=.true.)
    lt = lb !top will always correspond to bottom of previous
            !find input layer containing bottom output interface
    do while ((pi(lb+1) < zb) .and. (lb < ki))
      lb = lb+1
    enddo
    if ((zb-zt <= thin) .or. (zt >= zx)) then
      if (k /= 1) then
!
! ---       thin or bottomed layer, values taken from layer above
!
        do i=1,ks
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ks
          so(k,i) = si(k,i)
        enddo !i
      endif
    else
!
!         form layer averages.
!         use PPM-like logic (may not have minimum operation count)
!
!         if     (pi(lb) > zt) then
!           write(mesg,*) 'bad lb = ',lb
!           call MOM_error(FATAL, mesg)
!         endif
      xt = (zt-pi(lt)) / max(dpi(lt),thin)
      xb = (zb-pi(lb)) / max(dpi(lb),thin)
      if (lt /= lb) then  !multiple layers
        qt0 = (1.0-xt)
        qt1 = (1.0-xt**2)*0.5
        qb0 = xb
        qb1 = xb**2 *0.5
        do i=1,ks
          o  = si((lt+lb)/2,i)  !offset to reduce round-off
          c0 = si(lt,i) - o - 0.5*ci(lt,i)
          sz=  dpi(lt)*(c0*qt0 + ci(lt,i)*qt1)
          do l=lt+1,lb-1
            sz = sz+dpi(l)*(si(l,i) - o)
          enddo !l
          c0 = si(lb,i) - o - 0.5*ci(lb,i)
          sz = sz+dpi(lb)*(c0*qb0 + ci(lb,i)*qb1)
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      else  !single layer
        qt1 = (xb**2-xt**2 - (xb-xt))*0.5
        do i=1,ks
          sz = dpi(lt)*(ci(lt,i)*qt1)
          so(k,i) = si(lt,i) + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      endif
    endif !thin:std layer
  enddo !k

end subroutine hybgen_plm_remap

!> Set up the coefficients for PPM remapping of a set of scalars
subroutine hybgen_ppm_coefs(s, dp, lc, ci, kk, ks, thin)
  integer, intent(in)  :: kk        !< The number of input layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: s(kk,ks)  !< The input scalar fields [A]
  real,    intent(in)  :: dp(kk)    !< The input grid layer thicknesses [H ~> m or kg m-2]
  logical, intent(in)  :: lc(kk)    !< If true for a layer, use PCM remapping for that layer
  real,    intent(out) :: ci(kk,ks,3) !< The PPM coefficients describing the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic piecewise parabolic across each input cell
!
!     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       dp    - initial layer thicknesses (>=thin)
!       lc    - use PCM for selected layers
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       ci    - coefficients for hybgen_ppm_remap
!                profile(y) = ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  integer j, i
  real :: da, a6, slj, scj, srj
  real :: as(kk), al(kk), ar(kk)
  real :: dpjp(kk), dp2jp(kk), dpj2p(kk)
  real :: qdpjp(kk), qdp2jp(kk), qdpj2p(kk), dpq3(kk), qdp4(kk)
!
  !compute grid metrics
  do j=1,kk-1
    dpjp( j) = dp(j)   + dp(j+1)
    dp2jp(j) = dp(j)   + dpjp(j)
    dpj2p(j) = dpjp(j) + dp(j+1)
    qdpjp( j) = 1.0/dpjp( j)
    qdp2jp(j) = 1.0/dp2jp(j)
    qdpj2p(j) = 1.0/dpj2p(j)
  enddo !j
    dpq3(2) = dp(2)/(dp(1)+dpjp(2))
  do j=3,kk-1
    dpq3(j) = dp(j)/(dp(j-1)+dpjp(j)) !dp(j)/      (dp(j-1)+dp(j)+dp(j+1))
    qdp4(j) = 1.0/(dpjp(j-2)+dpjp(j)) !1.0/(dp(j-2)+dp(j-1)+dp(j)+dp(j+1))
  enddo !j
!
  do i=1,ks
    !Compute average slopes: Colella, Eq. (1.8)
    as(1) = 0.
    do j=2,kk-1
      if     (lc(j) .or. dp(j) <= thin) then  !use PCM
        as(j) = 0.0
      else
        slj = s(j,  i)-s(j-1,i)
        srj = s(j+1,i)-s(j,  i)
        if (slj*srj > 0.) then
          scj = dpq3(j)*( dp2jp(j-1)*srj*qdpjp(j) &
                       +dpj2p(j)  *slj*qdpjp(j-1) )
          as(j) = sign(min(abs(2.0*slj),abs(scj),abs(2.0*srj)),scj)
        else
          as(j) = 0.
        endif
      endif  !PCM:PPM
    enddo !j
    as(kk) = 0.
    !Compute "first guess" edge values: Colella, Eq. (1.6)
    al(1) = s(1,i)  !1st layer PCM
    ar(1) = s(1,i)  !1st layer PCM
    al(2) = s(1,i)  !1st layer PCM
    do j=3,kk-1
      al(j) = s(j-1,i)+dp(j-1)*(s(j,i)-s(j-1,i))*qdpjp(j-1) &
           +qdp4(j)*( &
              2.*dp(j)*dp(j-1)*qdpjp(j-1)*(s(j,i)-s(j-1,i))* &
              ( dpjp(j-2)*qdp2jp(j-1) &
               -dpjp(j)  *qdpj2p(j-1) ) &
              -dp(j-1)*as(j)  *dpjp(j-2)*qdp2jp(j-1) &
              +dp(j)  *as(j-1)*dpjp(j)  *qdpj2p(j-1) &
                )
      ar(j-1) = al(j)
    enddo !j
    ar(kk-1) = s(kk,i)  !last layer PCM
    al(kk)  = s(kk,i)  !last layer PCM
    ar(kk)  = s(kk,i)  !last layer PCM
    !Impose monotonicity: Colella, Eq. (1.10)
    do j=2,kk-1
      if (lc(j) .or. dp(j) <= thin) then  !use PCM
        al(j) = s(j,i)
        ar(j) = s(j,i)
      elseif ((s(j+1,i)-s(j,i))*(s(j,i)-s(j-1,i)) <= 0.) then !local extremum
        al(j) = s(j,i)
        ar(j) = s(j,i)
      else
        da=ar(j)-al(j)
        a6=6.0*s(j,i)-3.0*(al(j)+ar(j))
        if (da*a6 > da*da) then !peak in right half of zone
          al(j) = 3.0*s(j,i)-2.0*ar(j)
        elseif (da*a6  <  -da*da) then !peak in left half of zone
          ar(j) = 3.0*s(j,i)-2.0*al(j)
        endif
      endif
    enddo !j
    !Set coefficients
    do j=1,kk
      if (al(j) /= ar(j)) then
        ci(j,i,1) = al(j)
        ci(j,i,2) = ar(j)-al(j)
        ci(j,i,3) = 6.0*s(j,i)-3.0*(al(j)+ar(j))
      else !PCM
        ci(j,i,1) = al(j)
        ci(j,i,2) = 0.0
        ci(j,i,3) = 0.0
      endif
    enddo !j
  enddo !i

end subroutine hybgen_ppm_coefs

!> Do piecewise parabolic remapping for a set of scalars
subroutine hybgen_ppm_remap(si, pi, dpi, ci, so, po, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: ci(ki,ks,3) !< The PPM coefficients describing the scalar fields [A]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
!  real,    intent(in) :: dpo(ko+1) !< The output layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: monotonic piecewise parabolic across each input cell
!             the output is the average of the interpolation
!             profile across each output cell.
!     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       pi    - initial layer interface depths (non-negative)
!                  pi(   1) is the surface
!                  pi(ki+1) is the bathymetry
!                  pi(k+1) >= pi(k)
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       ci    - coefficients from hybgen_ppm_coefs
!                profile(y) = ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
!       ki    - number of  input layers
!       ko    - number of output layers
!       ks    - number of fields
!       po    - target interface depths (non-negative)
!                  po(   1) is the surface
!                  po(ko+1) is the bathymetry (== pi(ki+1))
!                  po(k+1) >= po(k)
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       so    - scalar fields in po-layer space
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  integer i,k,l,lb,lt
  real    qb0, qb1, qb2, qt0, qt1, qt2,xb,xt,zb,zt,zx,o
  real*8  sz
  real    si_min(ks),si_max(ks)
!
! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
!
  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i
!
  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while (pi(lb+1) < zb .and. lb < ki)
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
!       write(mesg,*) 'k,zt,zb = ',k,zt,zb
!       call MOM_mesg(mesg, all_print=.true.)
    lt=lb !top will always correspond to bottom of previous
          !find input layer containing bottom output interface
    do while (pi(lb+1) < zb .and. lb < ki)
      lb = lb+1
    enddo
    if     (zb-zt <= thin .or. zt >= zx) then
      if     (k /= 1) then
!
! ---       thin or bottomed layer, values taken from layer above
!
        do i=1,ks
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ks
          so(k,i) = si(k,i)
        enddo !i
      endif
    else
      ! Form layer averages.
!
!         if     (pi(lb) > zt) then
!           write(mesg,*) 'bad lb = ',lb
!           call MOM_error(FATAL, mesg)
!         endif
      xt=(zt-pi(lt))/max(dpi(lt),thin)
      xb = (zb-pi(lb))/max(dpi(lb),thin)
      if     (lt /= lb) then  !multiple layers
        qt0 = (1.0-xt)
        qt1 = (1.-xt**2)*0.5
        qt2 = (1.-xt**3)/3.0
        qb0 =     xb
        qb1 =     xb**2 *0.5
        qb2 =     xb**3 /3.0
        do i=1,ks
          o  = si((lt+lb)/2,i)  !offset to reduce round-off
          sz = dpi(lt)*(  (ci(lt,i,1)-o)*qt0 &
                         +(ci(lt,i,2) + ci(lt,i,3) ) *qt1 &
                          -ci(lt,i,3)   *qt2 )
          do l=lt+1,lb-1
            sz = sz+dpi(l)*(si(l,i)-o)
          enddo !l
          sz = sz+dpi(lb)*( (ci(lb,i,1)-o)*qb0 &
                           +(ci(lb,i,2) + ci(lb,i,3) ) *qb1 &
                            -ci(lb,i,3)   *qb2 )
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      else  !single layer
        qt0 = (xb-xt)
        qt1 = (xb**2-xt**2)*0.5
        qt2 = (xb**3-xt**3)/3.0
        do i=1,ks
          o  = si( lt,i)  !offset to reduce round-off
          sz = dpi(lt)*( (ci(lt,i,1)-o)*qt0 &
                        +(ci(lt,i,2) + ci(lt,i,3) ) *qt1 &
                         -ci(lt,i,3)   *qt2 )
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      endif
    endif !thin:std layer
  enddo !k

end subroutine hybgen_ppm_remap

!> Set up the coefficients for PPM remapping of a set of scalars
subroutine hybgen_weno_coefs(s, dp,lc,ci,kk,ks,thin)
  integer, intent(in)  :: kk        !< The number of input layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: s(kk,ks)  !< The input scalar fields [A]
  real,    intent(in)  :: dp(kk)    !< The input grid layer thicknesses [H ~> m or kg m-2]
  logical, intent(in)  :: lc(kk)    !< If true for a layer, use PCM remapping for that layer
  real,    intent(out) :: ci(kk,ks,2) !< The WENO coefficients describing the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the
!             interfacial values
!
!     REFERENCE?
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       dp    - initial layer thicknesses (>=thin)
!       lc    - use PCM for selected layers
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       ci    - coefficients for hybgen_weno_remap
!                ci.1 is value at interface above
!                ci.2 is value at interface below
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
!-----------------------------------------------------------------------
!
  real, parameter :: dsmll=1.0e-8
!
  integer j,i
  real    q, q01, q02, q001, q002
  real    qdpjm(kk), qdpjmjp(kk), dpjm2jp(kk)
  real    zw(kk+1,3)

  !compute grid metrics
  do j=2,kk-1
    qdpjm(  j) = 1.0 / (dp(j-1) + dp(j))
    qdpjmjp(j) = 1.0 / (dp(j-1) + dp(j) + dp(j+1))
    dpjm2jp(j) = dp(j-1) + 2.0*dp(j) + dp(j+1)
  enddo !j
  j = kk
  qdpjm(  j) = 1.0 / (dp(j-1) + dp(j))
!
  do i=1,ks
    do j=2,kk
      zw(j,3) = qdpjm(j) * (s(j,i)-s(j-1,i))
    enddo !j
      j = 1  !PCM first layer
        ci(j,i,1) = s(j,i)
        ci(j,i,2) = s(j,i)
        zw(j,  1) = 0.0
        zw(j,  2) = 0.0
    do j=2,kk-1
      if (lc(j) .or. dp(j) <= thin) then  !use PCM
        ci(j,i,1) = s(j,i)
        ci(j,i,2) = s(j,i)
        zw(j,  1) = 0.0
        zw(j,  2) = 0.0
      else
        q001 = dp(j)*zw(j+1,3)
        q002 = dp(j)*zw(j,  3)
        if (q001*q002 < 0.0) then
          q001 = 0.0
          q002 = 0.0
        endif
        q01 = dpjm2jp(j)*zw(j+1,3)
        q02 = dpjm2jp(j)*zw(j,  3)
        if (abs(q001) > abs(q02)) then
          q001 = q02
        endif
        if (abs(q002) > abs(q01)) then
          q002 = q01
        endif
        q    = (q001-q002)*qdpjmjp(j)
        q001 = q001-q*dp(j+1)
        q002 = q002+q*dp(j-1)

        ci(j,i,2) = s(j,i)+q001
        ci(j,i,1) = s(j,i)-q002
        zw(  j,1) = (2.0*q001-q002)**2
        zw(  j,2) = (2.0*q002-q001)**2
      endif  !PCM:WEND
    enddo !j
      j = kk  !PCM last layer
        ci(j,i,1) = s(j,i)
        ci(j,i,2) = s(j,i)
        zw(j,  1) = 0.0
        zw(j,  2) = 0.0

    do j=2,kk
      q002 = max(zw(j-1,2), dsmll)
      q001 = max(zw(j,  1), dsmll)
      zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
    enddo !j
      zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
      zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

    do j=2,kk-1
      if     (.not.(lc(j) .or. dp(j) <= thin)) then  !don't use PCM
        q01  = zw(j+1,3)-s(j,i)
        q02  = s(j,i)-zw(j,3)
        q001 = 2.0*q01
        q002 = 2.0*q02
        if     (q01*q02 < 0.0) then
          q01 = 0.0
          q02 = 0.0
        elseif (abs(q01) > abs(q002)) then
          q01 = q002
        elseif (abs(q02) > abs(q001)) then
          q02 = q001
        endif
        ci(j,i,1) = s(j,i)-q02
        ci(j,i,2) = s(j,i)+q01
      endif  !PCM:WEND
    enddo !j
  enddo !i

end subroutine hybgen_weno_coefs

!> Do WENO remapping for a set of scalars
subroutine hybgen_weno_remap(si, pi, dpi, ci, so, po, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: ci(ki,ks,2) !< The WENO coefficients describing the scalar fields [A]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
!  real,    intent(in) :: dpo(ko+1) !< The output layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the
!             interfacial values
!             the output is the average of the interpolation
!             profile across each output cell.
!
!     REFERENCE?
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       pi    - initial layer interface depths (non-negative)
!                  pi(   1) is the surface
!                  pi(ki+1) is the bathymetry
!                  pi(k+1) >= pi(k)
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       ci    - coefficients from hybgen_weno_coefs
!                ci.1 is value at interface above
!                ci.2 is value at interface below
!       ki    - number of  input layers
!       ko    - number of output layers
!       ks    - number of fields
!       po    - target interface depths (non-negative)
!                  po(   1) is the surface
!                  po(ko+1) is the bathymetry (== pi(ki+1))
!                  po(k+1) >= po(k)
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       so    - scalar fields in po-layer space
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  integer i,k,l,lb,lt
  real    dpb, dpt, qb0, qb1, qb2, qt0, qt1, qt2,xb,xt,zb,zt,zx,o
  real*8  sz
  real    si_min(ks),si_max(ks)
!
! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
!
  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i

  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while (pi(lb+1) < zb .and. lb < ki)
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
!       write(mesg,*) 'k,zt,zb = ',k,zt,zb
!       call MOM_mesg(mesg, all_print=.true.)
    lt = lb ! top will always correspond to bottom of previous
            ! find input layer containing bottom output interface
    do while (pi(lb+1) < zb .and. lb < ki)
      lb = lb+1
    enddo
    if (zb-zt <= thin .or. zt >= zx) then
      if (k /= 1) then
!
! ---       thin or bottomed layer, values taken from layer above
!
        do i=1,ks
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ks
          so(k,i) = si(k,i)
        enddo !i
      endif
    else
!
!         form layer averages.
!
!         if     (pi(lb) > zt) then
!           write(mesg,*) 'bad lb = ',lb
!           call MOM_error(FATAL, mesg)
!         endif
      xt = (zt-pi(lt)) / max(dpi(lt),thin)
      xb = (zb-pi(lb)) / max(dpi(lb),thin)
      if (lt /= lb) then  !multiple layers
        dpt = pi(lt+1)-zt
        dpb = zb-pi(lb)
        qt1 = xt*(xt-1.0)
        qt2 = qt1+xt
        qt0 = 1.0-qt1-qt2
        qb1 = (xb-1.0)**2
        qb2 = qb1-1.0+xb
        qb0 = 1.0-qb1-qb2
        do i=1,ks
          o = si((lt+lb)/2,i)  !offset to reduce round-off
          sz = dpt*(qt0*(si(lt,i)  -o) + &
                    qt1*(ci(lt,i,1)-o) + &
                    qt2*(ci(lt,i,2)-o)  )
          do l=lt+1,lb-1
            sz = sz+dpi(l)*(si(l,i) - o)
          enddo !l
          sz  = sz + dpb*(qb0*(si(lb,i)  -o) + &
                          qb1*(ci(lb,i,1)-o) + &
                          qb2*(ci(lb,i,2)-o)  )
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      else !single layer
        qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
        qt2 = qt1 - 1.0 + (xb+xt)
        do i=1,ks
          o =     si(lt,i)  !offset to reduce round-off
          sz = qt1*(ci(lt,i,1)-o) + qt2*(ci(lt,i,2)-o)
          so(k,i) = o + sz
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      endif !layers
    endif !thin:std layer
  enddo !k

end subroutine hybgen_weno_remap

!> Specify layer thicknesses at velocity points for remapping
subroutine dpudpv(dpu, dpv, dp, G, GV, OBC)
  type(ocean_grid_type),                        intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)  :: GV !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(out) :: dpu !< Thicknesses at u points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(out) :: dpv !< Thicknesses at v points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)  :: dp  !< Interface thicknesses [H ~> m or kg m-2]
  type(ocean_OBC_type),                         pointer     :: OBC !< Open boundary structure

  ! Local variables
  real :: p(SZI_(G),SZJ_(G),SZK_(GV)+1) ! Interface positions relative to the surface [H ~> m or kg m-2]
  real :: dp_tot_uv                     ! The minimum total thickness at a velocity point [H ~> m or kg m-2]
  integer :: i, j, k, nk

  nk = GV%ke

  ! Calculate the interface positions from the layer thicknesses
  ! --- p and dp have a valid 1-point halo
  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    p(i,j,1) = 0.0
  enddo ; enddo
  do k=1,nk ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    p(i,j,K+1) = p(i,j,K) + dp(i,j,k)
  enddo ; enddo ; enddo

  dpu(:,:,:) = 0.0

  do j=G%jsc,G%jec
    do k=1,nk
      do I=G%IscB,G%IecB
        if (G%mask2dCu(I,j) > 0.) then
          dp_tot_uv = min(p(i,j,nk+1), p(i+1,j,nk+1))
          dpu(I,j,k) = max(0., &
               min(dp_tot_uv, .5*(p(i,j,k+1)+p(i+1,j,k+1)))- &
               min(dp_tot_uv, .5*(p(i,j,k  )+p(i+1,j,k  ))))
          ! This variant does not restrict thicknesses using topography.
          ! dpu(I,j,k) = 0.5*(p(i,j,k+1)+p(i+1,j,k+1)) - &
          !              0.5*(p(i,j,k  )+p(i+1,j,k  ))
          if (associated(OBC)) then ; if (OBC%segnum_u(I,j) /= 0) then
            if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
              dpu(I,j,k) = dp(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              dpu(I,j,k) = dp(i+1,j,k)
            endif
          endif ; endif
        endif !iu
      enddo
    enddo
  enddo !j

  dpv(:,:,:) = 0.0
  do J=G%JscB,G%JecB
    do k=1,nk
      do i=G%isc,G%iec
        if (G%mask2dCv(i,J) > 0.) then
          dp_tot_uv = min(p(i,j,nk+1), p(i,j+1,nk+1))
          dpv(i,J,k) = max(0., &
               min(dp_tot_uv, .5*(p(i,j,k+1)+p(i,j+1,k+1)))- &
               min(dp_tot_uv, .5*(p(i,j,k  )+p(i,j+1,k  ))))
          ! This variant does not restrict thicknesses using topography.
          ! dpv(i,J,k) = 0.5*(p(i,j,k+1)+p(i,j+1,k+1)) - &
          !              0.5*(p(i,j,k  )+p(i,j+1,k  ))
          if (associated(OBC)) then ; if (OBC%segnum_v(i,J) /= 0) then
            if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
              dpu(I,j,k) = dp(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              dpu(I,j,k) = dp(i,j+1,k)
            endif
          endif ; endif
        endif !iv
      enddo !i
    enddo
  enddo !j

end subroutine dpudpv

end module MOM_hybgen_remap
