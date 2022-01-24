!> This module contains the hybgen regridding routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_regrid

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,             only : EOS_type, calculate_density
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING, assert
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_regrid_CS ; private

  !> Number of layers on the target grid
  integer :: nk

!  !> Minimum thickness allowed for layers, often in [H ~> m or kg m-2]
!  real :: min_thickness = 0.

  !> Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure

  !> Hybgen uses PCM if layer is within hybiso of target density [kg m-3]
  real :: hybiso
  !> Number of hybrid levels used by HYBGEN (0=all isopycnal)
  integer :: nhybrid
  !> Number of sigma levels used by HYBGEN (nhybrid-nsigma z-levels)
  integer :: nsigma
  !> Deep isopycnal spacing minimum thickness (m)
  real :: dp00i
  !> Hybgen relaxation coefficient (inverse baroclinic time steps) [s-1]
  real :: qhybrlx
  !> If true, Hybgen uses PCM to remap isopycnal layers
  logical :: isopcm

  !> Reference density for anomalies [R ~> kg m-3]
  real :: thbase

  real, allocatable, dimension(:) ::  &
    dp0k, & !< minimum deep    z-layer separation [H ~> m or kg m-2]
    ds0k    !< minimum shallow z-layer separation [H ~> m or kg m-2]

  real :: dpns  !< depth to start terrain following [H ~> m or kg m-2]
  real :: dsns  !< depth to stop terrain following [H ~> m or kg m-2]

  real :: thkbot !< Thickness of a bottom boundary layer, within which hybgen does
                 !! something different. [H ~> m or kg m-2]

  !> Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real :: topiso_const
  ! real, dimension(:,:), allocatable :: topiso

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  real :: onem       !< Nominally one m in thickness units [H ~> m or kg m-2]

  logical :: debug   !< If true, write verbose checksums for debugging purposes.

end type hybgen_regrid_CS


public hybgen_regrid, init_hybgen_regrid, end_hybgen_regrid
public hybgen_column_init, set_hybgen_regrid_params

contains

!> Initialise a hybgen_regrid_CS control structure and store its parameters
subroutine init_hybgen_regrid(CS, GV, US, param_file)
  type(hybgen_regrid_CS),  pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file

  character(len=40)               :: mdl = "MOM_hybgen" ! This module's name.
  integer :: k

  if (associated(CS)) call MOM_error(FATAL, "init_hybgen_regrid: CS already associated!")
  allocate(CS)

  CS%nk = GV%ke

  allocate(CS%target_density(CS%nk))
  allocate(CS%dp0k(CS%nk), source=0.0) ! minimum deep z-layer separation
  allocate(CS%ds0k(CS%nk), source=0.0) ! minimum shallow z-layer separation

  call get_param(param_file, mdl, "P_REF", CS%ref_pressure, &
                 "The pressure that is used for calculating the coordinate "//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) "//&
                 "This is only used if USE_EOS and ENABLE_THERMODYNAMICS are true.", &
                 units="Pa", default=2.0e7, scale=US%kg_m3_to_R*US%m_s_to_L_T**2)

  call get_param(param_file, mdl, "HYBGEN_N_HYBRID", CS%nhybrid, &
                 "The number of hybrid layers with Hybgen regridding, or 0 to use all "//&
                 "isopycnal layers.", default=0)
  call get_param(param_file, mdl, "HYBGEN_N_SIGMA", CS%nsigma, &
                 "The number of sigma-coordinate (terrain-following) layers with Hybgen regridding.", &
                 default=0)
  call get_param(param_file, mdl, "HYBGEN_DEEP_DZ_PR0FILE", CS%dp0k, &
                 "The layerwise list of deep z-level minimum thicknesses for Hybgen (dp0k in Hycom).", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_SHALLOW_DZ_PR0FILE", CS%ds0k, &
                 "The layerwise list of shallow z-level minimum thicknesses for Hybgen (ds0k in Hycom).", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_ISOPYCNAL_DZ_MIN", CS%dp00i, &
                 "The Hybgen deep isopycnal spacing minimum thickness (dp00i in Hycom)", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_MIN_ISO_DEPTH", CS%topiso_const, &
                 "The Hybgen shallowest depth for isopycnal layers (isotop in Hycom)", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_RELAX_PERIOD", CS%qhybrlx, &
                 "The Hybgen relaxation inteval in timesteps, or 1 for no relaxation (qhybrlx in Hycom)", &
                 units="timesteps", default=1.0)
  call get_param(param_file, mdl, "HYBGEN_BBL_THICKNESS", CS%thkbot, &
                 "A bottom boundary layer thickness within which Hybgen is able to move "//&
                 "overlying layers upward to match a target density.", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_REMAP_DENSITY_MATCH", CS%hybiso, &
                 "A tolerance between the layer densities and their target, within which "//&
                 "Hybgen determines that remapping uses PCM for a layer.", &
                 units="kg m-3", default=0.0, scale=US%kg_m3_to_R)
  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  CS%onem = 1.0 * GV%m_to_H

  do k=1,CS%nk ;  CS%target_density(k) = GV%Rlay(k) ; enddo

  ! reference density for anomalies [R ~> kg m-3]
  CS%thbase = 1000.0*US%kg_m3_to_R

  ! Determine the depth range over which to use a sigma (terrain-following) coordinate.
  ! --- terrain following starts at depth dpns and ends at depth dsns
  if (CS%nsigma == 0) then
    CS%dpns = CS%dp0k(1)
    CS%dsns = 0.0
  else
    CS%dpns = 0.0
    CS%dsns = 0.0
    do k=1,CS%nsigma
      CS%dpns = CS%dpns + CS%dp0k(k)
      CS%dsns = CS%dsns + CS%ds0k(k)
    enddo !k
  endif !nsigma

end subroutine init_hybgen_regrid

!> This subroutine deallocates memory in the control structure for the hybgen module
subroutine end_hybgen_regrid(CS)
  type(hybgen_regrid_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return

  deallocate(CS%target_density)
  deallocate(CS%dp0k, CS%ds0k)
  deallocate(CS)
end subroutine end_hybgen_regrid

!> This subroutine can be used to set the parameters for the hybgen module
subroutine set_hybgen_regrid_params(CS, min_thickness)
  type(hybgen_regrid_CS),  pointer    :: CS !< Coordinate regridding control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hybgen_params: CS not associated")

!  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_hybgen_regrid_params


!> Modify the input grid to give a new vertical grid based on the HYCOM hybgen code.
subroutine hybgen_regrid(G, GV, US, dp, tv, CS, dzInterface, PCM_cell)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp  !< Source grid layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv  !< Thermodynamics structure
  type(hybgen_regrid_CS),  intent(in)    :: CS  !< hybgen control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1), &
                           intent(inout) :: dzInterface !< The change in height of each interface,
                                                !! using a sign convention opposite to the change
                                                !! in pressure [H ~> m or kg m-2]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(inout) :: PCM_cell !< If true, PCM remapping should be used in a cell.
                                                !! This is effectively intent out, but values in wide
                                                !! halo regions and land points are reused.

!
! --- -------------------------------------
! --- hybrid grid generator from HYCOM
! --- -------------------------------------
! These notes on the parameters for the hybrid grid generator are inhereted from the
! Hycom source code for these algorithms.
!
! From blkdat.input (units may have changed from m to pressure):
!
! --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
! --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
! --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
! ---              k=1,kdm; dp0k must be zero for k>nhybrd
! --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
! ---              k=1,nsigma
! --- 'dp00i'  = deep isopycnal spacing minimum thickness (m)
! --- 'isotop' = shallowest depth for isopycnal layers     (m)
!                now in topiso(:,:)
! --- 'sigma ' = isopycnal layer target densities (sigma units)
! ---            now in theta(:,:,1:kdm)
!
! --- the above specifies a vertical coord. that is isopycnal or:
! ---  near surface z in    deep water, based on dp0k
! ---  near surface z in shallow water, based on ds0k and nsigma
! ---   terrain-following between them, based on ds0k and nsigma
!
! --- terrain following starts at depth dpns=sum(dp0k(k),k=1,nsigma) and
! --- ends at depth dsns=sum(ds0k(k),k=1,nsigma), and the depth of the
! --- k-th layer interface varies linearly with total depth between
! --- these two reference depths, i.e. a z-sigma-z fixed coordinate.
!
! --- near the surface (i.e. shallower than isotop), layers are always
! --- fixed depth (z or sigma).
! --  layer 1 is always fixed, so isotop=0.0 is not realizable.
! --- near surface layers can also be forced to be fixed depth
! --- by setting target densities (sigma(k)) very small.
!
! --- away from the surface, the minimum layer thickness is dp00i.
!
! --- for fixed depth targets to be:
! ---  z-only set nsigma=0,
! ---  sigma-z (shallow-deep) use a very small ds0k(:),
! ---  sigma-only set nsigma=kdm, dp0k large, and ds0k small.
!

  ! These arrays work with the input column
  real :: p_col(GV%ke)      ! A column of reference pressures [R L2 T-2 ~> Pa]
  real :: temp_in(GV%ke)    ! A column of input potential temperatures [degC]
  real :: saln_in(GV%ke)    ! A column of input layer salinities [ppt]
  real :: th3d_in(GV%ke)    ! An input column of coordinate potential density [R ~> kg m-3]
  real :: dp_in(GV%ke)      ! The input column of layer thicknesses [H ~> m or kg m-2]
  real :: pres_in(max(GV%ke,CS%nk)+1) ! original layer interface positions, padded with extra
                            ! massless layers if the input column has fewer layers than
                            ! the new grid. [H ~> m or kg m-2]
  logical :: PCM_lay(GV%ke) ! If true for a layer, use PCM remapping for that layer

  ! These arrays are on the target grid.
  real :: theta_i_j(CS%nk)  ! Target potential density [R ~> kg m-3]
  real :: th3d_i_j(CS%nk)   ! Initial values of coordinate potential density on the target grid [R ~> kg m-3]
  real :: dp_i_j(CS%nk)     ! A column of layer thicknesses [H ~> m or kg m-2]
  real :: p_i_j(CS%nk+1)    ! A column of interface depths [H ~> m or kg m-2]
  real :: dz_int(CS%nk+1)   ! The change in interface height due to remapping [H ~> m or kg m-2]
  real :: th3d_integral     ! Integrated coordinate potential density in a layer [R H ~> kg m-2 or kg2 m-5]

  real :: qdep              ! fraction not terrain following [nondim]
  real :: qhrlx( CS%nk+1)   ! relaxation coefficient [inverse timesteps?]
  real :: dp0ij( CS%nk)     ! minimum layer thickness [H ~> m or kg m-2]
  real :: dp0cum(CS%nk+1)   ! minimum interface depth [H ~> m or kg m-2]

  real :: h_tot             ! Total thickness of the water column [H ~> m or kg m-2]
  real :: dpthin            ! A very thin layer thickness, that is remapped differently [H ~> m or kg m-2]
  integer :: fixlay         ! Deepest fixed coordinate layer
  integer, dimension(0:CS%nk) :: k_end ! The index of the deepest source layer that contributes to
                            ! each target layer, in the unusual case where the the input grid is
                            ! larger than the new grid.  This situation only occurs during certain
                            ! types of initialization or when generating output diagnostics.
  integer :: i, j, k, kdm, m, k2, nk_in

  kdm = CS%nk

  p_col(:) = CS%ref_pressure
  dpthin = 1.0e-6*CS%onem

  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then

    ! --- store one-dimensional arrays of -p- for the 'old'  vertical grid before regridding
    pres_in(1) = 0.0
    do K=1,GV%ke
      temp_in(k) = tv%T(i,j,k)
      saln_in(k) = tv%S(i,j,k)
      dp_in(k) = dp(i,j,k)
      pres_in(K+1) = pres_in(K) + dp_in(k)
    enddo

    ! This sets the input column's coordinate potential density from T and S.
    call calculate_density(temp_in, saln_in, p_col, th3d_in, tv%eqn_of_state)

    ! Set the initial properties on the new grid from the old grid.
    p_i_j(1) = pres_in(1)
    nk_in = GV%ke
    if (GV%ke > CS%nk) then ; do k=GV%ke,CS%nk+1,-1
      ! Remove any excess massless layers from the bottom of the input column.
      if (dp_in(k) > 0.0) exit
      nk_in = k-1
    enddo ; endif

    if (CS%nk >= nk_in) then
      ! Simply copy over the common layers.  This is the usual case.
      do k=1,min(CS%nk,GV%ke)
        dp_i_j(k) = dp_in(k)
        p_i_j(K+1) = p_i_j(K) + dp_i_j(k)
        th3d_i_j(k) = th3d_in(k)
      enddo
      if (CS%nk > GV%ke) then
        ! Pad out the input column with additional massless layers with the bottom properties.
        ! This case only occurs during initialization or perhaps when writing diagnostics.
        do k=GV%ke+1,CS%nk
          pres_in(K+1) = pres_in(GV%ke+1)
          th3d_i_j(k) = th3d_in(GV%ke)
          dp_i_j(k) = 0.0
        enddo
      endif
    else ! (CS%nk < nk_in)
      ! The input column has more data than the output.  For now, combine layers to
      ! make them the same size, but there may be better approaches that should be taken.
      ! This case only occurs during initialization or perhaps when writing diagnostics.
      ! This case was not handled by the original Hycom code in hybgen.F90.
      do k=0,CS%nk ; k_end(k) = (k * nk_in) / CS%nk ; enddo
      do k=1,CS%nk
        dp_i_j(k) = 0.0 ; th3d_integral = 0.0
        do k2=k_end(k-1)+1,k_end(k)
          dp_i_j(k) = dp_i_j(k) + dp_in(k2)
          th3d_integral = th3d_integral + dp_in(k2)*th3d_in(k2)
        enddo
        if (dp_i_j(k) > GV%H_subroundoff) then
          ! Take the volume-weighted average properties.
          th3d_i_j(k) = th3d_integral / dp_i_j(k)
        else ! Take the properties of the topmost source layer that contributes.
          th3d_i_j(k) = th3d_in(k_end(k-1)+1)
        endif
      enddo
    endif

    ! Set the target densities for the new layers.
    do k=1,CS%nk
      ! theta_i_j(k) = theta(i,j,k)  ! If a 3-d target density were set up in theta, use that here.
      theta_i_j(k) = CS%target_density(k)  ! MOM6 does not yet support 3-d target densities.
    enddo

    h_tot = pres_in(GV%ke+1)

    call hybgen_column_init(kdm, CS%nhybrid, CS%nsigma, CS%dp0k, CS%ds0k, CS%dp00i, &
                            CS%topiso_const, CS%qhybrlx, CS%dpns, CS%dsns, h_tot, &
                            dp_i_j, fixlay, qdep, qhrlx, dp0ij, dp0cum, p_i_j)

    ! Determine whether to require the use of PCM remapping from each source layer.
    do k=1,GV%ke
      if (CS%hybiso > 0.0) then
        ! --- thin or isopycnal source layers are remapped with PCM.
        PCM_lay(k) = (k > fixlay) .and. &
            ((k > CS%nhybrid) .or. (pres_in(k+1)-pres_in(k) <= dpthin) .or. &
             (abs(th3d_i_j(k)-theta_i_j(k)) < CS%hybiso))
      else ! hybiso==0.0, so purely isopycnal layers use PCM
        PCM_lay(k) = (k > CS%nhybrid)
      endif ! hybiso
    enddo !k

    call hybgenaij_regrid(CS, kdm, CS%nhybrid, CS%thbase, CS%thkbot, &
                          CS%onem, 1.0e-11*US%kg_m3_to_R, &
                          theta_i_j, fixlay, qhrlx, dp0ij, dp0cum, &
                          th3d_i_j, dp_i_j, p_i_j, dz_int)


    ! Store the output from hybgenaij_regrid in 3-d arrays.
    if (present(PCM_cell)) then ; do k=1,GV%ke
      PCM_cell(i,j,k) = PCM_lay(k)
    enddo ; endif

    do K=1,kdm+1
      ! Note that dzInterface uses the opposite sign convention from the change in p.
      dzInterface(i,j,K) = -dz_int(K)
    enddo

    do K=1,kdm+1
      if (abs(dz_int(K) + (pres_in(K) - p_i_j(K))) > max(1.0e-13*p_i_j(kdm+1),kdm*GV%Angstrom_H,GV%H_subroundoff)) then
        call MOM_error(FATAL, "Mismatched interface height changes in hybgen_regrid.")
      endif
    enddo
  else
    if (present(PCM_cell)) then ; do k=1,GV%ke
      PCM_cell(i,j,k) = .false.
    enddo ; endif
    do k=1,CS%nk+1 ; dzInterface(i,j,k) = 0.0 ; enddo
  endif ; enddo ; enddo !i & j.

end subroutine hybgen_regrid

!> Initialize some of the variables that are used for regridding or unmixing, including the
!! previous interface heights and contraits on where the new interfaces can be.
subroutine hybgen_column_init(kdm, nhybrd, nsigma, dp0k, ds0k, dp00i, topiso_i_j, &
                          qhybrlx, dpns, dsns, h_tot, dp_i_j, &
                          fixlay, qdep, qhrlx, dp0ij, dp0cum, p_i_j)
  integer, intent(in)    :: kdm          !< The number of layers in the new grid
  integer, intent(in)    :: nhybrd       !< The number of hybrid layers (typically kdm)
  integer, intent(in)    :: nsigma       !< The number of sigma  levels (nhybrd-nsigma z-levels)
  real,    intent(in)    :: dp0k(kdm)    !< Layer deep z-level spacing minimum thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: ds0k(nsigma) !< Layer shallow z-level spacing minimum thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dp00i        !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real,    intent(in)    :: topiso_i_j   !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real,    intent(in)    :: qhybrlx      !< relaxation coefficient, 1/s?
  real,    intent(in)    :: h_tot        !< Total thickness of the water column [H ~> m or kg m-2]
  real,    intent(in)    :: dp_i_j(kdm)  !< Initial layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dpns         !< Vertical sum of dp0k [H ~> m or kg m-2]
  real,    intent(in)    :: dsns         !< Vertical sum of ds0k [H ~> m or kg m-2]
  integer, intent(out)   :: fixlay       !< Deepest fixed coordinate layer
  real,    intent(out)   :: qdep         !< fraction dp0k (vs ds0k) [nondim]
  real,    intent(out)   :: qhrlx( kdm+1) !< Interface relaxation coefficient [timesteps-1?]
  real,    intent(out)   :: dp0ij( kdm)   !< minimum layer thickness [H ~> m or kg m-2]
  real,    intent(out)   :: dp0cum(kdm+1) !< minimum interface depth [H ~> m or kg m-2]
  real,    intent(out)   :: p_i_j(kdm+1)  !< p(i,j,:), interface depths [H ~> m or kg m-2]
!
! --- --------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - initialization.
! --- --------------------------------------------------------------
!
  character(len=256) :: mesg  ! A string for output messages
  real :: hybrlx  ! The relaxation rate in the hybrid region [timestep-1]?
  real :: q       ! A portion of the thickness that contributes to the new cell [H ~> m or kg m-2]
  integer :: k, fixall
!
  hybrlx = 1.0 / qhybrlx
!
! --- dpns = sum(dp0k(k),k=1,nsigma)
! --- dsns = sum(ds0k(k),k=1,nsigma)
! --- terrain following starts (on the deep side) at depth dpns and ends (on the
! --- shallow side) at depth dsns and the depth of the k-th layer interface varies
! --- linearly with total depth between these two reference depths.
!
  if ((dpns <= dsns) .or. (h_tot >= dpns)) then
    qdep = 1.0  !not terrain following
  else
    qdep = max( 0.0, min( 1.0, (h_tot - dsns) / (dpns - dsns)) )
  endif

  if (qdep < 1.0) then
! ---   terrain following, qhrlx=1 and ignore dp00
    p_i_j( 1) = 0.0
    dp0cum(1) = 0.0
    qhrlx( 1) = 1.0
    dp0ij( 1) = qdep*dp0k(1) + (1.0-qdep)*ds0k(1)

    dp0cum(2) = dp0cum(1)+dp0ij(1)
    qhrlx( 2) = 1.0
    p_i_j( 2) = p_i_j(1)+dp_i_j(1)
    do k=2,kdm
      qhrlx( k+1) = 1.0
      dp0ij( k)   = qdep*dp0k(k) + (1.0-qdep)*ds0k(k)
      dp0cum(k+1) = dp0cum(k)+dp0ij(k)
      p_i_j( k+1) = p_i_j(k)+dp_i_j(k)
    enddo !k
  else
! ---   not terrain following
    p_i_j( 1) = 0.0
    dp0cum(1) = 0.0
    qhrlx( 1) = 1.0 !no relaxation in top layer
    dp0ij( 1) = dp0k(1)

    dp0cum(2) = dp0cum(1)+dp0ij(1)
    qhrlx( 2) = 1.0 !no relaxation in top layer
    p_i_j( 2) = p_i_j(1)+dp_i_j(1)
    do k=2,kdm
      if ((dp0k(k) <= dp00i) .or. (dp0k(k) >= p_i_j(k)-dp0cum(k))) then
        ! This layer is in fixed surface coordinates.
        dp0ij(k) = dp0k(k)
        qhrlx(k+1) = 1.0 ! 1 at  dp0k
      else
        q = dp0k(k) * (dp0k(k) / ( p_i_j(k)-dp0cum(k)) ) ! A fraction between 0 and 1 of dp0 to use here.
        if (dp00i >= q) then
          ! This layer is much deeper than the fixed surface coordinates.
          dp0ij(k) = dp00i
          qhrlx(k+1) = qhybrlx  ! 1 at  dp0k, qhybrlx at dp00i
        else
          ! This layer spans the margines of the fixed surface coordinates.
          ! In this case dp00i < q < dp0k.
          dp0ij(k) = q
          qhrlx(k+1) = qhybrlx * (dp0k(k) - dp00i) / &
                       ((dp0k(k)-q) + (q-dp00i)*qhybrlx)  !1 at  dp0k, qhybrlx at dp00i
        endif

        ! The old equivalent code is:
        ! hybrlx = 1.0 / qhybrlx
        ! q = max( dp00i, dp0k(k) * (dp0k(k) / max(dp0k( k), p_i_j(k)-dp0cum(k)) ) )
        ! qts = 1.0 - (q-dp00i) / (dp0k(k)-dp00i)  !0 at q = dp0k, 1 at q=dp00i
        ! qhrlx( k+1) = 1.0 / (1.0 + qts*(hybrlx-1.0))  !1 at dp0k, qhybrlx at dp00i
      endif
      dp0cum(k+1) = dp0cum(k) + dp0ij(k)
      p_i_j( k+1) = p_i_j(k) + dp_i_j(k)
    enddo !k
  endif !qdep<1:else

! --- identify the current fixed coordinate layers
  fixlay = 1  !layer 1 always fixed
  do k=2,nhybrd
    if (dp0cum(k) >= topiso_i_j) then
      exit  !layers k to nhybrd might be isopycnal
    endif
! ---   top of layer is above topiso, i.e. always fixed coordinate layer
    qhrlx(k+1) = 1.0  !no relaxation in fixed layers
    fixlay     = fixlay+1
  enddo !k

  fixall = fixlay
  do k=fixall+1,nhybrd
    if (p_i_j(k+1) > dp0cum(k+1)+0.1*dp0ij(k)) then
      if ( (fixlay > fixall) .and. (p_i_j(k) > dp0cum(k)) ) then
        ! --- The previous layer should remain fixed.
        fixlay = fixlay-1
      endif
      exit  !layers k to nhybrd might be isopycnal
    endif
! ---   sometimes fixed coordinate layer
    qhrlx(k) = 1.0  !no relaxation in fixed layers
    fixlay   = fixlay+1
  enddo !k

end subroutine hybgen_column_init

!> The cushion function from Bleck & Benjamin, 1992, which returns a smoothly varying
!! but limited value that goes between dp0 and delp
real function cushn(delp, dp0)
  real, intent(in) :: delp  ! A thickness change
  real, intent(in) :: dp0   ! A thickness change

  real :: qq  ! A limited ratio of delp/dp0 [nondim]

  ! These are the nondimensional parameters that define the cushion function.
  real, parameter :: qqmn=-4.0, qqmx=2.0  ! shifted range for cushn [nondim]
! real, parameter :: qqmn=-2.0, qqmx=4.0  ! traditional range for cushn [nondim]
! real, parameter :: qqmn=-4.0, qqmx=6.0  ! somewhat wider range for cushn [nondim]
  real, parameter :: cusha = qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2
  real, parameter :: cushb=1.0/qqmn

  ! --- if delp >= qqmx*dp0 >>  dp0, cushn returns delp.
  ! --- if delp <= qqmn*dp0 << -dp0, cushn returns dp0.

  qq = max(qqmn, min(qqmx, delp/dp0))
  cushn = dp0 * (1.0 + cusha * (1.0-cushb*qq)**2) * max(1.0, delp/(dp0*qqmx))

end function cushn

!> Create a new grid for a column of water using the Hybgen algorithm.
subroutine hybgenaij_regrid(CS, kdm, nhybrd, thbase, thkbot, &
                            onem, epsil, theta_i_j, &
                            fixlay, qhrlx, dp0ij, dp0cum, &
                            th3d_i_j, h_in, p_i_j, dp_int)
!
  type(hybgen_regrid_CS), intent(in)    :: CS  !< hybgen regridding control structure
  integer, intent(in)    :: kdm            !< number of layers
  integer, intent(in)    :: nhybrd         !< number of hybrid layers (typically kdm)
  real,    intent(in)    :: thbase         !< reference density [R ~> kg m-3]
  real,    intent(in)    :: thkbot         !< thickness of bottom boundary layer [H ~> m or kg m-2]
  real,    intent(in)    :: onem           !< one m in pressure units [H ~> m or kg m-2]
  real,    intent(in)    :: epsil          !< small nonzero density to prevent division by zero [R ~> kg m-3]
  real,    intent(in)    :: theta_i_j(kdm) !< theta(i,j,:) target density [R ~> kg m-3]
  integer, intent(in)    :: fixlay         !< deepest fixed coordinate layer
  real,    intent(in)    :: qhrlx( kdm+1)  !< relaxation coefficient per timestep [nondim]
  real,    intent(in)    :: dp0ij( kdm)    !< minimum layer thickness [H ~> m or kg m-2]
  real,    intent(in)    :: dp0cum(kdm+1)  !< minimum interface depth [H ~> m or kg m-2]
  real,    intent(in)    :: th3d_i_j(kdm)  !< Coordinate potential density [R ~> kg m-3]
  real,    intent(in)    :: h_in(kdm)      !< Layer thicknesses [H ~> m or kg m-2]
  real,    intent(inout) :: p_i_j(kdm+1)   !< layer interface positions [H ~> m or kg m-2]
  real,    intent(out)   :: dp_int(kdm+1)  !< The change in interface positions [H ~> m or kg m-2]

! --- ------------------------------------------------------
! --- hybrid grid generator, single column(part A) - regrid.
! --- ------------------------------------------------------
  real :: p_new  ! A new interface position [H ~> m or kg m-2]
!  real :: p_i_j(kdm+1) ! layer interface positions [H ~> m or kg m-2]
  real :: h_i_j(kdm)    ! layer thicknesses [H ~> m or kg m-2]
  real :: dp_rem        ! Remaining water to move [H ~> m or kg m-2]
  real :: q_frac ! A fraction of a layer to entrain [nondim]
  real :: h_hat3 ! Thickness movement upward across the interface between layers k-2 and k-3 [H ~> m or kg m-2]
  real :: h_hat2 ! Thickness movement upward across the interface between layers k-1 and k-2 [H ~> m or kg m-2]
  real :: h_hat  ! Thickness movement upward across the interface between layers k and k-1 [H ~> m or kg m-2]
  real :: h_hat0 ! A first guess at thickness movement upward across the interface
                 ! between layers k and k-1 [H ~> m or kg m-2]
  real :: dh_cor ! Thickness changes [H ~> m or kg m-2]
  real :: h1_tgt ! A target thickness for the top layer [H ~> m or kg m-2]
  real :: tenm   ! ten m  in pressure units [H ~> m or kg m-2]
  real :: onemm  ! one mm in pressure units [H ~> m or kg m-2]
  integer k, ka, ktr
  character(len=256) :: mesg  ! A string for output messages

  ! This line needs to be consistent with the parameters set in cushn().
  real, parameter :: qqmn=-4.0, qqmx=2.0  ! shifted range for cushn
! real, parameter :: qqmn=-2.0, qqmx=4.0  ! traditional range for cushn
! real, parameter :: qqmn=-4.0, qqmx=6.0  ! somewhat wider range for cushn

  tenm = 10.0*onem
  onemm = 0.001*onem

  do K=1,kdm+1 ; dp_int(K) = 0.0 ; enddo

  p_i_j(1) = 0.0
  do k=1,kdm
    h_i_j(k) = h_in(k)
    p_i_j(K+1) = p_i_j(K) + h_i_j(k)
  enddo

! --- try to restore isopycnic conditions by moving layer interfaces
! --- qhrlx(k) are relaxation coefficients (inverse baroclinic time steps)

  if (fixlay >= 1) then
! ---   maintain constant thickness in layer k = 1, entraining from layers below as necessary.
    k = 1
    h1_tgt = min(dp0ij(k), p_i_j(kdm+1))
    dp_rem = h1_tgt - h_i_j(1)
    h_i_j(1) = h1_tgt
    p_i_j(2) = h1_tgt
    dp_int(2) = dp_rem
    if (dp_rem <= h_i_j(2)) then
      h_i_j(2) = h_i_j(2) - dp_rem
      dp_rem = 0.0
    else
      h_i_j(2) = 0.0
      dp_rem = dp_rem - h_i_j(2)
      do K=2,kdm
        dp_int(K+1) = dp_int(K+1) + (h1_tgt - p_i_j(K+1))
        p_i_j(K+1) = h1_tgt
        if (dp_rem <= h_i_j(k+1)) then
          h_i_j(k+1) = h_i_j(k+1) - dp_rem
          dp_rem = 0.0
          exit  ! usually get here quickly
        else
          dp_rem = dp_rem - h_i_j(k+1)
          h_i_j(k+1) = 0.0
        endif
      enddo
    endif
  endif

  do k=2,nhybrd

    if (k <= fixlay) then
! ---     maintain constant thickness, k <= fixlay
      if (k < kdm) then  !p.kdm+1 not changed
        p_new = min(dp0cum(k+1), p_i_j(kdm+1))
        h_i_j(k)   = h_i_j(k)   + (p_new - p_i_j(K+1))
        h_i_j(k+1) = h_i_j(k+1) - (p_new - p_i_j(K+1))
        dp_int(K+1) = dp_int(K+1) + (p_new - p_i_j(K+1))
        p_i_j(k+1) = p_new
        if (k == fixlay) then
! ---         enforce interface order (may not be necessary).
          do ka=k+2,kdm
            if (p_i_j(Ka) >= p_i_j(K+1)) then
              exit  ! usually get here quickly
            else
              h_i_j(ka-1) = h_i_j(ka-1) + (p_i_j(K+1) - p_i_j(Ka))
              h_i_j(ka)   = h_i_j(ka)   - (p_i_j(K+1) - p_i_j(Ka))
              dp_int(Ka) = dp_int(Ka) + (p_i_j(K+1) - p_i_j(Ka))
              p_i_j(Ka) = p_i_j(K+1)
            endif
          enddo !ka
        endif !k == fixlay
      endif !k < kdm

    else
! ---     do not maintain constant thickness, k > fixlay

      if ((th3d_i_j(k) > theta_i_j(k)+epsil) .and.  (k > fixlay+1)) then
! ---       water in layer k is too dense
! ---       try to dilute with water from layer k-1
! ---       do not move interface if k = fixlay + 1

        if (th3d_i_j(k-1) >= theta_i_j(k-1) .or. &
            p_i_j(k) <= dp0cum(k)+onem .or. &
            h_i_j(k) <= h_i_j(k-1)) then
! ---         if layer k-1 is too light, thicken the thinner of the two,
! ---         i.e. skip this layer if it is thicker.

          if     ((theta_i_j(k)-th3d_i_j(k-1)) <= epsil) then
!               layer k-1 much too dense, take entire layer
            h_hat0 = 0.0  ! This line was not in the Hycom version of hybgen.F90.
            h_hat = dp0ij(k-1) - h_i_j(k-1)
          else
            q_frac = (theta_i_j(k)-th3d_i_j(k)) / &
                     (theta_i_j(k)-th3d_i_j(k-1))    ! -1 <= q_frac < 0
            h_hat0 = q_frac*h_i_j(k)  ! -h_i_j(k-1) <= h_hat0 < 0
            if (k == fixlay+2) then
! ---             treat layer k-1 as fixed.
              h_hat = max(h_hat0, dp0ij(k-1) - h_i_j(k-1))
            else
! ---             maintain minimum thickess of layer k-1.
              h_hat = cushn(h_hat0 + h_i_j(k-1), dp0ij(k-1)) - h_i_j(k-1)
            endif !fixlay+2:else
          end if
          ! h_hat is negative, so this check is unnecessary.
          h_hat = min(h_hat, p_i_j(kdm+1) - p_i_j(k))

! ---         if isopycnic conditions cannot be achieved because of a blocking
! ---         layer in the interior ocean, move interface k-1 (and k-2 if
! ---         necessary) upward
          if (k <= fixlay+2) then
! ---           do nothing.
          else if ((h_hat >= 0.0) .and. &
                   p_i_j(k-1) > dp0cum(k-1)+tenm .and. &
                  (p_i_j(kdm+1)-p_i_j(k-1) < thkbot .or. &
                   h_i_j(k-2) > qqmx*dp0ij(k-2))) then ! k > 2
            if (k == fixlay+3) then
! ---             treat layer k-2 as fixed.
              h_hat2 = max(h_hat0 - h_hat, dp0ij(k-2) - h_i_j(k-2))
            else
! ---             maintain minimum thickness of layer k-2.
              h_hat2 = cushn(h_i_j(k-2) + (h_hat0 - h_hat), dp0ij(k-2)) - h_i_j(k-2)
            endif !fixlay+3:else
            if (h_hat2 < -onemm) then
              dh_cor = qhrlx(k-1) * max(h_hat2, -h_hat - h_i_j(k-1))
              h_i_j(k-2)  = h_i_j(k-2) + dh_cor
              h_i_j(k-1)  = h_i_j(k-1) - dh_cor
              dp_int(k-1) = dp_int(k-1) + dh_cor
              p_i_j(k-1)  = p_i_j(k-1) + dh_cor
              h_hat = cushn(h_hat0 + h_i_j(k-1), dp0ij(k-1)) - h_i_j(k-1)
            elseif (k <= fixlay+3) then
! ---             do nothing.
            elseif (p_i_j(k-2) > dp0cum(k-2)+tenm .and. &
                   (p_i_j(kdm+1)-p_i_j(k-2) < thkbot .or. &
                    h_i_j(k-3) > qqmx*dp0ij(k-3))) then
              if (k == fixlay+4) then
! ---               treat layer k-3 as fixed.
                h_hat3 = max(h_hat0 - h_hat, dp0ij(k-3) - h_i_j(k-3))
              else
! ---               maintain minimum thickess of layer k-3.
                h_hat3 = cushn(h_i_j(k-3) + (h_hat0 - h_hat), dp0ij(k-3)) - h_i_j(k-3)
              endif !fixlay+4:else
              if (h_hat3 < -onemm) then
                dh_cor = qhrlx(k-2) * max(h_hat3, -h_i_j(k-2))
                h_i_j(k-3) = h_i_j(k-3) + dh_cor
                h_i_j(k-2) = h_i_j(k-2) - dh_cor
                dp_int(k-2) = dp_int(k-2) + dh_cor
                p_i_j(k-2) = p_i_j(k-2) + dh_cor

                h_hat2 = cushn(h_i_j(k-2) + (h_hat0 - h_hat), dp0ij(k-2)) - h_i_j(k-2)
                if (h_hat2 < -onemm) then
                  dh_cor = qhrlx(k-1) * (max(h_hat2, -h_hat - h_i_j(k-1)) )
                  h_i_j(k-2) = h_i_j(k-2) + dh_cor
                  h_i_j(k-1) = h_i_j(k-1) - dh_cor
                  dp_int(k-1) = dp_int(k-1) + dh_cor
                  p_i_j(k-1) = p_i_j(k-1) + dh_cor
                  h_hat = cushn(h_hat0 + h_i_j(k-1), dp0ij(k-1)) - h_i_j(k-1)
                endif !h_hat2
              endif !h_hat3
            endif !h_hat2:blocking
          endif !blocking
!
          if (h_hat < 0.0) then
! ---           entrain layer k-1 water into layer k, move interface up.
            dh_cor = qhrlx(k) * h_hat
            h_i_j(k-1) = h_i_j(k-1) + dh_cor
            h_i_j(k)   = h_i_j(k)   - dh_cor
            dp_int(k) = dp_int(k) + dh_cor
            p_i_j(k) = p_i_j(k) + dh_cor
          endif !entrain

        endif  !too-dense adjustment
!
      elseif (th3d_i_j(k) < theta_i_j(k)-epsil) then   ! layer too light
!
! ---       water in layer k is too light
! ---       try to dilute with water from layer k+1
! ---       do not entrain if layer k touches bottom
!
        if (p_i_j(k+1) < p_i_j(kdm+1)) then  ! k<kdm
          if (th3d_i_j(k+1) <= theta_i_j(k+1) .or. &
              p_i_j(k+1) <= dp0cum(k+1)+onem  .or. &
              h_i_j(k) < h_i_j(k+1)) then
! ---           if layer k+1 is too dense, thicken the thinner of the
! ---           two, i.e. skip this layer (never get here) if it is not
! ---           thinner than the other.

            if     ((th3d_i_j(k+1)-theta_i_j(k)) <= epsil) then
!                 layer k-1 too light, take entire layer
              h_hat = h_i_j(k+1)
            else
              q_frac = (theta_i_j(k) - th3d_i_j(k)) / (th3d_i_j(k+1)-theta_i_j(k)) ! 0 < q_frac <= 1
              h_hat = q_frac*h_i_j(k)
            endif
!
! ---           if layer k+1, or layer k+2, does not touch the bottom
! ---           then maintain minimum thicknesses of layers k and k+1 as
! ---           much as possible. otherwise, permit layers to collapse
! ---           to zero thickness at the bottom.
            if (p_i_j(min(k+3,kdm+1)) < p_i_j(kdm+1)) then
              if (p_i_j(kdm+1)-p_i_j(k) >  dp0ij(k)+dp0ij(k+1)) then
                h_hat = h_i_j(k+1) - cushn(h_i_j(k+1)-h_hat, dp0ij(k+1))
              endif
              h_hat = max(h_hat, dp0ij(k) - h_i_j(k))
              h_hat = min(h_hat, max(0.5*h_i_j(k+1), h_i_j(k+1)-dp0ij(k+1)) )
            else
              h_hat = min(h_i_j(k+1), h_hat)
            endif !p.k+2<p.kdm+1

            if (h_hat > 0.0) then
! ---             entrain layer k+1 water into layer k.
              dh_cor = qhrlx(k+1) * h_hat
              h_i_j(k)   = h_i_j(k)   + dh_cor
              h_i_j(k+1) = h_i_j(k+1) - dh_cor
              dp_int(k+1) = dp_int(k+1) + dh_cor
              p_i_j(k+1) = p_i_j(k+1) + dh_cor
            endif !entrain

          endif !too-light adjustment
        endif !above bottom
      endif !too dense or too light
!
! ---     if layer above is still too thin, move interface down.
      dh_cor = min(qhrlx(k-1) * min(dp0ij(k-1) - h_i_j(k-1), p_i_j(kdm+1) - p_i_j(k)), h_i_j(k))
      if (dh_cor > 0.0) then
        h_i_j(k-1) = h_i_j(k-1) + dh_cor
        h_i_j(k)   = h_i_j(k)   - dh_cor
        dp_int(k) = dp_int(k) + dh_cor
        p_i_j(k) = p_i_j(k) + dh_cor
      endif

    endif !k <= fixlay:else

  enddo !k  vertical coordinate relocation

end subroutine hybgenaij_regrid

end module MOM_hybgen_regrid

! This code was translated in 2022 from the HYCOM hybgen code, which was primarily developed
! between 2000 and 2015, with some minor subsequent changes.
