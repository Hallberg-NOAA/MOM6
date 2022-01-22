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

  !> Global i-index of a point where detailed diagnostics of hybgen are desired
  integer :: itest = -1
  !> Global j-index of a point where detailed diagnostics of hybgen are desired
  integer :: jtest = -1

  real :: thkbot !< Thickness of a bottom boundary layer, within which hybgen does
                 !! something different. (m)

  !> Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real :: topiso_const
  ! real, dimension(:,:), allocatable :: topiso

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  real :: onem       !< Nominally one m in thickness units [H ~> m or kg m-2]

  logical :: debug   !< If true, write verbose checksums for debugging purposes.

end type hybgen_regrid_CS

public hybgen_regrid, init_hybgen_regrid, end_hybgen_regrid
public set_hybgen_regrid_params

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
                 "The Hybgen relaxation inteval in timesteps, or 1 for no relaxation (qhbrlx in Hycom)", &
                 units="timesteps", default=1.0)
  call get_param(param_file, mdl, "HYBGEN_BBL_THICKNESS", CS%thkbot, &
                 "A bottom boundary layer thickness within which Hybgen is able to move "//&
                 "overlying layers upward to match a target density.", &
                 units="m", default=0.0, scale=1.0)
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
subroutine hybgen_regrid(G, GV, US, CS, dp, tv, h_new, dzInterface, PCM_cell)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(hybgen_regrid_CS),  intent(in)    :: CS  !< hybgen control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp  !< Source grid layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv  !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk), &
                           intent(inout) :: h_new !< Destination grid layer thicknesses [H ~> m or kg m-2]
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
  real :: th3d_integral     ! Integrated coordinate potential density in a layer [R H ~> kg m-2 or kg2 m-5]

  real :: qdep              ! fraction not terrain following [nondim]
  real :: qhrlx( CS%nk+1)   ! relaxation coefficient [inverse timesteps?]
  real :: dp0ij( CS%nk)     ! minimum layer thickness [H ~> m or kg m-2]
  real :: dp0cum(CS%nk+1)   ! minimum interface depth [H ~> m or kg m-2]

  real :: depths_i_j        ! Bottom depth in thickness units [H ~> m or kg m-2]
  real :: dpthin            ! A very thin layer thickness, that is remapped differently [H ~> m or kg m-2]
  integer :: fixlay         ! Deepest fixed coordinate layer
  logical :: test_col       ! If true, write verbose debugging for this column.
  integer, dimension(0:CS%nk) :: k_end ! The index of the deepest source layer that contributes to
                            ! each target layer, in the unusual case where the the input grid is
                            ! larger than the new grid.  This situation only occurs during certain
                            ! types of initialization or when generating output diagnostics.
  integer :: i, j, k, kdm, m, k2, nk_in

  kdm = CS%nk

  p_col(:) = CS%ref_pressure
  dpthin = 1.0e-6*CS%onem

  test_col = .false.
  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then
!diag   test_col = ((i == CS%itest + G%isd - G%isd_global) .and. &
!diag               (j == CS%jtest + G%jsd - G%jsd_global))

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

    depths_i_j = GV%Z_to_H * G%bathyT(i,j)
    !### depths_i_j = pres_in(GV%ke+1)

    call hybgenaij_init(   kdm, CS%nhybrid, CS%nsigma, &
                           CS%dp0k, CS%ds0k, CS%dp00i, CS%topiso_const, CS%qhybrlx, &
                           CS%dpns, CS%dsns, &
                           depths_i_j, dp_i_j, &
                           fixlay, qdep, qhrlx, dp0ij, dp0cum, &
                           p_i_j, test_col)

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
                          th3d_i_j, p_i_j, pres_in, test_col)


    ! Store the output from hybgenaij_regrid in 3-d arrays.
    if (present(PCM_cell)) then ; do k=1,GV%ke
      PCM_cell(i,j,k) = PCM_lay(k)
    enddo ; endif
    do k=1,kdm
      h_new(i,j,k) = max(p_i_j(K+1)-p_i_j(K), 0.0)
      !### dp(i,j,k) = max(dp_i_j(k), 0.0)
    enddo

    do k=1,kdm
      ! Note that dzInterface uses the opposite sign convention from the change in p.
      dzInterface(i,j,K) = pres_in(K) - p_i_j(K)
    enddo
    dzInterface(i,j,1) = 0.0
    dzInterface(i,j,kdm+1) = 0.0
  else
    if (present(PCM_cell)) then ; do k=1,GV%ke
      PCM_cell(i,j,k) = .false.
    enddo ; endif
    do k=1,CS%nk+1 ; dzInterface(i,j,k) = 0.0 ; enddo
    do k=1,min(CS%nk,GV%ke) ; h_new(i,j,k) = dp(i,j,k) ; enddo
    do k=GV%ke+1,CS%nk ; h_new(i,j,k) = 0.0 ; enddo
  endif ; enddo ; enddo !i & j.

end subroutine hybgen_regrid

subroutine hybgenaij_init(kdm, nhybrd, nsigma, &
                                dp0k, ds0k, dp00i,topiso_i_j, qhybrlx, &
                                dpns, dsns, &
                                depths_i_j, &
                                dp_i_j, &
                                fixlay, qdep, qhrlx, dp0ij, dp0cum, &
                                p_i_j, test_col)
!
  integer, intent(in)    :: &
            kdm, &          ! The number of layers
            nhybrd, &       !number of hybrid layers (typically kdm)
            nsigma          !number of sigma  levels (nhybrd-nsigma z-levels)
  real,    intent(in)    :: &
            dp0k(kdm),    & !layer deep    z-level spacing minimum thicknesses [H ~> m or kg m-2]
            ds0k(nsigma), & !layer shallow z-level spacing minimum thicknesses [H ~> m or kg m-2]
            dp00i, &        !deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
            topiso_i_j, &   !shallowest depth for isopycnal layers [H ~> m or kg m-2]
            qhybrlx, &      !relaxation coefficient, 1/s
            depths_i_j, &   ! Bottom depth om thickness units [H ~> m or kg m-2]
            dp_i_j(kdm)     ! dp(i,j,:,n), layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dpns ! Vertical sum of dp0k [H ~> m or kg m-2]
  real,    intent(in)    :: dsns ! Vertical sum of ds0k [H ~> m or kg m-2]
  integer, intent(out)   :: fixlay           !deepest fixed coordinate layer
  real,    intent(out)   :: qdep             !fraction dp0k (vs ds0k) [nondim]
  real,    intent(out)   :: qhrlx( kdm+1)    !relaxation coefficient [timesteps-1?]
  real,    intent(out)   :: dp0ij( kdm)      !minimum layer thickness [H ~> m or kg m-2]
  real,    intent(out)   :: dp0cum(kdm+1)    !minimum interface depth [H ~> m or kg m-2]
  real,    intent(out)   :: p_i_j(kdm+1)     !p(i,j,:), interface depths [H ~> m or kg m-2]
  logical, intent(in)    :: test_col !< If true, write verbose debugging for this column.
!
! --- --------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - initialization.
! --- depth_i_j is in m, everything else is in pressure units
! --- --------------------------------------------------------------
!
  character(len=256) :: mesg  ! A string for output messages
  real :: hybrlx  ! The relaxation rate in the hybrid region [timestep-1]?
  real :: q       ! A portion of the thickness that contributes to the new cell [H ~> m or kg m-2]
  real :: qts  !### A temporary variable that should be refactored out for accuracy.
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
  if ((dpns <= dsns) .or. (depths_i_j >= dpns)) then
    qdep = 1.0  !not terrain following
  else
    qdep = max( 0.0, min( 1.0, (depths_i_j - dsns) / (dpns - dsns)) )
  endif
!
  if (qdep < 1.0) then
! ---   terrain following, qhrlx=1 and ignore dp00
    p_i_j( 1) = 0.0
    dp0cum(1) = 0.0
    qhrlx( 1) = 1.0
    dp0ij( 1) = qdep*dp0k(1) + (1.0-qdep)*ds0k(1)
!diag       if (test_col) then
!diag         k=1
!diag         write (mesg,*) 'qdep = ', qdep
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag         write (mesg,'(a/i6,1x,4f9.3/a)') &
!diag         '     k     dp0ij     ds0k     dp0k        p', &
!diag         k, dp0ij(k)*GV%H_to_m, ds0k(1)*GV%H_to_m, dp0k(1)*GV%H_to_m, &
!diag                                  p_i_j(k)*GV%H_to_m, &
!diag         '     k     dp0ij    p-cum        p   dp0cum'
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag       endif !debug
    dp0cum(2) = dp0cum(1)+dp0ij(1)
    qhrlx( 2) = 1.0
    p_i_j( 2) = p_i_j(1)+dp_i_j(1)
    do k=2,kdm
      qhrlx( k+1) = 1.0
      dp0ij( k)   = qdep*dp0k(k) + (1.0-qdep)*ds0k(k)
      dp0cum(k+1) = dp0cum(k)+dp0ij(k)
      p_i_j( k+1) = p_i_j(k)+dp_i_j(k)
!diag         if (test_col) then
!diag           write (mesg,'(i6,1x,4f9.3)') &
!diag           k, dp0ij(k)*GV%H_to_m, p_i_j(k)*GV%H_to_m-dp0cum(k)*GV%H_to_m, &
!diag                            p_i_j(k)*GV%H_to_m, dp0cum(k)*GV%H_to_m
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag         endif !debug
    enddo !k
  else
! ---   not terrain following
    p_i_j( 1) = 0.0
    dp0cum(1) = 0.0
    qhrlx( 1) = 1.0 !no relaxation in top layer
    dp0ij( 1) = dp0k(1)
!diag       if (test_col) then
!diag         k=1
!diag         write (mesg,*) 'qdep = ', qdep
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag         write (mesg,'(a/i6,1x,f9.3)') &
!diag       '     k     dp0ij     dp0k        q    p-cum        p   dp0cum', &
!diag             k, dp0ij(k)*GV%H_to_m
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag       endif !debug
    dp0cum(2) = dp0cum(1)+dp0ij(1)
    qhrlx( 2) = 1.0 !no relaxation in top layer
    p_i_j( 2) = p_i_j(1)+dp_i_j(1)
    do k=2,kdm
! ---     q is dp0k(k) when in surface fixed coordinates
! ---     q is dp00i   when much deeper than surface fixed coordinates
      if     (dp0k(k) <= dp00i) then
        q  = dp0k(k)
        qts = 0.0     !0 at dp0k
      else
        q  = max( dp00i, dp0k(k) * dp0k(k)/ &
                         max( dp0k( k), p_i_j(k)-dp0cum(k) ) )
        qts = 1.0 - (q-dp00i) / (dp0k(k)-dp00i)  !0 at dp0k, 1 at dp00i
      endif
      qhrlx( k+1) = 1.0 / (1.0 + qts*(hybrlx-1.0))  !1 at  dp0k, qhybrlx at dp00i
      dp0ij( k)   = min( q, dp0k(k) )
      dp0cum(k+1) = dp0cum(k)+dp0ij(k)
      p_i_j( k+1) = p_i_j(k)+dp_i_j(k)
!diag         if (test_col) then
!diag           write (mesg,'(i6,1x,6f9.3)') &
!diag             k, dp0ij(k)*GV%H_to_m, dp0k(k)*GV%H_to_m, q*GV%H_to_m, &
!diag             p_i_j(k)*GV%H_to_m-dp0cum(k)*GV%H_to_m, &
!diag             p_i_j(k)*GV%H_to_m, dp0cum(k)*GV%H_to_m
!diag           call MOM_mesg(mesg, all_print=.true.)
!diag         endif !debug
    enddo !k
  endif !qdep<1:else
!
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
!diag      if (test_col) then
!diag        write(mesg,'(a,i3)') &
!diag              'hybgen, always-fixed coordinate layers: 1 to ', fixlay
!diag        call MOM_mesg(mesg, all_print=.true.)
!diag      endif !debug
!
  fixall = fixlay
  do k=fixall+1,nhybrd
!diag        if (test_col) then
!diag          write (mesg,'(i6,1x,2f9.3)') &
!diag            k, p_i_j(k+1)*GV%H_to_m, dp0cum(k+1)*GV%H_to_m
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        endif !debug
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
!diag      if (test_col) then
!diag        write(mesg, '(a,i3)') &
!diag              'hybgen,        fixed coordinate layers: 1 to ', fixlay
!diag        call MOM_mesg(mesg, all_print=.true.)

!diag        call MOM_mesg('hybgen:   thkns  minthk     dpth  mindpth   hybrlx', all_print=.true.)
!diag        do k=1,kdm
!diag          write (mesg,'(i6,1x,2f8.3,2f9.3,f9.3)') &
!diag              k, dp_i_j(k)*GV%H_to_m,   dp0ij(k)*GV%H_to_m, &
!diag              p_i_j(k+1)*GV%H_to_m, dp0cum(k+1)*GV%H_to_m, 1.0/qhrlx(k+1)
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        enddo
!diag      endif !debug

end subroutine hybgenaij_init


subroutine hybgenaij_regrid(CS, kdm, nhybrd, thbase, thkbot, &
                            onem, epsil, theta_i_j, &
                            fixlay, qhrlx, dp0ij, dp0cum, &
                            th3d_i_j, p_i_j, pres, test_col)
!
  type(hybgen_regrid_CS), intent(in)    :: CS  !< hybgen regridding control structure
  integer, intent(in)    :: kdm            !< number of layers
  integer, intent(in)    :: nhybrd         !< number of hybrid layers (typically kdm)
  real,    intent(in)    :: thbase         !< reference density (sigma units)
  real,    intent(in)    :: thkbot         !< thickness of bottom boundary layer (m)
  real,    intent(in)    :: onem           !< one m in pressure units [H ~> m or kg m-2]
  real,    intent(in)    :: epsil          !< small nonzero density to prevent division by zero
  real,    intent(in)    :: theta_i_j(kdm) !< theta(i,j,:) target density
  integer, intent(in)    :: fixlay         !< deepest fixed coordinate layer
  real,    intent(in)    :: qhrlx( kdm+1)  !< relaxation coefficient
  real,    intent(in)    :: dp0ij( kdm)    !< minimum layer thickness
  real,    intent(in)    :: dp0cum(kdm+1)  !< minimum interface depth
  real,    intent(in)    :: th3d_i_j(kdm)  !< Coordinate potential density [R ~> kg m-3]
  real,    intent(inout) :: p_i_j(kdm+1)   !< layer interface positions [H ~> m or kg m-2]
  real,    intent(in)    :: pres(kdm+1)    !< original layer interfaces [H ~> m or kg m-2]
  logical, intent(in)    :: test_col !< If true, write verbose debugging for this column.

! --- ------------------------------------------------------
! --- hybrid grid generator, single column(part A) - regrid.
! --- ------------------------------------------------------
  real :: p_hat, p_hat0, p_hat2, p_hat3
  real :: q, qtr, thkbop
  real :: zthk, dpthin
  real :: tenm  ! ten m  in pressure units
  real :: onemm ! one mm in pressure units
  integer k, ka, ktr
  character(len=256) :: mesg  ! A string for output messages
!diag character*12 cinfo

  double precision, parameter ::   zp5=0.5    !for sign function
! --- c u s h i o n   function (from Bleck & Benjamin, 1992):
! --- if delp >= qqmx*dp0 >>  dp0, -cushn- returns -delp-
! --- if delp <= qqmn*dp0 << -dp0, -cushn- returns  -dp0-
  real :: qqmn, qqmx, cusha, cushb
  parameter (qqmn=-4.0, qqmx=2.0)  ! shifted range
!     parameter (qqmn=-2.0, qqmx=4.0)  ! traditional range
!     parameter (qqmn=-4.0, qqmx=6.0)  ! somewhat wider range
  parameter (cusha=qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2)
  parameter (cushb=1.0/qqmn)

  real :: qq, cushn, delp, dp0

  qq(   delp, dp0) = max(qqmn, min(qqmx, delp/dp0))
  cushn(delp, dp0) = dp0* &
                  (1.0+cusha*(1.0-cushb*qq(delp, dp0))**2)* &
                  max(1.0, delp/(dp0*qqmx))

  tenm = 10.0*onem
  onemm = 0.001*onem
  thkbop = thkbot*onem

! --- try to restore isopycnic conditions by moving layer interfaces
! --- qhrlx(k) are relaxation coefficients (inverse baroclinic time steps)

  if (fixlay >= 1) then
! ---   maintain constant thickness, layer k = 1
    k = 1
    p_hat = p_i_j(k)+dp0ij(k)
    p_i_j(k+1) = p_hat
    do k=2,kdm
      if     (p_i_j(k+1) >= p_hat) then
        exit  ! usually get here quickly
      endif
      p_i_j(k+1) = p_hat
    enddo !k
  endif

  do k=2,nhybrd
!diag   if (test_col) then
!diag     write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
!diag     write(mesg,'(i9,2i5,a,a)') nstep, CS%itest, CS%jtest, &
!diag       cinfo,':     othkns    odpth    nthkns    ndpth'
!diag     call MOM_mesg(mesg, all_print=.true.)
!diag     do ka=1,kdm
!diag       if     (pres(ka+1) == p_i_j(ka+1) .and. &
!diag               pres(ka  ) == p_i_j(ka  )      ) then
!diag         write(mesg,'(i9,8x,a,a,i3,f10.3,f9.3)') &
!diag              nstep,cinfo,':',ka, &
!diag             (pres(ka+1) -pres(ka))*GV%H_to_m, pres(ka+1)*GV%H_to_m
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag       else
!diag         write(mesg,'(i9,8x,a,a,i3,f10.3,f9.3,f10.3,f9.3)') &
!diag             nstep,cinfo,':',ka, &
!diag             (pres(ka+1) - pres(ka))*GV%H_to_m, pres(ka+1)*GV%H_to_m, &
!diag             (p_i_j(ka+1) - p_i_j(ka))*GV%H_to_m, p_i_j(ka+1)*GV%H_to_m
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag       endif
!diag     enddo !ka
!diag   endif !debug

    if (k <= fixlay) then
! ---     maintain constant thickness, k <= fixlay
      if (k < kdm) then  !p.kdm+1 not changed
        p_i_j(k+1) = min(dp0cum(k+1), p_i_j(kdm+1))
        if (k == fixlay) then
! ---         enforce interface order (may not be necessary).
          do ka=k+2,kdm
            if     (p_i_j(ka) >= p_i_j(k+1)) then
              exit  ! usually get here quickly
            else
              p_i_j(ka) = p_i_j(k+1)
            endif
          enddo !ka
        endif !k == fixlay
      endif !k < kdm

!diag     if (test_col) then
!diag       write(mesg,'(a,i3.2,f8.2)') 'hybgen, fixlay :', &
!diag                                 k+1, p_i_j(k+1)*GV%H_to_m
!diag       call MOM_mesg(mesg, all_print=.true.)
!diag     endif !debug

    else
! ---     do not maintain constant thickness, k > fixlay

      if ((th3d_i_j(k) > theta_i_j(k)+epsil) .and.  (k > fixlay+1)) then
! ---       water in layer k is too dense
! ---       try to dilute with water from layer k-1
! ---       do not move interface if k = fixlay + 1

        if (th3d_i_j(k-1) >= theta_i_j(k-1) .or. &
            p_i_j(k) <= dp0cum(k)+onem .or. &
            p_i_j(k+1)-p_i_j(k) <= p_i_j(k)-p_i_j(k-1)) then
! ---         if layer k-1 is too light, thicken the thinner of the two,
! ---         i.e. skip this layer if it is thicker.

!diag         if (test_col) then
!diag           write(mesg,'(a,3x,i2.2,1pe13.5)') &
!diag                 'hybgen, too dense:',k,th3d_i_j(k)-theta_i_j(k)
!diag           call MOM_mesg(mesg, all_print=.true.)
!diag         endif !debug

          if     ((theta_i_j(k)-th3d_i_j(k-1)) <= epsil) then
!               layer k-1 much too dense, take entire layer
            p_hat = p_i_j(k-1)+dp0ij(k-1)
          else
            q = (theta_i_j(k)-th3d_i_j(k)) / &
                (theta_i_j(k)-th3d_i_j(k-1))         ! -1 <= q < 0
            p_hat0 = p_i_j(k) + q*(p_i_j(k+1)-p_i_j(k))  ! <p_i_j(k)
            if     (k == fixlay+2) then
! ---             treat layer k-1 as fixed.
              p_hat = p_i_j(k-1) +  max(p_hat0-p_i_j(k-1), dp0ij(k-1))
            else
! ---             maintain minimum thickess of layer k-1.
              p_hat = p_i_j(k-1) + cushn(p_hat0-p_i_j(k-1), dp0ij(k-1))
            endif !fixlay+2:else
          end if
          p_hat = min(p_hat, p_i_j(kdm+1))

! ---         if isopycnic conditions cannot be achieved because of a blocking
! ---         layer in the interior ocean, move interface k-1 (and k-2 if
! ---         necessary) upward
          if (k <= fixlay+2) then
! ---           do nothing.
          else if (p_hat >= p_i_j(k) .and. &
                   p_i_j(k-1) > dp0cum(k-1)+tenm .and. &
                  (p_i_j(kdm+1)-p_i_j(k-1) < thkbop .or. &
                   p_i_j(k-1)  -p_i_j(k-2) > qqmx*dp0ij(k-2))) then ! k > 2
            if (k == fixlay+3) then
! ---             treat layer k-2 as fixed.
              p_hat2 = p_i_j(k-2) + &
                       max(p_i_j(k-1)-p_hat+p_hat0-p_i_j(k-2), dp0ij(k-2))
            else
! ---             maintain minimum thickness of layer k-2.
              p_hat2 = p_i_j(k-2) + &
                       cushn(p_i_j(k-1)-p_hat+p_hat0-p_i_j(k-2), dp0ij(k-2))
            endif !fixlay+3:else
            if (p_hat2 < p_i_j(k-1)-onemm) then
              p_i_j(k-1) = (1.0-qhrlx(k-1))*p_i_j(k-1) + &
                           qhrlx(k-1) * max(p_hat2, 2.0*p_i_j(k-1)-p_hat)
!diag             if (test_col) then
!diag               write(mesg,'(a,i3.2,f8.2)') 'hybgen,  1blocking :', &
!diag                     k-1, p_i_j(k-1)*GV%H_to_m
!diag               call MOM_mesg(mesg, all_print=.true.)
!diag             endif !debug
              p_hat = p_i_j(k-1) + cushn(p_hat0-p_i_j(k-1), dp0ij(k-1))
            elseif (k <= fixlay+3) then
! ---             do nothing.
            elseif (p_i_j(k-2) > dp0cum(k-2)+tenm .and. &
                   (p_i_j(kdm+1)-p_i_j(k-2) < thkbop .or. &
                    p_i_j(k-2)  -p_i_j(k-3) > qqmx*dp0ij(k-3))) then
              if (k == fixlay+4) then
! ---               treat layer k-3 as fixed.
                p_hat3 = p_i_j(k-3) + max(p_i_j(k-2)-p_hat + p_hat0-p_i_j(k-3), dp0ij(k-3))
              else
! ---               maintain minimum thickess of layer k-3.
                p_hat3 = p_i_j(k-3) + cushn(p_i_j(k-2)-p_hat+ p_hat0-p_i_j(k-3), dp0ij(k-3))
              endif !fixlay+4:else
              if (p_hat3 < p_i_j(k-2)-onemm) then
                p_i_j(k-2) = (1.0-qhrlx(k-2))*p_i_j(k-2) + &
                             qhrlx(k-2)*max(p_hat3, 2.0*p_i_j(k-2)-p_i_j(k-1))
!diag               if (test_col) then
!diag                 write(mesg,'(a,i3.2,f8.2)') 'hybgen,  2blocking :', &
!diag                       k-2, p_i_j(k-2)*GV%H_to_m
!diag                 call MOM_mesg(mesg, all_print=.true.)
!diag               endif !debug
                p_hat2 = p_i_j(k-2) + cushn(p_i_j(k-1)-p_hat + p_hat0-p_i_j(k-2), dp0ij(k-2))
                if (p_hat2 < p_i_j(k-1)-onemm) then
                  p_i_j(k-1) = (1.0-qhrlx(k-1)) * p_i_j(k-1) + &
                               qhrlx(k-1) * max(p_hat2, 2.0*p_i_j(k-1)-p_hat)
!diag                 if (test_col) then
!diag                   write(mesg,'(a,i3.2,f8.2)') 'hybgen,  3blocking :', &
!diag                              k-1, p_i_j(k-1)*GV%H_to_m
!diag                   call MOM_mesg(mesg, all_print=.true.)
!diag                 endif !debug
                  p_hat = p_i_j(k-1) + cushn(p_hat0-p_i_j(k-1), dp0ij(k-1))
                endif !p_hat2
              endif !p_hat3
            endif !p_hat2:blocking
          endif !blocking
!
          if (p_hat < p_i_j(k)) then
! ---           entrain layer k-1 water into layer k, move interface up.
            p_i_j(k) = (1.0-qhrlx(k)) * p_i_j(k) + qhrlx(k) * p_hat
          endif !entrain
!
!diag         if (test_col) then
!diag           write(mesg,'(a,i3.2,f8.2)') 'hybgen, entrain(k) :', &
!diag                                     k, p_i_j(k)*GV%H_to_m
!diag           call MOM_mesg(mesg, all_print=.true.)
!diag         endif !debug
!
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
              p_i_j(k+1)-p_i_j(k) < p_i_j(k+2)-p_i_j(k+1)) then
! ---           if layer k+1 is too dense, thicken the thinner of the
! ---           two, i.e. skip this layer (never get here) if it is not
! ---           thinner than the other.

!diag           if (test_col) then
!diag             write(mesg,'(a,3x,i2.2,1pe13.5)') &
!diag                  'hybgen, too light:', k, theta_i_j(k)-th3d_i_j(k)
!diag             call MOM_mesg(mesg, all_print=.true.)
!diag           endif !debug

            if     ((th3d_i_j(k+1)-theta_i_j(k)) <= epsil) then
!                 layer k-1 too light, take entire layer
              p_hat = p_i_j(k+2)
            else
              q = (th3d_i_j(k)  -theta_i_j(k)) / &
                  (th3d_i_j(k+1)-theta_i_j(k))          !-1 <= q < 0
              p_hat = p_i_j(k+1) + q*(p_i_j(k)-p_i_j(k+1))  !>p_i_j(k+1)
            endif
!
! ---           if layer k+1, or layer k+2, does not touch the bottom
! ---           then maintain minimum thicknesses of layers k and k+1 as
! ---           much as possible. otherwise, permit layers to collapse
! ---           to zero thickness at the bottom.
!
            if (p_i_j(min(k+3,kdm+1)) < p_i_j(kdm+1)) then
              if (p_i_j(kdm+1)-p_i_j(k) >  dp0ij(k)+dp0ij(k+1)) then
                p_hat = p_i_j(k+2) - cushn(p_i_j(k+2)-p_hat, dp0ij(k+1))
              endif
              p_hat = p_i_j(k) + max(p_hat-p_i_j(k), dp0ij(k))
              p_hat = min(p_hat, max(0.5*(p_i_j(k+1)+p_i_j(k+2)), &
                                     p_i_j(k+2)-dp0ij(k+1)) )
            else
              p_hat = min(p_i_j(k+2), p_hat)
            endif !p.k+2<p.kdm+1
            if (p_hat > p_i_j(k+1)) then
! ---             entrain layer k+1 water into layer k.
              p_i_j(k+1) = (1.0-qhrlx(k+1)) * p_i_j(k+1) + qhrlx(k+1) * p_hat
            endif !entrain

!diag           if (test_col) then
!diag             write(mesg,'(a,i3.2,f8.2)') &
!diag                  'hybgen, entrain(k+):',k, p_i_j(k+1)*GV%H_to_m
!diag             call MOM_mesg(mesg, all_print=.true.)
!diag           endif !debug

          endif !too-light adjustment
        endif !above bottom
      endif !too dense or too light
!
! ---     if layer above is still too thin, move interface down.
      p_hat0 = min(p_i_j(k-1)+dp0ij(k-1), p_i_j(kdm+1))
      if (p_hat0 > p_i_j(k)) then
        p_hat = (1.0-qhrlx(k-1)) * p_i_j(k) + qhrlx(k-1) * p_hat0
        p_i_j(k) = min(p_hat, p_i_j(k+1))
!
!diag       if (test_col) then
!diag         write(mesg,'(a,i3.2,f8.2)') &
!diag              'hybgen, min. thknss (k+):',k-1, p_i_j(k)*GV%H_to_m
!diag         call MOM_mesg(mesg, all_print=.true.)
!diag       endif !debug
      endif

    endif !k <= fixlay:else

  enddo !k  vertical coordinate relocation

end subroutine hybgenaij_regrid

end module MOM_hybgen_regrid

!
!> HYCOM hybgen Revision history:
!>
!> Feb. 2000 -- total rewrite to convert to 'newzp' approach
!> Jul. 2000 -- added hybgenj for OpenMP parallelization
!> Oct. 2000 -- added hybgenbj to simplify OpenMP logic
!> Nov. 2000 -- fill massless layers on sea floor with salinity from above
!> Nov. 2000 -- unmixing of deepest inflated layer uses th&T&S from above
!> Nov. 2000 -- ignored isopycnic variance is now 0.002
!> Nov. 2000 -- iterate to correct for cabbeling
!> Nov. 2000 -- allow for "blocking" interior layers
!> Nov. 2000 -- hybflg selects conserved fields (any two of T/S/th)
!> Nov. 2002 -- replace PCM remapping with PLM when non-isopycnal
!> Apr. 2003 -- added dp00i for thinner minimum layers away from the surface
!> Dec. 2003 -- fixed tracer bug when deepest inflated layer is too light
!> Dec. 2003 -- improved water column conservation
!> Dec. 2003 -- compile time option for explicit water column conservation
!> Dec. 2003 -- ignored isopycnic variance is now 0.0001
!> Jan. 2004 -- shifted qqmn, qqmx range now used in cushion function
!> Mar. 2004 -- minimum thickness no longer enforced in near-bottom layers
!> Mar. 2004 -- ignored isopycnic variance is now epsil (i.e. very small)
!> Mar. 2004 -- relaxation to isopycnic layers controled via hybrlx
!> Mar. 2004 -- relaxation removes the need to correct for cabbeling
!> Mar. 2004 -- modified unmixing selection criteria
!> Mar. 2004 -- added isotop (topiso) for isopycnal layer minimum depths
!> Jun. 2005 -- hybrlx (qhybrlx) now input via blkdat.input
!> Jan. 2007 -- hybrlx now only active below "fixed coordinate" surface layers
!> Aug. 2007 -- removed mxlkta logic
!> Sep. 2007 -- added hybmap and hybiso for PCM,PLM,PPM remaper selection
!> Jan. 2008 -- updated logic for two layers (one too dense, other too light)
!> Jul. 2008 -- Added WENO-like, and bugfix to PPM for lcm==.true.
!> Aug. 2008 -- Use WENO-like (vs PPM) for most velocity remapping
!> Aug. 2008 -- Switch more thin near-isopycnal layers to PCM remapping
!> May  2009 -- New action when deepest inflated layer is very light
!> Oct  2010 -- updated test for deepest inflated layer too light
!> Oct  2010 -- updated sanity check on deltm
!> July 2010 -- Maintain vertical minval and maxval in remapping
!> Aug  2011 -- Option to apply Robert-Asselin filter to hybgen's updated dp
!> Mar  2012 -- Replaced dssk with dpns and dsns, see blkdat.F for info
!> July 2012 -- Bugfix for tracer in too-light deepest inflated layer
!> Sep  2012 -- Added ndebug_tracer for tracer debuging
!> Sep  2012 -- Don't unmix in terrain-following regime
!> Apr. 2013 -- Detect all fixed coordinate layers
!> Apr. 2013 -- bugfix for constant thickness layer k = 1
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Feb. 2015 - when too light, include layer k+2 in "near bottom" logic
!> Apr. 2015 - fixed a k+3 for k=kk-1 (k+3=kk+2) bug
!> Aug. 2015 - changed k-2 to k-N in lowest layer "runaway" unmixing test
!> Aug. 2015 - overturn with the layer above if bottom layer is very light
!> Aug. 2015 - entrain into too dense layer (only move upper interface up)
!> Aug. 2015 - allow entrainment to increase fixlay by 1
!> Nov. 2019 - avoid overflow in calculation of qdep
!> May  2021 - removed unneeded dpmixl halo update
