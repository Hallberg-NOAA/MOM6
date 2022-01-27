!> This module contains the hybgen unmixing routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_unmix

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,             only : EOS_type, calculate_density, calculate_density_derivs
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_hybgen_regrid,   only : hybgen_column_init
use MOM_hybgen_regrid,   only : hybgen_regrid_CS, get_hybgen_regrid_params
use MOM_tracer_registry, only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_unmix_CS ; private

  integer :: nhybrid !< Number of hybrid levels used by HYBGEN (0=all isopycnal)
  integer :: nsigma  !< Number of sigma levels used by HYBGEN (nhybrid-nsigma z-levels)
  real :: hybiso     !< Hybgen uses PCM if layer is within hybiso of target density [kg m-3]

  real :: dp00i   !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real :: qhybrlx !< Hybgen relaxation coefficient (inverse baroclinic time steps) [s-1]

  real, allocatable, dimension(:) ::  &
    dp0k, &     !< minimum deep    z-layer separation [H ~> m or kg m-2]
    ds0k        !< minimum shallow z-layer separation [H ~> m or kg m-2]

  real :: dpns  !< depth to start terrain following [H ~> m or kg m-2]
  real :: dsns  !< depth to stop terrain following [H ~> m or kg m-2]
  real :: min_dilate !< The minimum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when wetting occurs.
  real :: max_dilate !< The maximum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when drying occurs.

  real :: topiso_const !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  ! real, dimension(:,:), allocatable :: topiso

  real :: ref_pressure !< Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real :: thbase  !< Reference density for anomalies [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density !< Nominal density of interfaces [R ~> kg m-3]

end type hybgen_unmix_CS

public hybgen_unmix, init_hybgen_unmix, end_hybgen_unmix
public set_hybgen_unmix_params

contains

!> Initialise a hybgen_unmix_CS control structure and store its parameters
subroutine init_hybgen_unmix(CS, GV, US, param_file, hybgen_regridCS)
  type(hybgen_unmix_CS),   pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(hybgen_regrid_CS),  pointer    :: hybgen_regridCS !< Control structure for hybgen
                                             !! regridding for sharing parameters.

  character(len=40)               :: mdl = "MOM_hybgen" ! This module's name.
  integer :: k

  if (associated(CS)) call MOM_error(FATAL, "init_hybgen_unmix: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(GV%ke))

  allocate(CS%dp0k(GV%ke), source=0.0) ! minimum deep z-layer separation
  allocate(CS%ds0k(GV%ke), source=0.0) ! minimum shallow z-layer separation

  ! Set the parameters for the hybgen unmixing from a hybgen regridding control structure.
  call get_hybgen_regrid_params(hybgen_regridCS, ref_pressure=CS%ref_pressure, &
                nhybrid=CS%nhybrid, nsigma=CS%nsigma, dp0k=CS%dp0k, ds0k=CS%ds0k, &
                dp00i=CS%dp00i, topiso_const=CS%topiso_const, qhybrlx=CS%qhybrlx, &
                hybiso=CS%hybiso, min_dilate=CS%min_dilate, max_dilate=CS%max_dilate, &
                target_density=CS%target_density)

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

end subroutine init_hybgen_unmix

!> This subroutine deallocates memory in the control structure for the hybgen unmixing module
subroutine end_hybgen_unmix(CS)
  type(hybgen_unmix_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return

  deallocate(CS%target_density)
  deallocate(CS%dp0k, CS%ds0k)
  deallocate(CS)
end subroutine end_hybgen_unmix

!> This subroutine can be used to set the parameters for the hybgen module
subroutine set_hybgen_unmix_params(CS, min_thickness)
  type(hybgen_unmix_CS),  pointer    :: CS !< Coordinate unmixing control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hybgen_params: CS not associated")

!  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_hybgen_unmix_params


!> Unmix the properties in the lowest layer with mass if it is too light, and make
!! any other changes to the water column to prepare for regridding.
subroutine hybgen_unmix(G, GV, US, CS, tv, Reg, ntracer, dp)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(hybgen_unmix_CS),   intent(in)    :: CS  !< hybgen control structure
  type(thermo_var_ptrs),   intent(inout) :: tv  !< Thermodynamics structure
  type(tracer_registry_type), pointer    :: Reg !< Tracer registry structure
  integer,                 intent(in)    :: ntracer !< The number of tracers in the registry, or
                                                !! 0 if the registry is not in use.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: dp  !< Layer thicknesses [H ~> m or kg m-2]

! --- --------------------------------------------
! --- hybrid grid generator, single j-row (part A).
! --- --------------------------------------------

  integer :: fixlay         ! deepest fixed coordinate layer
  real :: qhrlx( GV%ke+1)   ! relaxation coefficient
  real :: dp0ij( GV%ke)     ! minimum layer thickness
  real :: dp0cum(GV%ke+1)   ! minimum interface depth

  real :: theta_i_j(GV%ke)  ! Target potential density [R ~> kg m-3]
  real ::  temp_i_j(GV%ke)  ! A column of potential temperature [degC]
  real ::  saln_i_j(GV%ke)  ! A column of salinity [ppt]
  real ::  th3d_i_j(GV%ke)  ! A column of coordinate potential density
  real ::    dp_i_j(GV%ke)  ! A column of layer thicknesses
  real :: p_col(GV%ke)      ! A column of reference pressures [R L2 T-2 ~> Pa]
  real :: tracer_i_j(GV%ke,max(ntracer,1))  ! Columns of each tracer [Conc]
  real :: h_tot             ! Total thickness of the water column [H ~> m or kg m-2]
  real :: nominalDepth      ! Depth of ocean bottom (positive downward) [H ~> m or kg m-2]
  real :: dpthin            ! A negligibly small thickness to identify essentially
                            ! vanished layers [H ~> m or kg m-2]
  real :: dilate            ! A factor by which to dilate the target positions from z to z* [nondim]
  logical :: terrain_following  ! True if this column is terrain following.
  integer :: trcflg(max(ntracer,1))  ! Hycom tracer type flag for each tracer
  integer :: i, j, k, kdm, m

  kdm = GV%ke

  ! Set all tracers to be passive.  Setting this to 2 treats a tracer like temperature.
  trcflg(:) = 3

  dpthin = 1e-6*GV%m_to_H

  p_col(:) = CS%ref_pressure

  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then

    h_tot = 0.0
    do k=1,kdm
      ! theta_i_j(k) = theta(i,j,k)  ! If a 3-d target density were set up in theta, use that here.
      theta_i_j(k) = CS%target_density(k)  ! MOM6 does not yet support 3-d target densities.
      dp_i_j(k) = dp(i,j,k)
      h_tot = h_tot + dp_i_j(k)
      temp_i_j(k) = tv%T(i,j,k)
      saln_i_j(k) = tv%S(i,j,k)
    enddo

    ! This sets the potential density from T and S.
    call calculate_density(temp_i_j, saln_i_j, p_col, th3d_i_j, tv%eqn_of_state)

    do m=1,ntracer ; do k=1,kdm
      tracer_i_j(k,m) = Reg%Tr(m)%t(i,j,k)
    enddo ; enddo

    ! The following block of code is used to trigger z* stretching of the targets heights.
    nominalDepth = (G%bathyT(i,j)+G%Z_ref)*GV%Z_to_H
    if (h_tot <= CS%min_dilate*nominalDepth) then
      dilate = CS%min_dilate
    elseif (h_tot >= CS%max_dilate*nominalDepth) then
      dilate = CS%max_dilate
    else
      dilate = h_tot / nominalDepth
    endif

    terrain_following = (h_tot < dilate*CS%dpns) .and. (CS%dpns >= CS%dsns)

    ! Convert the regridding parameters into specific constraints for this column.
    call hybgen_column_init(kdm, CS%nhybrid, CS%nsigma, CS%dp0k, CS%ds0k, CS%dp00i, &
                            CS%topiso_const, CS%qhybrlx, CS%dpns, CS%dsns, h_tot, dilate, &
                            dp_i_j, fixlay, qhrlx, dp0ij, dp0cum)

    ! Do any unmixing of the column that is needed to move the layer properties toward their targets.
    call hybgenaij_unmix(CS, kdm, theta_i_j, temp_i_j, saln_i_j, th3d_i_j, tv%eqn_of_state, &
                         ntracer, tracer_i_j, trcflg, fixlay, qhrlx, &
                         dp_i_j, terrain_following, dpthin, 1.0e-11*US%kg_m3_to_R)

    ! Store the output from hybgen_unmix in the 3-d arrays.
    do k=1,kdm
      dp(i,j,k) = dp_i_j(k)
    enddo
    ! Note that temperature and salinity are among the tracers unmixed here.
    do m=1,ntracer ; do k=1,kdm
      Reg%Tr(m)%t(i,j,k) = tracer_i_j(k,m)
    enddo ; enddo
  endif ; enddo ; enddo !i & j.

end subroutine hybgen_unmix


!> Unmix the properties in the lowest layer if it is too light.
subroutine hybgenaij_unmix(CS, kdm, theta_i_j, temp_i_j, saln_i_j, th3d_i_j, eqn_of_state, &
                           ntracr, trac_i_j, trcflg, fixlay, qhrlx, dp_i_j, &
                           terrain_following, dpthin, epsil)
  type(hybgen_unmix_CS), intent(in) :: CS  !< hybgen unmixing control structure
  integer,        intent(in)    :: kdm           !< The number of layers
  integer,        intent(in)    :: fixlay        !< deepest fixed coordinate layer
  real,           intent(in)    :: qhrlx( kdm+1) !< relaxation coefficient [s-1]
  real,           intent(in)    :: theta_i_j(kdm) !< Target density [R ~> kg m-3]
  real,           intent(inout) :: temp_i_j(kdm) !< A column of potential temperature [degC]
  real,           intent(inout) :: saln_i_j(kdm) !< A column of salinity [ppt]
  real,           intent(inout) :: th3d_i_j(kdm) !< Coordinate potential density [R ~> kg m-3]
  type(EOS_type), intent(in)    :: eqn_of_state  !< Equation of state structure
  integer,        intent(in)    :: ntracr        !< The number of registered passive tracers
  real,           intent(inout) :: trac_i_j(kdm, max(ntracr,1)) !< Columns of the passive tracers [Conc]
  integer,        intent(in)    :: trcflg(max(ntracr,1)) !< Hycom tracer type flag for each tracer
  real,           intent(inout) :: dp_i_j(kdm+1) !< Layer thicknesses [H ~> m or kg m-2]
  logical,        intent(in)    :: terrain_following !< True if this column is terrain following
  real,           intent(in)    :: dpthin        !< A negligibly small thickness to identify
                                                 !! essentially vanished layers [H ~> m or kg m-2]
  real,           intent(in)    :: epsil         !< small nonzero density difference to prevent
                                                 !! division by zero [R ~> kg m-3]
!
! --- ------------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - ummix lowest layer.
! --- ------------------------------------------------------------------
!
  character(len=256) :: mesg  ! A string for output messages
  real :: p_hat       ! A portion of a layer to move across an interface [H ~> m or kg m-2]
  real :: delt, deltm ! Temperature differences between successive layers [degC]
  real :: dels, delsm ! Salinity differences between successive layers [ppt]
  real :: q, qtr, qts ! Nondimensional fractions [nondim]
  real :: s1d(kdm, ntracr+4) !original scalar fields
  logical, parameter :: lunmix=.true.     ! unmix a too light deepest layer
  integer :: k, ka, kk, kp, ktr, fixall

  kk = kdm

  ! --- identify the deepest layer kp with significant thickness (> dpthin)
  kp = 2  !minimum allowed value
  do k=kk,3,-1
    if (dp_i_j(k) >= dpthin) then
      kp = k
      exit
    endif
  enddo !k

  k  = kp  !at least 2
  ka = max(k-2,1)  !k might be 2
!
  if ( ((k > fixlay+1) .and. (.not.terrain_following)) .and. & ! layer not fixed depth
       (dp_i_j(k-1) >= dpthin)              .and. & ! layer above not too thin
       (theta_i_j(k)-epsil > th3d_i_j(k))   .and. & ! layer is lighter than its target
       ((th3d_i_j(k-1) > th3d_i_j(k)) .and. (th3d_i_j(ka) > th3d_i_j(k))) ) then
!
! ---   water in the deepest inflated layer with significant thickness
! ---   (kp) is too light, and it is lighter than the two layers above.
! ---
! ---   this should only occur when relaxing or nudging layer thickness
! ---   and is a bug (bad interaction with tsadvc) even in those cases
! ---
! ---   entrain the entire layer into the one above
!---    note the double negative in T=T-q*(T-T'), equiv. to T=T+q*(T'-T)
    q = dp_i_j(k) / (dp_i_j(k) + dp_i_j(k-1))
    temp_i_j(k-1) = temp_i_j(k-1) - q*(temp_i_j(k-1) - temp_i_j(k))
    saln_i_j(k-1) = saln_i_j(k-1) - q*(saln_i_j(k-1) - saln_i_j(k))
    call calculate_density(temp_i_j(k-1), saln_i_j(k-1), CS%ref_pressure, &
                           th3d_i_j(k-1), eqn_of_state, CS%thbase)

    do ktr= 1,ntracr
      trac_i_j(k-1,ktr) = trac_i_j(k-1,ktr) - &
                          q*(trac_i_j(k-1,ktr) - trac_i_j(k,ktr) )
    enddo !ktr
! ---   entrained the entire layer into the one above, so now kp=kp-1
    dp_i_j(k-1) = dp_i_j(k-1) + dp_i_j(k)
    dp_i_j(k) = 0.0
    kp = k-1
  elseif ( ((k > fixlay+1) .and. (.not.terrain_following)) .and. & ! layer not fixed depth
           (dp_i_j(k-1) >= dpthin)              .and. & ! layer above not too thin
           (theta_i_j(k)-epsil > th3d_i_j(k))   .and. & ! layer is lighter than its target
           (th3d_i_j(k-1) > th3d_i_j(k)) ) then
! ---   water in the deepest inflated layer with significant thickness
! ---   (kp) is too light, and it is lighter than the layer above, but not the layer two above.
! ---
! ---   swap the entire layer with the one above.
    if (dp_i_j(k) <= dp_i_j(k-1)) then
! ---     bottom layer is thinner, take entire bottom layer
!---      note the double negative in T=T-q*(T-T'), equiv. to T=T+q*(T'-T)
      s1d(k-1,1) = temp_i_j(k-1)
      s1d(k-1,2) = saln_i_j(k-1)
      s1d(k-1,3) = th3d_i_j(k-1)
      q = dp_i_j(k) / dp_i_j(k-1)  !<=1.0

      temp_i_j(k-1) = temp_i_j(k-1) - q*(temp_i_j(k-1) - temp_i_j(k))
      saln_i_j(k-1) = saln_i_j(k-1) - q*(saln_i_j(k-1) - saln_i_j(k))
      call calculate_density(temp_i_j(k-1), saln_i_j(k-1), CS%ref_pressure, &
                             th3d_i_j(k-1), eqn_of_state, CS%thbase)
        ! th3d_i_j(k-1) = sig(temp_i_j(k-1),saln_i_j(k-1))-CS%thbase

      temp_i_j(k) = s1d(k-1,1)
      saln_i_j(k) = s1d(k-1,2)
      th3d_i_j(k) = s1d(k-1,3)
      do ktr= 1,ntracr
        s1d(k-1,2+ktr)    = trac_i_j(k-1,ktr)
        trac_i_j(k-1,ktr) = trac_i_j(k-1,ktr) - &
                              q * (trac_i_j(k-1,ktr) - trac_i_j(k,ktr))
        trac_i_j(k,  ktr) = s1d(k-1,2+ktr)
      enddo !ktr
    else
! ---     bottom layer is thicker, take entire layer above
      s1d(k,1) = temp_i_j(k)
      s1d(k,2) = saln_i_j(k)
      s1d(k,3) = th3d_i_j(k)
      q = dp_i_j(k-1) / dp_i_j(k)  !<1.0

      temp_i_j(k) = temp_i_j(k) + q*(temp_i_j(k-1) - temp_i_j(k))
      saln_i_j(k) = saln_i_j(k) + q*(saln_i_j(k-1) - saln_i_j(k))
      ! th3d_i_j(k) = sig(temp_i_j(k),saln_i_j(k))-CS%thbase
      call calculate_density(temp_i_j(k), saln_i_j(k), CS%ref_pressure, &
                             th3d_i_j(k), eqn_of_state, CS%thbase)

      temp_i_j(k-1) = s1d(k,1)
      saln_i_j(k-1) = s1d(k,2)
      th3d_i_j(k-1) = s1d(k,3)
      do ktr= 1,ntracr
        s1d(k,2+ktr)      = trac_i_j(k,ktr)
        trac_i_j(k,  ktr) = trac_i_j(k,ktr) + &
                              q * (trac_i_j(k-1,ktr) - trac_i_j(k,ktr))
        trac_i_j(k-1,ktr) = s1d(k,2+ktr)
      enddo !ktr
    endif !bottom too light
  endif
!
  k  = kp  !at least 2
  ka = max(k-2,1)  !k might be 2
!
  if ( lunmix .and.  & !usually .true.
       ((k > fixlay+1) .and. (.not.terrain_following)) .and. & ! layer not fixed depth
       (dp_i_j(k-1) >= dpthin)      .and. & ! layer above not too thin
       (theta_i_j(k)-epsil > th3d_i_j(k))   .and. & ! layer is lighter than its target
       (theta_i_j(k-1) < th3d_i_j(k))       .and. & ! layer is denser than the target above
       (abs(theta_i_j(k-1) - th3d_i_j(k-1)) < CS%hybiso) .and. & ! layer above is near its target
       (th3d_i_j(k) - th3d_i_j(k-1) > 0.001*(theta_i_j(k) - theta_i_j(k-1))) ) then
!
! ---   water in the deepest inflated layer with significant thickness (kp) is too
! ---   light but denser than the layer above, with the layer above near-isopycnal
! ---
! ---   split layer into 2 sublayers, one near the desired density
! ---   and one exactly matching the T&S properties of layer k-1.
! ---   To prevent "runaway" T or S, the result satisfies either
! ---     abs(T.k - T.k-1) <= abs(T.k-N - T.k-1) or
! ---     abs(S.k - S.k-1) <= abs(S.k-N - S.k-1) where
! ---     th3d.k-1 - th3d.k-N is at least theta(k-1) - theta(k-2)
! ---   It is also limited to a 50% change in layer thickness.

    ka = 1
    do ktr=k-2,2,-1
      if ( th3d_i_j(k-1) - th3d_i_j(ktr) >= theta_i_j(k-1) - theta_i_j(k-2) ) then
        ka = ktr  !usually k-2
        exit
      endif
    enddo !ktr
!
    delsm = abs(saln_i_j(ka) - saln_i_j(k-1))
    dels = abs(saln_i_j(k-1) - saln_i_j(k))
    deltm = abs(temp_i_j(ka) - temp_i_j(k-1))
    delt = abs(temp_i_j(k-1) - temp_i_j(k))
! ---   sanity check on deltm and delsm
    q = min(temp_i_j(ka), temp_i_j(k-1), temp_i_j(k))
    if   (q > 6.0) then
      deltm = min( deltm,  6.0*(theta_i_j(k)-theta_i_j(k-1)) )
    else  !(q <= 6.0)
      deltm = min( deltm, 10.0*(theta_i_j(k)-theta_i_j(k-1)) )
    endif
    delsm = min( delsm, 1.3*(theta_i_j(k)-theta_i_j(k-1)) )
    qts = 0.0
    if (delt > 1.0d-11) then  ! The hard-coded limit here is copied over from cb_arrays.F90
      qts = max(qts, (min(deltm, 2.0*delt)-delt)/delt)  ! qts<=1.0
    endif
    if (dels > 1.0d-11) then  ! The hard-coded limit here is copied over from cb_arrays.F90
      qts = max(qts, (min(delsm, 2.0*dels)-dels)/dels)  ! qts<=1.0
    endif
    q = (theta_i_j(k)-th3d_i_j(k)) / (theta_i_j(k)-th3d_i_j(k-1))
    q = min(q, qts/(1.0+qts))  ! upper sublayer <= 50% of total
    q = qhrlx(k)*q
! ---   qhrlx is relaxation coefficient (inverse baroclinic time steps)
    p_hat = q * dp_i_j(k)
    dp_i_j(k-1) = dp_i_j(k-1) + p_hat
    dp_i_j(k) = dp_i_j(k) - p_hat

    temp_i_j(k) = temp_i_j(k) + (q/(1.0-q)) * (temp_i_j(k) - temp_i_j(k-1))
    saln_i_j(k) = saln_i_j(k) + (q/(1.0-q)) * (saln_i_j(k) - saln_i_j(k-1))
    ! th3d_i_j(k) = sig(temp_i_j(k),saln_i_j(k))-CS%thbase
    call calculate_density(temp_i_j(k), saln_i_j(k), CS%ref_pressure, &
                           th3d_i_j(k), eqn_of_state, CS%thbase)

    if ((ntracr > 0) .and. (p_hat /= 0.0)) then
! ---     fraction of new upper layer from old lower layer
      qtr = p_hat / max(p_hat, dp_i_j(k))  !between 0 and 1
      do ktr= 1,ntracr
        if (trcflg(ktr) == 2) then !temperature tracer
          trac_i_j(k,ktr) = trac_i_j(k,ktr) + &
                            (q/(1.0-q)) * (trac_i_j(k,ktr) - trac_i_j(k-1,ktr))
        else !standard tracer - not split into two sub-layers
          trac_i_j(k-1,ktr) = trac_i_j(k-1,ktr) + &
                              qtr * (trac_i_j(k,ktr) - trac_i_j(k-1,ktr))
        endif !trcflg
      enddo !ktr
    endif !tracers
  endif !too light
!
! --- massless or near-massless (thickness < dpthin) layers
!
  do k=kp+1,kk
    if (k <= CS%nhybrid) then
      ! --- fill thin and massless layers on sea floor with fluid from above
      th3d_i_j(k) = th3d_i_j(k-1)
      saln_i_j(k) = saln_i_j(k-1)
      temp_i_j(k) = temp_i_j(k-1)
    elseif (th3d_i_j(k) /= theta_i_j(k)) then
      ! ! This is the code used in Hycom, but TofSig is a much more dangerous function
      ! ! than SofSig because of the nonlinearities of the equation of state.
      ! ! --- fill with saln from above
      ! th3d_i_j(k) = max(theta_i_j(k), th3d_i_j(k-1))
      ! saln_i_j(k) = saln_i_j(k-1)
      ! temp_i_j(k) = Tofsig(th3d_i_j(k), temp_i_j(k), saln_i_j(k), CS%ref_pressure, eqn_of_state, CS%thbase)
      ! saln_i_j(k) = Sofsig(th3d_i_j(k), temp_i_j(k), saln_i_j(k), CS%ref_pressure, eqn_of_state, CS%thbase)

      ! --- fill with temperature from above
      th3d_i_j(k) = max(theta_i_j(k), th3d_i_j(k-1))
      temp_i_j(k) = temp_i_j(k-1)
      saln_i_j(k) = Sofsig(th3d_i_j(k), temp_i_j(k), saln_i_j(k), CS%ref_pressure, eqn_of_state, CS%thbase)
    endif
    do ktr= 1,ntracr
      trac_i_j(k,ktr) = trac_i_j(k-1,ktr)
    enddo !ktr
  enddo !k

end subroutine hybgenaij_unmix


!> Determine the potential temperature that is consistent with a given density anomaly,
!! salinity and reference pressure, using Newton's method.
function Tofsig(sigma, T_input, salin, ref_pres, eqn_of_state, sigma_ref)
  real, intent(in) :: sigma     !< The density anomaly [R ~> kg m-3]
  real, intent(in) :: T_input   !< An input first guess at the temperature [degC]
  real, intent(in) :: salin     !< Salinity [ppt]
  real, intent(in) :: ref_pres  !< The reference pressure for sigma [R L2 T-2 ~> Pa]
  type(EOS_type), intent(in) :: eqn_of_state !< Equation of state structure
  real, intent(in) :: sigma_ref  !< The difference between sigma and density [R ~> kg m-3]
  real :: tofsig !< The potential temperature that is consistent with sigma and salin [degC]

  real :: T_guess   ! A guess at the potential temperature [degC]
  real :: sigma_guess ! The density anomaly with the salinity and guessed temperature [R ~> kg m-3]
  real :: dRho_dT   ! The partial derivative of density with potential temperature [R degC-1 ~> kg m-3 degC-1]
  real :: dRho_dS   ! The partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1]
  integer :: itt, max_itt

  call MOM_error(FATAL, "Tofsig still needs to be written with error bounds.")
  T_guess = T_input
  max_itt = 5
  call calculate_density(T_guess, salin, ref_pres, sigma_guess, eqn_of_state, sigma_ref)
  do itt=1,max_itt
    call calculate_density_derivs(T_guess, salin, ref_pres, dRho_dT, dRho_dS, eqn_of_state)
    if (abs(sigma_guess - sigma) < 1e-15*dRho_dT*T_guess) exit
    T_guess = T_guess + (sigma - sigma_guess) / dRho_dT
    call calculate_density(T_guess, salin, ref_pres, sigma_guess, eqn_of_state, sigma_ref)
  enddo
  tofsig = T_guess
end function Tofsig

!> Determine the salinity that is consistent with a given density anomaly,
!! potential temperature and reference pressure, using Newton's method.
function Sofsig(sigma, Temp, S_input, ref_pres, eqn_of_state, sigma_ref)
  real, intent(in) :: sigma     !< The density anomaly [R ~> kg m-3]
  real, intent(in) :: Temp      !< Potential temperature [degC]
  real, intent(in) :: S_input   !< An input first guess at the salinity [ppt]
  real, intent(in) :: ref_pres  !< The reference pressure for sigma [R L2 T-2 ~> Pa]
  type(EOS_type), intent(in) :: eqn_of_state !< Equation of state structure
  real, intent(in) :: sigma_ref  !< The difference between sigma and density [R ~> kg m-3]
  real :: Sofsig !< The salinity that is consistent with sigma and temp [ppt]

  real :: S_guess   ! A guess at the salinity [ppt]
  real :: sigma_guess ! The density anomaly with the temperature and guessed salinity [R ~> kg m-3]
  real :: dRho_dT   ! The partial derivative of density with potential temperature [R degC-1 ~> kg m-3 degC-1]
  real :: dRho_dS   ! The partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1]
  integer :: itt, max_itt

  ! call MOM_error(FATAL, "Sofsig still needs to be written with error bounds.")
  S_guess = S_input
  max_itt = 5
  call calculate_density(Temp, S_guess, ref_pres, sigma_guess, eqn_of_state, sigma_ref)
  do itt=1,max_itt
    call calculate_density_derivs(Temp, S_guess, ref_pres, dRho_dT, dRho_dS, eqn_of_state)
    if (abs(sigma_guess - sigma) < 1e-15*dRho_dS*S_guess) exit
    S_guess = S_guess + (sigma - sigma_guess) / dRho_dS
    call calculate_density(Temp, S_guess, ref_pres, sigma_guess, eqn_of_state, sigma_ref)
  enddo
  Sofsig = S_guess
end function Sofsig

end module MOM_hybgen_unmix
