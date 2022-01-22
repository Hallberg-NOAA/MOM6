!> This module contains the hybgen unmixing routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_unmix

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,             only : EOS_type, calculate_density, calculate_density_derivs
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_tracer_registry, only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_unmix_CS ; private

!  !> Number of layers
!  integer :: nk

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
  !> Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real :: dp00i
  !> Hybgen relaxation coefficient (inverse baroclinic time steps) [s-1]
  real :: qhybrlx

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

  !> Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real :: topiso_const
  ! real, dimension(:,:), allocatable :: topiso

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  real :: onem       !< Nominally one m in thickness units [H ~> m or kg m-2]

  logical :: debug   !< If true, write verbose checksums for debugging purposes.

end type hybgen_unmix_CS

public hybgen_unmix, init_hybgen_unmix, end_hybgen_unmix
public set_hybgen_unmix_params

contains

!> Initialise a hybgen_unmix_CS control structure and store its parameters
subroutine init_hybgen_unmix(CS, GV, US, param_file)
  type(hybgen_unmix_CS),   pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file

  character(len=40)               :: mdl = "MOM_hybgen" ! This module's name.
  integer :: k

  if (associated(CS)) call MOM_error(FATAL, "init_hybgen_unmix: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(GV%ke))

  allocate(CS%dp0k(GV%ke), source=0.0) ! minimum deep z-layer separation
  allocate(CS%ds0k(GV%ke), source=0.0) ! minimum shallow z-layer separation

!  CS%nk = GV%ke
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
  call get_param(param_file, mdl, "HYBGEN_REMAP_DENSITY_MATCH", CS%hybiso, &
                 "A tolerance between the layer densities and their target, within which "//&
                 "Hybgen determines that remapping uses PCM for a layer.", &
                 units="kg m-3", default=0.0, scale=US%kg_m3_to_R)
  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  CS%onem = 1.0 * GV%m_to_H

  do k=1,GV%ke ;  CS%target_density(k) = GV%Rlay(k) ; enddo

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
!   integer,                 intent(in)    :: j   !< The j-slice to work on

!
! --- --------------------------------------------
! --- hybrid grid generator, single j-row (part A).
! --- --------------------------------------------
!
!

  integer :: fixlay            !deepest fixed coordinate layer
  real :: qdep              !fraction not terrain following
  real :: qhrlx( GV%ke+1)   !relaxation coefficient
  real :: dp0ij( GV%ke)     !minimum layer thickness
  real :: dp0cum(GV%ke+1)   !minimum interface depth
  real :: pres(GV%ke+1)     !original layer interfaces
!
  real :: theta_i_j(GV%ke)   ! Target potential density [R ~> kg m-3]
  real ::  temp_i_j(GV%ke)   !  temp(i,j,:) potential temperature [degC]
  real ::  saln_i_j(GV%ke)   !  saln(i,j,:) salinity [ppt]
  real ::  th3d_i_j(GV%ke)   !  Coordinate potential density
  real ::    dp_i_j(GV%ke)   !    dp(i,j,:,n) layer thicknesses
  real ::     p_i_j(GV%ke+1) !     p(i,j,:)   interface depths
  real :: p_col(GV%ke)       ! A column of reference pressures [R L2 T-2 ~> Pa]
  real :: tracer_i_j(GV%ke,max(ntracer,1))  !  Columns of each tracer [Conc]
  real :: depths_i_j         ! Bottom depth in thickness units [H ~> m or kg m-2]
  real :: onemm ! one mm in thickness units [H ~> m or kg m-2]
  logical :: test_col    ! If true, write verbose debugging for this column.
  integer :: trcflg(max(ntracer,1))  ! Hycom tracer type flag for each tracer
  integer :: i, j, k, kdm, m
!
  kdm = GV%ke
  onemm = 0.001*CS%onem

  ! Set all tracers to be passive.  Setting this to 2 treats a tracer like temperature.
  trcflg(:) = 3

          ! theta_i_j(k) = theta(i,j,k)

  p_col(:) = CS%ref_pressure

  test_col = .false.
  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then
!diag   test_col = ((i == CS%itest + G%isd - G%isd_global) .and. &
!diag               (j == CS%jtest + G%jsd - G%jsd_global))

    !### p_i_j(1) = 0.0
    do k=1,kdm
      ! theta_i_j(k) = theta(i,j,k)  ! If a 3-d target density were set up in theta, use that here.
      theta_i_j(k) = CS%target_density(k)  ! MOM6 does not yet support 3-d target densities.
      dp_i_j(k) = dp(i,j,k)
      !### p_i_j(K+1) = p_i_j(K) + dp_i_j(k)
      temp_i_j(k) = tv%T(i,j,k)
      saln_i_j(k) = tv%S(i,j,k)
    enddo
    depths_i_j = GV%Z_to_H * G%bathyT(i,j)
    !### depths_i_j = p_i_j(kdm+1)

    ! This sets the potential density from T and S.
    call calculate_density(temp_i_j, saln_i_j, p_col, th3d_i_j, tv%eqn_of_state)

    do m=1,ntracer ; do k=1,kdm
      tracer_i_j(k,m) = Reg%Tr(m)%t(i,j,k)
    enddo ; enddo

    call hybgenaij_init(   kdm, CS%nhybrid, CS%nsigma, &
                           CS%dp0k, CS%ds0k, CS%dp00i, CS%topiso_const, CS%qhybrlx, &
                           CS%dpns, CS%dsns, &
                           depths_i_j, dp_i_j, &
                           fixlay, qdep, qhrlx, dp0ij, dp0cum, &
                           p_i_j, test_col)
    call hybgenaij_unmix(  CS, kdm, theta_i_j, &
                           temp_i_j, saln_i_j, th3d_i_j, tv%eqn_of_state, &
                           ntracer, tracer_i_j, trcflg, fixlay, qdep, qhrlx, &
                           dp_i_j, onemm, 1.0e-11*US%kg_m3_to_R, test_col)

    ! Store the output from hybgen_unmix
    do k=1,kdm
      dp(i,j,k) = dp_i_j(k)
    enddo
    ! Note that temperature and salinity are among the tracers unmixed here.
    do m=1,ntracer ; do k=1,kdm
      Reg%Tr(m)%t(i,j,k) = tracer_i_j(k,m)
    enddo ; enddo
  endif ; enddo ; enddo !i & j.
end subroutine hybgen_unmix


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

!> Unmix the properties in the lowest layer if it is too light.
subroutine hybgenaij_unmix(CS, kdm, theta_i_j, temp_i_j, saln_i_j, th3d_i_j, eqn_of_state, &
                           ntracr, trac_i_j, trcflg, fixlay, qdep, qhrlx, &
                           dp_i_j, onemm, epsil, test_col)
  type(hybgen_unmix_CS), intent(in) :: CS  !< hybgen unmixing control structure
  integer,        intent(in)    :: kdm           !< The number of layers
  integer,        intent(in)    :: fixlay        !< deepest fixed coordinate layer
  real,           intent(in)    :: qdep          !< fraction dp0k (vs ds0k)
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
  real,           intent(in)    :: onemm         !< one mm in pressure units
  real,           intent(in)    :: epsil         !< small nonzero density difference to prevent
                                                 !! division by zero [R ~> kg m-3]
  logical,        intent(in)    :: test_col      !< If true, write verbose debugging for this column.
!
! --- ------------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - ummix lowest layer.
! --- ------------------------------------------------------------------
!
      logical, parameter :: lunmix=.true.     !unmix a too light deepest layer
      integer, parameter :: ndebug_tracer=0   !tracer to debug, usually 0 (off)
!
  integer :: k, ka, kk, kp, ktr, fixall
  character(len=256) :: mesg  ! A string for output messages
  real :: p_hat, dpthin
  real :: delt, deltm, dels, delsm, q, qtr, qts
  real :: s1d(kdm, ntracr+4) !original scalar fields

  kk = kdm
  dpthin = 1e-6*CS%onem

  ! --- identify the deepest layer kp with significant thickness (> dpthin)
  kp = 2  !minimum allowed value
  do k=kk,3,-1
    if (dp_i_j(k) >= dpthin) then
      kp = k
      exit
    endif
  enddo !k

!diag      if (test_col) then
!diag        write(mesg,'(a,i3)') &
!diag              'hybgen, deepest inflated layer:',kp
!diag        call MOM_mesg(mesg, all_print=.true.)
!diag      endif !debug

  k  = kp  !at least 2
  ka = max(k-2,1)  !k might be 2
!
  if ( ((k > fixlay+1) .and. (qdep == 1.0)) .and. & ! layer not fixed depth
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

          if (ndebug_tracer > 0 .and. ndebug_tracer <= ntracr .and. test_col) then
            ktr = ndebug_tracer
            write(mesg,'(a,i3,f6.3,f9.4)') &
              'hybgen, 11(+):', k-1, 0.0, trac_i_j(k-1,ktr)
            call MOM_mesg(mesg, all_print=.true.)
          endif !debug_tracer
    do ktr= 1,ntracr
      trac_i_j(k-1,ktr) = trac_i_j(k-1,ktr) - &
                          q*(trac_i_j(k-1,ktr) - trac_i_j(k,ktr) )
    enddo !ktr
          if (ndebug_tracer > 0 .and. ndebug_tracer <= ntracr .and. test_col) then
            ktr = ndebug_tracer
            write(mesg,'(a,i3,f6.3,f9.4)') &
              'hybgen, 11(+):', k-1, q, trac_i_j(k-1,ktr)
            call MOM_mesg(mesg, all_print=.true.)
            write(mesg,'(a,i3,f6.3,f9.4)') &
              'hybgen, 11(+):', k, q, trac_i_j(k,ktr)
            call MOM_mesg(mesg, all_print=.true.)
            write(mesg,'(a,i3)') &
                  'hybgen, deepest inflated layer:', kp
            call MOM_mesg(mesg, all_print=.true.)
          endif !debug_tracer
! ---   entrained the entire layer into the one above, so now kp=kp-1
    dp_i_j(k-1) = dp_i_j(k-1) + dp_i_j(k)
    dp_i_j(k) = 0.0
    kp = k-1
!diag        if (test_col) then
!diag          write(mesg,'(a,i3,f6.3,5f8.3)') &
!diag              'hybgen, 11(+):', k-1, q, temp_i_j(k-1),saln_i_j(k-1), &
!diag              th3d_i_j(k-1)+CS%thbase,theta_i_j(k-1)+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag          write(mesg,'(a,i3)') &
!diag                'hybgen, deepest inflated layer:',kp
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        endif !debug
  elseif ( ((k > fixlay+1) .and. (qdep == 1.0)) .and. & ! layer not fixed depth
           (dp_i_j(k-1) >= dpthin)              .and. & ! layer above not too thin
           (theta_i_j(k)-epsil > th3d_i_j(k))   .and. & ! layer is lighter than its target
           (th3d_i_j(k-1) > th3d_i_j(k)) ) then
!
! ---   water in the deepest inflated layer with significant thickness
! ---   (kp) is too light, and it is lighter than the layer above, but not the layer two above.
! ---
! ---   swap the entire layer with the one above.
!diag        if (test_col) then
!diag          write(mesg,'(a,i3,f8.5,5f10.5)') &
!diag            'hybgen, original:', &
!diag            k-1,0.0,temp_i_j(k-1),saln_i_j(k-1), &
!diag                th3d_i_j(k-1)+CS%thbase,theta_i_j(k-1)+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag          write(mesg,'(a,i3,f8.5,5f10.5)') &
!diag            'hybgen, original:', &
!diag            k,0.0,temp_i_j(k),saln_i_j(k), &
!diag                th3d_i_j(k)+CS%thbase,theta_i_j(k  )+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        endif !debug
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
!diag        if (test_col) then
!diag          write(mesg,'(a,i3,f8.5,5f10.5)') &
!diag            'hybgen, overturn:', k-1, q,temp_i_j(k-1),saln_i_j(k-1), &
!diag                th3d_i_j(k-1)+CS%thbase,theta_i_j(k-1)+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag          write(mesg,'(a,i3,f8.5,5f10.5)') &
!diag            'hybgen, overturn:', k,  q,temp_i_j(k),saln_i_j(k), &
!diag                th3d_i_j(k)+CS%thbase,theta_i_j(k  )+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        endif !debug
  endif
!
  k  = kp  !at least 2
  ka = max(k-2,1)  !k might be 2
!
  if ( lunmix .and.  & !usually .true.
       ((k > fixlay+1) .and. (qdep == 1.0)) .and. & ! layer not fixed depth
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
!
!diag   if (test_col) then
!diag     write(mesg,'(a,i3)') &
!diag       'hybgen, deepest inflated layer too light   (stable):',k
!diag     call MOM_mesg(mesg, all_print=.true.)
!diag   endif !debug
!
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
            if (ndebug_tracer > 0 .and. ndebug_tracer <= ntracr .and. &
                test_col) then
              ktr = ndebug_tracer
              write(mesg,'(a,i3,f6.3,f9.4)') &
                'hybgen, 10(+):', k-1, 0.0, trac_i_j(k-1,ktr)
              call MOM_mesg(mesg, all_print=.true.)
            endif !debug_tracer
! ---     fraction of new upper layer from old lower layer
      qtr = p_hat / max(p_hat, dp_i_j(k))  !between 0 and 1
      do ktr= 1,ntracr
        if (trcflg(ktr) == 2) then !temperature tracer
          trac_i_j(k,ktr) = trac_i_j(k,ktr) + &
                            (q/(1.0-q)) * (trac_i_j(k,ktr) - trac_i_j(k-1,ktr))
        else !standard tracer - not split into two sub-layers
          trac_i_j(k-1,ktr) = trac_i_j(k-1,ktr) + &
                              qtr * (trac_i_j(k,ktr) - trac_i_j(k-1,ktr))
!diag              if (test_col) then
!diag                write(mesg,'(a,i4,i3,5e12.3)') &
!diag                    'hybgen, 10(+):', k,ktr, p_hat, &
!diag                    p_i_j(k), p_i_j(k-1), qtr, trac_i_j(k-1,ktr)
!diag                call MOM_mesg(mesg, all_print=.true.)
!diag              endif !debug
        endif !trcflg
      enddo !ktr
            if (ndebug_tracer > 0 .and. ndebug_tracer <= ntracr .and. test_col) then
              ktr = ndebug_tracer
              write(mesg,'(a,i3,f6.3,f9.4)') &
                'hybgen, 10(+):', k-1, qtr, trac_i_j(k-1,ktr)
              call MOM_mesg(mesg, all_print=.true.)
              write(mesg,'(a,i3,f6.3,f9.4)') &
                'hybgen, 10(+):', k, qtr, trac_i_j(k,ktr)
              call MOM_mesg(mesg, all_print=.true.)
              write(mesg,'(a,i3)') &
                    'hybgen, deepest inflated layer:', kp
              call MOM_mesg(mesg, all_print=.true.)
            endif !debug_tracer
    endif !tracers
!diag        if (test_col) then
!diag          write(mesg,'(a,i3,f6.3,5f8.3)') &
!diag              'hybgen, 10(+):', k, q, temp_i_j(k), saln_i_j(k), &
!diag                th3d_i_j(k)+CS%thbase, theta_i_j(k)+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        endif !debug
!diag        if (test_col) then
!diag          write(mesg,'(a,i3,f6.3,5f8.3)') &
!diag              'hybgen, 10(-):', k, 0.0, temp_i_j(k), saln_i_j(k), &
!diag              th3d_i_j(k)+CS%thbase, theta_i_j(k)+CS%thbase
!diag          call MOM_mesg(mesg, all_print=.true.)
!diag        endif !debug
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
    if (ndebug_tracer > 0 .and. ndebug_tracer <= ntracr .and. test_col) then
      ktr = ndebug_tracer
      write(mesg,'(a,i3,f9.4)') 'hybgen, massless:', k, trac_i_j(k,ktr)
      call MOM_mesg(mesg, all_print=.true.)
    endif !debug_tracer
  enddo !k

end subroutine hybgenaij_unmix


!> Determine the potential temperature that is consistent with a given density anomaly,
!! salinity and reference pressure, using Newton's method.
function Tofsig(sigma, T_input, salin, ref_pres, eqn_of_state, sigma_ref)
  real, intent(in) :: sigma     ! The density anomaly [R ~> kg m-3]
  real, intent(in) :: T_input   ! An input first guess at the temperature [degC]
  real, intent(in) :: salin     ! Salinity [ppt]
  real, intent(in) :: ref_pres  ! The reference pressure for sigma [R L2 T-2 ~> Pa]
  type(EOS_type), intent(in) :: eqn_of_state !< Equation of state structure
  real, intent(in) :: sigma_ref  ! The difference between sigma and density [R ~> kg m-3]
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
  real, intent(in) :: sigma     ! The density anomaly [R ~> kg m-3]
  real, intent(in) :: Temp      ! Potential temperature [degC]
  real, intent(in) :: S_input   ! An input first guess at the salinity [ppt]
  real, intent(in) :: ref_pres  ! The reference pressure for sigma [R L2 T-2 ~> Pa]
  type(EOS_type), intent(in) :: eqn_of_state !< Equation of state structure
  real, intent(in) :: sigma_ref  ! The difference between sigma and density [R ~> kg m-3]
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
