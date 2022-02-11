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

  !> Hybgen remapping flag (0=PCM,1=PLM,2=PPM,3=WENO-like), by default 3
  integer :: hybmap

  !> Hybgen velocity remapping flag (0=PCM,1=PLM,3=WENO-like), by default 3
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

public init_hybgen_remap, hybgen_remap_column, get_hybgen_remap_params

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

!> This subroutine can be used to retrieve the parameters for the hybgen_remap module
subroutine get_hybgen_remap_params(CS, hybmap, hybmap_vel)
  type(hybgen_remap_CS), intent(in)  :: CS     !< hybgen remapping control structure
  integer,     optional, intent(out) :: hybmap !< A flag indicating the scheme to use for
                                               !! hybgen remapping of tracers
  integer,     optional, intent(out) :: hybmap_vel  !< A flag indicating the scheme to use
                                               !! for hybgen remapping of velocities
  if (present(hybmap))     hybmap     = CS%hybmap
  if (present(hybmap_vel)) hybmap_vel = CS%hybmap_vel
end subroutine get_hybgen_remap_params


!> Vertically remap a column of scalars (such as tracers or velocity components) to the new grid
subroutine hybgen_remap_column(remap_scheme, nfld, nks, h_src, fld_src, nkt, h_tgt, fld_tgt, dpthin, PCM_cell)
  integer,           intent(in)    :: remap_scheme !< A coded integer indicating the remapping scheme to use
  integer,           intent(in)    :: nfld         !< The number of fields to remap
  integer,           intent(in)    :: nks          !< Number of layers in the source column
  real,              intent(in)    :: h_src(nks)   !< Source grid layer thicknesses [H ~> m or kg m-2]
  real,              intent(in)    :: fld_src(nks,nfld) !< Columns of the fields on the source grid
                                                   !! in arbitrary units [A]
  integer,           intent(in)    :: nkt          !< Number of layers in the target column
  real,              intent(in)    :: h_tgt(nkt)    !< Target grid layer thicknesses [H ~> m or kg m-2]
  real,              intent(out)   :: fld_tgt(nkt,nfld) !< Columns of the fields on the target grid
                                                   !! in the same arbitrary units as fld_src [A]
  real,              intent(in)    :: dpthin       !< A negligibly small thickness for the purpose of cell
                                                   !! reconstructions [H ~> m or kg m-2].
  logical, optional, intent(in)    :: PCM_cell(nks) !< If true for a layer, use PCM remapping for that layer

! --- -------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - remap scalars.
! --- -------------------------------------------------------------
  real :: c1d(nks,nfld,3)   ! Interpolation coefficients

  ! --- remap scalar field profiles from the source vertical grid
  ! --- grid onto the new target vertical grid.
  if     (remap_scheme == REMAPPING_PCM) then !PCM
    call hybgen_pcm_remap(fld_src, h_src, fld_tgt, h_tgt, nks, nkt, nfld, dpthin)
  elseif (remap_scheme == REMAPPING_PLM) then !PLM (as in 2.1.08)
    call hybgen_plm_coefs(fld_src, h_src, c1d, nks, nfld, dpthin, PCM_cell)
    call hybgen_plm_remap(fld_src, h_src, c1d, fld_tgt, h_tgt, nks, nkt, nfld, dpthin)
  elseif (remap_scheme == REMAPPING_PPM_H4) then !PPM
    call hybgen_ppm_coefs(fld_src, h_src, c1d, nks, nfld, dpthin, PCM_cell)
    call hybgen_ppm_remap(fld_src, h_src, c1d, fld_tgt, h_tgt, nks, nkt, nfld, dpthin)
  elseif (remap_scheme == REMAPPING_WENO) then !WENO-like
    call hybgen_weno_coefs(fld_src, h_src, c1d, nks, nfld, dpthin, PCM_cell)
    call hybgen_weno_remap(fld_src, h_src, c1d, fld_tgt, h_tgt, nks, nkt, nfld, dpthin)
  endif

end subroutine hybgen_remap_column


!> Do piecewise constant remapping for a set of scalars
subroutine hybgen_pcm_remap(si, dpi, so, dpo, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: dpo(ko)   !< The output layer thicknesses [H ~> m or kg m-2]
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
!       dpo   - target layer thicknesses (dpo(k) = po(k+1)-po(k))
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
  real  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
  integer :: i, k, l, lb, lt
  real :: dpb, dpt, xb, xt, zb, zt, zx, o
  real*8 :: sz
  real :: si_min(ks), si_max(ks)

! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i

  pi(1) = 0.0
  do k=1,ki ; pi(k+1) = pi(K) + dpi(k) ; enddo
  po(1) = 0.0
  do k=1,ko ; po(k+1) = po(K) + dpo(k) ; enddo

  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while (pi(lb+1) < zb .and. lb < ki)
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
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
subroutine hybgen_plm_coefs(si, dpi, ci, kk, ks, thin, PCM_lay)
  integer, intent(in)  :: kk        !< The number of input layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(kk,ks) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(kk)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: ci(kk,ks) !< The PLM coefficients (slopes) describing the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(kk) !< If true for a layer, use PCM remapping for that layer

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
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
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
    if (dpi(k) <= thin) then  !use PCM
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

  if (present(PCM_lay)) then
    do k=1,kk ; if (PCM_lay(k)) then
      do i=1,ks ; ci(k,i) = 0.0 ; enddo
    endif ; enddo
  endif

end subroutine hybgen_plm_coefs

!> Do piecewise linear remapping for a set of scalars
subroutine hybgen_plm_remap(si, dpi, ci, so, dpo, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: ci(ki,ks) !< The PLM coefficients (slopes) describing the scalar fields [A]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: dpo(ko)   !< The output layer thicknesses [H ~> m or kg m-2]
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
!       dpo   - target layer thicknesses (dpo(k) = po(k+1)-po(k))
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
  real  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
  integer :: i, k, l, lb, lt
  real :: c0, qb0, qb1, qt0, qt1, xb, xt, zb, zt, zx, o
  real*8 :: sz
  real :: si_min(ks),si_max(ks)

! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
!
  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i

  pi(1) = 0.0
  do k=1,ki ; pi(k+1) = pi(K) + dpi(k) ; enddo
  po(1) = 0.0
  do k=1,ko ; po(k+1) = po(K) + dpo(k) ; enddo

  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while ((pi(lb+1) < zb) .and. (lb < ki))
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
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
subroutine hybgen_ppm_coefs(s, h_src, ci, kk, ks, thin, PCM_lay)
  integer, intent(in)  :: kk        !< The number of input layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: s(kk,ks)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(kk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: ci(kk,ks,3) !< The PPM coefficients describing the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(kk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic piecewise parabolic across each input cell
!
!     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
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
  real :: dp(kk) ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  logical :: PCM_layer(kk) ! True for layers that should use PCM remapping, either because they are
                           ! very thin, or because this is specified by PCM_lay.
  real :: da, a6, slj, scj, srj
  real :: as(kk), al(kk), ar(kk)
  real :: dpjp(kk), dp2jp(kk), dpj2p(kk)
  real :: qdpjp(kk), qdp2jp(kk), qdpj2p(kk), dpq3(kk), qdp4(kk)

  ! This PPM remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
  do j=1,kk ; dp(j) = max(h_src(j), thin) ; enddo

  ! Specify the layers that will use PCM remapping.
  if (present(PCM_lay)) then
    do j=1,kk ; PCM_layer(j) = (PCM_lay(j) .or. dp(j) <= thin) ; enddo
  else
    do j=1,kk ; PCM_layer(j) = (dp(j) <= thin) ; enddo
  endif

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
      if (PCM_layer(j)) then  !use PCM
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
    ar(kk-1) = s(kk,i) !last layer PCM
    al(kk)  = s(kk,i)  !last layer PCM
    ar(kk)  = s(kk,i)  !last layer PCM
    !Impose monotonicity: Colella, Eq. (1.10)
    do j=2,kk-1
      if ((PCM_layer(j)) .or. ((s(j+1,i)-s(j,i))*(s(j,i)-s(j-1,i)) <= 0.)) then !local extremum
        al(j) = s(j,i)
        ar(j) = s(j,i)
      else
        da = ar(j)-al(j)
        a6 = 6.0*s(j,i) - 3.0*(al(j)+ar(j))
        if (da*a6 > da*da) then !peak in right half of zone
          al(j) = 3.0*s(j,i) - 2.0*ar(j)
        elseif (da*a6 < -da*da) then !peak in left half of zone
          ar(j) = 3.0*s(j,i) - 2.0*al(j)
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
subroutine hybgen_ppm_remap(si, dpi, ci, so, dpo, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: ci(ki,ks,3) !< The PPM coefficients describing the scalar fields [A]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: dpo(ko)   !< The output layer thicknesses [H ~> m or kg m-2]
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
!       dpo   - target layer thicknesses (dpo(k) = po(k+1)-po(k))
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
  real  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
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

  pi(1) = 0.0
  do k=1,ki ; pi(k+1) = pi(K) + dpi(k) ; enddo
  po(1) = 0.0
  do k=1,ko ; po(k+1) = po(K) + dpo(k) ; enddo

  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while (pi(lb+1) < zb .and. lb < ki)
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
    lt = lb !top will always correspond to bottom of previous
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
subroutine hybgen_weno_coefs(s, h_src, ci, kk, ks, thin, PCM_lay)
  integer, intent(in)  :: kk        !< The number of input layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: s(kk,ks)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(kk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: ci(kk,ks,2) !< The WENO coefficients describing the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(kk) !< If true for a layer, use PCM remapping for that layer

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
!       h_src - initial layer thicknesses (>=0)
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
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
  logical :: PCM_layer(kk) ! True for layers that should use PCM remapping, either because they are
                           ! very thin, or because this is specified by PCM_lay.
  real :: dp(kk) ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  real    qdpjm(kk), qdpjmjp(kk), dpjm2jp(kk)
  real    zw(kk+1,3)

  ! The WENO remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
  do j=1,kk ; dp(j) = max(h_src(j), thin) ; enddo

  ! Specify the layers that will use PCM remapping.
  if (present(PCM_lay)) then
    do j=1,kk ; PCM_layer(j) = (PCM_lay(j) .or. dp(j) <= thin) ; enddo
  else
    do j=1,kk ; PCM_layer(j) = (dp(j) <= thin) ; enddo
  endif

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
      if (PCM_layer(j)) then  !use PCM
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
      if (.not.PCM_layer(j)) then  !don't use PCM
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
subroutine hybgen_weno_remap(si, dpi, ci, so, dpo, ki, ko, ks, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ks        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ks) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: ci(ki,ks,2) !< The WENO coefficients describing the scalar fields [A]
  real,    intent(out) :: so(ko,ks) !< The output scalar fields [A]
  real,    intent(in)  :: dpo(ko)   !< The output layer thicknesses [H ~> m or kg m-2]
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
!       dpo   - target layer thicknesses (dpo(k) = po(k+1)-po(k))
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
  real  :: pi(ki+1)  !< The input interface positions relative to the surface [H ~> m or kg m-2]
  real  :: po(ko+1)  !< The output interface positions relative to the surface [H ~> m or kg m-2]
  integer i,k,l,lb,lt
  real    dpb, dpt, qb0, qb1, qb2, qt0, qt1, qt2,xb,xt,zb,zt,zx,o
  real*8  sz
  real    si_min(ks),si_max(ks)

! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision

  do i=1,ks
    si_min(i) = minval(si(:,i))
    si_max(i) = maxval(si(:,i))
  enddo !i

  pi(1) = 0.0
  do k=1,ki ; pi(k+1) = pi(K) + dpi(k) ; enddo
  po(1) = 0.0
  do k=1,ko ; po(k+1) = po(K) + dpo(k) ; enddo

  zx = pi(ki+1) !maximum depth
  zb = max(po(1), pi(1))
  lb = 1
  do while (pi(lb+1) < zb .and. lb < ki)
    lb = lb+1
  enddo
  do k= 1,ko  !output layers
    zt = zb
    zb = min(po(k+1),zx)
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

end module MOM_hybgen_remap
