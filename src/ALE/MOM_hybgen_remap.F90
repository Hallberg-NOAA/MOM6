!> This module contains the hybgen remapping routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_remap

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, param_file_type, log_param
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,    only : unit_scale_type
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
public hybgen_pcm_remap, hybgen_plm_coefs, hybgen_plm_remap
public hybgen_ppm_coefs, hybgen_weno_coefs, hybgen_ppm_remap

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
subroutine hybgen_remap_column(remap_scheme, nfld, nks, h_src, fld_src, nkt, h_tgt, fld_tgt, h_thin, PCM_cell)
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
  real,              intent(in)    :: h_thin       !< A negligibly small thickness for the purpose of cell
                                                   !! reconstructions [H ~> m or kg m-2].
  logical, optional, intent(in)    :: PCM_cell(nks) !< If true for a layer, use PCM remapping for that layer

! --- -------------------------------------------------------------
! --- hybrid grid generator, single column(part A) - remap scalars.
! --- -------------------------------------------------------------
  real :: edges(nks,2,nfld) ! Interpolation edge values [A]
  real :: slope(nks,nfld)   ! Interpolation slopes per cell width [A]

  ! --- remap scalar field profiles from the source vertical grid
  ! --- grid onto the new target vertical grid.
  if     (remap_scheme == REMAPPING_PCM) then !PCM
    call hybgen_pcm_remap(fld_src, h_src, fld_tgt, h_tgt, nks, nkt, nfld, h_thin)
  elseif (remap_scheme == REMAPPING_PLM) then !PLM (as in 2.1.08)
    call hybgen_plm_coefs(fld_src, h_src, slope, nks, nfld, h_thin, PCM_cell)
    call hybgen_plm_remap(fld_src, h_src, slope, fld_tgt, h_tgt, nks, nkt, nfld, h_thin)
  elseif (remap_scheme == REMAPPING_PPM_H4) then !PPM
    call hybgen_ppm_coefs(fld_src, h_src, edges, nks, nfld, h_thin, PCM_cell)
    call hybgen_ppm_remap(fld_src, h_src, edges, fld_tgt, h_tgt, nks, nkt, nfld, h_thin)
  elseif (remap_scheme == REMAPPING_WENO) then !WENO-like
    call hybgen_weno_coefs(fld_src, h_src, edges, nks, nfld, h_thin, PCM_cell)
    call hybgen_ppm_remap(fld_src, h_src, edges, fld_tgt, h_tgt, nks, nkt, nfld, h_thin)
  endif

end subroutine hybgen_remap_column


!> Do piecewise constant remapping for a set of scalars
subroutine hybgen_pcm_remap(si, dpi, so, dpo, ki, ko, ns, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: si(ki,ns) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: so(ko,ns) !< The output scalar fields [A]
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
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       ki    - number of  input layers
!       ko    - number of output layers
!       ns    - number of fields
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
  real :: si_min(ns), si_max(ns)

! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ns
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
  do i=1,ns
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
        do i=1,ns
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ns
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
        do i=1,ns
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
        do i=1,ns
          so(k,i) = si(lt,i)
        enddo !i
      endif
    endif !thin:std layer
  enddo !k

end subroutine hybgen_pcm_remap

!> Set up the coefficients for PLM remapping of a set of scalars
subroutine hybgen_plm_coefs(si, dpi, slope, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: si(nk,ns) !< The cell-averaged input scalar fields [A]
  real,    intent(in)  :: dpi(nk)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: slope(nk,ns) !< The PLM slope times cell width [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

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
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       slope - coefficients for hybgen_plm_remap
!                profile(y) = si+slope*(y-1),  -0.5 <= y <= 0.5
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  real :: qcen   ! A layer's thickness divided by the distance between the centers
                 ! of the adjacent cells, usually ~0.5, but always <= 1 [nondim]
  real :: zbot, zcen, ztop ! Tracer slopes times the layer thickness [A]
  integer :: i, k

  do i=1,ns
    slope(1, i) = 0.0
    slope(nk,i) = 0.0
  enddo !i
  do k= 2,nk-1
    if (dpi(k) <= thin) then  !use PCM
      do i=1,ns ; slope(k,i) = 0.0 ; enddo
    else
! ---     use qcen in place of 0.5 to allow for non-uniform grid
      qcen = dpi(k) / (dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))  !dpi(k)>thin
      do i=1,ns
! ---       PLM (non-zero slope, but no new extrema)
! ---       layer value is si-0.5*slope at top    interface,
! ---                  and si+0.5*slope at bottom interface.
!
! ---       monotonized central-difference limiter (van Leer, 1977,
! ---       JCP 23 pp 276-299).  For a discussion of PLM limiters, see
! ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
        ztop = 2.0*(si(k,  i)-si(k-1,i))
        zbot = 2.0*(si(k+1,i)-si(k,  i))
        zcen = qcen*(si(k+1,i)-si(k-1,i))
        if     (ztop*zbot > 0.0) then !ztop,zbot are the same sign
          slope(k,i) = sign(min(abs(zcen),abs(zbot),abs(ztop)), zbot)
        else
          slope(k,i) = 0.0  !local extrema, so no slope
        endif
      enddo !i
    endif  !PCM:PLM
  enddo !k

  if (present(PCM_lay)) then
    do k=1,nk ; if (PCM_lay(k)) then
      do i=1,ns ; slope(k,i) = 0.0 ; enddo
    endif ; enddo
  endif

end subroutine hybgen_plm_coefs

!> Do piecewise linear remapping for a set of scalars
subroutine hybgen_plm_remap(si, dpi, slope, so, dpo, ki, ko, ns, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: si(ki,ns) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: slope(ki,ns) !< The PLM slope times cell width [A]
  real,    intent(out) :: so(ko,ns) !< The output scalar fields [A]
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
!       slope - slopes from hybgen_plm_coefs
!                profile(y) = si+slope*(y-1),  -0.5 <= y <= 0.5
!       ki    - number of  input layers
!       ko    - number of output layers
!       ns    - number of fields
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
  real :: si_min(ns),si_max(ns)

! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ns
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision
!
  do i=1,ns
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
  do k=1,ko  !output layers
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
        do i=1,ns
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ns
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
        do i=1,ns
          o  = si((lt+lb)/2,i)  !offset to reduce round-off
          c0 = si(lt,i) - o - 0.5*slope(lt,i)
          sz=  dpi(lt)*(c0*qt0 + slope(lt,i)*qt1)
          do l=lt+1,lb-1
            sz = sz+dpi(l)*(si(l,i) - o)
          enddo !l
          c0 = si(lb,i) - o - 0.5*slope(lb,i)
          sz = sz+dpi(lb)*(c0*qb0 + slope(lb,i)*qb1)
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      else  !single layer
        qt1 = (xb**2-xt**2 - (xb-xt))*0.5
        do i=1,ns
          sz = dpi(lt)*(slope(lt,i)*qt1)
          so(k,i) = si(lt,i) + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      endif
    endif !thin:std layer
  enddo !k

end subroutine hybgen_plm_remap

!> Set up the coefficients for PPM remapping of a set of scalars
subroutine hybgen_ppm_coefs(s, h_src, edges, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The scalar fields to work on
  real,    intent(in)  :: s(nk,ns)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2,ns) !< The PPM interpolation edge values of the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic piecewise parabolic across each input cell
!
!     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       edges - cell edge scalar values for the PPM reconstruction
!                edges.1 is value at interface above
!                edges.2 is value at interface below
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  real :: dp(nk) ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  logical :: PCM_layer(nk) ! True for layers that should use PCM remapping, either because they are
                           ! very thin, or because this is specified by PCM_lay.
  real :: da        ! Difference between the unlimited scalar edge value estimates [A]
  real :: a6        ! Scalar field differences that are proportional to the curvature [A]
  real :: slk, srk  ! Differences between adjacent cell averages of scalars [A]
  real :: sck       ! Scalar differences across a cell.
  real :: as(nk)    ! Scalar field difference across each cell [A]
  real :: al(nk), ar(nk)   ! Scalar field at the left and right edges of a cell [A]
  real :: h112(nk+1), h122(nk+1)  ! Combinations of thicknesses [H ~> m or kg m-2]
  real :: I_h12(nk+1) ! Inverses of combinations of thickesses [H-1 ~> m-1 or m2 kg-1]
  real :: h2_h123(nk)  ! A ratio of a layer thickness of the sum of 3 adjacent thicknesses [nondim]
  real :: I_h0123(nk)     ! Inverse of the sum of 4 adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: h01_h112(nk+1) ! A ratio of sums of adjacent thicknesses [nondim], 2/3 in the limit of uniform thicknesses.
  real :: h23_h122(nk+1) ! A ratio of sums of adjacent thicknesses [nondim], 2/3 in the limit of uniform thicknesses.
  integer :: k, i

  ! This PPM remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
  do k=1,nk ; dp(k) = max(h_src(k), thin) ; enddo

  ! Specify the layers that will use PCM remapping.
  if (present(PCM_lay)) then
    do k=1,nk ; PCM_layer(k) = (PCM_lay(k) .or. dp(k) <= thin) ; enddo
  else
    do k=1,nk ; PCM_layer(k) = (dp(k) <= thin) ; enddo
  endif

  !compute grid metrics
  do k=2,nk
    h112(K) = 2.*dp(k-1) + dp(k)
    h122(K) = dp(k-1) + 2.*dp(k)
    I_h12(K) = 1.0 / (dp(k-1) + dp(k))
  enddo !k
  do k=2,nk-1
    h2_h123(k) = dp(k) / (dp(k) + (dp(k-1)+dp(k+1)))
  enddo
  do K=3,nk-1
    I_h0123(K) = 1.0 / ((dp(k-2) + dp(k-1)) + (dp(k) + dp(k+1)))

    h01_h112(K) = (dp(k-2) + dp(k-1)) / (2.0*dp(k-1) + dp(k))
    h23_h122(K) = (dp(k) + dp(k+1))   / (dp(k-1) + 2.0*dp(k))
  enddo

  do i=1,ns
    !Compute average slopes: Colella, Eq. (1.8)
    as(1) = 0.
    do k=2,nk-1
      if (PCM_layer(k)) then  !use PCM
        as(k) = 0.0
      else
        slk = s(k,  i)-s(k-1,i)
        srk = s(k+1,i)-s(k,  i)
        if (slk*srk > 0.) then
          sck = h2_h123(k)*( h112(K)*srk*I_h12(K+1) + h122(K+1)*slk*I_h12(K) )
          as(k) = sign(min(abs(2.0*slk), abs(sck), abs(2.0*srk)), sck)
        else
          as(k) = 0.
        endif
      endif  !PCM:PPM
    enddo !k
    as(nk) = 0.
    !Compute "first guess" edge values: Colella, Eq. (1.6)
    al(1) = s(1,i)  ! 1st layer PCM
    ar(1) = s(1,i)  ! 1st layer PCM
    al(2) = s(1,i)  ! 1st layer PCM
    do K=3,nk-1
      ! This is a 4th order explicit edge value estimate.
      al(k) = (dp(k)*s(k-1,i) + dp(k-1)*s(k,i)) * I_h12(K) &
            + I_h0123(K)*( 2.*dp(k)*dp(k-1)*I_h12(K)*(s(k,i)-s(k-1,i)) * &
                           ( h01_h112(K) - h23_h122(K) ) &
                    + (dp(k)*as(k-1)*h23_h122(K) - dp(k-1)*as(k)*h01_h112(K)) )
     ar(k-1) = al(k)
    enddo !k
    ar(nk-1) = s(nk,i) ! last layer PCM
    al(nk)  = s(nk,i)  ! last layer PCM
    ar(nk)  = s(nk,i)  ! last layer PCM
    !Impose monotonicity: Colella, Eq. (1.10)
    do k=2,nk-1
      if ((PCM_layer(k)) .or. ((s(k+1,i)-s(k,i))*(s(k,i)-s(k-1,i)) <= 0.)) then !local extremum
        al(k) = s(k,i)
        ar(k) = s(k,i)
      else
        da = ar(k)-al(k)
        a6 = 6.0*s(k,i) - 3.0*(al(k)+ar(k))
        if (da*a6 > da*da) then !peak in right half of zone
          al(k) = 3.0*s(k,i) - 2.0*ar(k)
        elseif (da*a6 < -da*da) then !peak in left half of zone
          ar(k) = 3.0*s(k,i) - 2.0*al(k)
        endif
      endif
    enddo !k
    !Set coefficients
    do k=1,nk
      edges(k,1,i) = al(k)
      edges(k,2,i) = ar(k)
    enddo !k
  enddo !i

end subroutine hybgen_ppm_coefs


!> Set up the coefficients for PPM remapping of a set of scalars
subroutine hybgen_weno_coefs(s, h_src, edges, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: s(nk,ns)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2,ns) !< The WENO interpolation edge values of the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the
!             interfacial values
!
!     This scheme might have ben developed by Shchepetkin. A.F., personal communication.
!     See also Engwirda, D., and M. Kelley, A WENO-type slope-limiter for a family of piecewise
!       polynomial methods, arXive:1606.08188v1, 27 June 2016.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       edges - cell edge scalar values for the WENO reconstruction
!                edges.1 is value at interface above
!                edges.2 is value at interface below
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
!-----------------------------------------------------------------------
!
!  real, parameter :: dsmll=1.0e-8  ! This has units of [A2], and hence can not be a parameter.
!
  real :: curv_cell   ! An estimate of the tracer curvature centered on a cell times the grid
                      ! spacing [A H-1 ~> A m-1 or A kg m-2]
  real :: seh1, seh2  ! Tracer slopes at the cell edges times the cell grid spacing [A]
  real :: q01, q02    ! Various tracer differences between a cell average and the edge values [A]
  real :: q001, q002  ! Tracer slopes at the cell edges times the cell grid spacing [A]
  real :: ds2a, ds2b  ! Squared tracer differences between a cell average and the edge values [A2]
  logical :: PCM_layer(nk) ! True for layers that should use PCM remapping, either because they are
                      ! very thin, or because this is specified by PCM_lay.
  real :: dp(nk)      ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  real :: qdpkm(nk)   ! Inverse of the sum of two adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: qdpkmkp(nk) ! Inverse of the sum of three adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: dpkm2kp(nk) ! Twice the distance between the centers of the layers two apart [H ~> m or kg m-2]
  real :: zw(nk,2)    ! Squared combinations of the differences between the the cell average tracer
                      ! concentrations and the left and right edges [A2]
  real :: min_ratio   ! The minimum ratio of the values of zw used to interpolate the edge values [nondim]
  real :: wt1         ! The weight of the upper layer in the interpolated shared edge value [nondim]
  real :: slope_edge(nk+1)  ! Tracer slopes at the edges [A H-1 ~> A m-1 or A kg m-2]
  real :: val_edge(nk+1)    ! A weighted average edge concentration [A]
  integer :: i, k

  min_ratio = 1.0e-8

  ! The WENO remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
  do k=1,nk ; dp(k) = max(h_src(k), thin) ; enddo

  ! Specify the layers that will use PCM remapping.
  if (present(PCM_lay)) then
    do k=1,nk ; PCM_layer(k) = (PCM_lay(k) .or. dp(k) <= thin) ; enddo
  else
    do k=1,nk ; PCM_layer(k) = (dp(k) <= thin) ; enddo
  endif

  !compute grid metrics
  do k=2,nk-1
    qdpkm(  K) = 1.0 / (dp(k-1) + dp(k))
    qdpkmkp(k) = 1.0 / (dp(k-1) + dp(k) + dp(k+1))
    dpkm2kp(k) = dp(k-1) + 2.0*dp(k) + dp(k+1)
  enddo !k
  qdpkm(nk) = 1.0 / (dp(nk-1) + dp(nk))

  do i=1,ns
    do K=2,nk
      slope_edge(K) = qdpkm(K) * (s(k,i)-s(k-1,i))
    enddo !k
    k = 1  !PCM first layer
    edges(k,1,i) = s(k,i)
    edges(k,2,i) = s(k,i)
    zw(k,1) = 0.0
    zw(k,2) = 0.0
    do k=2,nk-1
      if ((slope_edge(K)*slope_edge(K+1) < 0.0) .or. PCM_layer(k)) then  !use PCM
        edges(k,1,i) = s(k,i)
        edges(k,2,i) = s(k,i)
        zw(k,1) = 0.0
        zw(k,2) = 0.0
      else
        seh1 = dp(k)*slope_edge(K+1)
        seh2 = dp(k)*slope_edge(K)
        q01 = dpkm2kp(k)*slope_edge(K+1)
        q02 = dpkm2kp(k)*slope_edge(K)
        if (abs(seh1) > abs(q02)) then
          seh1 = q02
        endif
        if (abs(seh2) > abs(q01)) then
          seh2 = q01
        endif
        curv_cell = (seh1 - seh2) * qdpkmkp(k)
        q001 = seh1 - curv_cell*dp(k+1)
        q002 = seh2 + curv_cell*dp(k-1)
        ! q001 = (seh1 * (dp(k-1) + dp(k)) + seh2 * dp(k+1)) * qdpkmkp(k)
        ! q002 = (seh2 * (dp(k+1) + dp(k)) + seh1 * dp(k-1)) * qdpkmkp(k)

        edges(k,2,i) = s(k,i) + q001
        edges(k,1,i) = s(k,i) - q002
        zw(k,1) = (2.0*q001 - q002)**2
        zw(k,2) = (2.0*q002 - q001)**2
      endif  !PCM:WENO
    enddo !k
    k = nk  !PCM last layer
    edges(k,1,i) = s(k,i)
    edges(k,2,i) = s(k,i)
    zw(k,  1) = 0.0
    zw(k,  2) = 0.0

    do k=2,nk
      ! This was the original code based on that in Hycom, but because zw has
      ! dimensions of [A2], it can not use a constant (hard coded) value of dsmll.
      !   ds2a = max(zw(k-1,2), dsmll)
      !   ds2b = max(zw(k,  1), dsmll)
      !   val_edge(K) = (ds2b*edges(k-1,2,i)+ds2a*edges(k,1,i)) / (ds2b+ds2a)
      ! Use a weighted average of the two layers' estimated edge values as the actual edge value.
      if (zw(k,1) + zw(k-1,2) <= 0.0) then
        wt1 = 0.5
      elseif (zw(k,1) <= min_ratio * (zw(k,1) + zw(k-1,2))) then
        wt1 = min_ratio
      elseif (zw(k-1,2) <= min_ratio * (zw(k,1) + zw(k-1,2))) then
        wt1 = (1.0 - min_ratio)
      else
        wt1 = zw(k,1) / (zw(k,1) + zw(k-1,2))
      endif
      val_edge(k) = wt1*edges(k-1,2,i) + (1.0-wt1)*edges(k,1,i)
    enddo !k
    val_edge(   1) = 2.0*s( 1,i)-val_edge( 2)  !not used?
    val_edge(nk+1) = 2.0*s(nk,i)-val_edge(nk)  !not used?

    do k=2,nk-1
      if (.not.PCM_layer(k)) then  !don't use PCM
        q01 = val_edge(K+1) - s(k,i)
        q02 = s(k,i) - val_edge(K)
        if (q01*q02 < 0.0) then
          q01 = 0.0
          q02 = 0.0
        elseif (abs(q01) > abs(2.0*q02)) then
          q01 = 2.0*q02
        elseif (abs(q02) > abs(2.0*q01)) then
          q02 = 2.0*q01
        endif
        edges(k,1,i) = s(k,i) - q02
        edges(k,2,i) = s(k,i) + q01
      endif  ! PCM:WENO
    enddo !k
  enddo !i

end subroutine hybgen_weno_coefs

!> Do piecewise parabolic remapping for a set of scalars
subroutine hybgen_ppm_remap(si, dpi, edges, so, dpo, ki, ko, ns, thin)
  integer, intent(in)  :: ki        !< The number of input layers
  integer, intent(in)  :: ko        !< The number of output layers
  integer, intent(in)  :: ns        !< The scalar fields to work on
  real,    intent(in)  :: si(ki,ns) !< The input scalar fields [A]
  real,    intent(in)  :: dpi(ki)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: edges(ki,2,ns) !< The interpolation edge values of the scalar fields [A]
  real,    intent(out) :: so(ko,ns) !< The output scalar fields [A]
  real,    intent(in)  :: dpo(ko)   !< The output layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]

!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: A monotonic piecewise parabolic subgrid cell distribution within each
!             input cell is uniquely determined from the cell average and edge values.
!             The output is the average of the interpolation
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
!       edges - cell edge scalar values for the PPM reconstruction
!                  edges.1 is value at interface above
!                  edges.2 is value at interface below
!       ki    - number of  input layers
!       ko    - number of output layers
!       ns    - number of fields
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
  real    si_min(ns),si_max(ns)

! --- enforce minval(si(:,i)) <= minval(so(:,i)) and
! ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ns
! --- in particular this enforces non-negativity, e.g. of tracers
! --- only required due to finite precision

  do i=1,ns
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
  do k=1,ko  !output layers
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
        do i=1,ns
          so(k,i) = so(k-1,i)
        enddo !i
      else !thin surface layer
        do i=1,ns
          so(k,i) = si(k,i)
        enddo !i
      endif
    else
      ! Form layer averages.
!         if (pi(lb) > zt) then
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
        do i=1,ns
          o = si((lt+lb)/2,i)  !offset to reduce round-off
          sz = dpt*(qt0*(si(lt,i)  -o) + &
                    qt1*(edges(lt,1,i)-o) + &
                    qt2*(edges(lt,2,i)-o)  )
          do l=lt+1,lb-1
            sz = sz+dpi(l)*(si(l,i) - o)
          enddo !l
          sz  = sz + dpb*(qb0*(si(lb,i)  -o) + &
                          qb1*(edges(lb,1,i)-o) + &
                          qb2*(edges(lb,2,i)-o)  )
          so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      else !single layer
        qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
        qt2 = qt1 - 1.0 + (xb+xt)
        do i=1,ns
          o =     si(lt,i)  !offset to reduce round-off
          sz = qt1*(edges(lt,1,i)-o) + qt2*(edges(lt,2,i)-o)
          so(k,i) = o + sz
          so(k,i) = max( si_min(i), so(k,i) )
          so(k,i) = min( si_max(i), so(k,i) )
        enddo !i
      endif !layers
    endif !thin:std layer
  enddo !k

end subroutine hybgen_ppm_remap

end module MOM_hybgen_remap
