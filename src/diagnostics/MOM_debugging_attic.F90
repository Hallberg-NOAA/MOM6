!> Provides archaic checksumming functions for debugging
!!
!! This module contains subroutines that had previously been used perform various
!! debugging functions for MOM6 but are no longer called from anywhere in the
!! MOM6 code.  They are retained because they might be useful as temporary additions
!! to the code for future debugging exercises.
module MOM_debugging_attic

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums, only : hchksum, Bchksum, uvchksum
use MOM_coms, only : reproducing_sum
use MOM_debugging, only : check_redundant_C, check_redundant_B, check_redundant_T, check_redundant
use MOM_domains, only : BGRID_NE, AGRID, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : log_version, param_file_type, get_param
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type

implicit none ; private

public :: vec_chksum, vec_chksum_C, vec_chksum_B, vec_chksum_A
public :: MOM_debugging_attic_init, totalStuff, totalTandS

!> Do checksums on the components of a C-grid vector
interface vec_chksum
  module procedure chksum_vec_C3d, chksum_vec_C2d
end interface vec_chksum
!> Do checksums on the components of a C-grid vector
interface vec_chksum_C
  module procedure chksum_vec_C3d, chksum_vec_C2d
end interface vec_chksum_C
!> Do checksums on the components of a B-grid vector
interface vec_chksum_B
  module procedure chksum_vec_B3d, chksum_vec_B2d
end interface vec_chksum_B
!> Do checksums on the components of an A-grid vector
interface vec_chksum_A
  module procedure chksum_vec_A3d, chksum_vec_A2d
end interface vec_chksum_A

! Note: these parameters are module data but ONLY used when debugging and
!       so can violate the thread-safe requirement of no module/global data.
logical :: debug = .false. !< Write out verbose debugging data
logical :: debug_chksums = .true. !< Perform checksums on arrays
logical :: debug_redundant = .true. !< Check redundant values on PE boundaries

contains

!> MOM_debugging_init initializes the MOM_debugging module, and sets
!! the parameterts that control which checks are active for MOM6.
subroutine MOM_debugging_attic_init(param_file)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_debugging_attic" ! This module's name.

  call log_version(param_file, mdl, version, debugging=.true.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_CHKSUMS", debug_chksums, &
                 "If true, checksums are performed on arrays in the "//&
                 "various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_REDUNDANT", debug_redundant, &
                 "If true, debug redundant data points during calls to "//&
                 "the various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)

end subroutine MOM_debugging_attic_init

!> Do a checksum and redundant point check on a 3d C-grid vector.
subroutine chksum_vec_C3d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                  intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),             intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector [a]
  real, dimension(G%isd:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector [a]
  integer,                 optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                 optional, intent(in)    :: scalars !< If true this is a pair of
                                                             !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call uvchksum(mesg, u_comp, v_comp, G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_C3d

!> Do a checksum and redundant point check on a 2d C-grid vector.
subroutine chksum_vec_C2d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector [a]
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector [a]
  integer,               optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,               optional, intent(in)    :: scalars !< If true this is a pair of
                                                           !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call uvchksum(mesg, u_comp, v_comp, G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_C2d

!> Do a checksum and redundant point check on a 3d B-grid vector.
subroutine chksum_vec_B3d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                   intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),              intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: u_comp !< The u-component of the vector [a]
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector [a]
  integer,                  optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                  optional, intent(in)    :: scalars !< If true this is a pair of
                                                              !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call Bchksum(u_comp, mesg//"(u)", G%HI, halos)
    call Bchksum(v_comp, mesg//"(v)", G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_B3d

! Do a checksum and redundant point check on a 2d B-grid vector.
subroutine chksum_vec_B2d(mesg, u_comp, v_comp, G, halos, scalars, symmetric)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: u_comp !< The u-component of the vector [a]
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector [a]
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  logical,                optional, intent(in)    :: symmetric !< If true, do the checksums on the
                                                            !! full symmetric computational domain.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call Bchksum(u_comp, mesg//"(u)", G%HI, halos, symmetric=symmetric)
    call Bchksum(v_comp, mesg//"(v)", G%HI, halos, symmetric=symmetric)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_B2d

!> Do a checksum and redundant point check on a 3d C-grid vector.
subroutine chksum_vec_A3d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector [a]
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: v_comp !< The v-component of the vector [a]
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call hchksum(u_comp, mesg//"(u)", G%HI, halos)
    call hchksum(v_comp, mesg//"(v)", G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_A3d

!> Do a checksum and redundant point check on a 2d C-grid vector.
subroutine chksum_vec_A2d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),               intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),          intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector [a]
  real, dimension(G%isd:,G%jsd:), intent(in)    :: v_comp !< The v-component of the vector [a]
  integer,              optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,              optional, intent(in)    :: scalars !< If true this is a pair of
                                                          !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call hchksum(u_comp, mesg//"(u)", G%HI, halos)
    call hchksum(v_comp, mesg//"(v)", G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_A2d

!> This function returns the sum over computational domain of all
!! processors of hThick*stuff, where stuff is a 3-d array at tracer points.
function totalStuff(HI, hThick, areaT, stuff)
  type(hor_index_type),               intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights [m]
  real, dimension(HI%isd:,HI%jsd:),   intent(in) :: areaT  !< The array of cell areas [m2]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: stuff  !< The array of stuff to be summed [a]
  real                                         :: totalStuff !< the globally integrated amount of stuff [a m3]
  ! Local variables
  real, dimension(HI%isc:HI%iec, HI%jsc:HI%jec) :: tmp_for_sum ! The column integrated amout of stuff in a cell [a m3]
  integer :: i, j, k, nz

  nz = size(hThick,3)
  tmp_for_sum(:,:) = 0.0
  do k=1,nz ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    tmp_for_sum(i,j) = tmp_for_sum(i,j) + hThick(i,j,k) * stuff(i,j,k) * areaT(i,j)
  enddo ; enddo ; enddo
  totalStuff = reproducing_sum(tmp_for_sum)

end function totalStuff

!> This subroutine display the total thickness, temperature and salinity
!! as well as the change since the last call.
subroutine totalTandS(HI, hThick, areaT, temperature, salinity, mesg)
  type(hor_index_type),               intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights [m]
  real, dimension(HI%isd:,HI%jsd:),   intent(in) :: areaT  !< The array of cell areas [m2]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: temperature !< The temperature field to sum [degC]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: salinity    !< The salinity field to sum [ppt]
  character(len=*),                   intent(in) :: mesg        !< An identifying message
  ! NOTE: This subroutine uses "save" data which is not thread safe and is purely for
  ! extreme debugging without a proper debugger.
  real, save :: totalH = 0.   ! The total ocean volume, saved for the next call [m3]
  real, save :: totalT = 0.   ! The total volume integrated ocean temperature, saved for the next call [degC m3]
  real, save :: totalS = 0.   ! The total volume integrated ocean salinity, saved for the next call [ppt m3]
  ! Local variables
  logical, save :: firstCall = .true.
  real, dimension(HI%isc:HI%iec, HI%jsc:HI%jec) :: tmp_for_sum ! The volume of each column [m3]
  real :: thisH, delH  ! The total ocean volume and the change from the last call [m3]
  real :: thisT, delT  ! The current total volume integrated temperature and the change from the last call [degC m3]
  real :: thisS, delS  ! The current total volume integrated salinity and the change from the last call [ppt m3]
  integer :: i, j, k, nz

  nz = size(hThick,3)
  tmp_for_sum(:,:) = 0.0
  do k=1,nz ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    tmp_for_sum(i,j) = tmp_for_sum(i,j) + hThick(i,j,k) * areaT(i,j)
  enddo ; enddo ; enddo
  thisH = reproducing_sum(tmp_for_sum)
  thisT = totalStuff(HI, hThick, areaT, temperature)
  thisS = totalStuff(HI, hThick, areaT, salinity)

  if (is_root_pe()) then
    if (firstCall) then
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Totals H,T,S:',thisH,thisT,thisS,' ',mesg
      firstCall = .false.
    else
      delH = thisH - totalH
      delT = thisT - totalT
      delS = thisS - totalS
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Tot/del H,T,S:',thisH,thisT,thisS,delH,delT,delS,' ',mesg
    endif
  endif

end subroutine totalTandS

end module MOM_debugging_attic
