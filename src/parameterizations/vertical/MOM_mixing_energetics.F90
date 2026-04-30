!> Energetically consistent planetary boundary layer parameterization
module MOM_mixing_energetics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_coms,           only : EFP_type, real_to_EFP, EFP_to_real, operator(+), assignment(=), EFP_sum_across_PEs
use MOM_debugging,      only : hchksum, is_NaN
use MOM_diag_mediator,  only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator,  only : time_type, diag_ctrl
use MOM_domains,        only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler,  only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser,    only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,   only : forcing
use MOM_grid,           only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_intrinsic_functions, only : cuberoot
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid,   only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public energetic_mixing, energetic_mixing_init, energetic_mixing_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters for the MOM_energetic_mixing module
type, public :: energetic_mixing_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.

  logical :: tridiagonal_w2 = .false.  !< If true, use a linearized tridiagonal solver for w**2
                             !! after the first iteration.
  logical :: use_prior_diffusivity  !< If true, include a background diffusivity from other
                             !! processes in the energetic_mixing diffusivity calculation.

  real    :: L_interior      !< The mixing distance far from the edges, in thickness units [H ~> m or kg m-2]
  real    :: Bdry_fac        !< A scaling factor for the mixing length distance near the top and bottom
                             !! boundaries, perhaps the Von Karman constant. [nondim]
  real    :: w_tol           !< The tolerance for convergence of the iterations for the turbulent velocity
                             !! in the turbulent kinetic energy equation [Z T-1 ~> m s-1]
  real    :: decay_rate      !< A background rate of energy decay in the turbulent kinetic energy
                             !! equation [T-1 ~> s-1]
  real    :: decay_w_L_scale !< A coefficient relating the turbulent velocity divided by the local
                             !! mixing distance to a turbulent kinetic energy decay rate [nondim]

  !/ Constants
  real    :: VonKar          !< The von Karman coefficient as used in the eMix module [nondim]


  !/ Mixing Length terms
  integer :: max_mixing_its  !< The maximum number of iterations that can be used to find a
                             !! self-consistent set of diffusivities with Use_MLD_iteration.


  !/ Options for documenting differences from parameter choices
  integer :: options_diff    !< If positive, this is a coded integer indicating a pair of
                             !! settings whose differences are diagnosed in a passive diagnostic mode
                             !! via extra calls to eMix_column.  If this is 0 or negative no extra
                             !! calls occur.

  !/ Others
  type(time_type), pointer :: Time=>NULL() !< A pointer to the ocean model's clock.

  logical :: TKE_diagnostics = .false. !< If true, diagnostics of the TKE budget are being calculated.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.

  type(EFP_type), dimension(2) :: sum_its !< The total number of iterations and columns worked on

  !>@{ Diagnostic IDs
  integer :: id_TKE_mixing = -1
  integer :: id_TKE_forcing = -1
  integer :: id_frac_en_diff = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1

  ! The next options are used when passively diagnosing sensitivities from parameter choices
  integer :: id_opt_diff_Kd_eMix = -1, id_opt_maxdiff_Kd_eMix = -1
  !>@}
end type energetic_mixing_CS

logical :: report_avg_its = .false.  !< Report the average number of eMix iterations for debugging.

!> A type for conveniently passing around eMix diagnostics for a column.
type, public :: eMix_column_diags ; private
  !>@{ Local column copies of energy change diagnostics, all in [R Z3 T-3 ~> W m-2].
  real :: dTKE_forcing, dTKE_mixing ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  !>@}
  integer :: eMix_its !< The number of iterations used to find self-consistent mixing
end type eMix_column_diags

contains

!>   This subroutine determines interior ocean diffusivities from implicit energetics, using a
!! simple fixed turbulent length scale and a supplied source of turbulent kinetic energy from
!! another parameterization.  All calculations are done implicitly, and there is no stability
!! limit on the time step.
subroutine energetic_mixing(h_3d, tv, dSV_dT, dSV_dS, TKE_forcing, dt, Kd_int, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_3d   !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in)    :: TKE_forcing !< The rate of TKE forcing that is being
                                                   !! applied to the water surrounding each interface
                                                   !! using a finite-volume integral
                                                   !! [R Z3 T-3 ~> W m-2].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: Kd_int !< The total diffusivities at interfaces
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(energetic_mixing_CS), pointer     :: CS     !< Energetic mixing control structure

  !   This subroutine determines interior ocean diffusivities from implicit energetics,
  ! using a simple fixed turbulent length scale and a supplied source of turbulent kinetic energy
  ! from another parameterization.  This source could be considered \Gamma \epsilon from other
  ! simple closures.  All calculations are done implicitly, and there
  ! is no stability limit on the time step.  In the limit where the boudaries are far away and
  ! the water is well stratified and the energy source is slowly varying, the returned diffusivities
  ! approximate the Osborn relation results (\kappa = \Gamma \epsilon / N^2), but this version works
  ! equally well when the stratification vanishes and it avoids overly large or strongly varying
  ! diffusivities.  The two parameters are the peak interior length scale and a TKE decay rate, but
  ! the well-stratified interior ocean diffusivities vary only weakly with these two parameters.

  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    h_2d, &         ! A 2-d slice of the layer thickness [H ~> m or kg m-2].
    dz_2d, &        ! A 2-d slice of the vertical distance across layers [Z ~> m].
    T_2d, &         ! A 2-d slice of the layer temperatures [C ~> degC].
    S_2d, &         ! A 2-d slice of the layer salinities [S ~> ppt].
    dSV_dT_2d, &    ! A 2-d slice of dSV_dT [R-1 C-1 ~> m3 kg-1 degC-1].
    dSV_dS_2d       ! A 2-d slice of dSV_dS [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    Kd_2d, &        ! A 2-d version of the added diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    Kd_in_2d, &     ! A 2-d version of the input diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    TKE_forcing_2d  ! A 2-d slice of TKE_forcing [R Z3 T-3 ~> W m-2].
  real, dimension(SZK_(GV)) :: &
    h, &            ! The layer thickness [H ~> m or kg m-2].
    dz, &           ! The vertical distance across layers [Z ~> m].
    T0, &           ! The initial layer temperatures [C ~> degC].
    S0, &           ! The initial layer salinities [S ~> ppt].
    dSV_dT_1d, &    ! The partial derivatives of specific volume with temperature [R-1 C-1 ~> m3 kg-1 degC-1].
    dSV_dS_1d       ! The partial derivatives of specific volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)+1) :: &
    Kd, &           ! The diapycnal diffusivity due to eMix [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
    Kd_other, &     ! A diapycnal diffusivity due to other processes whose energetics are dealt
                    ! with elsewhere or a molecular diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    TKE_source, &   ! Forcing of the TKE in the layer coming from TKE_forcing integrated
                    ! through a timestep [R Z3 T-2 ~> J m-2].
    mixvel, &       ! A turbulent mixing velocity [Z T-1 ~> m s-1].
    mixlen          ! A turbulent mixing length [Z ~> m].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: I_rho     ! The inverse of the Boussinesq reference density [R-1 ~> m3 kg-1]
  real :: I_dt      ! The Adcroft reciprocal of the timestep [T-1 ~> s-1]
  real :: I_rho0dt  ! The inverse of the Boussinesq reference density times the time
                    ! step [R-1 T-1 ~> m3 kg-1 s-1]

  type(eMix_column_diags) :: eCD ! A container for passing around diagnostics.

  ! The following variables are used for diagnostics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    diag_Velocity_Scale, & ! The velocity scale used in getting Kd [Z T-1 ~> m s-1]
    diag_Mixing_Length     ! The length scale used in getting Kd [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    ! The next 7 diagnostics are terms in the mixed layer TKE budget, all in [R Z3 T-3 ~> W m-2].
    diag_TKE_forcing, & ! The TKE sink required to mix surface penetrating shortwave heating [R Z3 T-3 ~> W m-2]
    diag_TKE_mixing, &  ! The work done by TKE to deepen the mixed layer [R Z3 T-3 ~> W m-2]
    frac_en_diff        ! The fractional difference between the energy used and the energy input [nondim]

  ! The following variables are used for debugging.
  real :: pres(SZK_(GV)+1)    ! Interface pressures [R L2 T-2 ~> Pa].
  real :: Kddt_h(SZK_(GV)+1)  ! The diapycnal diffusivity times a timestep divided by the
                              ! average thicknesses around a layer [H ~> m or kg m-2].
  real :: Tf(SZK_(GV))        ! Updated values of the temperatures after mixing [C ~> degC]
  real :: Sf(SZK_(GV))        ! Updated values of the salinities after mixing [S ~> ppt].
  real :: dMass               ! The mass per unit area within a layer [R Z ~> kg m-2].
  real :: dPres               ! The hydrostatic pressure change across a layer [R L2 T-2 ~> Pa].
  real :: dT_to_dPE(SZK_(GV)) ! Partial derivative of column potential energy with the temperature
                              ! changes within a layer [R Z L2 T-2 C-1 ~> J m-2 degC-1]
  real :: dS_to_dPE(SZK_(GV)) ! Partial derivative of column potential energy with the salinity
                              ! changes within a layer [R Z L2 T-2 S-1 ~> J m-2 ppt-1]
  real :: PE_chg_tot1D        ! Changes in column potential energy [R Z L2 T-2 ~> J m-2]
  real :: TKE_force_tot1D     ! The time-integrated column-integrated energy driving mixing
                              ! within a timestep [R Z L2 T-2 ~> J m-2]

  ! The following variables are only used for diagnosing sensitivities to eMix settings
  real, dimension(SZK_(GV)+1) :: &
    Kd_1, Kd_2      ! Diapycnal diffusivities found with different eMix options [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: diff_Kd(SZI_(G),SZJ_(G),SZK_(GV)+1) ! The change in diapycnal diffusivities found with different
                        ! eMix options [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: max_abs_diff_Kd(SZI_(G),SZJ_(G))  ! The column maximum magnitude of the change in diapycnal
                        ! diffusivities found with different eMix options [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(eMix_column_diags) :: eCD_tmp   ! A container for not passing around diagnostics.
  type(energetic_mixing_CS)  :: CS_tmp1, CS_tmp2 ! Copies of the energetic PBL control structure that
                                       ! can be modified to test for sensitivities
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "energetic_mixing: "// &
         "Module must be initialized before it is used.")
  if (.not. CS%initialized) call MOM_error(FATAL, "energetic_mixing: "//&
         "Module must be initialized before it is used.")
  if (.not. associated(tv%eqn_of_state)) call MOM_error(FATAL, &
      "energetic_mixing: Temperature, salinity and an equation of state "//&
      "must now be used.")

  h_neglect = GV%H_subroundoff
  I_rho = GV%H_to_Z * GV%RZ_to_H ! == 1.0 / GV%Rho0 ! This is not used when fully non-Boussinesq.
  I_dt = 0.0 ; if (dt > 0.0) I_dt = 1.0 / dt
  I_rho0dt = 1.0 / (GV%Rho0 * dt)  ! This is not used when fully non-Boussinesq.

  ! Zero out diagnostics before accumulation.
  if (CS%TKE_diagnostics) then
    !!OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      diag_TKE_forcing(i,j) = 0.0
      diag_TKE_mixing(i,j) = 0.0
      !; diag_TKE_unbalanced(i,j) = 0.0
    enddo ; enddo
  endif
  if (CS%debug .or. (CS%id_Mixing_Length>0)) diag_Mixing_Length(:,:,:) = 0.0
  if (CS%debug .or. (CS%id_Velocity_Scale>0)) diag_Velocity_Scale(:,:,:) = 0.0

  ! CS_tmp is used to test sensitivity to parameter setting changes.
  ! if (CS%options_diff > 0) then
  !   CS_tmp1 = CS ; CS_tmp2 = CS

  !   elseif (CS%options_diff == 4) then
  !   elseif (CS%options_diff == 5) then
  !   endif

  !   if (CS%id_opt_diff_Kd_eMix > 0)    diff_Kd(:,:,:) = 0.0
  !   if (CS%id_opt_maxdiff_Kd_eMix > 0) max_abs_diff_Kd(:,:) = 0.0
  ! endif

  !!OMP parallel do default(private) shared(js,je,nz,is,ie,h_3d,tv,dt,I_dt, &
  !!OMP                                  CS,G,GV,US,fluxes,TKE_forcing,dSV_dT,dSV_dS,Kd_int)
  do j=js,je
    ! Copy the thicknesses and other fields to 2-d arrays.
    do k=1,nz ; do i=is,ie
      h_2d(i,k) = h_3d(i,j,k) ; T_2d(i,k) = tv%T(i,j,k) ; S_2d(i,k) = tv%S(i,j,k)
      TKE_forcing_2d(i,k) = TKE_forcing(i,j,k)
      dSV_dT_2d(i,k) = dSV_dT(i,j,k) ; dSV_dS_2d(i,k) = dSV_dS(i,j,k)
    enddo ; enddo
    do K=1,nz+1 ; do i=is,ie
      TKE_forcing_2d(i,K) = TKE_forcing(i,j,K)
      Kd_in_2d(i,K) = Kd_int(i,j,K)
    enddo ; enddo
    call thickness_to_dz(h_3d, tv, dz_2d, j, G, GV)

    !   Determine the mixing energy available integrated over each layer.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then

      ! Copy the thicknesses and other fields to 1-d arrays.
      do k=1,nz
        h(k) = h_2d(i,k) + GV%H_subroundoff ; dz(k) = dz_2d(i,k) + GV%dZ_subroundoff
        T0(k) = T_2d(i,k) ; S0(k) = S_2d(i,k)
        dSV_dT_1d(k) = dSV_dT_2d(i,k) ; dSV_dS_1d(k) = dSV_dS_2d(i,k)
      enddo
      ! do K=1,nz+1 ; Kd(K) = 0.0 ; enddo
      do K=1,nz+1
        TKE_source(K) =  GV%RZ_to_H * TKE_forcing_2d(i,K)*dt
        Kd_other(K) = Kd_in_2d(i,K)
      enddo

      call eMix_column(h, dz, T0, S0, dSV_dT_1d, dSV_dS_1d, TKE_source, &
                       dt, Kd, mixvel, mixlen, GV, US, CS, eCD, Kd_other)

      ! Copy the diffusivities to a 2-d array.
      do K=1,nz+1
        Kd_2d(i,K) = Kd(K)
      enddo

      if (CS%TKE_diagnostics) then
        diag_TKE_forcing(i,j) = diag_TKE_forcing(i,j) + eCD%dTKE_forcing
        diag_TKE_mixing(i,j) = diag_TKE_mixing(i,j) + eCD%dTKE_mixing
       ! diag_TKE_unbalanced(i,j) = diag_TKE_unbalanced(i,j) + eCD%dTKE_unbalanced
      endif
      ! Write mixing length and velocity scale to 3-D arrays for diagnostic output
      if (CS%debug .or. (CS%id_Mixing_Length > 0)) then ; do K=1,nz+1
        diag_Mixing_Length(i,j,K) = mixlen(K)
      enddo ; endif
      if (CS%debug .or. (CS%id_Velocity_Scale > 0)) then ; do K=1,nz+1
        if (.not. is_NaN(mixvel(K))) &
          diag_Velocity_Scale(i,j,K) = mixvel(K)
      enddo ; endif
      if (report_avg_its) then
        CS%sum_its(1) = CS%sum_its(1) + real_to_EFP(real(eCD%eMix_its))
        CS%sum_its(2) = CS%sum_its(2) + real_to_EFP(1.0)
      endif

      if (CS%options_diff > 0) then
        ! Call eMix_column with different parameter settings to diagnose sensitivities.
        ! These do not change the model state, and are only used for diagnostic purposes.
        if (CS%options_diff < 4) then
          call eMix_column(h, dz, T0, S0, dSV_dT_1d, dSV_dS_1d, TKE_source, dt, Kd_1, &
                           mixvel, mixlen, GV, US, CS_tmp1, eCD_tmp, Kd_other)
          call eMix_column(h, dz, T0, S0, dSV_dT_1d, dSV_dS_1d, TKE_source, dt, Kd_2, &
                           mixvel, mixlen, GV, US, CS_tmp2, eCD_tmp, Kd_other)
        endif

        if (CS%id_opt_diff_Kd_eMix > 0) then
          do K=1,nz+1 ; diff_Kd(i,j,K) = Kd_1(K) - Kd_2(K) ; enddo
        endif
        if (CS%id_opt_maxdiff_Kd_eMix > 0) then
          max_abs_diff_Kd(i,j) = 0.0
          do K=1,nz+1 ; max_abs_diff_Kd(i,j) = max(max_abs_diff_Kd(i,j), abs(Kd_1(K) - Kd_2(K))) ; enddo
        endif
      endif

      if (CS%id_frac_en_diff > 0) then
        ! Recalculate the change in column integrated potential energy to verify
        ! the plausibility of the calculations.

        Kddt_h(1) = 0.0 ; Kddt_h(nz+1) = 0.0
        do K=2,nz
          Kddt_h(K) = dt * Kd(K) / max(0.5*(dz(k-1) + dz(k)), GV%dZ_subroundoff)
        enddo

        call tridiag_T_down(T0, Kddt_h, h, GV, Tf)
        call tridiag_T_down(S0, Kddt_h, h, GV, Sf)

        pres(1) = 0.0
        PE_chg_tot1D = 0.0
        do k=1,nz
          dMass = GV%H_to_RZ * h(k)
          dPres = (GV%g_Earth * GV%H_to_RZ) * h(k)
          dT_to_dPE(k) = (dMass * (pres(K) + 0.5*dPres)) * dSV_dT_1d(k)
          dS_to_dPE(k) = (dMass * (pres(K) + 0.5*dPres)) * dSV_dS_1d(k)
          ! dT_to_dColHt(k) = dMass * dSV_dT(k) * CS%ColHt_scaling
          ! dS_to_dColHt(k) = dMass * dSV_dS(k) * CS%ColHt_scaling
          pres(K+1) = pres(K) + (GV%g_Earth * GV%H_to_RZ) * h(k)

          PE_chg_tot1D = PE_chg_tot1D + (dT_to_dPE(k) * (Tf(k) - T0(k)) + &
                                         dS_to_dPE(k) * (Sf(k) - S0(k)))
        enddo
        do K=1,nz+1
          TKE_force_tot1D = TKE_force_tot1D + TKE_forcing_2d(i,K)*dt
        enddo
        frac_en_diff(i,j) = (PE_chg_tot1D - TKE_force_tot1D) / (0.5*(TKE_force_tot1D + PE_chg_tot1D))
      endif

    else ! End of the ocean-point part of the i-loop
      ! For masked points, Kd_int must still be set (to 0) because it has intent out.
      do K=1,nz+1 ; Kd_2d(i,K) = 0. ; enddo
    endif ; enddo ! Close of i-loop - Note the unusual loop order, with k-loops inside i-loops.

    do K=1,nz+1 ; do i=is,ie
      if (.not. is_NaN(Kd_2d(i,K))) then
        Kd_int(i,j,K) = Kd_int(i,j,K) + Kd_2d(i,K)
      endif
    enddo ; enddo

  enddo ! j-loop

  if (CS%debug) then
    call hchksum(diag_Mixing_Length, "energetic_mixing Mixing_Length", G%HI, unscale=US%Z_to_m)
    call hchksum(diag_Velocity_Scale, "energetic_mixing Vel scale", G%HI, unscale=US%Z_to_m*US%s_to_T)
  endif

  if (CS%id_TKE_forcing > 0) call post_data(CS%id_TKE_forcing, diag_TKE_forcing, CS%diag)
  if (CS%id_TKE_mixing > 0) call post_data(CS%id_TKE_mixing, diag_TKE_mixing, CS%diag)
  if (CS%id_Mixing_Length > 0) call post_data(CS%id_Mixing_Length, diag_Mixing_Length, CS%diag)
  if (CS%id_Velocity_Scale >0) call post_data(CS%id_Velocity_Scale, diag_Velocity_Scale, CS%diag)
  if (CS%id_frac_en_diff >0) call post_data(CS%id_frac_en_diff, frac_en_diff, CS%diag)

  if (CS%options_diff > 0) then
    ! These diagnostics are only for determining sensitivities to different eMix settings.
    if (CS%id_opt_diff_Kd_eMix > 0)    call post_data(CS%id_opt_diff_Kd_eMix, diff_Kd, CS%diag)
    if (CS%id_opt_maxdiff_Kd_eMix > 0) call post_data(CS%id_opt_maxdiff_Kd_eMix, max_abs_diff_Kd, CS%diag)
  endif

end subroutine energetic_mixing


!> This subroutine determines the diffusivities from the implicit energetics with a simple mixing
!! length prescription for the mixing of a single column of water, including spreading of mixing
!! by diffusion from other mixing sources and auto-diffusion.
subroutine eMix_column(h, dz, T0, S0, dSV_dT, dSV_dS, TKE_source, &
                       dt, Kd, mixvel, mixlen, GV, US, CS, eCD, Kd_other)
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZK_(GV)), intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m].
  real, dimension(SZK_(GV)), intent(in)  :: T0     !< The initial layer temperatures [C ~> degC].
  real, dimension(SZK_(GV)), intent(in)  :: S0     !< The initial layer salinities [S ~> ppt].

  real, dimension(SZK_(GV)), intent(in)  :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZK_(GV)), intent(in)  :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)+1), intent(in)  :: TKE_source !< The forcing that is applied to the water
                                                   !! around each interface integrated over a timestep
                                                   !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: Kd     !< The diagnosed diffusivities at interfaces
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixvel !< The mixing velocity scale used in Kd
                                                   !! [Z T-1 ~> m s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixlen !< The mixing length scale used in Kd [Z ~> m].
  type(energetic_mixing_CS), intent(in)  :: CS     !< Energetic PBL control structure
  type(eMix_column_diags), intent(inout) :: eCD    !< A container for passing around diagnostics.
  real, dimension(SZK_(GV)+1), intent(in) :: Kd_other !< A diapycnal diffusivity due to other processes
                                                   !! whose energetics are dealt with elsewhere or a
                                                   !! molecular diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  ! Local variables
  real, dimension(SZK_(GV)+1) :: &
    pres_Z, &       ! Interface pressures with a rescaling factor to convert interface height
                    ! movements into changes in column potential energy [R Z2 T-2 ~> kg m-1 s-2].
    L_mix, &        ! The mixing length at an interface [H ~> m or kg m-2]
    H_int, &        ! The finite volume thicknesses associated with the interfaces [H ~> m or kg m-2]
    dt_dz_int, &    ! The timestep divided by a bounded vertical distance between layer centers [T Z-1 ~> s m-1]
    Ldt_h, &        ! The product of L_mix and dt_h [H T Z-1 ~> s or kg s m-3]
    c1_a, &         ! c1_a is used in downward passes of the tridiagonal solver for temperature [nondim].
    c1_b, &         ! c1_b is used in upward passes of the tridiagonal solver for temperature [nondim].
    hp_w2_a, &      ! An effective pivot thickness of the interface including the effects
                    ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_w2_b, &      ! An effective pivot thickness of the interface including the effects
                    ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    w_prev, &       ! The value of w after the previous iteration [Z2 T-2 ~> m2 s-2]
    w, &            ! The turblent velocity magnitude at interfaces [Z T-1 ~> m s-1]
    w_in, &         ! The turblent velocity magnitude at interfaces at the end of the
                    ! previous timestep [Z T-1 ~> m s-1]
    w2_term, &      ! A term related to the interface mass and TKE decay that scales with the
                    ! turbulent velocity squared [H ~> m or kg m-2]
    w3_term, &      ! A term related to the TKE decay that scales with the turbulent velocity
                    ! cubed [H T Z-1 ~> s or s kg m-3]
    H_int_eff, &    ! The effective thickness associated with the interface in the w2 equation,
                    ! including terms from the partial derivative of energy sinks with changes
                    ! in w**2 [H ~> m or kg m-2]
    H_from_top, &   ! The summed thicknesses between an interface and the top of the water column [H ~> m or kg m-2]
    H_from_bot, &   ! The summed thicknesses between an interface and the bottom of the water column [H ~> m or kg m-2]
    net_TKE_input   ! All the source of mixing TKE associated with an interface, including the
                    ! TKE present at the start ot the timestep [H Z2 T-2 ~> m3 s-2 or J m-2]

  real, dimension(SZK_(GV)) :: &
    dT_to_dColHt, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer [Z C-1 ~> m degC-1].
    dS_to_dColHt, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer  [Z S-1 ~> m ppt-1].
    dT_to_dPE, &    ! Partial derivatives of column potential energy with the temperature
                    ! changes within a layer, in [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE, &    ! Partial derivatives of column potential energy with the salinity changes
                    ! within a layer, in [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    dT_to_dColHt_a, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z C-1 ~> m degC-1].
    dS_to_dColHt_a, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z S-1 ~> m ppt-1].
    dT_to_dColHt_b, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers lower
                    ! in the water column [Z C-1 ~> m degC-1].
    dS_to_dColHt_b, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers lower
                    ! in the water column [Z S-1 ~> m ppt-1].
    dT_to_dPE_a, &  ! Partial derivatives of column potential energy with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE_a, &  ! Partial derivative of column potential energy with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    dT_to_dPE_b, &  ! Partial derivative of column potential energy with the temperature changes
                    ! within a layer, including the implicit effects  of mixing with layers lower
                    ! in the water column, in units of [R Z L2 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE_b, &  ! Partial derivative of column potential energy with the salinity changes
                    ! within a layer, including the implicit effects  of mixing with layers lower
                    ! in the water column, in units of [R Z L2 T-2 S-1 ~> J m-2 ppt-1].
    ! Tf, &           ! Final values of T in the column [C ~> degC].
    ! Sf, &           ! Final values of S in the column [S ~> ppt].
    hp_a, &         ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_b            ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in an upward-oriented tridiagonal solver.
  real :: Th_a      ! An effective temperature times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [C H ~> degC m or degC kg m-2].
  real :: Sh_a      ! An effective salinity times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [S H ~> ppt m or ppt kg m-2].
  real :: Th_b      ! An effective temperature times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [C H ~> degC m or degC kg m-2].
  real :: Sh_b      ! An effective salinity times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [S H ~> ppt m or ppt kg m-2].

  ! Note that the following arrays have extra (ficticious) interfaces above or below the
  ! water column for code convenience
  real, dimension(0:GV%ke+1) :: &
    w2_a        ! Running incomplete estimates of the new values of w2 in a downward pass [Z2 T-2 ~> m2 s-2]
  real, dimension(GV%ke+2) :: &
    w2_b        ! Running incomplete estimates of the new values of w2 in an upward pass [Z2 T-2 ~> m2 s-2]
  ! Note that the following arrays have extra (ficticious) layers above or below the
  ! water column for code convenience
  real, dimension(0:GV%ke) :: &
    Te_a, &         ! Running incomplete estimates of the new temperatures in a downward pass [C ~> degC]
    Se_a            ! Running incomplete estimates of the new salinities in a downward pass [S ~> ppt]
  real, dimension(GV%ke+1) :: &
    Te_b, &         ! Running incomplete estimates of the new temperatures in an upward pass [C ~> degC]
    Se_b            ! Running incomplete estimates of the new salinities in an upward pass [S ~> ppt]

  real, dimension(SZK_(GV)+1) :: &
    Kd_so_far, &    ! The total diffusivity (including Kd_other) at an interface currently being
                    ! used times a timestep divided by the average thicknesses around an
                    ! interface [H ~> m or kg m-2].
    ! Kddt_h, &       ! The total diapycnal diffusivity at an interface times a timestep divided by the
    !                 ! average thicknesses around an interface [H ~> m or kg m-2].
    Kddt_dz_other   ! The total diapycnal diffusivity at an interface due to other processes whose
                    ! energetics have been dealt with previously, times a timestep divided by the
                    ! average thicknesses around an interface [H ~> m or kg m-2].
  real, dimension(SZK_(GV)) :: &
    dt_dz, &        ! dt_dz is used to set the diffusive exchange between interfaces [T Z-1 ~> s m-1]
    Kddt_other_lay, & ! The total diapycnal diffusivity centered in a layer due to other processes whose
                    ! energetics have been dealt with previously, times a timestep divided by the
                    ! average thicknesses around an interface [H ~> m or kg m-2].
    Ldt_dzlay, &    ! The product of the mixing length in the middle of a layer and
                    ! dt_dz [H T Z-1 ~> s or kg s m-3]
    Kddt_hlay       ! The diapycnal diffusivity centered on layers times a timestep divided by the
                    ! layer thickness [H ~> m or kg m-2]
  real :: b1        ! Inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dz_neglect ! A vertical distance that is so small it is usually lost
                    ! in roundoff and can be neglected [Z ~> m].
  real :: h2_neglect ! The square of h_neglect [H2 ~> m2 or kg2 m-4]
  real :: dMass     ! The mass per unit area within a layer [Z R ~> kg m-2].
  real :: dPres     ! The hydrostatic pressure change across a layer [R Z2 T-2 ~> Pa] or
                    ! equivalently [R Z2 T-2 ~> J m-3].
  real :: w2_dgnl   ! The diagonal element in the solution for w2 [H ~> m or kg m-2]
  real :: w2_RHS    ! The right hand side, including the energy source, in the w2 equation
                    ! in units of [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: w2_err    ! The imbalance in the present iteration in the w2 equation
                    ! in units of [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: dFn_dw    ! The partial derivative of the TKE equation with turbulent velocity
                    ! in units of [H Z T-1 ~> m2 s-1 or J s m-3]

  real :: h_tot     ! The total thickness in the water column [H ~> m or kg m-2]
  real :: dztot     ! The total depth of the layers above an interface or across layers [Z ~> m]
  real :: h_huge    ! A large value compared with the thickness of the water column [H ~> m or kg m-2]
  real :: h2_topbot ! The product of the total thickness from the top and bottom [H2 ~> m2 or kg2 m-4]
  real :: H4_tot_int ! The product of the square of the total diffusivity and the square of the interior
                    ! length scale, used to help set the mixing length [H4 ~> m2 or kg4 m-8]

  real :: PE_chg_core  ! The local diffusivity invariant coefficient out front of an expression
                       ! relating the change in column potential energy from applying
                       ! Kddt_h at the present interface to that diffusivity
                       ! [H3 R Z L2 T-2 ~> J m or J kg3 m-8]
                       ! [H4 L2 T-2 ~> m6 s-2 or J kg3 m-8]
  real :: PE_chg_w0    ! The derviative of the potential energy change with w when w is 0 [R Z2 T-1 ~> J s m-3]
                       ! [H Z T-1 ~> m2 s-1 or J s m-3]
  real :: PE_chg    ! The change in potential energy due to mixing at an interface [R Z3 T-2 ~> J m-2],
                    ! [H Z2 T-2 ~> m3 s-2 or J m-2]
                    ! positive for the column increasing in potential energy (i.e., consuming TKE).
  real :: dPE_dw2   ! The partial derivative of the potential energy change with w**2 [R Z ~> J s2] [H ~> m or J s2]
!  logical :: use_Newt  ! Use Newton's method for the next guess at Kddt_h(K).

  real :: hps       ! The sum of the two effective pivot thicknesses [H ~> m or kg m-2]
  real :: bdt1      ! A product of the two pivot thicknesses plus a diffusive term [H2 ~> m2 or kg2 m-4]
  real :: lam_dt    ! The timestep times the turbulence decay rate [nondim]

  ! The following is only used for diagnostics.
  real :: I_dtdiag  !  = 1.0 / dt [T-1 ~> s-1].

  integer :: eMix_it        ! Iteration counter

  logical :: calc_Te    ! If true calculate the expected final temperature and salinity values.
  logical :: debug      ! This is used as a hard-coded value for debugging.
  logical :: converged  ! If true, the iteration for each layer of w(K) has converged to within
                        ! the specified tolerance.

  !  The following arrays are used only for debugging purposes.
!  real :: dPE_debug     ! An estimate of the potential energy change [R Z3 T-2 ~> J m-2]
!  real :: mixing_debug  ! An estimate of the rate of change of potential energy due to mixing [R Z3 T-3 ~> W m-2]
!  real, dimension(20) :: PE_chg_itt     ! The value of PE_chg after each iteration [R Z3 T-2 ~> J m-2]
!  real, dimension(20) :: Kddt_h_itt     ! The value of Kddt_h_guess after each iteration [H ~> m or kg m-2]
!  real, dimension(SZK_(GV)) :: dT_expect ! Expected temperature changes [C ~> degC]
!  real, dimension(SZK_(GV)) :: dS_expect ! Expected salinity changes [S ~> ppt]
!  integer, dimension(SZK_(GV)) :: num_itts

  integer :: k, nz, local_itt, max_itt

  nz = GV%ke

  debug = .false.  ! Change this hard-coded value for debugging.
  calc_Te = .true.

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dZ_subroundoff

  I_dtdiag = 1.0 / dt
  max_itt = 20

  lam_dt = dt*CS%decay_rate

  do K=1,nz+1 ; Kd(K) = 0.0 ; enddo

  !### Need to add code to store and set w_in
  do K=1,nz+1 ; w_in(K) = 0.0 ; enddo

  !### Optionally merge sufficiently massless layers?

  ! Zero out the temperature and salinity estimates in the extra (ficticious) layers.
  ! The actual values set here are irrelevant (so long as they are not NaNs) because they
  ! are always multiplied by a zero diffusivity reflecting the no-flux boundary condition.
  Te_a(0) = 0.0 ; Se_a(0) = 0.0 ; Te_b(nz+1) = 0.0 ; Se_b(nz+1) = 0.0

  H2_neglect = GV%H_subroundoff**2

  pres_Z(1) = 0.0
  dztot = 0.0
  h_tot = 0.0
  h_from_top(1) = 0.0
  do k=1,nz
    dMass = GV%H_to_RZ * h(k)
    dPres = GV%g_Earth_Z_T2 * dMass
    dT_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dT(k)
    dS_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dS(k)
    dT_to_dColHt(k) = dMass * dSV_dT(k)
    dS_to_dColHt(k) = dMass * dSV_dS(k)
    dztot = dztot + dz(k)
    h_tot = h_tot + h(k)

    pres_Z(K+1) = pres_Z(K) + dPres
    h_from_top(K+1) = h_from_top(K) + h(k)
  enddo
  h_from_bot(nz+1) = 0.0
  do k=nz,1,-1
    h_from_bot(K) = h_from_bot(K+1) + h(k)
  enddo

  H4_tot_int = CS%Bdry_fac*(h_tot**2*CS%L_interior**2)
  H_int(1) = max(0.5*h(1), GV%H_subroundoff)
  L_mix(1) = 0.0 ; dt_dz_int(1) = 0.0 ; Ldt_h(1) = 0.0
  w3_term(1) = 0.0
  Kddt_dz_other(1) = 0.0
  do K=2,nz
    H_int(K) = max(0.5*(h(k-1)+h(k)), GV%H_subroundoff)
    ! The following line is equivalent to
    ! L_mix(K) = sqrt(1.0 / ((1.0/CS%L_interior**2) + CS%Bdry_fac*(h_tot**2 / (H_from_top(K)**2*H_from_bot(K)**2))))
    ! which satisfies the law-of-the-wall lengthscale and goes smoothly to L_interior far from the boundaries.
    ! L_mix(K) has units of [H ~> m or kg m-2].

    h2_topbot = H_from_top(K) * H_from_bot(K)
    L_mix(K) = (CS%L_interior * h2_topbot) / max(sqrt(h2_topbot**2 + H4_tot_int), H2_neglect)
    dt_dz_int(K) = dt / (max(0.5*(dz(k-1) + dz(k)), 1e-15*dztot, GV%dz_subroundoff))
    Kddt_dz_other(K) = Kd_other(K)*dt_dz_int(K)

    Ldt_h(K) = L_mix(K)*dt_dz_int(K)
    w2_term(K) = H_int(K) * (1.0 + lam_dt)
    w3_term(K) = dt*CS%decay_w_L_scale * H_int(K) / (US%H_to_Z*L_mix(K))
  enddo
  Kddt_dz_other(nz+1) = 0.0
  H_int(nz+1) = max(0.5*h(nz), GV%H_subroundoff)
  L_mix(nz+1) = 0.0 ; dt_dz_int(nz+1) = 0.0 ; Ldt_h(nz+1) = 0.0
  w3_term(nz+1) = 0.0

  ! Setting the thicknesses associated with the top and bottom boundary conditions is one way to
  ! enforce the Direchlet boundary conditions that w2 = 0.
  h_huge = 1.0e30 * (h_tot + GV%H_subroundoff)
  h_int(1) = h_huge ; h_int(nz+1) = h_huge

  ! dt_dz is used to set the diffusive exchange between interfaces.
  do k=1,nz
    dt_dz(k) = dt / (max(dz(k), 1e-15*dztot, GV%dz_subroundoff))

    h2_topbot = (H_from_top(K)+0.5*h(k)) * (H_from_bot(K+1)+0.5*h(k))
    Ldt_dzlay(k) = dt_dz(k) * (CS%L_interior * h2_topbot) / &
                              (max(sqrt(h2_topbot**2 + H4_tot_int), H2_neglect))
  enddo

  ! These ficticious values do not matter because Kddt_dz_other must be 0 at the insulating boundaries.
  Te_a(0) = 0.0 ; Se_a(0) = 0.0 ; Te_b(nz+1) = 0.0 ; Se_b(nz+1) = 0.0

  if (CS%use_prior_diffusivity) then
    Kddt_dz_other(1) = 0.0
    do K=2,nz ; Kddt_dz_other(K) = Kd_other(K)*dt_dz_int(K) ; enddo
    Kddt_dz_other(nz+1) = 0.0

    Kddt_other_lay(1) = 0.5*(Kd_other(2)) * dt_dz(1)
    do k=2,nz-1 ; Kddt_other_lay(k) = 0.5*(Kd_other(K) + Kd_other(K+1)) * dt_dz(k) ; enddo
    Kddt_other_lay(nz) = 0.5*(Kd_other(nz)) * dt_dz(nz)

    ! Set values that are appropriate for use with a previous diffusivity.
    hp_a(1) = h(1) ; hp_b(nz) = h(nz)
    dT_to_dPE_a(1) = dT_to_dPE(1) ; dS_to_dPE_a(1) = dS_to_dPE(1)
    dT_to_dColHt_a(1) = dT_to_dColHt(1) ; dS_to_dColHt_a(1) = dS_to_dColHt(1)
    dT_to_dPE_b(nz) = dT_to_dPE(nz) ; dS_to_dPE_b(nz) = dS_to_dPE(nz)
    dT_to_dColHt_b(nz) = dT_to_dColHt(nz) ; dS_to_dColHt_b(nz) = dS_to_dColHt(nz)
    do K=2,nz
      b1 = 1.0 / (hp_a(k-1) + Kddt_dz_other(K))
      Te_a(k-1) = b1 * (h(k-1)*T0(k-1) + Kddt_dz_other(K-1)*Te_a(k-2))
      Se_a(k-1) = b1 * (h(k-1)*S0(k-1) + Kddt_dz_other(K-1)*Se_a(k-2))

      c1_a(K) = Kddt_dz_other(K) * b1
      hp_a(k) = h(k) + (hp_a(k-1) * b1) * Kddt_dz_other(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1_a(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1_a(K)*dS_to_dPE_a(k-1)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1_a(K)*dT_to_dColHt_a(k-1)
      dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1_a(K)*dS_to_dColHt_a(k-1)
    enddo

    do K=nz,2,-1  ! Loop over interior interfaces.
      b1 = 1.0 / (hp_b(k) + Kddt_dz_other(K))
      Te_b(k) = b1 * (h(k) * T0(k) + Kddt_dz_other(K+1) * Te_b(k+1))
      Se_b(k) = b1 * (h(k) * S0(k) + Kddt_dz_other(K+1) * Se_b(k+1))

      c1_b(K) = Kddt_dz_other(K) * b1
      hp_b(k-1) = h(k-1) + (hp_b(k) * b1) * Kddt_dz_other(K)
      dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1_b(K)*dT_to_dPE_b(k)
      dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1_b(K)*dS_to_dPE_b(k)
      dT_to_dColHt_b(k-1) = dT_to_dColHt(k-1) + c1_b(K)*dT_to_dColHt_b(k)
      dS_to_dColHt_b(k-1) = dS_to_dColHt(k-1) + c1_b(K)*dS_to_dColHt_b(k)
    enddo

    if (CS%tridiagonal_w2) then
      ! Set up the coefficients for an inside-out tridiagonal solver at interfaces for w2,
      ! subject to a Dirichlet boundary condition of 0 on w2 at the top and bottom.
      ! Setting hp_w2_a(1) and hp_w2_b(nz+1) to huge values gives a de-facto Dirichlet boundary condition.
      hp_w2_a(1) = H_huge ; hp_w2_b(nz+1) = H_huge
      w2_a(1) = 0.0 ; w2_b(nz+1) = 0.0
    endif

  else
    do K=1,nz+1 ; Kddt_dz_other(K) = 0.0 ; enddo
    do k=1,nz ; Kddt_other_lay(k) = 0.0 ; enddo

    ! Set up values appropriate for no background diffusivity.
    do k=1,nz
      hp_a(k) = h(k) ; hp_b(k) = h(k)
      dT_to_dPE_a(k) = dT_to_dPE(k) ; dS_to_dPE_a(k) = dS_to_dPE(k)
      dT_to_dPE_b(k) = dT_to_dPE(k) ; dS_to_dPE_b(k) = dS_to_dPE(k)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) ; dS_to_dColHt_a(k) = dS_to_dColHt(k)
      dT_to_dColHt_b(k) = dT_to_dColHt(k) ; dS_to_dColHt_b(k) = dS_to_dColHt(k)
    enddo
  endif


  do K=1,nz+1
    Kd_so_far(K) = Kddt_dz_other(K) ; w(K) = 0.0
  enddo
  do k=1,nz
    Kddt_hlay(k) = Kddt_other_lay(k) + 0.5*(w_in(K)+w_in(K+1)) * Ldt_dzlay(k)
  enddo

  ! Solve for the TKE balance at each interface [R Z3 T-2 ~> J m-2] ! This is actually [H Z2 T-2 ~> m3 s-2 or J m-2]

  ! Iterate upward and downward to capture the connection between the diffusivities at adjacent
  ! layers.  The first iteration never uses a tridiagonal solver for w**2 because the
  ! linearization of the potential energy with w**2 is not robust until the answers are close.
  do eMix_it=1,CS%Max_mixing_its

    do K=2,nz

      Th_a = h(k-1) * T0(k-1) + Kd_so_far(K-1) * Te_a(k-2)
      Sh_a = h(k-1) * S0(k-1) + Kd_so_far(K-1) * Se_a(k-2)
      Th_b = h(k) * T0(k) + Kd_so_far(K+1) * Te_b(k+1)
      Sh_b = h(k) * S0(k) + Kd_so_far(K+1) * Se_b(k+1)

      PE_chg_core = GV%RZ_to_H * find_PE_chg_core(hp_a(k-1), hp_b(k), Th_a, Sh_a, Th_b, Sh_b, &
          dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
          pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), dT_to_dColHt_b(k), dS_to_dColHt_b(k))

      ! These terms are independent of w(K):
      hps = hp_a(k-1) + hp_b(k)
      bdt1 = hp_a(k-1) * hp_b(k) + hps * Kddt_dz_other(K)

      ! kappa_dt(K) = w(K)*Ldt_h(K)

      ! Find the increase in column potential energy due to the change in the
      ! diffusivity at this interface by w(K)*L_mix(K):
      !   PE_chg = PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps)))

      ! Solve for the value of w(K) such that Fn_w = w2_RHS.
      ! Fn_w is a monotonically increasing function of positive w(K) that is 0 at w=0,
      ! while the right hand side is always positive.
      ! Fn_w = w2_dgnl*w(K)**2 + PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps)))

      if ((eMix_it > 1) .and. (CS%tridiagonal_w2)) then
        ! This turns the expression into the pivot for a tridiagonal solver for w2.
        w2_dgnl = (w2_term(K) + (Kddt_hlay(k)*hp_w2_b(K+1) / (hp_w2_b(K+1) + Kddt_hlay(k)) + &
                                 Kddt_hlay(k-1)*hp_w2_a(K-1) / (hp_w2_a(K-1) + Kddt_hlay(k-1))) )
        w2_RHS = (TKE_source(K) + H_int(K)*(w_in(K)**2)) + (Kddt_hlay(k)*w2_b(K+1) + Kddt_hlay(k-1)*w2_a(K-1))
      else
        ! This is the pivot for a simple local solution.
        w2_dgnl = (w2_term(K) + (Kddt_hlay(k) + Kddt_hlay(k-1)))
        w2_RHS = (TKE_source(K) + H_int(K)*(w_in(K)**2)) + (Kddt_hlay(k) * w(K+1)**2 + Kddt_hlay(k-1)*w(K-1)**2)
      endif

      w_prev(K) = w(K)

      if (w(K) <= 0.0) then
        ! For the first guess, ignore the cubic terms, giving an underestimate.
        PE_chg_w0 = PE_chg_core * Ldt_h(K) / (bdt1 * bdt1)
        !  w2_dgnl*w(K)**2 + PE_chg_w0 * w(K) - w2_RHS = 0.0
        !  For accuracy with the relevant root, avoid subtraction by solving for 1/w(K) and inverting:
        w(K) = 2.0*w2_RHS / (PE_chg_w0 + sqrt(PE_chg_w0**2 + 4.0*w2_RHS*w2_dgnl))
        ! Note that this agrees with the Osborn relation limit, in which w2_dgnl is small and
        ! w2_RHS ~= TKE_source(K), in which case this expression is  w(K) = TKE_source(K) / PE_chg_w0.
      endif

      ! Iterate with Newton's method, either in w or w2, depending on which is more nearly linear.
      do local_itt=1,max_itt
        w2_err = w2_RHS - ( (w2_dgnl + w3_term(K)*w(K))*w(K)**2 + &
                            PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps))) )

        ! For positive w, dFn_dw is always positive:
        dFn_dw = (2.0*w2_dgnl + 3.0*w3_term(K)*w(K))*w(K) + &
                 PE_chg_core * Ldt_h(K) / (bdt1 + w(K)*Ldt_h(K) * hps)**2

        if (PE_chg_core * Ldt_h(K) / (bdt1 + w(K)*Ldt_h(K) * hps)**2 >= &
            (2.0*w2_dgnl + 3.0*w3_term(K)*w(K))*w(K) ) then

          ! Use Newton's method in w, because the linear term is dominant.
          w(K) = w(K) + w2_err / dFn_dw

          ! The sign of the second derivative determines whether Newton's method converges from above or below.
          ! d2Fn_dw_2 = (2.0*w2_dgnl + 6.0*w3_term(K)*w(K)) - &
          !             2.0*PE_chg_core * Ldt_h(K)**2 * hps / (bdt1 + w(K)*Ldt_h(K) * hps)**3
          ! if ((w2_dgnl+3.0*w3_term(K)*w(K))*(bdt1 + w(K)*Ldt_h(K) * hps)**3 >= PE_chg_core * Ldt_h(K)**2 * hps) then
          !   ! At this value of w(K), the Newton's method solution will converge monotonically from above
          ! else
          !   ! At this value of w(K), the Newton's method solution will converge monotonically from below
          ! endif
        else
          ! Use Newton's method in w**2, because the quadratic term is dominant.
          ! dFn_dw2 = dFn_dw / (2.0*w(K))
          !  w(K) = sqrt(w(K)**2 + w2_err / dFn_dw2)
          w(K) = sqrt(w(K)**2 + 2.0 * w(K) * w2_err / dFn_dw)
        endif

        if (abs(w2_err) < CS%w_tol * dFn_dw) exit  ! Stop after this iteration.
      enddo

      ! verify =  -H_int(K)*(w(K)**2 - w_in(K)**2) + TKE_source(K) - PE_chg - &
      !     (lam_dt*H_int(K)*(w(K)**2) + w3_term(K)*w(K)**3) + &
      !     ((L_dz(k)*(w(K+1) + w(K)) + Kddt_hlay(k)) * (w(K+1)**2 - w(K)**2) + &
      !      (L_dz(k-1)*(w(K-1) + w(K)) + Kddt_hlay(k-1))*(w(K-1)**2 - w(K)**2))
      ! Folding find_PE_chg into the expression above gives a quartic equation for w(K) that
      ! can be solved iteratively.

      !   If the layer-center diffusivities are already known from a previous iteration, solve for w(K):
      ! ((w2_term(K) + w3_term(K)*w(K)) + (Kddt_hlay(k) + Kddt_hlay(k-1)))*w(K)**2 + &
      !    PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps))) = &
      !       (TKE_source(K) + H_int(K)*(w_in(K)**2)) + (Kddt_hlay(k) * w(K+1)**2 + Kddt_hlay(k-1)*w(K-1)**2)

      if (CS%tridiagonal_w2) then ! Because of the change in w(K), hp_w2_a(K) and w2_a(K) need to be updated.
        PE_chg = PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps)))
        dPE_dw2 = 0.5*PE_chg_core * Ldt_h(K) / (w(K) * (bdt1 + w(K)*Ldt_h(K) * hps)**2)
        net_TKE_input(K) = h_int(K)*w_in(K)**2 + (TKE_source(K) - (PE_chg - dPE_dw2*w(K)**2))
        H_int_eff(K) = (w2_term(K) + 1.5*w3_term(K)*w(K)) + dPE_dw2

        b1 = 1.0 / (hp_w2_a(K-1) + Kddt_hlay(k-1))
        hp_w2_a(K) = H_int_eff(K) + (hp_w2_a(K-1) * b1)*Kddt_hlay(k-1)

        b1 = 1.0 / (hp_w2_a(K) + Kddt_hlay(k))
        ! w2_a includes a linearization of the potential energy change.
        w2_a(K) = b1 * (net_TKE_input(K) + Kddt_hlay(k-1)*w2_a(K-1))
      endif

      Kd_so_far(K) = Kddt_dz_other(K) + w(K)*Ldt_h(K)

      b1 = 1.0 / (hp_a(k-1) + Kd_so_far(K))
      c1_a(K) = Kd_so_far(K) * b1

      Te_a(k-1) = b1 * (h(k-1) * T0(k-1) + Kd_so_far(K-1) * Te_a(k-2))
      Se_a(k-1) = b1 * (h(k-1) * S0(k-1) + Kd_so_far(K-1) * Se_a(k-2))

      hp_a(k) = h(k) + (hp_a(k-1) * b1) * Kd_so_far(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1_a(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1_a(K)*dS_to_dPE_a(k-1)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1_a(K)*dT_to_dColHt_a(k-1)
      dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1_a(K)*dS_to_dColHt_a(k-1)
    enddo

    ! Update the layer-centered diffusivities used for the w2 equation.
    do k=1,nz
      Kddt_hlay(k) = 0.25*((w_in(K)+w_in(K+1)) + (w(K)+w(K+1))) * Ldt_dzlay(k) + Kddt_other_lay(k)
    enddo

    if (CS%tridiagonal_w2) then
      ! We need to reestimate w2_b, hp_w2_b, w2_a and hp_w2_a because they depend on Kddt_hlay.
      hp_w2_a(1) = H_huge
      w2_a(1) = 0.0
      do K=2,nz-1
        b1 = 1.0 / (hp_w2_a(K-1) + Kddt_hlay(k-1))
        hp_w2_a(K) = H_int_eff(K) + (hp_w2_a(K-1) * b1)*Kddt_hlay(k-1)

        b1 = 1.0 / (hp_w2_a(K) + Kddt_hlay(k))
        ! w2_a includes a linearization of the potential energy change.
        w2_a(K) = b1 * (net_TKE_input(K) + Kddt_hlay(k-1)*w2_a(K-1))
      enddo
      hp_w2_b(nz+1) = H_huge
      w2_b(nz+1) = 0.0
      do K=nz,3,-1
        b1 = 1.0 / (hp_w2_b(K+1) + Kddt_hlay(k))
        hp_w2_b(K) = H_int_eff(K) + (hp_w2_b(K+1) * b1)*Kddt_hlay(k)

        b1 = 1.0 / (hp_w2_b(K) + Kddt_hlay(k-1))
        ! w2_b includes a linearization of the potential energy change.
        w2_b(K) = b1 * (net_TKE_input(K) + Kddt_hlay(k)*w2_b(K+1))
      enddo
    endif

    converged = .true.  ! This may be reset to false in the loop below.

    ! Calculate the dependencies on layers below.
    do K=nz,2,-1  ! Loop over interior interfaces.
      ! This block repeats the calculation in the downward block above.
      ! First calculate some terms that are independent of w(K).

      Th_a = h(k-1) * T0(k-1) + Kd_so_far(K-1) * Te_a(k-2)
      Sh_a = h(k-1) * S0(k-1) + Kd_so_far(K-1) * Se_a(k-2)
      Th_b = h(k) * T0(k) + Kd_so_far(K+1) * Te_b(k+1)
      Sh_b = h(k) * S0(k) + Kd_so_far(K+1) * Se_b(k+1)

      PE_chg_core = GV%RZ_to_H * find_PE_chg_core(hp_a(k-1), hp_b(k), Th_a, Sh_a, Th_b, Sh_b, &
          dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
          pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), dT_to_dColHt_b(k), dS_to_dColHt_b(k))

      ! These terms are independent of w(K):
      hps = hp_a(k-1) + hp_b(k)
      bdt1 = hp_a(k-1) * hp_b(k) + hps * Kddt_dz_other(K)

      ! Store the current value of w to test for convergence.
      w_prev(K) = w(K)

      ! Solve for an updated value of w(K) such that Fn_w = w2_RHS.
      ! Fn_w = w2_dgnl*w(K)**2 + w3_term(K)*w(K)**3 + &
      !        PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps)))
      if ((eMix_it > 1) .and. (CS%tridiagonal_w2)) then
        ! This turns the expression into the pivot for a tridiagonal solver for w2.
        w2_dgnl = w2_term(K) + (Kddt_hlay(k)*hp_w2_b(K+1) / (hp_w2_b(K+1) + Kddt_hlay(k)) + &
                                Kddt_hlay(k-1)*hp_w2_a(K-1) / (hp_w2_a(K-1) + Kddt_hlay(k-1)))
        w2_RHS = (TKE_source(K) + H_int(K)*(w_in(K)**2)) + (Kddt_hlay(k) * w2_b(K+1) + Kddt_hlay(k-1)*w2_a(K-1))
      else
        ! This is the simpler local solution.
        w2_dgnl = w2_term(K) + (Kddt_hlay(k) + Kddt_hlay(k-1))
        w2_RHS = (TKE_source(K) + H_int(K)*(w_in(K)**2)) + (Kddt_hlay(k) * w(K+1)**2 + Kddt_hlay(k-1)*w(K-1)**2)
      endif

      ! Iterate with Newton's method, either in w or w2, depending on which is more nearly linear.
      do local_itt=1,max_itt
        w2_err = w2_RHS - ( (w2_dgnl + w3_term(K)*w(K)) * w(K)**2 + &
                 PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps))) )

        ! dFn_dw is always positive:
        dFn_dw = (2.0*w2_dgnl + 3.0*w3_term(K)*w(K))*w(K) + &
                 PE_chg_core * Ldt_h(K) / (bdt1 + w(K)*Ldt_h(K) * hps)**2

        if (PE_chg_core * Ldt_h(K) / (bdt1 + w(K)*Ldt_h(K) * hps)**2 >= &
            (2.0*w2_dgnl + 3.0*w3_term(K)*w(K))*w(K) ) then
          ! Use Newton's method in w, because the linear term is dominant.
          w(K) = w(K) + w2_err / dFn_dw
        else
          ! Use Newton's method in w**2, because the quadratic term is dominant.
          w(K) = sqrt(w(K)**2 + 2.0 * w(K) * w2_err / dFn_dw)
        endif

        if (abs(w2_err) < CS%w_tol * dFn_dw) exit  ! Stop after this iteration.
      enddo

      if (CS%tridiagonal_w2) then ! Because of the change in w(K), hp_w2_b(K) and w2_b(K) need to be updated.
        PE_chg = PE_chg_core * (w(K)*Ldt_h(K) / (bdt1 * (bdt1 + w(K)*Ldt_h(K) * hps)))
        dPE_dw2 = 0.5*PE_chg_core * Ldt_h(K) / (w(K) * (bdt1 + w(K)*Ldt_h(K) * hps)**2)
        net_TKE_input(K) = h_int(K)*w_in(K)**2 + (TKE_source(K) - (PE_chg - dPE_dw2*w(K)**2))
        H_int_eff(K) = (w2_term(K) + 1.5*w3_term(K)*w(K)) + dPE_dw2

        b1 = 1.0 / (hp_w2_b(K+1) + Kddt_hlay(k))
        hp_w2_b(K) = H_int_eff(K) + (hp_w2_b(K+1) * b1)*Kddt_hlay(k)

        b1 = 1.0 / (hp_w2_b(K) + Kddt_hlay(k-1))
        ! w2_b includes a linearization of the potential energy change.
        w2_b(K) = b1 * (net_TKE_input(K) + Kddt_hlay(k)*w2_b(K+1))
      endif

      Kd_so_far(K) = Kddt_dz_other(K) + w(K)*Ldt_h(K)

      ! Use the change in w(K) between iterations to determine convergence at this interface.
      ! Stop when all interfaces have converged.  This may be more aggressive than necessary.
      if (abs(w(K) - w_prev(K)) > CS%w_tol) converged = .false.

      b1 = 1.0 / (hp_b(k) + Kd_so_far(K))
      c1_b(K) = Kd_so_far(K) * b1

      Te_b(k) = b1 * (h(k) * T0(k) + Kd_so_far(K+1) * Te_b(k+1))
      Se_b(k) = b1 * (h(k) * S0(k) + Kd_so_far(k+1) * Se_b(k+1))

      hp_b(k-1) = h(k-1) + (hp_b(k) * b1) * Kd_so_far(K)
      dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1_b(K)*dT_to_dPE_b(k)
      dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1_b(K)*dS_to_dPE_b(k)
      dT_to_dColHt_b(k-1) = dT_to_dColHt(k-1) + c1_b(K)*dT_to_dColHt_b(k)
      dS_to_dColHt_b(k-1) = dS_to_dColHt(k-1) + c1_b(K)*dS_to_dColHt_b(k)

    enddo

    if (converged) exit

    ! Update the layer-centered diffusivities used for the w2 equation.
    do k=1,nz
      Kddt_hlay(k) = 0.25*((w_in(K)+w_in(K+1)) + (w(K)+w(K+1))) * Ldt_dzlay(k) + Kddt_other_lay(k)
    enddo

    if (CS%tridiagonal_w2) then
      ! We need to reestimate w2_b, hp_w2_b, w2_a and hp_w2_a because they depend on Kddt_hlay.
      hp_w2_a(1) = H_huge
      w2_a(1) = 0.0
      do K=2,nz-1
        b1 = 1.0 / (hp_w2_a(K-1) + Kddt_hlay(k-1))
        hp_w2_a(K) = H_int_eff(K) + (hp_w2_a(K-1) * b1)*Kddt_hlay(k-1)

        b1 = 1.0 / (hp_w2_a(K) + Kddt_hlay(k))
        ! w2_a includes a linearization of the potential energy change.
        w2_a(K) = b1 * (net_TKE_input(K) + Kddt_hlay(k-1)*w2_a(K-1))
      enddo
      hp_w2_b(nz+1) = H_huge
      w2_b(nz+1) = 0.0
      do K=nz,3,-1
        b1 = 1.0 / (hp_w2_b(K+1) + Kddt_hlay(k))
        hp_w2_b(K) = H_int_eff(K) + (hp_w2_b(K+1) * b1)*Kddt_hlay(k)

        b1 = 1.0 / (hp_w2_b(K) + Kddt_hlay(k-1))
        ! w2_b includes a linearization of the potential energy change.
        w2_b(K) = b1 * (net_TKE_input(K) + Kddt_hlay(k)*w2_b(K+1))
      enddo
    endif

    ! Now update the other layer working down from the top to determine the
    ! final temperatures and salinities.  This is only needed for debugging.
!    b1 = 1.0 / (hp_b(1))
!    Tf(1) = b1 * (h(1) * T0(1) + Kd_so_far(2) * Te_b(2))
!    Sf(1) = b1 * (h(1) * S0(1) + Kd_so_far(2) * Se_b(2))
!    do k=2,nz
!      Tf(k) = Te_b(k) + c1_b(K)*Tf(k-1)
!      Sf(k) = Se_b(k) + c1_b(K)*Sf(k-1)
!    enddo

  enddo ! eMix_it

  ! Optionally interpolate Kd back to the original layers, filling in massless values.

  ! The diffusivities go quadratically to zero at the Direchlet boundaries, as both
  ! w(z) and L_mix(Z) go to zero.  This can be accomplished by linearly interpolating w(z)
  ! between the interfaces of the coarsened grid to the interfaces of the new grid
  ! and using the known expressions for L_mix(z).

  ! Store the final diffusivities, stting the no-flux boundary conditions:
  Kd(1) = 0.0 ; mixvel(1) = 0.0 ; mixlen(1) = 0.0
  do K=2,nz
    ! The diagnosed added diffusivities at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
    Kd(K) = w(K)*L_mix(K)
    ! Change Kd to be the diagnosed total diffusivities at interfaces
    ! Kd(K) = Kd_other(K) + Kd(K)
    mixvel(K) = w(K)      ! w(K) is in [Z T-1 ~> m s-1]
    mixlen(K) = GV%H_to_Z*L_mix(K)  ! L_mix(K) is in [H ~> m or kg m-2]
  enddo
  Kd(nz+1) = 0.0 ; mixvel(nz+1) = 0.0 ; mixlen(nz+1) = 0.0

!  This would be better in non-Boussinesq mode.
!  if (GV%Boussinesq) then
!    do K=1,nz+1 ; mixlen(K) = GV%H_to_Z * L_mix(K) ; enddo
!  else
!    mixlen(1) = 0.0
!    do K=2,nz
!      mixlen(K) =  L_mix(K) * ((dz(k-1) + dz(k) + dz_neglect) / (h(k-1) + h(k) + h_neglect))
!    enddo
!    mixlen(nz+1) = 0.0
!  endif

!    if (debug) then
!      ! Complete the tridiagonal solve for Te.
!      b1 = 1.0 / hp_a(nz)
!      Te(nz) = b1 * (h(nz) * T0(nz) + Kd_so_far(nz) * Te(nz-1))
!      Se(nz) = b1 * (h(nz) * S0(nz) + Kd_so_far(nz) * Se(nz-1))
!      dT_expect(nz) = Te(nz) - T0(nz) ; dS_expect(nz) = Se(nz) - S0(nz)
!      do k=nz-1,1,-1
!        Te(k) = Te(k) + c1(K+1)*Te(k+1)
!        Se(k) = Se(k) + c1(K+1)*Se(k+1)
!        dT_expect(k) = Te(k) - T0(k) ; dS_expect(k) = Se(k) - S0(k)
!      enddo
!    endif

!    if (debug) then
!      dPE_debug = 0.0
!      do k=1,nz
!        dPE_debug = dPE_debug + (dT_to_dPE(k) * (Tf(k) - T0(k)) + &
!                                 dS_to_dPE(k) * (Sf(k) - S0(k)))
!      enddo
!      mixing_debug = dPE_debug * I_dtdiag
!    endif

  eCD%eMix_its = min(eMix_it, CS%max_mixing_its)

end subroutine eMix_column



!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interface's diapycnal diffusivity times a timestep.
function find_PE_chg_core(hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, pres_Z, &
                       dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b) result(PE_chg_core)
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface below [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The rescaled hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].

  real :: PE_chg_core           !< The local diffusivity invariant coefficient out front of an expression
                                !! relating the change in column potential energy from applying
                                !! Kddt_h at the present interface to that diffusivity
                                !! [H3 R Z L2 T-2 ~> J m or J kg3 m-8]

  ! Local variables
  real :: dT_c ! The core term in the expressions for the temperature changes [C H2 ~> degC m2 or degC kg2 m-4].
  real :: dS_c ! The core term in the expressions for the salinity changes [S H2 ~> ppt m2 or ppt kg2 m-4].
  real :: PEc_core ! The diffusivity-independent core term in the expressions
                   ! for the potential energy changes [H3 R Z L2 T-2 ~> J m or J kg3 m-8]
  real :: ColHt_core ! The diffusivity-independent core term in the expressions
                     ! for the column height changes [H3 Z ~> m4 or kg3 m-5].

  !   The expression for the change in potential energy used here is derived from the
  ! expressions for the final estimates of the changes in temperature and salinities, and then
  ! extensively manipulated to get it into its most succinct form. The derivation is not
  !  necessarily obvious, but it demonstrably works by comparison with separate calculations
  ! of the energy changes after the tridiagonal solver for the final changes in temperature
  ! and salinity are applied.  It also appears in find_PE_chg in MOM_diapyc_energy_req.F90.
  !   It is then restricted for use here to only consider stably stratified profiles.

  dT_c = hp_a * Th_b - hp_b * Th_a
  dS_c = hp_a * Sh_b - hp_b * Sh_a
  ! The max and min here treat convectively unstable profiles as though they are neutral.
  PEc_core = max((hp_b * (dT_to_dPE_a * dT_c + dS_to_dPE_a * dS_c) - &
                  hp_a * (dT_to_dPE_b * dT_c + dS_to_dPE_b * dS_c)), 0.0)
  ColHt_core = min((hp_b * (dT_to_dColHt_a * dT_c + dS_to_dColHt_a * dS_c) - &
                    hp_a * (dT_to_dColHt_b * dT_c + dS_to_dColHt_b * dS_c)), 0.0)
  PE_chg_core = (PEc_core - pres_Z * ColHt_core)

  ! Find the increase in column potential energy due to the change in the diffusivity at this
  ! interface (converted to dKddt_h [H ~> m or kg m-2]) with the following expressions:
  !  hps = hp_a + hp_b
  !  bdt1 = hp_a * hp_b + Kddt_h0 * hps
  !  PE_chg = PE_chg_core * (dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps)))

end function find_PE_chg_core


!> This subroutine initializes the energetic_mixing module
subroutine energetic_mixing_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(energetic_mixing_CS), pointer     :: CS   !< Energetic PBL control structure

  ! Local variables
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_energetic_mixing"  ! This module's name.
  character(len=120) :: diff_text ! A clause describing parameter setting that differ.
  logical :: debug
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) then ; allocate(CS)
  else ; return ; endif

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")

!/1. General eMix settings
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "MIXING_EN_DEBUG", CS%debug, &
                 "If true, write out verbose debugging data for energetic_mixing.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "MIXING_EN_TRIDIAGONAL_W2", CS%tridiagonal_w2, &
                 "If true, use a linearized tridiagonal solver for w**2 after the first "//&
                 "pair of iterations in the energetic_mixing diffusivity calculation.", &
                 default=.false.)
  call get_param(param_file, mdl, "MIXING_EN_PRIOR_DIFF", CS%use_prior_diffusivity, &
                 "If true, include a background diffusivity from other processes in the "//&
                 "energetic_mixing diffusivity calculation.", default=.false.)
  call get_param(param_file, mdl, "MIXING_EN_L_INTERIOR", CS%L_interior, &
                 "The mixing distance far from the edges in the "//&
                 "energetic_mixing diffusivity calculation.", &
                 units="m", default=10.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%vonKar, &
                 'The value the von Karman constant as used for mixed layer viscosity.', &
                 units='nondim', default=0.41)
  call get_param(param_file, mdl, "MIXING_EN_BOUNDARY_SCALE_FAC", CS%Bdry_fac, &
                 "A scaling factor for the mixing length distance near the top and bottom "//&
                 "boundaries, perhaps the Von Karman constant.", &
                 units="nondim", default=CS%vonKar)
  call get_param(param_file, mdl, "MIXING_EN_W_TOLERANCE", CS%w_tol, &
                 "The tolerance for convergence of the iterations for the turbulent velocity "//&
                 "in the turbulent kinetic energy equation", &
                 units="m s-1", default=1.0e-7, scale=US%m_to_Z*US%T_to_s)
  call get_param(param_file, mdl, "MIXING_EN_DECAY_RATE", CS%decay_rate, &
                 "A background rate of energy decay in the turbulent kinetic energy equation", &
                 units="s-1", default=1.0e-5, scale=US%T_to_s)
  call get_param(param_file, mdl, "MIXING_EN_W_L_DECAY", CS%decay_w_L_scale, &
                 "A coefficient relating the turbulent velocity divided by the local mixing "//&
                 "distance to a turbulent kinetic energy decay rate", &
                  units="nondim", default=0.0)
  call get_param(param_file, mdl, "EMIX_MAX_ITS", CS%max_mixing_its, &
                 "The maximum number of iterations that can be used to find a self-consistent "//&
                 "diffusivity profile with the energetics-based interior mixing.", &
                 default=20)

  !/ Options for documenting differences from parameter choices
  call get_param(param_file, mdl, "EMIX_OPTIONS_DIFF", CS%options_diff, &
                 "If positive, this is a coded integer indicating a pair of settings whose "//&
                 "differences are diagnosed in a passive diagnostic mode via extra calls to "//&
                 "eMix_column.  If this is 0 or negative no extra calls occur.", &
                 default=0)
  if (CS%options_diff > 0) then
    if (CS%options_diff == 1) then
      diff_text = "eMix_ORIGINAL_PE_CALC settings"
    elseif (CS%options_diff == 2) then
      diff_text = "eMix_ANSWER_DATE settings"
    elseif (CS%options_diff == 3) then
      diff_text = "DIRECT_eMix_MIXING_CALC settings"
    else
      diff_text = "unchanged settings"
    endif
  endif


!/ Checking output flags
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'eMix_TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_Mixing_Length = register_diag_field('ocean_model', 'eMix_Mixing_Length', diag%axesTi, &
      Time, 'Mixing Length that is used', units='m', conversion=US%Z_to_m)
  CS%id_Velocity_Scale = register_diag_field('ocean_model', 'eMix_Velocity_Scale', diag%axesTi, &
      Time, 'Velocity Scale that is used.', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_frac_en_diff = register_diag_field('ocean_model', 'eMix_frac_en_diff', diag%axesT1, &
      Time, 'The fractional difference between the energy used and the energy input', units='nondim')

  if (CS%options_diff > 0) then
    CS%id_opt_diff_Kd_eMix = register_diag_field('ocean_model', 'eMix_opt_diff_Kd_eMix', diag%axesTi, &
        Time, 'Change in eMix diapycnal diffusivity at interfaces due to '//trim(diff_text), &
        units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_opt_maxdiff_Kd_eMix = register_diag_field('ocean_model', 'eMix_opt_maxdiff_Kd_eMix', diag%axesT1, &
        Time, 'Column maximum change in eMix diapycnal diffusivity at interfaces due to '//trim(diff_text), &
        units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
  endif

  if (report_avg_its) then
    CS%sum_its(1) = real_to_EFP(0.0) ; CS%sum_its(2) = real_to_EFP(0.0)
  endif

  CS%TKE_diagnostics = (max(CS%id_TKE_mixing, CS%id_TKE_forcing) > 0)

end subroutine energetic_mixing_init

!> Clean up and deallocate memory associated with the energetic_mixing module.
subroutine energetic_mixing_end(CS)
  type(energetic_mixing_CS), pointer :: CS !< energetic_mixing control structure

  character(len=256) :: mesg
  real :: avg_its ! The averaged number of iterations used by eMix [nondim]

  if (.not.associated(CS)) return

  if (report_avg_its) then
    call EFP_sum_across_PEs(CS%sum_its, 2)
    avg_its = EFP_to_real(CS%sum_its(1)) / EFP_to_real(CS%sum_its(2))
    write (mesg,*) "Average eMix iterations = ", avg_its
    call MOM_mesg(mesg)
  endif
  deallocate(CS)
end subroutine energetic_mixing_end

! The following extra code includes various working tridiagonal solvers that were useful in
! developing or debugging this code.

!> A simple downward-first tridiagonal solver for Tf
subroutine tridiag_T_down(T, Kddt_h, h, GV, Tf)
  type(verticalGrid_type),     intent(in)  :: GV     !< The ocean's vertical grid structure
  real, dimension(SZK_(GV)),   intent(in)  :: T      !< The initial layer temperatures [C ~> degC]
  real, dimension(SZK_(GV)+1), intent(in)  :: Kddt_h !< The diapycnal diffusivity times a timestep
                                                     !! divided by the average thicknesses around a
                                                     !! layer [H ~> m or kg m-2]
  real, dimension(SZK_(GV)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)),   intent(out) :: Tf     !< The updated layer temperatures [C ~> degC]

  ! Local variables
  real, dimension(GV%ke) :: &
    hp_a, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    c1_a, &     ! c1_a is used by a downward-oriented tridiagonal solver [nondim].
    h_tr        ! h_tr is h at tracer points with a h_neglect added to
                ! ensure positive definiteness [H ~> m or kg m-2].
  ! Note that the following array has extra (ficticious) layers above the
  ! water column for code convenience.
  real, dimension(0:GV%ke) :: &
    Te_a        ! Running incomplete estimates of the new temperatures in a downward pass [C ~> degC]
  real :: b1    ! Inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  integer :: K, nz

  ! This subroutine is verifiably equivalent to
  !   call tridiag_T_cent(T, Kddt_h, h, GV, GV%ke, Tf)

  nz = GV%ke
  do k=1,nz ; h_tr(k) = max(h(k), GV%H_subroundoff) ; enddo

  ! This value does not matter because Kddt_h must be 0 at the insulating boundaries.
  Te_a(0) = 0.0

  hp_a(1) = h_tr(1)
  do K=2,nz
    b1 = 1.0 / (hp_a(k-1) + Kddt_h(K))
    Te_a(k-1) = b1 * (h_tr(k-1)*T(k-1) + Kddt_h(K-1)*Te_a(k-2))

    c1_a(K) = Kddt_h(K) * b1
    hp_a(k) = h_tr(k) + (hp_a(k-1) * b1)*Kddt_h(K)
  enddo
  b1 = 1.0 / (hp_a(nz)) ! (+ Kddt_h(nz+1) == 0.)
  Tf(nz) = b1 * (h_tr(nz)*T(nz) + Kddt_h(nz)*Te_a(nz-1))
  do k=nz-1,1,-1
    Tf(k) = Te_a(k) + c1_a(K+1)*Tf(k+1)
  enddo

end subroutine tridiag_T_down

!> A simple upward-first tridiagonal solver for Tf
subroutine tridiag_T_up(T, Kddt_h, h, GV, Tf)
  type(verticalGrid_type),     intent(in)  :: GV     !< The ocean's vertical grid structure
  real, dimension(SZK_(GV)),   intent(in)  :: T      !< The initial layer temperatures [C ~> degC]
  real, dimension(SZK_(GV)+1), intent(in)  :: Kddt_h !< The diapycnal diffusivity times a timestep
                                                     !! divided by the average thicknesses around a
                                                     !! layer [H ~> m or kg m-2]
  real, dimension(SZK_(GV)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)),   intent(out) :: Tf     !< The updated layer temperatures [C ~> degC]

  ! Local variables
  real, dimension(GV%ke) :: &
    hp_b, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    c1_b, &     ! c1_b is used by an upward-oriented tridiagonal solver [nondim].
    h_tr        ! h_tr is h at tracer points with a h_neglect added to
                ! ensure positive definiteness [H ~> m or kg m-2].
  ! Note that the following array has extra (ficticious) layers below the
  ! water column for code convenience.
  real, dimension(GV%ke+1) :: &
    Te_b        ! Running incomplete estimates of the new temperatures in an upward pass [C ~> degC]
  real :: b1    ! Inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  integer :: K, nz

  ! This subroutine is verifiably equivalent to
  !   call tridiag_T_cent(T, Kddt_h, h, GV, 1, Tf)

  nz = GV%ke
  do k=1,nz ; h_tr(k) = max(h(k), GV%H_subroundoff) ; enddo

  ! This value does not matter because Kddt_h must be 0 at the insulating boundaries.
  Te_b(nz+1) = 0.0

  hp_b(nz) = h_tr(nz)
  do K=nz,2,-1
    b1 = 1.0 / (hp_b(k) + Kddt_h(K))
    Te_b(k) = b1 * (h_tr(k)*T(k) + Kddt_h(K+1)*Te_b(k+1))

    c1_b(K) = Kddt_h(K) * b1
    hp_b(k-1) = h_tr(k-1) + (hp_b(k) * b1)*Kddt_h(K)
  enddo
  b1 = 1.0 / (hp_b(1)) ! (+ Kddt_h(1) == 0.)
  Tf(1) = b1 * (h_tr(1)*T(1) + Kddt_h(2)*Te_b(2))
  do k=2,nz
    Tf(k) = Te_b(k) + c1_b(K)*Tf(k-1)
  enddo

end subroutine tridiag_T_up

!> This is a tridiagonal solver centered on interface kl_cent for updating
!! temperatures due to diffusion.
subroutine tridiag_T_cent(T, Kddt_h, h, GV, k_cent, Tf)
  type(verticalGrid_type),     intent(in)  :: GV     !< The ocean's vertical grid structure
  real, dimension(SZK_(GV)),   intent(in)  :: T      !< The initial layer temperatures [C ~> degC]
  real, dimension(SZK_(GV)+1), intent(in)  :: Kddt_h !< The diapycnal diffusivity times a timestep
                                                     !! divided by the average thicknesses around a
                                                     !! layer [H ~> m or kg m-2]
  real, dimension(SZK_(GV)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  integer,                     intent(in)  :: k_cent !< The layer about which the tridiagonal calculation is centered.
  real, dimension(SZK_(GV)),   intent(out) :: Tf     !< The updated layer temperatures [C ~> degC]

  ! Local variables
  real, dimension(GV%ke) :: &
    hp_a, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_b, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    c1_a, &     ! c1_a is used by a downward-oriented tridiagonal solver [nondim].
    c1_b, &     ! c1_b is used by an upward-oriented tridiagonal solver [nondim].
    h_tr        ! h_tr is h at tracer points with a h_neglect added to
                ! ensure positive definiteness [H ~> m or kg m-2].
  ! Note that the following arrays have extra (ficticious) layers above or below the
  ! water column for code convenience
  real, dimension(0:GV%ke) :: &
    Te_a        ! Running incomplete estimates of the new temperatures in a downward pass [C ~> degC]
  real, dimension(GV%ke+1) :: &
    Te_b        ! Running incomplete estimates of the new temperatures in an upward pass [C ~> degC]
  real :: b1    ! Inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  real :: hp_ab ! A limited version of the product of pivot thicknessness plus the diffusive
                ! coupling distance for the layers and interfaces above and below layer
                ! kl_cent with the layer-centered tridiagonal solver [H2 ~> m2 or kg2 m-4].
  integer :: K, kl_cent, nz

  nz = GV%ke
  kl_cent = max(min(k_cent, nz), 1)
  do k=1,nz ; h_tr(k) = max(h(k), GV%H_subroundoff) ; enddo

  ! These values do not matter because Kddt_h must be 0 at the insulating boundaries.
  Te_a(0) = 0.0 ; Te_b(nz+1) = 0.0

  ! Inside out tridiagonal solver for Tf, centered on layer kl_cent
  hp_a(1) = h_tr(1)
  do K=2,kl_cent
    b1 = 1.0 / (hp_a(k-1) + Kddt_h(K))
    Te_a(k-1) = b1 * (h_tr(k-1)*T(k-1) + Kddt_h(K-1)*Te_a(k-2))

    c1_a(K) = Kddt_h(K) * b1
    hp_a(k) = h_tr(k) + (hp_a(k-1) * b1)*Kddt_h(K)
  enddo
  hp_b(nz) = h_tr(nz)
  do K=nz,kl_cent+1,-1
    b1 = 1.0 / (hp_b(k) + Kddt_h(K))
    Te_b(k) = b1 * (h_tr(k)*T(k) + Kddt_h(K+1)*Te_b(k+1))

    c1_b(K) = Kddt_h(K) * b1
    hp_b(k-1) = h_tr(k-1) + (hp_b(k) * b1)*Kddt_h(K)
  enddo
  k = kl_cent
  if ((Kddt_h(K) == 0.0) .and. (Kddt_h(K+1) == 0.0)) then
    b1 = 1.0 / h_tr(k)
    ! Equivalent to Tf(k) = T(k)
  elseif (Kddt_h(K) == 0.0) then
    b1 = (hp_b(k+1) + Kddt_h(K+1)) / ((hp_b(k+1) + Kddt_h(K+1)) * h_tr(k) + ( Kddt_h(K+1)*hp_b(k+1) ) )
  elseif (Kddt_h(K+1) == 0.0) then
    b1 = (hp_a(k-1) + Kddt_h(K)) / ((hp_a(k-1) + Kddt_h(K)) * h_tr(k) + ( Kddt_h(K)*hp_a(k-1) ) )
  else
    ! The derivation of the general expression for b1 at the central layer can be found
    ! in diapyc_energy_req_calc.
    hp_ab = max((hp_b(k+1) + Kddt_h(K+1)) * (hp_a(k-1) + Kddt_h(K)), GV%H_subroundoff**2)
    b1 = hp_ab / (hp_ab * h_tr(k) + ( Kddt_h(K+1)*hp_b(k+1) * (hp_a(k-1) + Kddt_h(K)) + &
                                      Kddt_h(K)*hp_a(k-1) * (hp_b(k+1) + Kddt_h(K+1)) ) )
  endif
  ! This is the back-substitute to find the final value of Tf.
  Tf(k) = b1 * (h_tr(k)*T(k) + (Kddt_h(K+1)*Te_b(k+1) + Kddt_h(K)*Te_a(k-1)))
  do k=kl_cent-1,1,-1
    Tf(k) = Te_a(k) + c1_a(K+1)*Tf(k+1)
  enddo
  do k=kl_cent+1,nz
    Tf(k) = Te_b(k) + c1_b(K)*Tf(k-1)
  enddo

end subroutine tridiag_T_cent

!> This is a tridiagonal solver centered on interface K_int_cent for updating
!! temperatures due to diffusion.
subroutine tridiag_T_interface_cent(T0, Kddt_h, h, GV, K_int_cent, Tf)
  type(verticalGrid_type),     intent(in)  :: GV     !< The ocean's vertical grid structure
  real, dimension(SZK_(GV)),   intent(in)  :: T0     !< The initial layer temperatures [C ~> degC]
  real, dimension(SZK_(GV)+1), intent(in)  :: Kddt_h !< The diapycnal diffusivity times a timestep
                                                     !! divided by the average thicknesses around a
                                                     !! layer [H ~> m or kg m-2]
  real, dimension(SZK_(GV)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  integer,                     intent(in)  :: K_int_cent !< The interface about which the tridiagonal
                                                     !! calculation is centered.
  real, dimension(SZK_(GV)),   intent(out) :: Tf     !< The updated layer temperatures [C ~> degC]

  ! Local variables
  real, dimension(GV%ke) :: &
    hp_a, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_b, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    c1_a, &     ! c1_a is used by a downward-oriented tridiagonal solver [nondim].
    c1_b, &     ! c1_b is used by an upward-oriented tridiagonal solver [nondim].
    h_tr        ! h_tr is h at tracer points with a h_neglect added to
                ! ensure positive definiteness [H ~> m or kg m-2].
  ! Note that the following arrays have extra (ficticious) layers above or below the
  ! water column for code convenience
  real, dimension(0:GV%ke) :: &
    Te_a        ! Running incomplete estimates of the new temperatures in a downward pass [C ~> degC]
  real, dimension(GV%ke+1) :: &
    Te_b        ! Running incomplete estimates of the new temperatures in an upward pass [C ~> degC]
  real :: Th_a  ! An effective temperature times a thickness in the layer above, including implicit
                ! mixing effects with other yet higher layers [C H ~> degC m or degC kg m-2].
  real :: Th_b  ! An effective temperature times a thickness in the layer below, including implicit
                ! mixing effects with other yet lower layers [C H ~> degC m or degC kg m-2].
  real :: b1    ! Inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  integer :: K, K_cent, nz

  nz = GV%ke
  K_cent = max(min(K_int_cent, nz), 2)
  do k=1,nz ; h_tr(k) = max(h(k), GV%H_subroundoff) ; enddo

  ! These values do not matter because Kddt_h must be 0 at the insulating boundaries.
  Te_a(0) = 0.0 ; Te_b(nz+1) = 0.0

  ! This is the update for the temperatures due to a change in diffusivity at interface K_cent:
  Te_a(0) = 0.0 ; Te_b(nz+1) = 0.0 ! We must have ; Kddt_h(1) = 0.0 ; Kddt_h(nz+1) = 0.0
  do K=2,K_cent-1  ! Loop downward over interior interfaces.
    b1 = 1.0 / (hp_a(k-1) + Kddt_h(K))
    c1_a(K) = Kddt_h(K) * b1
    Te_a(k-1) = b1 * (h(k-1) * T0(k-1) + Kddt_h(K-1) * Te_a(k-2))
  enddo

  do K=nz,K_cent+1,-1  ! Loop over interior interfaces.
    b1 = 1.0 / (hp_b(k) + Kddt_h(K))
    c1_b(K) = Kddt_h(K) * b1
    Te_b(k) = b1 * (h(k) * T0(k) + Kddt_h(K+1) * Te_b(k+1))
  enddo

  K = K_cent
  Th_a = h(k-1) * T0(k-1) + Kddt_h(K-1) * Te_a(k-2)
  Th_b = h(k) * T0(k) + Kddt_h(K+1) * Te_b(k+1)

  b1 = 1.0 / (hp_a(k-1)*hp_b(k) + Kddt_h(K)*(hp_a(k-1) + hp_b(k)))
  Tf(k-1) = ((hp_b(k) + Kddt_h(K)) * Th_a + Kddt_h(K) * Th_b ) * b1
  Tf(k) = (Kddt_h(K) * Th_a + (hp_a(k-1) + Kddt_h(K)) * Th_b ) * b1

  c1_a(K) = Kddt_h(K) / (hp_a(k-1) + Kddt_h(K))
  c1_b(K) = Kddt_h(K) / (hp_b(k) + Kddt_h(K))

  ! Now update the other layer working outward from k_cent to determine the final
  ! final temperatures.
  do k=K_cent-2,1,-1
    Tf(k) = Te_a(k) + c1_a(K+1)*Tf(k+1)
  enddo
  do k=K_cent+1,nz
    Tf(k) = Te_b(k) + c1_b(K)*Tf(k-1)
  enddo

end subroutine tridiag_T_interface_cent


!> This is a tridiagonal solver centered on interface K_int_cent for diffusively updating
!! a variable w2 with Dirichlet boundary conditions of 0 at the top and bottom.
subroutine tridiag_w2_cent(w2_in, Kddt_hlay, h, GV, K_int_cent, w2f, debug)
  type(verticalGrid_type),     intent(in)  :: GV     !< The ocean's vertical grid structure
  real, dimension(SZK_(GV)+1), intent(in)  :: w2_in  !< The initial value of w2 in arbitrary units [A ~> a]
  real, dimension(SZK_(GV)),   intent(in)  :: Kddt_hlay !< The diapycnal diffusivity centered on layers times
                                                     !! a timestep divided by the layer thickness [H ~> m or kg m-2]
  real, dimension(SZK_(GV)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  integer,                     intent(in)  :: K_int_cent !< The interface about which the tridiagonal
                                                     !! calculation is centered.
  real, dimension(SZK_(GV)+1), intent(out) :: w2f    !< The updated layer values of w2 [A ~> a]
  logical,           optional, intent(in)  :: debug  !< If present and true, double check the solutions for correctness.

  ! Local variables
  real, dimension(GV%ke+1) :: &
    h_int, &    ! The finite volume thicknesses associated with the interfaces [H ~> m or kg m-2]
    hp_w2_a, &  ! An effective pivot thickness of the interface including the effects
                ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_w2_b, &  ! An effective pivot thickness of the interface including the effects
                ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    c_w2_a, &   ! c_w2_a is used by a downward-oriented tridiagonal solver [nondim].
    c_w2_b, &   ! c_w2_b is used by an upward-oriented tridiagonal solver [nondim].
    error       ! The error at each interface of the final solution [H A ~> a m or a kg m-2]
  ! Note that the following arrays have extra (ficticious) interfaces above or below the
  ! water column for code convenience
  real, dimension(0:GV%ke+1) :: &
    w2_a        ! Running incomplete estimates of the new values of w2 in a downward pass [A ~> a]
  real, dimension(GV%ke+2) :: &
    w2_b        ! Running incomplete estimates of the new values of w2 in an upward pass [A ~> a]
  real :: sum_abs_terms ! The sum of the absolute values of all the terms in the equation for
                ! determining whether the errors are at the level of roundoff [H A ~> a m or a kg m-2]
  real :: b1    ! Inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  real :: hp_w2_ab ! A limited version of the product of pivot thicknessness plus the diffusive
                ! coupling distance for the layers and interfaces above and below layer
                ! kl_cent with the layer-centered tridiagonal solver [H2 ~> m2 or kg2 m-4].
  integer :: K, KI_cent, nz

  nz = GV%ke
  ! KI_cent = max(min(K_int_cent, nz+1), 1)
  ! With spacified top and bottom values, the permitted range for KI_cent is reduced.
  KI_cent = max(min(K_int_cent, nz), 2)

  ! Set the effective finite volume thicknesses associated with the interfaces.
  h_int(1) = max(0.5*h(1), GV%H_subroundoff) ; h_int(nz+1) = max(0.5*h(nz), GV%H_subroundoff)
  do K=2,nz ; h_int(K) = max(0.5*(h(k-1) + h(k)), GV%H_subroundoff) ; enddo

  ! This block of code is an inside-out tridiagonal solver for w2 at interfaces,
  ! subject to a Dirichlet boundary condition of 0 on w2 at the top and bottom.
  ! h_int(1), hp_w2_a(1), h_int(nz+1), hp_w2_b(nz+1), c_w2_b(nz) and c_w2_a(1) are not used.
  ! Alternately, h_int(1), hp_w2_a(1), h_int(nz+1) and hp_w2_b(nz+1) could be set to huge values.
  w2_a(1) = 0.0
  hp_w2_a(2) = h_int(2) + Kddt_hlay(1)
  do k=2,KI_cent-1
    b1 = 1.0 / (hp_w2_a(K) + Kddt_hlay(k))
    w2_a(K) = b1 * (h_int(K)*w2_in(K) + Kddt_hlay(k-1)*w2_a(K-1))

    c_w2_a(k) = Kddt_hlay(k) * b1
    hp_w2_a(K+1) = h_int(K+1) + (hp_w2_a(K) * b1)*Kddt_hlay(k)
  enddo
  ! hp_w2_b(nz+1) = h_int(nz+1) = HUGE
  w2_b(nz+1) = 0.0
  hp_w2_b(nz) = h_int(nz) + Kddt_hlay(nz)
  do k=nz-1,KI_cent,-1
    b1 = 1.0 / (hp_w2_b(K+1) + Kddt_hlay(k))
    w2_b(K+1) = b1 * (h_int(K+1)*w2_in(K+1) + Kddt_hlay(k+1)*w2_b(K+2))

    c_w2_b(k) = Kddt_hlay(k) * b1
    hp_w2_b(K) = h_int(K) + (hp_w2_b(K+1) * b1)*Kddt_hlay(k)
  enddo
  if ((KI_cent >= 2) .and. (KI_cent <= nz)) then
    K = KI_cent
    if (((Kddt_hlay(k-1) == 0.0) .and. (Kddt_hlay(k) == 0.0)) .or. (nz == 2)) then
      ! Equivalent to w2f(k) = w2_in(k)
      ! or K == 2 and hp_w2_a(K-1) = inf and hp_w2_b(K+1) = inf
      b1 = 1.0 / (h_int(K) + (Kddt_hlay(k-1) + Kddt_hlay(k)) )
    elseif (K == 2) then   ! hp_w2_a(K-1) = inf
      b1 = (hp_w2_b(K+1) + Kddt_hlay(k)) / &
           ((hp_w2_b(K+1) + Kddt_hlay(k)) * (h_int(K) + Kddt_hlay(k-1)) + (Kddt_hlay(k)*hp_w2_b(K+1)) )
    elseif (K == nz) then   ! hp_w2_b(K+1) = inf
      b1 = (hp_w2_a(K-1) + Kddt_hlay(k-1)) / &
           ((hp_w2_a(K-1) + Kddt_hlay(k-1)) * (h_int(K) + Kddt_hlay(k)) + (Kddt_hlay(k-1)*hp_w2_a(K-1)) )
    elseif (Kddt_hlay(k-1) == 0.0) then
      b1 = (hp_w2_b(K+1) + Kddt_hlay(k)) / &
           ((hp_w2_b(K+1) + Kddt_hlay(k)) * h_int(K) + ( Kddt_hlay(k)*hp_w2_b(K+1) ) )
    elseif (Kddt_hlay(k) == 0.0) then
      b1 = (hp_w2_a(K-1) + Kddt_hlay(k-1)) / &
           ((hp_w2_a(K-1) + Kddt_hlay(k-1)) * h_int(K) + ( Kddt_hlay(k-1)*hp_w2_a(K-1) ) )
    else  ! This is the general case.  The other cases above can be obtained by setting various terms
          ! to be zero or infinite and cancelling out terms accordingly.
      ! b1 = 1.0 / (h_int(K) + (Kddt_hlay(k)*hp_w2_b(K+1) / (hp_w2_b(K+1) + Kddt_hlay(k)) + &
      !                         Kddt_hlay(k-1)*hp_w2_a(K-1) / (hp_w2_a(K-1) + Kddt_hlay(k-1))) )
      hp_w2_ab = (hp_w2_b(K+1) + Kddt_hlay(k)) * (hp_w2_a(K-1) + Kddt_hlay(k-1))
      b1 = hp_w2_ab / (hp_w2_ab * h_int(K) + ( Kddt_hlay(k)*hp_w2_b(K+1) * (hp_w2_a(K-1) + Kddt_hlay(k-1)) + &
                                               Kddt_hlay(k-1)*hp_w2_a(K-1) * (hp_w2_b(K+1) + Kddt_hlay(k)) ) )
    endif
    w2f(K) = b1 * (h_int(K)*w2_in(K) + (Kddt_hlay(k)*w2_b(K+1) + Kddt_hlay(k-1)*w2_a(K-1)))
  endif

  ! Re-impose Direchlet boundary conditions at the top and bottom.
  w2f(1) = 0.0
  w2f(nz+1) = 0.0
  do K=KI_cent-1,2,-1
    w2f(K) = w2_a(K) + c_w2_a(k)*w2f(K+1)
  enddo
  do k=KI_cent+1,nz
    w2f(K) = w2_b(K) + c_w2_b(k-1)*w2f(K-1)
  enddo

  ! For verification:
  if (present(debug)) then ; if (debug) then ; do K=2,nz
    error(K) = h_int(K)*(w2f(K) - w2_in(K)) - Kddt_hlay(k-1) * (w2f(K-1) - w2f(K)) + &
               Kddt_hlay(k) * (w2f(K) - w2f(K+1))
    sum_abs_terms = abs(h_int(K)*w2f(K)) + abs(h_int(K)*w2_in(K)) + &
                    abs(Kddt_hlay(k-1)*w2f(K-1)) + abs(Kddt_hlay(k-1)*w2f(K)) + &
                    abs(Kddt_hlay(k)*w2f(K)) + abs(Kddt_hlay(k)*w2f(K+1))
    if (error(K) > 1.0e-14*sum_abs_terms) &
      call MOM_error(FATAL, "Bad tridiagonal solver for w2f.")
  enddo ; endif ; endif

end subroutine tridiag_w2_cent

!> \namespace MOM_energetic_mixing
!!
!! By Robert Hallberg, 2026.
!!
!!   This file contains the subroutine (energetic_mixing) that uses the implicit
!! mixing calculation to convert the tubulent kinetic energy input into a diapycal
!! diffusivity.

end module MOM_mixing_energetics
