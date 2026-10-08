! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  BOX module for coupling with UKCA containing a wrapper subroutine for
!  setting up a specified configuration of the UKCA model.
!    Branched from atmos_ukca_setup_mod from UM version 12.1
!
! Method:
!
!  1. Call 'ukca_setup' to configure UKCA according to UM namelist inputs
!     and other configuration data.
!  2. Check that the UKCA configuration created is consistent with the
!     UM inputs. This protects against the possible omission of
!     configuration settings from the 'ukca_setup' call and/or any
!     unexpected interpretation of these settings by UKCA.
!  3. Retrieve values of any internal UKCA configuration variables that
!     are required by the UM.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE box_ukca_setup_mod

USE ukca_api_mod, ONLY:                                                        &
  ukca_setup, ukca_get_config,                                                 &
  ukca_maxlen_message, ukca_maxlen_procname,                                   &
  ukca_chem_off,                                                               &
  ukca_chem_trop,                                                              &
  ukca_chem_raq,                                                               &
  ukca_chem_offline_be,                                                        &
  ukca_chem_tropisop,                                                          &
  ukca_chem_strattrop,                                                         &
  ukca_chem_strat,                                                             &
  ukca_chem_cristrat,                                                          &
  ukca_int_method_be_explicit,                                                 &
  ukca_int_method_nr,                                                          &
  ukca_age_reset_by_level,                                                     &
  ukca_age_reset_by_height,                                                    &
  ukca_strat_lbc_off,                                                          &
  ukca_strat_lbc_wmoa1,                                                        &
  ukca_strat_lbc_env,                                                          &
  ukca_activation_arg,                                                         &
  ukca_activation_jones

USE ukca_option_mod, ONLY:                                                     &
  i_ukca_chem,                                                                 &
  i_ageair_reset_method,                                                       &
  max_ageair_reset_level,                                                      &
  max_ageair_reset_height,                                                     &
  l_ukca_chem_aero,                                                            &
  l_ukca_mode,                                                                 &
  l_ukca_ageair,                                                               &
  l_ukca_emissions_off,                                                        &
  l_fix_tropopause_level,                                                      &
  l_ukca_drydep_off,                                                           &
  l_ukca_wetdep_off,                                                           &
  l_tracer_lumping,                                                            &
  i_ukca_chem_version,                                                         &
  nrsteps,                                                                     &
  chem_timestep,                                                               &
  dts0,                                                                        &
  l_ukca_asad_columns,                                                         &
  l_ukca_asad_full,                                                            &
  l_ukca_debug_asad,                                                           &
  l_ukca_intdd,                                                                &
  l_ukca_ddepo3_ocean,                                                         &
  l_ukca_ddep_lev1,                                                            &
  l_ukca_dry_dep_so2wet,                                                       &
  nit,                                                                         &
  l_ukca_quasinewton,                                                          &
  i_ukca_quasinewton_start,                                                    &
  i_ukca_quasinewton_end,                                                      &
  max_z_for_offline_chem,                                                      &
  i_ukca_topboundary,                                                          &
  l_ukca_ro2_ntp,                                                              &
  l_ukca_ro2_perm,                                                             &
  l_ukca_persist_off,                                                          &
  l_ukca_het_psc,                                                              &
  i_ukca_hetconfig,                                                            &
  l_ukca_limit_nat,                                                            &
  l_ukca_sa_clim,                                                              &
  l_ukca_trophet,                                                              &
  l_ukca_classic_hetchem,                                                      &
  i_ukca_photol,                                                               &
  fastjx_mode,                                                                 &
  fastjx_numwl,                                                                &
  fastjx_prescutoff,                                                           &
  i_ukca_solcyc,                                                               &
  i_ukca_solcyc_start_year,                                                    &
! l_environ_jo2,                                                               &
! l_environ_jo2b,                                                              &
  l_ukca_ibvoc,                                                                &
  l_ukca_inferno,                                                              &
  l_ukca_inferno_ch4,                                                          &
  i_inferno_emi,                                                               &
  l_ukca_so2ems_expvolc,                                                       &
  l_ukca_so2ems_plumeria,                                                      &
  l_ukca_qch4inter,                                                            &
  l_ukca_emsdrvn_ch4,                                                          &
  mode_parfrac,                                                                &
  i_ukca_dms_flux,                                                             &
  l_ukca_scale_seadms_ems,                                                     &
  seadms_ems_scaling,                                                          &
  l_ukca_scale_sea_salt_ems,                                                   &
  sea_salt_ems_scaling,                                                        &
  l_ukca_scale_marine_pom_ems,                                                 &
  marine_pom_ems_scaling,                                                      &
  l_ukca_linox_scaling,                                                        &
  lightnox_scale_fac,                                                          &
  i_ukca_light_param,                                                          &
  l_ukca_scale_soa_yield_mt,                                                   &
  soa_yield_scaling_mt,                                                        &
  l_ukca_scale_soa_yield_isop,                                                 &
  soa_yield_scaling_isop,                                                      &
  l_ukca_h2o_feedback,                                                         &
  l_ukca_conserve_h,                                                           &
  l_ukca_set_trace_gases,                                                      &
  l_ukca_prescribech4,                                                         &
  l_ukca_chem,                                                                 &
  l_ukca_trop,                                                                 &
  l_ukca_raq,                                                                  &
  l_ukca_raqaero,                                                              &
  l_ukca_offline_be,                                                           &
  l_ukca_tropisop,                                                             &
  l_ukca_strattrop,                                                            &
  l_ukca_strat,                                                                &
  l_ukca_offline,                                                              &
  l_ukca_stratcfc,                                                             &
  ukca_int_method,                                                             &
  i_mode_nzts,                                                                 &
  i_mode_setup,                                                                &
  l_mode_bhn_on,                                                               &
  l_mode_bln_on,                                                               &
  i_mode_bln_param_method,                                                     &
  i_mode_nucscav,                                                              &
  mode_activation_dryr,                                                        &
  mode_incld_so2_rfrac,                                                        &
  l_ukca_plume_scav,                                                           &
  l_ukca_primss,                                                               &
  l_ukca_primsu,                                                               &
  l_ukca_primdu,                                                               &
  l_ukca_primbcoc,                                                             &
  l_ukca_prim_moc,                                                             &
  l_bcoc_bf,                                                                   &
  l_bcoc_bm,                                                                   &
  l_bcoc_ff,                                                                   &
  l_ukca_scale_biom_aer_ems,                                                   &
  biom_aer_ems_scaling,                                                        &
  l_ukca_fine_no3_prod,                                                        &
  l_ukca_coarse_no3_prod,                                                      &
  l_no3_prod_in_aero_step,                                                     &
  l_ukca_radaer,                                                               &
  i_ukca_tune_bc,                                                              &
  i_ukca_activation_scheme,                                                    &
  i_ukca_nwbins,                                                               &
  l_ukca_sfix,                                                                 &
  i_ukca_scenario,                                                             &
  i_ukca_scenario_um,                                                          &
  i_ukca_scenario_rcp,                                                         &
  i_ukca_scenario_wmoa1

USE photol_api_mod,  ONLY: photol_setup, photol_get_config, photol_off,        &
                           photol_strat_only, photol_2d, photol_fastjx

USE submodel_mod, ONLY: atmos_sm
USE nlsizes_namelist_mod, ONLY:  global_row_length,                            &
    row_length, rows, model_levels, bl_levels, ntiles, n_cca_lev
USE nlstgen_mod, ONLY: secs_per_periodim, steps_per_periodim
USE nlstcall_mod, ONLY: lcal360, ltimer
USE jules_surface_types_mod, ONLY: ntype, npft,                                &
    brd_leaf, brd_leaf_dec, brd_leaf_eg_trop, brd_leaf_eg_temp,                &
    ndl_leaf, ndl_leaf_dec, ndl_leaf_eg,                                       &
    c3_grass, c3_crop, c3_pasture,                                             &
    c4_grass, c4_crop, c4_pasture,                                             &
    shrub, shrub_dec, shrub_eg,                                                &
    urban, lake, soil, ice, elev_ice
USE jules_soil_mod, ONLY: dzsoil_io
USE jules_sea_seaice_mod, ONLY: l_ctile
USE run_aerosol_mod, ONLY:                                                     &
    l_sulpc_so2, l_soot, l_ocff, l_use_seasalt_sulpc, l_sulpc_dms
USE rad_input_mod, ONLY:                                                       &
    l_use_arclsulp, l_use_biogenic, l_use_seasalt_indirect,                    &
    l_use_seasalt_direct
USE mphys_inputs_mod, ONLY: l_use_seasalt_autoconv
USE bl_option_mod, ONLY: l_conv_tke
USE carbon_options_mod, ONLY: l_co2_interactive
USE cv_run_mod, ONLY: l_param_conv, l_3d_cca
USE ozone_inputs_mod, ONLY: zon_av_ozone
USE hybrid_control_mod, ONLY: i_activate_strip_ukca, ip_activate_run_in_snr,   &
    l_junior, l_senior
USE tuning_segments_mod, ONLY: ukca_chem_seg_size, ukca_mode_seg_size
USE science_fixes_mod, ONLY: l_fix_improve_drydep,                             &
    l_fix_ukca_h2dd_x, l_fix_drydep_so2_water, l_fix_ukca_offox_h2o_fac,       &
    l_fix_ukca_h2so4_ystore, l_fix_neg_pvol_wat,                               &
    l_fix_ukca_impscav, l_fix_nacl_density, l_fix_ukca_activate_pdf,           &
    l_fix_ukca_activate_vert_rep, l_improve_aero_drydep,                       &
    l_fix_ukca_hygroscopicities, l_fix_ukca_water_content

! Parameters for initialising photolysis
USE rad_pcf,  ONLY: ip_aerosol_param_moist, ip_accum_sulphate,                 &
                    ip_aitken_sulphate
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2

!++SAN add print statements
USE umPrintMgr, ONLY: Printstatus, ummessage, umprint, PrStatus_Diag
!--SAN

! check if using components used with UM or not - should be F for box model
USE ukca_um_legacy_mod, ONLY: l_um_infrastructure

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='BOX_UKCA_SETUP_MOD'

PUBLIC :: box_ukca_setup

! Configuration variables without direct UM equivalents

INTEGER, PARAMETER :: nlev_ent_tr_mix = 3
                                      ! Number of grid levels for
                                      ! entrainment-related fields used by
                                      ! tr_mix
INTEGER :: i_strat_lbc_source         ! Source for gas MMR values for lower
                                      ! boundary conditions in stratospheric
                                      ! chemistry schemes
INTEGER, ALLOCATABLE :: i_elev_ice(:) ! Indices for elevated ice surface type

REAL :: timestep                      ! Model time step (seconds)
REAL :: sigwmin                       ! Lower limit for std. dev. of updraft
                                      ! p.d.f. in Activate

LOGICAL, PARAMETER :: l_support_ems_vertprof = .TRUE.
                                      ! True to support full range of vertical
                                      ! scaling options for 2D emission fields.
LOGICAL, PARAMETER :: l_support_ems_gridbox_units = .TRUE.
                                      ! True to support the provision of
                                      ! offline emissions in grid-box units
!s LOGICAL, PARAMETER :: l_enable_diag_um = .TRUE.
!++SAN STASH diagnostics always off in box model
LOGICAL, PARAMETER :: l_enable_diag_um = .FALSE.
                                      ! True to enable UM diagnostics via STASH
LOGICAL :: l_use_classic_so4          ! True to use CLASSIC SO4 for
                                      ! heterogeneous chemistry
LOGICAL :: l_use_classic_soot         ! True to use CLASSIC black carbon for
                                      ! het. chem.
LOGICAL :: l_use_classic_ocff         ! True to use CLASSIC organic carbon from
                                      ! fossil fuels for het. chem.
LOGICAL :: l_use_classic_biogenic     ! True to use CLASSIC biogenic secondary
                                      ! organic aerosol for het. chem.
LOGICAL :: l_use_classic_seasalt      ! True to use CLASSIC sea salt for het.
                                      ! chem.
LOGICAL :: l_bug_repro_tke_index      ! True to reproduce results of an old bug
                                      ! in Activate related to incorrect TKE
                                      ! indexing
LOGICAL :: l_ntpreq_n_activ_sum       ! True to include total number
                                      ! concentration of active aerosol in
                                      ! non-transported prognostics
LOGICAL :: l_ntpreq_dryd_nuc_sol      ! True to include the dry diameter
                                      ! for nucleation soluble mode.
LOGICAL :: l_deposition_jules         ! True to use JULES-based interactive
                                      ! dry deposition routines
                                      ! ALWAYS FALSE FOR BOX MODEL
LOGICAL  :: l_ukca_strat_chem = .FALSE. ! for Photolysis, determine if UKCA
                                        ! is using a Stratospheric scheme

! Generic interface for subroutines to check an individual UKCA configuration
! value against a UM equivalent - overloaded according to the data type
INTERFACE check_config_value
  MODULE PROCEDURE check_config_integer
  MODULE PROCEDURE check_config_integer_vec
  MODULE PROCEDURE check_config_real
  MODULE PROCEDURE check_config_logical
END INTERFACE check_config_value

INTERFACE print_setting 
  MODULE PROCEDURE print_setting_integer
  MODULE PROCEDURE print_setting_integer_vec
  MODULE PROCEDURE print_setting_real
  MODULE PROCEDURE print_setting_logical
END INTERFACE print_setting 

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE box_ukca_setup()
! ----------------------------------------------------------------------
! Description:
!
! Set up the UKCA configuration based on the values of UM variables,
! check the consistency of the configuration created against the input
! values and retrieve any internal UKCA configuration variables
! required by the UM. Also set some UM variables that are dependent on
! the UKCA configuration.
! ----------------------------------------------------------------------

USE ukca_scavenging_mod,     ONLY: ukca_mode_scavcoeff
USE ukca_cdnc_mod,           ONLY: ukca_cdnc_init

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Local variables

INTEGER, PARAMETER :: n_elev_ice = 10  ! No. of elevated ice types in the 27
                                       ! surface type configuration

INTEGER :: i_photol_scheme             ! Photolysis scheme in use

! UKCA error reporting variables
CHARACTER(LEN=ukca_maxlen_message)  :: ukca_errmsg    ! Error return message
CHARACTER(LEN=ukca_maxlen_procname) :: ukca_errproc   ! Routine in which error
                                                      ! was trapped
! UM error reporting variables
INTEGER                             :: errcode=0      ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength)   :: cmessage       ! Error return message
CHARACTER(LEN=*), PARAMETER :: errproc_suffix = ' in UKCA'
                                                      ! Text to append to
                                                      ! routine name

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='BOX_UKCA_SETUP'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!++SAN
CALL umPrint('SAN: in BOX_UKCA_SETUP',src=RoutineName)
! Hardcode ntype and ntiles - will be needed for the dummy interactive dry dep scheme
ntype = 9
ntiles = 9
npft = 5

brd_leaf   = 1
ndl_leaf   = 2
c3_grass   = 3
c4_grass   = 4
shrub      = 5
urban      = 6
lake       = 7
soil       = 8
ice        = 9

!++SAN If using column call in box model, need to hardcode seg size
IF (l_ukca_asad_columns) THEN
  ukca_chem_seg_size = 1
END IF
!--SAN

! Check consistency of dry deposition scheme and number of tiles
IF ((ntiles /= ntype) .AND. l_ukca_intdd) THEN
  errcode = 1
  WRITE(cmessage,'(A,A,I4,A,I4)') 'Cannot use interactive dry dep.',           &
    ' You have ', ntiles, ' surface tiles but should have ', ntype
  CALL ereport(RoutineName,errcode,cmessage)
END IF

!++SAN
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: ntiles = ', ntiles
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: ntype = ', ntype
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: npft = ', npft
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: l_ukca_intdd = ', l_ukca_intdd
CALL umPrint(umMessage,src=RoutineName)
!-SAN

! Set elevated ice surface type array for 27 type configuration only
! (An allocatable array is required for UKCA API calls)
IF (ntype == 27) THEN
  ALLOCATE(i_elev_ice(n_elev_ice))
  i_elev_ice = elev_ice(1:n_elev_ice)
END IF

! Set model timestep
! (equivalent to UM variable 'secs_per_stepim' set later in 'inittime' routine)
timestep = REAL(secs_per_periodim(atmos_sm)) /                                 &
           REAL(steps_per_periodim(atmos_sm))

! Set source for gas MMR values for lower boundary conditions in stratospheric
! chemistry schemes. For scenario=um and =rcp values are read from namelists
! and external files respectively, and further collated/ lumped to form the
! actual lower boundary conditions. The scenario=wmoa1 uses values currently
! hard-wired in UKCA.

SELECT CASE(i_ukca_scenario)
CASE (i_ukca_scenario_um, i_ukca_scenario_rcp)
  i_strat_lbc_source = ukca_strat_lbc_env
CASE (i_ukca_scenario_wmoa1)
  i_strat_lbc_source = ukca_strat_lbc_wmoa1
CASE DEFAULT
  i_strat_lbc_source = ukca_strat_lbc_off
END SELECT


!++SAN - checking this bit, due to error
WRITE(cmessage,*) 'SAN: BOX_UKCA_SETUP - LBC source variables'
CALL umPrint(cmessage, src=RoutineName)
WRITE(cmessage,*) 'SAN: i_strat_lbc_source ', i_strat_lbc_source
CALL umPrint(cmessage, src=RoutineName)
WRITE(cmessage,*) 'SAN: i_ukca_scenario_rcp ', i_ukca_scenario_rcp
CALL umPrint(cmessage, src=RoutineName)
WRITE(cmessage,*) 'SAN: i_ukca_scenario_wmoa1 ', i_ukca_scenario_wmoa1
CALL umPrint(cmessage, src=RoutineName)
WRITE(cmessage,*) 'SAN: ukca_strat_lbc_wmoa1 ', ukca_strat_lbc_wmoa1
CALL umPrint(cmessage, src=RoutineName)
WRITE(cmessage,*) 'SAN: ukca_strat_lbc_off ', ukca_strat_lbc_off
CALL umPrint(cmessage, src=RoutineName)
WRITE(cmessage,*) 'SAN: i_ukca_scenario_um ', i_ukca_scenario_um
CALL umPrint(cmessage, src=RoutineName)
!--SAN


! Set environment controls for use of CLASSIC data.
! Data are used if available when the relevant heterogeneous chemistry option is
! selected.
l_use_classic_so4 = (l_ukca_het_psc .OR. l_ukca_classic_hetchem) .AND.         &
                    (l_sulpc_so2 .OR. l_use_arclsulp)
l_use_classic_soot = l_ukca_classic_hetchem .AND. l_soot
l_use_classic_ocff = l_ukca_classic_hetchem .AND. l_ocff
l_use_classic_biogenic = l_ukca_classic_hetchem .AND. l_use_biogenic
l_use_classic_seasalt = l_ukca_classic_hetchem .AND.                           &
                        (l_use_seasalt_autoconv .OR. l_use_seasalt_sulpc .OR.  &
                         l_use_seasalt_indirect .OR. l_use_seasalt_direct)

! Control the calculation of the std. dev. of the updraft p.d.f. in Activate
! which is based on the TKE diagnostic (stashcode 03473).
IF (l_conv_tke) THEN
   ! New calculation, used when the inclusion of convective effects in the TKE
   ! diagnostic is turned on.
  l_bug_repro_tke_index = .FALSE.
  sigwmin = 0.01
ELSE
   ! Preserve effects of TKE-related bug in updraft p.d.f. spread and
   ! also use the old value for the std. dev. lower limit.
   ! The bug was caused by incorrect indexing of the TKE input passed to UKCA
   ! which has now been fixed but its effects can still be simulated to
   ! preserve results.
  l_bug_repro_tke_index = .TRUE.
  sigwmin = 0.1
END IF

! Include total number concentration of active aerosol in NTPs if running
! hybrid resolution model
l_ntpreq_n_activ_sum = l_senior .OR. l_junior

! Include the dry diameter for nucleation solube mode if calculating
! activated aerosol in Senior component of hybrid resolution model
l_ntpreq_dryd_nuc_sol = (i_activate_strip_ukca == ip_activate_run_in_snr)

! Set l_deposition_jules for JULES with atmospheric deposition
! For coupled UM_JULES applications, the JULES-based atmospheric deposition
! routines can only be called from the UKCA routine ukca_chemistry_ctl
! (or its equivalents, ukca_chemistry_ctl_BE and ukca_chemistry_ctl_col).
! The switches l_deposition and l_deposition_from_ukca both need to be true.
! This is to demonstrate that the JULES-based routines give outputs that
! bit compare with UKCA known good outputs.
! Currently always .FALSE. for the box model as we do not use JULES
l_deposition_jules = .FALSE.

! Determine if UKCA is using one of the 'Stratospheric' Chemistry schemes
l_ukca_strat_chem = ( i_ukca_chem == ukca_chem_strat     .OR.                  &
                      i_ukca_chem == ukca_chem_strattrop .OR.                  &
                      i_ukca_chem == ukca_chem_cristrat )

! Initialise the photolysis configuration -done before UKCA since some
! environ fields and photolysis diagnostics are still controlled by UKCA based
! on the selected photolysis scheme
CALL photol_setup(i_photol_scheme=i_ukca_photol,                               &
                  error_code=errcode,                                          &
                  global_row_length=global_row_length,                         &
                  l_cal360=lcal360,                                            &
                  n_cca_lev=n_cca_lev,                                         &
                  timestep=timestep,                                           &
                  chem_timestep=chem_timestep,                                 &
                  fastjx_numwl=fastjx_numwl,                                   &
                  fastjx_mode=fastjx_mode,                                     &
                  fastjx_prescutoff=fastjx_prescutoff,                         &
                  i_solcylc_type=i_ukca_solcyc,                                &
                  solcylc_start_year=i_ukca_solcyc_start_year,                 &
!                 l_environ_jo2=l_environ_jo2,                                 &
!                 l_environ_jo2b=l_environ_jo2b,                               &
                  l_strat_chem=l_ukca_strat_chem,                              &
                  ip_aerosol_param_moist=ip_aerosol_param_moist,               &
                  ip_accum_sulphate=ip_accum_sulphate,                         &
                  ip_aitken_sulphate=ip_aitken_sulphate,                       &
                  l_cloud_pc2=(i_cld_vn == i_cld_pc2),                         &
                  l_3d_cca=l_3d_cca,                                           &
                  l_enable_diag_um=.TRUE.,                                     &
                  error_message=ukca_errmsg, error_routine=ukca_errproc)

IF (errcode > 0) THEN
  cmessage = ukca_errmsg
  CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
END IF

! Set up the UKCA configuration based on configuration variables from the UM
! and the working control variables set above
!!!! Also, pass details of the photolysis scheme to allow UKCA to manage
!!!! environmental driver requirements on behalf of the separate UKCA
!!!! Photolysis code for now. Note that this UKCA functionality is intended to
!!!! be removed when possible.

CALL photol_get_config(i_photol_scheme=i_photol_scheme)

CALL ukca_setup(errcode,                                                       &
  ! Context information
  row_length=row_length,                                                       &
  rows=rows,                                                                   &
  model_levels=model_levels,                                                   &
  bl_levels=bl_levels,                                                         &
  nlev_ent_tr_mix=nlev_ent_tr_mix,                                             &
  ntype=ntype,                                                                 &
  npft=npft,                                                                   &
  i_brd_leaf=brd_leaf,                                                         &
  i_brd_leaf_dec=brd_leaf_dec,                                                 &
  i_brd_leaf_eg_trop=brd_leaf_eg_trop,                                         &
  i_brd_leaf_eg_temp=brd_leaf_eg_temp,                                         &
  i_ndl_leaf=ndl_leaf,                                                         &
  i_ndl_leaf_dec=ndl_leaf_dec,                                                 &
  i_ndl_leaf_eg=ndl_leaf_eg,                                                   &
  i_c3_grass=c3_grass,                                                         &
  i_c3_crop=c3_crop,                                                           &
  i_c3_pasture=c3_pasture,                                                     &
  i_c4_grass=c4_grass,                                                         &
  i_c4_crop=c4_crop,                                                           &
  i_c4_pasture=c4_pasture,                                                     &
  i_shrub=shrub,                                                               &
  i_shrub_dec=shrub_dec,                                                       &
  i_shrub_eg=shrub_eg,                                                         &
  i_urban=urban,                                                               &
  i_lake=lake,                                                                 &
  i_soil=soil,                                                                 &
  i_ice=ice,                                                                   &
  i_elev_ice=i_elev_ice,                                                       &
  dzsoil_layer1=dzsoil_io(1),                                                  &
  l_cal360=lcal360,                                                            &
  timestep=timestep,                                                           &
  ! General UKCA configuration options
  i_ukca_chem=i_ukca_chem,                                                     &
  l_ukca_chem_aero=l_ukca_chem_aero,                                           &
  l_ukca_mode=l_ukca_mode,                                                     &
  l_fix_tropopause_level=l_fix_tropopause_level,                               &
  l_ukca_ageair=l_ukca_ageair,                                                 &
  l_ukca_emissions_off=l_ukca_emissions_off,                                   &
  i_ageair_reset_method=i_ageair_reset_method,                                 &
  max_ageair_reset_level=max_ageair_reset_level,                               &
  max_ageair_reset_height=max_ageair_reset_height,                             &
  l_enable_diag_um=l_enable_diag_um,                                           &
  l_ukca_persist_off = l_ukca_persist_off,                                     &
  l_timer=ltimer,                                                              &
  ! Chemistry configuration options
  i_ukca_chem_version = i_ukca_chem_version,                                   &
  nrsteps = nrsteps,                                                           &
  chem_timestep = chem_timestep,                                               &
  dts0 = dts0,                                                                 &
  l_ukca_asad_columns = l_ukca_asad_columns,                                   &
  l_ukca_asad_full = l_ukca_asad_full,                                         &
  l_ukca_debug_asad = l_ukca_debug_asad,                                       &
  l_ukca_intdd = l_ukca_intdd,                                                 &
  l_ukca_ddepo3_ocean = l_ukca_ddepo3_ocean,                                   &
  l_ukca_ddep_lev1 = l_ukca_ddep_lev1,                                         &
  l_ukca_dry_dep_so2wet = l_ukca_dry_dep_so2wet,                               &
  l_deposition_jules = l_deposition_jules,                                     &
  nit = nit,                                                                   &
  l_ukca_quasinewton = l_ukca_quasinewton,                                     &
  i_ukca_quasinewton_start = i_ukca_quasinewton_start,                         &
  i_ukca_quasinewton_end = i_ukca_quasinewton_end,                             &
  ukca_chem_seg_size = ukca_chem_seg_size,                                     &
  max_z_for_offline_chem = max_z_for_offline_chem,                             &
  i_ukca_topboundary = i_ukca_topboundary,                                     &
  l_ukca_ro2_ntp = l_ukca_ro2_ntp,                                             &
  l_ukca_ro2_perm = l_ukca_ro2_perm,                                           &
  ! Chemistry - Heterogeneous chemistry
  l_ukca_het_psc = l_ukca_het_psc,                                             &
  i_ukca_hetconfig = i_ukca_hetconfig,                                         &
  l_ukca_limit_nat = l_ukca_limit_nat,                                         &
  l_ukca_sa_clim = l_ukca_sa_clim,                                             &
  l_ukca_trophet = l_ukca_trophet,                                             &
  l_ukca_classic_hetchem = l_ukca_classic_hetchem,                             &
  l_use_photolysis = (i_ukca_photol /= photol_off ),                           &
  ! UKCA emissions configuration options
  l_ukca_ibvoc = l_ukca_ibvoc,                                                 &
  l_ukca_inferno = l_ukca_inferno,                                             &
  l_ukca_inferno_ch4 = l_ukca_inferno_ch4,                                     &
  i_inferno_emi = i_inferno_emi,                                               &
  l_ukca_so2ems_expvolc = l_ukca_so2ems_expvolc,                               &
  l_ukca_so2ems_plumeria = l_ukca_so2ems_plumeria,                             &
  l_ukca_qch4inter = l_ukca_qch4inter,                                         &
  l_ukca_emsdrvn_ch4 = l_ukca_emsdrvn_ch4,                                     &
  mode_parfrac = mode_parfrac,                                                 &
  l_ukca_enable_seadms_ems = .NOT. l_sulpc_dms,                                &
  i_ukca_dms_flux = i_ukca_dms_flux,                                           &
  l_ukca_scale_seadms_ems = l_ukca_scale_seadms_ems,                           &
  seadms_ems_scaling = seadms_ems_scaling,                                     &
  l_ukca_linox_scaling = l_ukca_linox_scaling,                                 &
  lightnox_scale_fac = lightnox_scale_fac,                                     &
  i_ukca_light_param = i_ukca_light_param,                                     &
  l_ukca_scale_soa_yield_mt = l_ukca_scale_soa_yield_mt,                       &
  soa_yield_scaling_mt = soa_yield_scaling_mt,                                 &
  l_ukca_scale_soa_yield_isop = l_ukca_scale_soa_yield_isop,                   &
  soa_yield_scaling_isop = soa_yield_scaling_isop,                             &
  l_support_ems_vertprof = l_support_ems_vertprof,                             &
  l_support_ems_gridbox_units = l_support_ems_gridbox_units,                   &
  ! UKCA feedback configuration options
  l_ukca_h2o_feedback = l_ukca_h2o_feedback,                                   &
  l_ukca_conserve_h = l_ukca_conserve_h,                                       &
  ! UKCA environmental driver configuration options
  l_param_conv = l_param_conv,                                                 &
  l_ctile = l_ctile,                                                           &
  l_zon_av_ozone = zon_av_ozone,                                               &
  i_strat_lbc_source = i_strat_lbc_source,                                     &
  l_chem_environ_gas_scalars = l_ukca_set_trace_gases,                         &
  l_chem_environ_co2_fld = l_co2_interactive,                                  &
  l_ukca_prescribech4 = l_ukca_prescribech4,                                   &
  l_use_classic_so4 = l_use_classic_so4,                                       &
  l_use_classic_soot = l_use_classic_soot,                                     &
  l_use_classic_ocff = l_use_classic_ocff,                                     &
  l_use_classic_biogenic = l_use_classic_biogenic,                             &
  l_use_classic_seasalt = l_use_classic_seasalt,                               &
  l_use_gridbox_volume = .FALSE.,                                              &
  l_use_gridbox_mass = .FALSE.,                                                &
!s  l_environ_z_top = .TRUE.,                                                    &
  ! z top false in box model
  l_environ_z_top = .FALSE.,                                                   &
  ! UKCA temporary logicals
  l_fix_improve_drydep = l_fix_improve_drydep,                                 &
  l_fix_ukca_h2dd_x = l_fix_ukca_h2dd_x,                                       &
  l_fix_drydep_so2_water = l_fix_drydep_so2_water,                             &
  l_fix_ukca_offox_h2o_fac = l_fix_ukca_offox_h2o_fac,                         &
  l_fix_ukca_h2so4_ystore = l_fix_ukca_h2so4_ystore,                           &
  l_improve_aero_drydep = l_improve_aero_drydep,                               &
  !!!! Temporary settings for managing photolysis environmental driver
  !!!! requirements while these are determined by UKCA
  i_photol_scheme = i_photol_scheme,                                           &
  i_photol_scheme_off = photol_off,                                            &
  i_photol_scheme_strat_only = photol_strat_only,                              &
  i_photol_scheme_2d = photol_2d,                                              &
  i_photol_scheme_fastjx = photol_fastjx,                                      &
  ! General GLOMAP configuration options
  i_mode_nzts = i_mode_nzts,                                                   &
  ukca_mode_seg_size = ukca_mode_seg_size,                                     &
  i_mode_setup = i_mode_setup,                                                 &
  l_mode_bhn_on = l_mode_bhn_on,                                               &
  l_mode_bln_on = l_mode_bln_on,                                               &
  i_mode_bln_param_method = i_mode_bln_param_method,                           &
  i_mode_nucscav = i_mode_nucscav,                                             &
  mode_activation_dryr = mode_activation_dryr,                                 &
  mode_incld_so2_rfrac = mode_incld_so2_rfrac,                                 &
  l_cv_rainout = .NOT. l_ukca_plume_scav,                                      &
  l_ddepaer = .FALSE.,                                                         &
  ! GLOMAP emissions configuration options
  l_ukca_primss = l_ukca_primss,                                               &
  l_ukca_primsu = l_ukca_primsu,                                               &
  l_ukca_primdu = l_ukca_primdu,                                               &
  l_ukca_primbcoc = l_ukca_primbcoc,                                           &
  l_ukca_prim_moc = l_ukca_prim_moc,                                           &
  l_bcoc_bf = l_bcoc_bf,                                                       &
  l_bcoc_bm = l_bcoc_bm,                                                       &
  l_bcoc_ff = l_bcoc_ff,                                                       &
  l_ukca_scale_biom_aer_ems = l_ukca_scale_biom_aer_ems,                       &
  biom_aer_ems_scaling = biom_aer_ems_scaling,                                 &
  l_ukca_fine_no3_prod = l_ukca_fine_no3_prod,                                 &
  l_ukca_coarse_no3_prod = l_ukca_coarse_no3_prod,                             &
  l_no3_prod_in_aero_step = l_no3_prod_in_aero_step,                           &
  l_ukca_scale_sea_salt_ems = l_ukca_scale_sea_salt_ems,                       &
  sea_salt_ems_scaling = sea_salt_ems_scaling,                                 &
  l_ukca_scale_marine_pom_ems = l_ukca_scale_marine_pom_ems,                   &
  marine_pom_ems_scaling = marine_pom_ems_scaling,                             &
  ! GLOMAP feedback configuration options
  l_ukca_radaer = l_ukca_radaer,                                               &
  i_ukca_tune_bc = i_ukca_tune_bc,                                             &
  i_ukca_activation_scheme = i_ukca_activation_scheme,                         &
  i_ukca_nwbins = i_ukca_nwbins,                                               &
  sigwmin = sigwmin,                                                           &
  l_ntpreq_n_activ_sum = l_ntpreq_n_activ_sum,                                 &
  l_ntpreq_dryd_nuc_sol = l_ntpreq_dryd_nuc_sol,                               &
  l_ukca_sfix = l_ukca_sfix,                                                   &
  ! GLOMAP temporary logicals
  l_fix_neg_pvol_wat = l_fix_neg_pvol_wat,                                     &
  l_fix_ukca_impscav = l_fix_ukca_impscav,                                     &
  l_fix_nacl_density = l_fix_nacl_density,                                     &
  l_fix_ukca_activate_pdf = l_fix_ukca_activate_pdf,                           &
  l_fix_ukca_activate_vert_rep = l_fix_ukca_activate_vert_rep,                 &
  l_bug_repro_tke_index = l_bug_repro_tke_index,                               &
  l_fix_ukca_hygroscopicities = l_fix_ukca_hygroscopicities,                   &
  l_fix_ukca_water_content = l_fix_ukca_water_content,                         &
  l_ukca_drydep_off = l_ukca_drydep_off,                                       &
  l_ukca_wetdep_off = l_ukca_wetdep_off,                                       &
  l_tracer_lumping = l_tracer_lumping,                                         &
  ! Error reporting variables
  error_message=ukca_errmsg, error_routine=ukca_errproc)

CALL print_config_settings(                                                    &
  ! Context information
  row_length=row_length,                                                       &
  rows=rows,                                                                   &
  model_levels=model_levels,                                                   &
  bl_levels=bl_levels,                                                         &
  nlev_ent_tr_mix=nlev_ent_tr_mix,                                             &
  ntype=ntype,                                                                 &
  npft=npft,                                                                   &
  i_brd_leaf=brd_leaf,                                                         &
  i_brd_leaf_dec=brd_leaf_dec,                                                 &
  i_brd_leaf_eg_trop=brd_leaf_eg_trop,                                         &
  i_brd_leaf_eg_temp=brd_leaf_eg_temp,                                         &
  i_ndl_leaf=ndl_leaf,                                                         &
  i_ndl_leaf_dec=ndl_leaf_dec,                                                 &
  i_ndl_leaf_eg=ndl_leaf_eg,                                                   &
  i_c3_grass=c3_grass,                                                         &
  i_c3_crop=c3_crop,                                                           &
  i_c3_pasture=c3_pasture,                                                     &
  i_c4_grass=c4_grass,                                                         &
  i_c4_crop=c4_crop,                                                           &
  i_c4_pasture=c4_pasture,                                                     &
  i_shrub=shrub,                                                               &
  i_shrub_dec=shrub_dec,                                                       &
  i_shrub_eg=shrub_eg,                                                         &
  i_urban=urban,                                                               &
  i_lake=lake,                                                                 &
  i_soil=soil,                                                                 &
  i_ice=ice,                                                                   &
  i_elev_ice=i_elev_ice,                                                       &
  dzsoil_layer1=dzsoil_io(1),                                                  &
  l_cal360=lcal360,                                                            &
  timestep=timestep,                                                           &
  ! General UKCA configuration options
  i_ukca_chem=i_ukca_chem,                                                     &
  l_ukca_chem_aero=l_ukca_chem_aero,                                           &
  l_ukca_mode=l_ukca_mode,                                                     &
  l_ukca_ageair=l_ukca_ageair,                                                 &
  l_ukca_emissions_off=l_ukca_emissions_off,                                   &
  i_ageair_reset_method=i_ageair_reset_method,                                 &
  max_ageair_reset_level=max_ageair_reset_level,                               &
  max_ageair_reset_height=max_ageair_reset_height,                             &
  l_enable_diag_um=l_enable_diag_um,                                           &
  l_ukca_persist_off = l_ukca_persist_off,                                     &
  l_timer=ltimer,                                                              &
  ! Chemistry configuration options
  i_ukca_chem_version = i_ukca_chem_version,                                   &
  nrsteps = nrsteps,                                                           &
  chem_timestep = chem_timestep,                                               &
  dts0 = dts0,                                                                 &
  l_ukca_asad_columns = l_ukca_asad_columns,                                   &
  l_ukca_asad_full = l_ukca_asad_full,                                         &
  l_ukca_debug_asad = l_ukca_debug_asad,                                       &
  l_ukca_intdd = l_ukca_intdd,                                                 &
  l_ukca_ddepo3_ocean = l_ukca_ddepo3_ocean,                                   &
  l_ukca_ddep_lev1 = l_ukca_ddep_lev1,                                         &
  l_ukca_dry_dep_so2wet = l_ukca_dry_dep_so2wet,                               &
  l_deposition_jules = l_deposition_jules,                                     &
  nit = nit,                                                                   &
  l_ukca_quasinewton = l_ukca_quasinewton,                                     &
  i_ukca_quasinewton_start = i_ukca_quasinewton_start,                         &
  i_ukca_quasinewton_end = i_ukca_quasinewton_end,                             &
  ukca_chem_seg_size = ukca_chem_seg_size,                                     &
  max_z_for_offline_chem = max_z_for_offline_chem,                             &
  i_ukca_topboundary = i_ukca_topboundary,                                     &
  l_ukca_ro2_ntp = l_ukca_ro2_ntp,                                             &
  l_ukca_ro2_perm = l_ukca_ro2_perm,                                           &
  ! Chemistry - Heterogeneous chemistry
  l_ukca_het_psc = l_ukca_het_psc,                                             &
  i_ukca_hetconfig = i_ukca_hetconfig,                                         &
  l_ukca_limit_nat = l_ukca_limit_nat,                                         &
  l_ukca_sa_clim = l_ukca_sa_clim,                                             &
  l_ukca_trophet = l_ukca_trophet,                                             &
  l_ukca_classic_hetchem = l_ukca_classic_hetchem,                             &
  ! Chemistry - Photolysis
  i_ukca_photol = i_ukca_photol,                                               &
  fastjx_mode = fastjx_mode,                                                   &
  fastjx_prescutoff = fastjx_prescutoff,                                       &
  i_ukca_solcyc = i_ukca_solcyc,                                               &
  i_ukca_solcyc_start_year = i_ukca_solcyc_start_year,                         &
  ! UKCA emissions configuration options
  l_ukca_ibvoc = l_ukca_ibvoc,                                                 &
  l_ukca_inferno = l_ukca_inferno,                                             &
  l_ukca_inferno_ch4 = l_ukca_inferno_ch4,                                     &
  i_inferno_emi = i_inferno_emi,                                               &
  l_ukca_so2ems_expvolc = l_ukca_so2ems_expvolc,                               &
  l_ukca_so2ems_plumeria = l_ukca_so2ems_plumeria,                             &
  l_ukca_qch4inter = l_ukca_qch4inter,                                         &
  l_ukca_emsdrvn_ch4 = l_ukca_emsdrvn_ch4,                                     &
  mode_parfrac = mode_parfrac,                                                 &
  l_ukca_enable_seadms_ems = .NOT. l_sulpc_dms,                                &
  i_ukca_dms_flux = i_ukca_dms_flux,                                           &
  l_ukca_scale_seadms_ems = l_ukca_scale_seadms_ems,                           &
  seadms_ems_scaling = seadms_ems_scaling,                                     &
  l_ukca_linox_scaling = l_ukca_linox_scaling,                                 &
  lightnox_scale_fac = lightnox_scale_fac,                                     &
  i_ukca_light_param = i_ukca_light_param,                                     &
  l_ukca_scale_soa_yield_mt = l_ukca_scale_soa_yield_mt,                       &
  soa_yield_scaling_mt = soa_yield_scaling_mt,                                 &
  l_ukca_scale_soa_yield_isop = l_ukca_scale_soa_yield_isop,                   &
  soa_yield_scaling_isop = soa_yield_scaling_isop,                             &
  l_support_ems_vertprof = l_support_ems_vertprof,                             &
  l_support_ems_gridbox_units = l_support_ems_gridbox_units,                   &
  ! UKCA feedback configuration options
  l_ukca_h2o_feedback = l_ukca_h2o_feedback,                                   &
  l_ukca_conserve_h = l_ukca_conserve_h,                                       &
  ! UKCA environmental driver configuration options
  l_param_conv = l_param_conv,                                                 &
  l_ctile = l_ctile,                                                           &
  l_zon_av_ozone = zon_av_ozone,                                               &
  i_strat_lbc_source = i_strat_lbc_source,                                     &
  l_chem_environ_gas_scalars = l_ukca_set_trace_gases,                         &
  l_chem_environ_co2_fld = l_co2_interactive,                                  &
  l_ukca_prescribech4 = l_ukca_prescribech4,                                   &
  l_use_classic_so4 = l_use_classic_so4,                                       &
  l_use_classic_soot = l_use_classic_soot,                                     &
  l_use_classic_ocff = l_use_classic_ocff,                                     &
  l_use_classic_biogenic = l_use_classic_biogenic,                             &
  l_use_classic_seasalt = l_use_classic_seasalt,                               &
  l_use_gridbox_volume = .FALSE.,                                              &
  l_use_gridbox_mass = .FALSE.,                                                &
  l_environ_z_top = .FALSE.,                                                   &
  ! UKCA temporary logicals
  l_fix_improve_drydep = l_fix_improve_drydep,                                 &
  l_fix_ukca_h2dd_x = l_fix_ukca_h2dd_x,                                       &
  l_fix_drydep_so2_water = l_fix_drydep_so2_water,                             &
  l_fix_ukca_offox_h2o_fac = l_fix_ukca_offox_h2o_fac,                         &
  l_fix_ukca_h2so4_ystore = l_fix_ukca_h2so4_ystore,                           &
  l_improve_aero_drydep = l_improve_aero_drydep,                               &
  ! General GLOMAP configuration options
  i_mode_nzts = i_mode_nzts,                                                   &
  ukca_mode_seg_size = ukca_mode_seg_size,                                     &
  i_mode_setup = i_mode_setup,                                                 &
  l_mode_bhn_on = l_mode_bhn_on,                                               &
  l_mode_bln_on = l_mode_bln_on,                                               &
  i_mode_bln_param_method = i_mode_bln_param_method,                           &
  i_mode_nucscav = i_mode_nucscav,                                             &
  mode_activation_dryr = mode_activation_dryr,                                 &
  mode_incld_so2_rfrac = mode_incld_so2_rfrac,                                 &
  l_cv_rainout = .TRUE.,                                                       &
  l_ddepaer = .TRUE.,                                                          &
  ! GLOMAP emissions configuration options
  l_ukca_primss = l_ukca_primss,                                               &
  l_ukca_primsu = l_ukca_primsu,                                               &
  l_ukca_primdu = l_ukca_primdu,                                               &
  l_ukca_primbcoc = l_ukca_primbcoc,                                           &
  l_ukca_prim_moc = l_ukca_prim_moc,                                           &
  l_bcoc_bf = l_bcoc_bf,                                                       &
  l_bcoc_bm = l_bcoc_bm,                                                       &
  l_bcoc_ff = l_bcoc_ff,                                                       &
  l_ukca_scale_biom_aer_ems = l_ukca_scale_biom_aer_ems,                       &
  biom_aer_ems_scaling = biom_aer_ems_scaling,                                 &
  l_ukca_fine_no3_prod = l_ukca_fine_no3_prod,                                 &
  l_ukca_coarse_no3_prod = l_ukca_coarse_no3_prod,                             &
  l_no3_prod_in_aero_step = l_no3_prod_in_aero_step,                           &
  l_ukca_scale_sea_salt_ems = l_ukca_scale_sea_salt_ems,                       &
  sea_salt_ems_scaling = sea_salt_ems_scaling,                                 &
  l_ukca_scale_marine_pom_ems = l_ukca_scale_marine_pom_ems,                   &
  marine_pom_ems_scaling = marine_pom_ems_scaling,                             &
  ! GLOMAP feedback configuration options
  l_ukca_radaer = l_ukca_radaer,                                               &
  i_ukca_tune_bc = i_ukca_tune_bc,                                             &
  i_ukca_activation_scheme = i_ukca_activation_scheme,                         &
  i_ukca_nwbins = i_ukca_nwbins,                                               &
  sigwmin = sigwmin,                                                           &
  l_ntpreq_n_activ_sum = l_ntpreq_n_activ_sum,                                 &
  l_ntpreq_dryd_nuc_sol = l_ntpreq_dryd_nuc_sol,                               &
  l_ukca_sfix = l_ukca_sfix,                                                   &
  ! GLOMAP temporary logicals
  l_fix_neg_pvol_wat = l_fix_neg_pvol_wat,                                     &
  l_fix_ukca_impscav = .TRUE.,                                                 &
  l_fix_nacl_density = l_fix_nacl_density,                                     &
  l_fix_ukca_activate_pdf = l_fix_ukca_activate_pdf,                           &
  l_fix_ukca_activate_vert_rep = l_fix_ukca_activate_vert_rep,                 &
  l_bug_repro_tke_index = l_bug_repro_tke_index,                               &
  l_fix_ukca_hygroscopicities = l_fix_ukca_hygroscopicities,                   &
  l_fix_ukca_water_content = l_fix_ukca_water_content,                         &
  ! New settings in Ticket #25
  l_ukca_drydep_off = l_ukca_drydep_off,                                       &
  l_ukca_wetdep_off = l_ukca_wetdep_off,                                       &
  nlev_above_trop_o3_env = 2,                                                  &
  nlev_ch4_stratloss = 1,                                                      &
  l_tracer_lumping = l_tracer_lumping,                                         &
  env_log_step = 12,                                                           &
  l_aero_rainout = .TRUE.,                                                          &
  l_impc_scav = .TRUE.)

!++SAN
CALL umPrint('SAN: checking GLOMAP mode values',src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: ukca_mode_seg_size = ', ukca_mode_seg_size
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: i_mode_setup = ', i_mode_setup
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: l_mode_bhn_on = ', l_mode_bhn_on
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: i_mode_nucscav = ', i_mode_nucscav
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: mode_activation_dryr = ', mode_activation_dryr
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,*) 'SAN BOX_UKCA_SETUP: biom_aer_ems_scaling = ', biom_aer_ems_scaling
CALL umPrint(umMessage,src=RoutineName)



IF (errcode > 0) THEN
  cmessage = ukca_errmsg
  CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
END IF

! Check that the UKCA configuration created is consistent with the UM inputs.
CALL check_ukca_configuration()

! Deallocate working array used in API calls
IF (ALLOCATED(i_elev_ice)) DEALLOCATE(i_elev_ice)

! Retrieve internal UKCA configuration variables that are required by the UM
CALL ukca_get_config(l_ukca_chem=l_ukca_chem,                                  &
                     l_ukca_trop=l_ukca_trop,                                  &
                     l_ukca_raq=l_ukca_raq,                                    &
                     l_ukca_raqaero=l_ukca_raqaero,                            &
                     l_ukca_offline_be=l_ukca_offline_be,                      &
                     l_ukca_tropisop=l_ukca_tropisop,                          &
                     l_ukca_strattrop=l_ukca_strattrop,                        &
                     l_ukca_strat=l_ukca_strat,                                &
                     l_ukca_offline=l_ukca_offline,                            &
                     l_ukca_stratcfc=l_ukca_stratcfc,                          &
                     ukca_int_method=ukca_int_method)

! Set up the scavenging coefficients for plume scavenging
IF (l_ukca_mode) CALL ukca_mode_scavcoeff()

! Set item numbers for CDNC etc
IF (i_ukca_activation_scheme == ukca_activation_jones .OR.                     &
    i_ukca_activation_scheme == ukca_activation_arg) THEN
  CALL ukca_cdnc_init()
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE box_ukca_setup

! ----------------------------------------------------------------------
SUBROUTINE check_ukca_configuration()
! ----------------------------------------------------------------------
! Description:
!
! Check that the UKCA configuration is consistent with the UM inputs.
! For each active configuration variable, there should be a direct match
! with the UM equivalent.
! ----------------------------------------------------------------------

USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Local variables

! Values of configuration options set in UKCA for checking against inputs

! -- Context information ---
INTEGER :: row_length_from_ukca
INTEGER :: rows_from_ukca
INTEGER :: model_levels_from_ukca
INTEGER :: ntype_from_ukca
INTEGER :: npft_from_ukca
INTEGER :: i_brd_leaf_from_ukca
INTEGER :: i_brd_leaf_dec_from_ukca
INTEGER :: i_brd_leaf_eg_trop_from_ukca
INTEGER :: i_brd_leaf_eg_temp_from_ukca
INTEGER :: i_ndl_leaf_from_ukca
INTEGER :: i_ndl_leaf_dec_from_ukca
INTEGER :: i_ndl_leaf_eg_from_ukca
INTEGER :: i_c3_grass_from_ukca
INTEGER :: i_c3_crop_from_ukca
INTEGER :: i_c3_pasture_from_ukca
INTEGER :: i_c4_grass_from_ukca
INTEGER :: i_c4_crop_from_ukca
INTEGER :: i_c4_pasture_from_ukca
INTEGER :: i_shrub_from_ukca
INTEGER :: i_shrub_dec_from_ukca
INTEGER :: i_shrub_eg_from_ukca
INTEGER :: i_urban_from_ukca
INTEGER :: i_lake_from_ukca
INTEGER :: i_soil_from_ukca
INTEGER :: i_ice_from_ukca
INTEGER, ALLOCATABLE :: i_elev_ice_from_ukca(:)
REAL :: dzsoil_layer1_from_ukca
REAL :: timestep_from_ukca

! -- General UKCA configuration options --
INTEGER :: i_ukca_chem_from_ukca
INTEGER :: i_ageair_reset_method_from_ukca
INTEGER :: max_ageair_reset_level_from_ukca
REAL :: max_ageair_reset_height_from_ukca
LOGICAL :: l_ukca_chem_aero_from_ukca
LOGICAL :: l_ukca_mode_from_ukca
LOGICAL :: l_ukca_ageair_from_ukca
LOGICAL :: l_ukca_emissions_off_from_ukca
LOGICAL :: l_fix_tropopause_level_from_ukca
LOGICAL :: l_ukca_drydep_off_from_ukca
LOGICAL :: l_ukca_wetdep_off_from_ukca
LOGICAL :: l_tracer_lumping_from_ukca
LOGICAL :: l_enable_diag_um_from_ukca
LOGICAL :: l_ukca_persist_off_from_ukca
LOGICAL :: l_timer_from_ukca

! -- Chemistry configuration options --
INTEGER :: i_ukca_chem_version_from_ukca
INTEGER :: nrsteps_from_ukca
INTEGER :: chem_timestep_from_ukca
INTEGER :: dts0_from_ukca
INTEGER :: nit_from_ukca
INTEGER :: i_ukca_quasinewton_start_from_ukca
INTEGER :: i_ukca_quasinewton_end_from_ukca
INTEGER :: ukca_chem_seg_size_from_ukca
INTEGER :: i_ukca_topboundary_from_ukca
REAL :: max_z_for_offline_chem_from_ukca
LOGICAL :: l_ukca_asad_columns_from_ukca
LOGICAL :: l_ukca_asad_full_from_ukca
LOGICAL :: l_ukca_debug_asad_from_ukca
LOGICAL :: l_ukca_ddepo3_ocean_from_ukca
LOGICAL :: l_ukca_ddep_lev1_from_ukca
LOGICAL :: l_ukca_dry_dep_so2wet_from_ukca
LOGICAL :: l_deposition_jules_from_ukca
LOGICAL :: l_ukca_quasinewton_from_ukca
LOGICAL :: l_ukca_ro2_ntp_from_ukca
LOGICAL :: l_ukca_ro2_perm_from_ukca

! -- Chemistry - Heterogeneous chemistry --
INTEGER :: i_ukca_hetconfig_from_ukca
LOGICAL :: l_ukca_het_psc_from_ukca
LOGICAL :: l_ukca_limit_nat_from_ukca
LOGICAL :: l_ukca_sa_clim_from_ukca
LOGICAL :: l_ukca_trophet_from_ukca
LOGICAL :: l_ukca_classic_hetchem_from_ukca

! -- Chemistry - Photolysis --
INTEGER :: i_scheme_from_photol
INTEGER :: fastjx_mode_from_photol
INTEGER :: i_solcyc_from_photol
INTEGER :: i_solcyc_start_year_from_photol
REAL :: fastjx_prescutoff_from_photol
!LOGICAL :: l_environ_jo2_from_ukca
!LOGICAL :: l_environ_jo2b_from_ukca

! -- UKCA feedback configuration options --
LOGICAL :: l_ukca_h2o_feedback_from_ukca
LOGICAL :: l_ukca_conserve_h_from_ukca

! -- UKCA environmental driver configuration options --
LOGICAL :: l_param_conv_from_ukca
LOGICAL :: l_ctile_from_ukca
LOGICAL :: l_zon_av_ozone_from_ukca
LOGICAL :: l_chem_environ_gas_scalars_from_ukca
LOGICAL :: l_chem_environ_co2_fld_from_ukca
LOGICAL :: l_use_classic_so4_from_ukca
LOGICAL :: l_use_classic_soot_from_ukca
LOGICAL :: l_use_classic_ocff_from_ukca
LOGICAL :: l_use_classic_biogenic_from_ukca
LOGICAL :: l_use_classic_seasalt_from_ukca
LOGICAL :: l_use_gridbox_volume_from_ukca
LOGICAL :: l_use_gridbox_mass_from_ukca
LOGICAL :: l_environ_z_top_from_ukca

! -- UKCA temporary logicals --
LOGICAL :: l_fix_drydep_so2_water_from_ukca
LOGICAL :: l_fix_ukca_offox_h2o_fac_from_ukca
LOGICAL :: l_fix_ukca_h2so4_ystore_from_ukca
LOGICAL :: l_improve_aero_drydep_from_ukca

! -- General GLOMAP configuration options --
INTEGER :: i_mode_nzts_from_ukca
INTEGER :: ukca_mode_seg_size_from_ukca
INTEGER :: i_mode_setup_from_ukca
INTEGER :: i_mode_bln_param_method_from_ukca
REAL :: mode_activation_dryr_from_ukca
LOGICAL :: l_mode_bhn_on_from_ukca
LOGICAL :: l_mode_bln_on_from_ukca
LOGICAL :: l_ddepaer_from_ukca

! -- GLOMAP feedback configuration options --
INTEGER :: i_ukca_activation_scheme_from_ukca
INTEGER :: i_ukca_nwbins_from_ukca
REAL :: sigwmin_from_ukca
LOGICAL :: l_ukca_radaer_from_ukca
INTEGER :: i_ukca_tune_bc_from_ukca
LOGICAL :: l_ntpreq_n_activ_sum_from_ukca
LOGICAL :: l_ntpreq_dryd_nuc_sol_from_ukca
LOGICAL :: l_ukca_sfix_from_ukca

! -- GLOMAP temporary logicals --
LOGICAL :: l_fix_neg_pvol_wat_from_ukca
LOGICAL :: l_fix_ukca_impscav_from_ukca
LOGICAL :: l_fix_nacl_density_from_ukca
LOGICAL :: l_fix_ukca_activate_pdf_from_ukca
LOGICAL :: l_fix_ukca_activate_vert_rep_from_ukca
LOGICAL :: l_bug_repro_tke_index_from_ukca
LOGICAL :: l_fix_ukca_hygroscopicities_from_ukca
LOGICAL :: l_fix_ukca_water_content_from_ukca

! Dr Hook
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_UKCA_CONFIGURATION'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! -- Check context information required for all configurations --

CALL ukca_get_config(                                                          &
  row_length=row_length_from_ukca,                                             &
  rows=rows_from_ukca,                                                         &
  model_levels=model_levels_from_ukca,                                         &
  timestep=timestep_from_ukca)

CALL check_config_value(row_length_from_ukca, row_length, 'row_length')
CALL check_config_value(rows_from_ukca, rows, 'rows')
CALL check_config_value(model_levels_from_ukca, model_levels, 'model_levels')
CALL check_config_value(timestep_from_ukca, timestep, 'timestep')

! -- Check general UKCA configuration options --

CALL ukca_get_config(                                                          &
  i_ukca_chem=i_ukca_chem_from_ukca,                                           &
  l_ukca_chem_aero=l_ukca_chem_aero_from_ukca,                                 &
  l_ukca_mode=l_ukca_mode_from_ukca,                                           &
  l_ukca_ageair=l_ukca_ageair_from_ukca,                                       &
  l_ukca_emissions_off=l_ukca_emissions_off_from_ukca,                         &
  l_fix_tropopause_level=l_fix_tropopause_level_from_ukca,                     &
  l_ukca_drydep_off=l_ukca_drydep_off_from_ukca,                               &
  l_ukca_wetdep_off=l_ukca_wetdep_off_from_ukca,                               &
  l_tracer_lumping=l_tracer_lumping_from_ukca,                                 &
  i_ageair_reset_method=i_ageair_reset_method_from_ukca,                       &
  max_ageair_reset_level=max_ageair_reset_level_from_ukca,                     &
  max_ageair_reset_height=max_ageair_reset_height_from_ukca,                   &
  l_enable_diag_um=l_enable_diag_um_from_ukca,                                 &
  l_ukca_persist_off=l_ukca_persist_off_from_ukca,                             &
  l_timer=l_timer_from_ukca)

CALL check_config_value(i_ukca_chem_from_ukca,                                 &
                        i_ukca_chem, 'i_ukca_chem')
CALL check_config_value(l_ukca_chem_aero_from_ukca,                            &
                        l_ukca_chem_aero, 'l_ukca_chem_aero')
CALL check_config_value(l_ukca_mode_from_ukca,                                 &
                        l_ukca_mode, 'l_ukca_mode')
CALL check_config_value(l_ukca_ageair_from_ukca,                               &
                        l_ukca_ageair, 'l_ukca_ageair')

CALL check_config_value(l_ukca_emissions_off_from_ukca,l_ukca_emissions_off,   &
                        'l_ukca_emissions_off')
CALL check_config_value(l_fix_tropopause_level_from_ukca,                      &
                        l_fix_tropopause_level, 'l_fix_tropopause_level')
CALL check_config_value(l_ukca_drydep_off_from_ukca,l_ukca_drydep_off,         &
                        'l_ukca_drydep_off')
CALL check_config_value(l_ukca_wetdep_off_from_ukca,l_ukca_wetdep_off,         &
                        'l_ukca_wetdep_off')
CALL check_config_value(l_tracer_lumping_from_ukca,l_tracer_lumping,           &
                        'l_tracer_lumping')

IF (l_ukca_ageair) THEN
  CALL check_config_value(i_ageair_reset_method_from_ukca,                     &
                          i_ageair_reset_method, 'i_ageair_reset_method')
  IF (i_ageair_reset_method == ukca_age_reset_by_level) THEN
    CALL check_config_value(max_ageair_reset_level_from_ukca,                  &
                            max_ageair_reset_level, 'max_ageair_reset_level')
  ELSE IF (i_ageair_reset_method == ukca_age_reset_by_height) THEN
    CALL check_config_value(max_ageair_reset_height_from_ukca,                 &
                            max_ageair_reset_height, 'max_ageair_reset_height')
  END IF
END IF

CALL check_config_value(l_enable_diag_um_from_ukca,                            &
                        l_enable_diag_um, 'l_enable_diag_um')
CALL check_config_value(l_ukca_persist_off_from_ukca, l_ukca_persist_off,      &
                        'l_ukca_persist_off')
CALL check_config_value(l_timer_from_ukca, ltimer, 'l_timer')

! Check chemistry configuration variables only if chemistry is on

IF (i_ukca_chem /= ukca_chem_off) THEN

  ! -- Check chemistry configuration options --

  CALL ukca_get_config(                                                        &
    i_ukca_chem_version=i_ukca_chem_version_from_ukca,                         &
    nrsteps=nrsteps_from_ukca,                                                 &
    chem_timestep=chem_timestep_from_ukca,                                     &
    dts0=dts0_from_ukca,                                                       &
    l_ukca_asad_columns=l_ukca_asad_columns_from_ukca,                         &
    l_ukca_asad_full=l_ukca_asad_full_from_ukca,                               &
    l_ukca_debug_asad=l_ukca_debug_asad_from_ukca,                             &
    nit=nit_from_ukca,                                                         &
    l_ukca_quasinewton=l_ukca_quasinewton_from_ukca,                           &
    i_ukca_quasinewton_start=i_ukca_quasinewton_start_from_ukca,               &
    i_ukca_quasinewton_end=i_ukca_quasinewton_end_from_ukca,                   &
    ukca_chem_seg_size=ukca_chem_seg_size_from_ukca,                           &
    max_z_for_offline_chem=max_z_for_offline_chem_from_ukca,                   &
    i_ukca_topboundary=i_ukca_topboundary_from_ukca,                           &
    l_ukca_ro2_ntp=l_ukca_ro2_ntp_from_ukca,                                   &
    l_ukca_ro2_perm=l_ukca_ro2_perm_from_ukca)

  CALL check_config_value(chem_timestep_from_ukca,                             &
                          chem_timestep, 'chem_timestep')
  CALL check_config_value(l_ukca_asad_columns_from_ukca,                       &
                          l_ukca_asad_columns, 'l_ukca_asad_columns')
  CALL check_config_value(l_ukca_asad_full_from_ukca,                          &
                          l_ukca_asad_full, 'l_ukca_asad_full')

  !!!! CHECK l_ukca_ro2_ntp_from_ukca & l_ukca_ro2_perm_from_ukca HERE
  !!!! CONDITIONAL ON i_ukca_chem VALUE

  ! Check solver configuration

  SELECT CASE (ukca_int_method)

  CASE (ukca_int_method_be_explicit)

    ! Check configuration values specific to B-E schemes

    CALL check_config_value(dts0_from_ukca, dts0, 'dts0')
    CALL check_config_value(nit_from_ukca, nit, 'nit')

  CASE (ukca_int_method_nr)

    ! Check configuration values specific to N-R schemes

    CALL check_config_value(i_ukca_chem_version_from_ukca,                     &
                            i_ukca_chem_version, 'i_ukca_chem_version')
    CALL check_config_value(nrsteps_from_ukca, nrsteps, 'nrsteps')
    CALL check_config_value(l_ukca_debug_asad_from_ukca,                       &
                            l_ukca_debug_asad, 'l_ukca_debug_asad')
    CALL check_config_value(l_ukca_quasinewton_from_ukca,                      &
                            l_ukca_quasinewton, 'l_ukca_quasinewton')

    IF (l_ukca_quasinewton) THEN
      CALL check_config_value(i_ukca_quasinewton_start_from_ukca,              &
                              i_ukca_quasinewton_start,                        &
                              'i_ukca_quasinewton_start')
      CALL check_config_value(i_ukca_quasinewton_end_from_ukca,                &
                              i_ukca_quasinewton_end, 'i_ukca_quasinewton_end')
    END IF

  END SELECT

  ! Check configuration for column-based processing
  IF (l_ukca_asad_columns) THEN
    CALL check_config_value(ukca_chem_seg_size_from_ukca,                      &
                            ukca_chem_seg_size, 'ukca_chem_seg_size')
  END IF

  ! Check other chemistry-scheme specific values

  IF (i_ukca_chem == ukca_chem_offline_be) THEN
    CALL check_config_value(max_z_for_offline_chem_from_ukca,                  &
                            max_z_for_offline_chem, 'max_z_for_offline_chem')
  END IF

  IF (i_ukca_chem == ukca_chem_strat .OR.                                      &
      i_ukca_chem == ukca_chem_strattrop .OR.                                  &
      i_ukca_chem == ukca_chem_cristrat) THEN
    CALL check_config_value(i_ukca_topboundary_from_ukca,                      &
                            i_ukca_topboundary, 'i_ukca_topboundary')
  END IF

  ! -- Chemistry - Heterogeneous chemistry --

  CALL ukca_get_config(                                                        &
    l_ukca_het_psc=l_ukca_het_psc_from_ukca,                                   &
    i_ukca_hetconfig=i_ukca_hetconfig_from_ukca,                               &
    l_ukca_limit_nat=l_ukca_limit_nat_from_ukca,                               &
    l_ukca_sa_clim=l_ukca_sa_clim_from_ukca,                                   &
    l_ukca_trophet=l_ukca_trophet_from_ukca,                                   &
    l_ukca_classic_hetchem=l_ukca_classic_hetchem_from_ukca)

  CALL check_config_value(l_ukca_het_psc_from_ukca,                            &
                          l_ukca_het_psc, 'l_ukca_het_psc')

  IF (l_ukca_het_psc) THEN
    CALL check_config_value(i_ukca_hetconfig_from_ukca,                        &
                            i_ukca_hetconfig, 'i_ukca_hetconfig')
    CALL check_config_value(l_ukca_limit_nat_from_ukca,                        &
                            l_ukca_limit_nat, 'l_ukca_limit_nat')
    CALL check_config_value(l_ukca_sa_clim_from_ukca,                          &
                            l_ukca_sa_clim, 'l_ukca_sa_clim')
  END IF

  CALL check_config_value(l_ukca_trophet_from_ukca,                            &
                          l_ukca_trophet, 'l_ukca_trophet')

  IF (i_ukca_chem == ukca_chem_raq) THEN
    CALL check_config_value(l_ukca_classic_hetchem_from_ukca,                  &
                            l_ukca_classic_hetchem, 'l_ukca_classic_hetchem')
  END IF

  ! -- Chemistry - Photolysis --

  CALL photol_get_config(                                                      &
    i_photol_scheme=i_scheme_from_photol,                                      &
    fastjx_mode=fastjx_mode_from_photol,                                       &
    fastjx_prescutoff=fastjx_prescutoff_from_photol,                           &
    i_solcylc_type=i_solcyc_from_photol,                                       &
    solcylc_start_year=i_solcyc_start_year_from_photol) !,                     &
!   l_environ_jo2=l_environ_jo2_from_ukca,                                     &
!   l_environ_jo2b=l_environ_jo2b_from_ukca)

  CALL check_config_value(i_scheme_from_photol,                                &
                          i_ukca_photol, 'i_ukca_photol')

  IF (i_scheme_from_photol == photol_fastjx) THEN
    CALL check_config_value(fastjx_mode_from_photol, fastjx_mode, 'fastjx_mode')
    CALL check_config_value(fastjx_prescutoff_from_photol,                     &
                            fastjx_prescutoff, 'fastjx_prescutoff')
    CALL check_config_value(i_solcyc_from_photol,                              &
                            i_ukca_solcyc, 'i_ukca_solcyc')
    CALL check_config_value(i_solcyc_start_year_from_photol,                   &
                            i_ukca_solcyc_start_year,                          &
                            'i_ukca_solcyc_start_year')
  END IF
!  CALL check_config_value(l_environ_jo2_from_ukca,                             &
!                          l_environ_jo2, 'l_environ_jo2')
!  CALL check_config_value(l_environ_jo2b_from_ukca,                            &
!                          l_environ_jo2b, 'l_environ_jo2b')

  ! -- Check UKCA feedback configuration options --

  CALL ukca_get_config(l_ukca_h2o_feedback=l_ukca_h2o_feedback_from_ukca,      &
                       l_ukca_conserve_h=l_ukca_conserve_h_from_ukca)

  CALL check_config_value(l_ukca_h2o_feedback_from_ukca,                       &
                          l_ukca_h2o_feedback, 'l_ukca_h2o_feedback')
  CALL check_config_value(l_ukca_conserve_h_from_ukca,                         &
                          l_ukca_conserve_h, 'l_ukca_conserve_h')

  ! -- Check UKCA environmental driver configuration options --

  CALL ukca_get_config(                                                        &
    l_param_conv=l_param_conv_from_ukca,                                       &
    l_ctile=l_ctile_from_ukca,                                                 &
    l_zon_av_ozone=l_zon_av_ozone_from_ukca,                                   &
    l_chem_environ_gas_scalars=l_chem_environ_gas_scalars_from_ukca,           &
    l_chem_environ_co2_fld=l_chem_environ_co2_fld_from_ukca,                   &
    l_use_classic_so4=l_use_classic_so4_from_ukca,                             &
    l_use_classic_soot=l_use_classic_soot_from_ukca,                           &
    l_use_classic_ocff=l_use_classic_ocff_from_ukca,                           &
    l_use_classic_biogenic=l_use_classic_biogenic_from_ukca,                   &
    l_use_classic_seasalt=l_use_classic_seasalt_from_ukca,                     &
    l_use_gridbox_volume=l_use_gridbox_volume_from_ukca,                       &
    l_use_gridbox_mass=l_use_gridbox_mass_from_ukca,                           &
    l_environ_z_top=l_environ_z_top_from_ukca)

  CALL check_config_value(l_param_conv_from_ukca,                              &
                          l_param_conv, 'l_param_conv')
  CALL check_config_value(l_ctile_from_ukca,                                   &
                          l_ctile, 'l_ctile')

  IF (i_ukca_chem == ukca_chem_trop .OR.                                       &
      i_ukca_chem == ukca_chem_raq .OR.                                        &
      i_ukca_chem == ukca_chem_tropisop) THEN
    CALL check_config_value(l_zon_av_ozone_from_ukca,                          &
                            zon_av_ozone, 'l_zon_av_ozone')
  END IF

  CALL check_config_value(l_chem_environ_gas_scalars_from_ukca,                &
                          l_ukca_set_trace_gases,                              &
                          'l_chem_environ_gas_scalars')
  CALL check_config_value(l_chem_environ_co2_fld_from_ukca,                    &
                          l_co2_interactive,                                   &
                          'l_chem_environ_co2_fld')

  IF (l_ukca_classic_hetchem .OR. l_ukca_sa_clim) THEN
    CALL check_config_value(l_use_classic_so4_from_ukca,                       &
                            l_use_classic_so4, 'l_use_classic_so4')
  END IF

  IF (l_ukca_classic_hetchem) THEN
    CALL check_config_value(l_use_classic_soot_from_ukca,                      &
                            l_use_classic_soot, 'l_use_classic_soot')
    CALL check_config_value(l_use_classic_ocff_from_ukca,                      &
                            l_use_classic_ocff, 'l_use_classic_ocff')
    CALL check_config_value(l_use_classic_biogenic_from_ukca,                  &
                            l_use_classic_biogenic, 'l_use_classic_biogenic')
    CALL check_config_value(l_use_classic_seasalt_from_ukca,                   &
                            l_use_classic_seasalt, 'l_use_classic_seasalt')
  END IF

  ! don't use gridbox mass calulcation as do not have solid angle etc. info
  CALL check_config_value(l_use_gridbox_mass_from_ukca,                        &
                          .FALSE., 'l_use_gridbox_mass')
  ! don't use gridbox volume calulcation
  CALL check_config_value(l_use_gridbox_volume_from_ukca,                        &
                          .FALSE., 'l_use_gridbox_volume')
!s  CALL check_config_value(l_environ_z_top_from_ukca, .TRUE., 'l_environ_z_top')
  ! l_environ_z_top is false in box model
  CALL check_config_value(l_environ_z_top_from_ukca, .FALSE., 'l_environ_z_top')

  ! -- Check UKCA temporary logicals --

  CALL ukca_get_config(                                                        &
    l_fix_drydep_so2_water=l_fix_drydep_so2_water_from_ukca,                   &
    l_fix_ukca_offox_h2o_fac=l_fix_ukca_offox_h2o_fac_from_ukca,               &
    l_fix_ukca_h2so4_ystore=l_fix_ukca_h2so4_ystore_from_ukca,                 &
    l_improve_aero_drydep=l_improve_aero_drydep_from_ukca)

  ! Check temporary logicals for specific chemistry schemes

  IF (i_ukca_chem == ukca_chem_offline_be) THEN
    CALL check_config_value(l_fix_ukca_offox_h2o_fac_from_ukca,                &
                            l_fix_ukca_offox_h2o_fac,                          &
                            'l_fix_ukca_offox_h2o_fac')
  END IF

  IF (l_ukca_mode .AND. (ukca_int_method == ukca_int_method_nr)) THEN
    CALL check_config_value(l_fix_ukca_h2so4_ystore_from_ukca,                 &
                            l_fix_ukca_h2so4_ystore, 'l_fix_ukca_h2so4_ystore')
  END IF

END IF ! i_ukca_chem /= ukca_chem_off

! Check GLOMAP-specific configuration variables only if GLOMAP is selected

IF (l_ukca_mode) THEN

  ! -- Check general GLOMAP configuration options --

  CALL ukca_get_config(                                                        &
    i_mode_nzts=i_mode_nzts_from_ukca,                                         &
    ukca_mode_seg_size=ukca_mode_seg_size_from_ukca,                           &
    i_mode_setup=i_mode_setup_from_ukca,                                       &
    l_mode_bhn_on=l_mode_bhn_on_from_ukca,                                     &
    l_mode_bln_on=l_mode_bln_on_from_ukca,                                     &
    i_mode_bln_param_method=i_mode_bln_param_method_from_ukca,                 &
    mode_activation_dryr=mode_activation_dryr_from_ukca,                       &
    l_ddepaer=l_ddepaer_from_ukca)

  CALL check_config_value(i_mode_nzts_from_ukca, i_mode_nzts, 'i_mode_nzts')
  CALL check_config_value(ukca_mode_seg_size_from_ukca,                        &
                          ukca_mode_seg_size, 'ukca_mode_seg_size')
  CALL check_config_value(i_mode_setup_from_ukca, i_mode_setup, 'i_mode_setup')
  CALL check_config_value(l_mode_bhn_on_from_ukca,                             &
                          l_mode_bhn_on, 'l_mode_bhn_on')
  CALL check_config_value(l_mode_bln_on_from_ukca,                             &
                          l_mode_bln_on, 'l_mode_bln_on')

  IF (l_mode_bln_on) THEN
    CALL check_config_value(i_mode_bln_param_method_from_ukca,                 &
                            i_mode_bln_param_method, 'i_mode_bln_param_method')
  END IF

  CALL check_config_value(mode_activation_dryr_from_ukca,                      &
                          mode_activation_dryr, 'mode_activation_dryr')
  ! Make sure aerosol dry deposition is off in box model
  CALL check_config_value(l_ddepaer_from_ukca,                                 &
                          .FALSE., 'l_ddepaer')

  ! -- Check GLOMAP feedback configuration options --

  CALL ukca_get_config(                                                        &
    l_ukca_radaer=l_ukca_radaer_from_ukca,                                     &
    i_ukca_tune_bc=i_ukca_tune_bc_from_ukca,                                   &
    i_ukca_activation_scheme=i_ukca_activation_scheme_from_ukca,               &
    i_ukca_nwbins=i_ukca_nwbins_from_ukca,                                     &
    sigwmin=sigwmin_from_ukca,                                                 &
    l_ntpreq_n_activ_sum=l_ntpreq_n_activ_sum_from_ukca,                       &
    l_ntpreq_dryd_nuc_sol=l_ntpreq_dryd_nuc_sol_from_ukca,                     &
    l_ukca_sfix=l_ukca_sfix_from_ukca)

  CALL check_config_value(l_ukca_radaer_from_ukca,                             &
                          l_ukca_radaer, 'l_ukca_radaer')
  IF (l_ukca_radaer) THEN
    CALL check_config_value(i_ukca_tune_bc_from_ukca, i_ukca_tune_bc,          &
                            'i_ukca_tune_bc')
  END IF
  CALL check_config_value(i_ukca_activation_scheme_from_ukca,                  &
                          i_ukca_activation_scheme, 'i_ukca_activation_scheme')

  ! Check activation-scheme specific values

  SELECT CASE (i_ukca_activation_scheme)

  CASE (ukca_activation_arg)
    CALL check_config_value(i_ukca_nwbins_from_ukca,                           &
                            i_ukca_nwbins,'i_ukca_nwbins')
    CALL check_config_value(sigwmin_from_ukca, sigwmin, 'sigwmin')
    CALL check_config_value(l_ntpreq_n_activ_sum_from_ukca,                    &
                            l_ntpreq_n_activ_sum, 'l_ntpreq_n_activ_sum')
    CALL check_config_value(l_ntpreq_dryd_nuc_sol_from_ukca,                   &
                            l_ntpreq_dryd_nuc_sol, 'l_ntpreq_dryd_nuc_sol')
    CALL check_config_value(l_ukca_sfix_from_ukca, l_ukca_sfix, 'l_ukca_sfix')

  CASE (ukca_activation_jones)
    CALL check_config_value(l_ntpreq_n_activ_sum_from_ukca,                    &
                            l_ntpreq_n_activ_sum, 'l_ntpreq_n_activ_sum')
    CALL check_config_value(l_ntpreq_dryd_nuc_sol_from_ukca,                   &
                            l_ntpreq_dryd_nuc_sol, 'l_ntpreq_dryd_nuc_sol')

  END SELECT

  ! -- Check GLOMAP temporary logicals --

  CALL ukca_get_config(                                                        &
    l_fix_neg_pvol_wat=l_fix_neg_pvol_wat_from_ukca,                           &
    l_fix_ukca_impscav=l_fix_ukca_impscav_from_ukca,                           &
    l_fix_nacl_density=l_fix_nacl_density_from_ukca,                           &
    l_fix_ukca_activate_pdf=l_fix_ukca_activate_pdf_from_ukca,                 &
    l_fix_ukca_activate_vert_rep=l_fix_ukca_activate_vert_rep_from_ukca,       &
    l_bug_repro_tke_index=l_bug_repro_tke_index_from_ukca,                     &
    l_improve_aero_drydep=l_improve_aero_drydep_from_ukca,                     &
    l_fix_ukca_hygroscopicities=l_fix_ukca_hygroscopicities_from_ukca,         &
    l_fix_ukca_water_content=l_fix_ukca_water_content_from_ukca)

  CALL check_config_value(l_fix_neg_pvol_wat_from_ukca,                        &
                          l_fix_neg_pvol_wat, 'l_fix_neg_pvol_wat')
  CALL check_config_value(l_fix_ukca_impscav_from_ukca,                        &
                          l_fix_ukca_impscav, 'l_fix_ukca_impscav')
  CALL check_config_value(l_fix_nacl_density_from_ukca,                        &
                          l_fix_nacl_density, 'l_fix_nacl_density')
  CALL check_config_value(l_improve_aero_drydep_from_ukca,                     &
                          l_improve_aero_drydep, 'l_improve_aero_drydep')

  ! Check temporary logicals for Activate scheme
  IF (i_ukca_activation_scheme == ukca_activation_arg) THEN
    CALL check_config_value(l_fix_ukca_activate_pdf_from_ukca,                 &
                            l_fix_ukca_activate_pdf, 'l_fix_ukca_activate_pdf')
    CALL check_config_value(l_fix_ukca_activate_vert_rep_from_ukca,            &
                            l_fix_ukca_activate_vert_rep,                      &
                            'l_fix_ukca_activate_vert_rep')
    CALL check_config_value(l_bug_repro_tke_index_from_ukca,                   &
                            l_bug_repro_tke_index,                             &
                            'l_bug_repro_tke_index')
  END IF

  CALL check_config_value(l_fix_ukca_hygroscopicities_from_ukca,               &
                            l_fix_ukca_hygroscopicities,                       &
                            'l_fix_ukca_hygroscopicities')

  CALL check_config_value(l_fix_ukca_water_content_from_ukca,                  &
                            l_fix_ukca_water_content,                          &
                            'l_fix_ukca_water_content')
  
END IF ! l_ukca_mode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE check_ukca_configuration

! ----------------------------------------------------------------------
SUBROUTINE check_config_integer(ukca_value, um_value, ukca_varname)
! ----------------------------------------------------------------------
! Description:
!
! Check that the given UKCA value for the integer configuration variable
! with the specified name matches the given UM value.
! ----------------------------------------------------------------------

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: ukca_value   ! Value of UKCA configuration variable
INTEGER, INTENT(IN) :: um_value     ! Value of UM equivalent
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname
                                    ! Name of UKCA configuration variable

! Local variables
INTEGER :: errcode                              ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength) :: cmessage   ! Error return message

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_CONFIG_INTEGER'

errcode = 0
IF (ukca_value /= um_value) THEN
  errcode = 1
  WRITE(cmessage,'(A,A,A,I0,A)') 'Value of UKCA configuration variable ',      &
    ukca_varname, ' = ', ukca_value, ' does not match the UM requirement'
  CALL ereport(RoutineName, errcode, cmessage)
END IF

RETURN

END SUBROUTINE check_config_integer

! ----------------------------------------------------------------------
SUBROUTINE check_config_integer_vec(ukca_value, um_value, ukca_varname)
! ----------------------------------------------------------------------
! Description:
!
! Check that the given UKCA value for the integer vector configuration
! variable with the specified name matches the given UM value.
! ----------------------------------------------------------------------

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, ALLOCATABLE, INTENT(IN) :: ukca_value(:) ! Value of UKCA config. var.
INTEGER, ALLOCATABLE, INTENT(IN) :: um_value(:)   ! Value of UM equivalent
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname      ! Name of UKCA config. var.

! Local variables
INTEGER :: errcode                                ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength) :: cmessage     ! Error return message

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_CONFIG_INTEGER_VEC'

errcode = 0
IF (ALLOCATED(ukca_value) .OR. ALLOCATED(um_value)) THEN
  IF (ALLOCATED(ukca_value) .AND. ALLOCATED(um_value)) THEN
    ! Both allocated - further checks required
    IF (SIZE(ukca_value) /= SIZE(um_value)) THEN
      errcode = 1
      WRITE(cmessage,'(A,A,A,I0,A)') 'Size of UKCA configuration variable ',   &
        ukca_varname, ' = ', SIZE(ukca_value),                                 &
        ' does not match the UM requirement'
    ELSE IF (ANY(ukca_value(:) /= um_value(:))) THEN
      errcode = 1
      WRITE(cmessage,'(A,A,A)')                                                &
        'One or more values of UKCA configuration variable ',                  &
        ukca_varname, ' do not match the UM requirement'
    END IF
  ELSE IF (ALLOCATED(um_value)) THEN
    ! Only UM array is allocated
    errcode = 1
    WRITE(cmessage,'(A,A,A)') 'Expected UKCA configuration variable ',         &
      ukca_varname, ' is not allocated'
  ELSE
    ! Only UKCA array is allocated
    errcode = 1
    WRITE(cmessage,'(A,A,A,A)') ' UKCA configuration variable ',               &
      ukca_varname, ' is allocated but was not expected'
  END IF
END IF

IF (errcode /= 0) CALL ereport(RoutineName, errcode, cmessage)

RETURN

END SUBROUTINE check_config_integer_vec

! ----------------------------------------------------------------------
SUBROUTINE check_config_real(ukca_value, um_value, ukca_varname)
! ----------------------------------------------------------------------
! Description:
!
! Check that the given UKCA value for the real configuration variable
! with the specified name matches the given UM value.
! ----------------------------------------------------------------------

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
REAL, INTENT(IN) :: ukca_value      ! Value of UKCA configuration variable
REAL, INTENT(IN) :: um_value        ! Value of UM equivalent
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname
                                    ! Name of UKCA configuration variable

! Local variables
INTEGER :: errcode                              ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength) :: cmessage   ! Error return message

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_CONFIG_REAL'

errcode = 0
IF (ukca_value /= um_value) THEN
  errcode = 1
  WRITE(cmessage,'(A,A,A,E12.3,A)') 'Value of UKCA configuration variable ',   &
    ukca_varname, ' = ', ukca_value, ' does not match the UM requirement'
  CALL ereport(RoutineName, errcode, cmessage)
END IF

RETURN

END SUBROUTINE check_config_real

! ----------------------------------------------------------------------
SUBROUTINE check_config_logical(ukca_value, um_value, ukca_varname)
! ----------------------------------------------------------------------
! Description:
!
! Check that the given UKCA value for the configuration variable with
! the specified name matches the given UM value.
! ----------------------------------------------------------------------

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
LOGICAL, INTENT(IN) :: ukca_value      ! Value of UKCA configuration variable
LOGICAL, INTENT(IN) :: um_value        ! Value of UM equivalent
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname
                                       ! Name of UKCA configuration variable

! Local variables
INTEGER :: errcode                              ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength) :: cmessage   ! Error return message

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_CONFIG_LOGICAL'

errcode = 0
IF (ukca_value .NEQV. um_value) THEN
  errcode = 1
  WRITE(cmessage,'(A,A,A,L1,A)') 'Value of UKCA configuration variable ',      &
    ukca_varname, ' = ', ukca_value, ' does not match the UM requirement'
  CALL ereport(RoutineName, errcode, cmessage)
END IF

RETURN

END SUBROUTINE check_config_logical

! ----------------------------------------------------------------------
SUBROUTINE print_config_settings(                                              &
  ! Context information
  row_length,                                                                  &
  rows,                                                                        &
  model_levels,                                                                &
  bl_levels,                                                                   &
  nlev_ent_tr_mix,                                                             &
  ntype,                                                                       &
  npft,                                                                        &
  i_brd_leaf,                                                                  &
  i_brd_leaf_dec,                                                              &
  i_brd_leaf_eg_trop,                                                          &
  i_brd_leaf_eg_temp,                                                          &
  i_ndl_leaf,                                                                  &
  i_ndl_leaf_dec,                                                              &
  i_ndl_leaf_eg,                                                               &
  i_c3_grass,                                                                  &
  i_c3_crop,                                                                   &
  i_c3_pasture,                                                                &
  i_c4_grass,                                                                  &
  i_c4_crop,                                                                   &
  i_c4_pasture,                                                                &
  i_shrub,                                                                     &
  i_shrub_dec,                                                                 &
  i_shrub_eg,                                                                  &
  i_urban,                                                                     &
  i_lake,                                                                      &
  i_soil,                                                                      &
  i_ice,                                                                       &
  i_elev_ice,                                                                  &
  dzsoil_layer1,                                                               &
  l_cal360,                                                                    &
  timestep,                                                                    &
  ! General UKCA configuration options
  i_ukca_chem,                                                                 &
  l_ukca_chem_aero,                                                            &
  l_ukca_mode,                                                                 &
  l_ukca_ageair,                                                               &
  l_ukca_emissions_off,                                                        &
  l_fix_tropopause_level,                                                      &
  i_ageair_reset_method,                                                       &
  max_ageair_reset_level,                                                      &
  max_ageair_reset_height,                                                     &
  l_enable_diag_um,                                                            &
  l_ukca_persist_off,                                                          &
  l_timer,                                                                     &
  ! Chemistry configuration options
  i_ukca_chem_version,                                                         &
  nrsteps,                                                                     &
  chem_timestep,                                                               &
  dts0,                                                                        &
  l_ukca_asad_columns,                                                         &
  l_ukca_asad_full,                                                            &
  l_ukca_debug_asad,                                                           &
  l_ukca_intdd,                                                                &
  l_ukca_ddepo3_ocean,                                                         &
  l_ukca_ddep_lev1,                                                            &
  l_ukca_dry_dep_so2wet,                                                       &
  l_deposition_jules,                                                          &
  nit,                                                                         &
  l_ukca_quasinewton,                                                          &
  i_ukca_quasinewton_start,                                                    &
  i_ukca_quasinewton_end,                                                      &
  ukca_chem_seg_size,                                                          &
  max_z_for_offline_chem,                                                      &
  i_ukca_topboundary,                                                          &
  l_ukca_ro2_ntp,                                                              &
  l_ukca_ro2_perm,                                                             &
  ! Chemistry - Heterogeneous chemistry
  l_ukca_het_psc,                                                              &
  i_ukca_hetconfig,                                                            &
  l_ukca_limit_nat,                                                            &
  l_ukca_sa_clim,                                                              &
  l_ukca_trophet,                                                              &
  l_ukca_classic_hetchem,                                                      &
  ! Chemistry - Photolysis
  i_ukca_photol,                                                               &
  fastjx_mode,                                                                 &
  fastjx_prescutoff,                                                           &
  i_ukca_solcyc,                                                               &
  i_ukca_solcyc_start_year,                                                    &
  ! UKCA emissions configuration options
  l_ukca_ibvoc,                                                                &
  l_ukca_inferno,                                                              &
  l_ukca_inferno_ch4,                                                          &
  i_inferno_emi,                                                               &
  l_ukca_so2ems_expvolc,                                                       &
  l_ukca_so2ems_plumeria,                                                      &
  l_ukca_qch4inter,                                                            &
  l_ukca_emsdrvn_ch4,                                                          &
  mode_parfrac,                                                                &
  l_ukca_enable_seadms_ems,                                                    &
  i_ukca_dms_flux,                                                             &
  l_ukca_scale_seadms_ems,                                                     &
  seadms_ems_scaling,                                                          &
  l_ukca_linox_scaling,                                                        &
  lightnox_scale_fac,                                                          &
  i_ukca_light_param,                                                          &
  l_ukca_scale_soa_yield_mt,                                                   &
  soa_yield_scaling_mt,                                                        &
  l_ukca_scale_soa_yield_isop,                                                 &
  soa_yield_scaling_isop,                                                      &
  l_support_ems_vertprof,                                                      &
  l_support_ems_gridbox_units,                                                 &
  ! UKCA feedback configuration options
  l_ukca_h2o_feedback,                                                         &
  l_ukca_conserve_h,                                                           &
  ! UKCA environmental driver configuration options
  l_param_conv,                                                                &
  l_ctile,                                                                     &
  l_zon_av_ozone,                                                              &
  i_strat_lbc_source,                                                          &
  l_chem_environ_gas_scalars,                                                  &
  l_chem_environ_co2_fld,                                                      &
  l_ukca_prescribech4,                                                         &
  l_use_classic_so4,                                                           &
  l_use_classic_soot,                                                          &
  l_use_classic_ocff,                                                          &
  l_use_classic_biogenic,                                                      &
  l_use_classic_seasalt,                                                       &
  l_use_gridbox_volume,                                                        &
  l_use_gridbox_mass,                                                          &
  l_environ_z_top,                                                             &
  ! UKCA temporary logicals
  l_fix_improve_drydep,                                                        &
  l_fix_ukca_h2dd_x,                                                           &
  l_fix_drydep_so2_water,                                                      &
  l_fix_ukca_offox_h2o_fac,                                                    &
  l_fix_ukca_h2so4_ystore,                                                     &
  l_improve_aero_drydep,                                                       &
  ! General GLOMAP configuration options
  i_mode_nzts,                                                                 &
  ukca_mode_seg_size,                                                          &
  i_mode_setup,                                                                &
  l_mode_bhn_on,                                                               &
  l_mode_bln_on,                                                               &
  i_mode_bln_param_method,                                                     &
  i_mode_nucscav,                                                              &
  mode_activation_dryr,                                                        &
  mode_incld_so2_rfrac,                                                        &
  l_cv_rainout,                                                                &
  l_ddepaer,                                                                   &
  l_dust_mp_slinn_impc_scav,                                                   &
  ! GLOMAP emissions configuration options
  l_ukca_primss,                                                               &
  l_ukca_primsu,                                                               &
  l_ukca_primdu,                                                               &
  l_ukca_primbcoc,                                                             &
  l_ukca_prim_moc,                                                             &
  l_bcoc_bf,                                                                   &
  l_bcoc_bm,                                                                   &
  l_bcoc_ff,                                                                   &
  l_ukca_scale_biom_aer_ems,                                                   &
  biom_aer_ems_scaling,                                                        &
  l_ukca_fine_no3_prod,                                                        &
  l_ukca_coarse_no3_prod,                                                      &
  l_no3_prod_in_aero_step,                                                     &
  l_ukca_scale_sea_salt_ems,                                                   &
  sea_salt_ems_scaling,                                                        &
  l_ukca_scale_marine_pom_ems,                                                 &
  marine_pom_ems_scaling,                                                      &
  ! GLOMAP feedback configuration options
  l_ukca_radaer,                                                               &
  i_ukca_tune_bc,                                                              &
  i_ukca_activation_scheme,                                                    &
  i_ukca_nwbins,                                                               &
  sigwmin,                                                                     &
  l_ntpreq_n_activ_sum,                                                        &
  l_ntpreq_dryd_nuc_sol,                                                       &
  l_ukca_sfix,                                                                 &
  ! GLOMAP temporary logicals
  l_fix_neg_pvol_wat,                                                          &
  l_fix_ukca_impscav,                                                          &
  l_fix_nacl_density,                                                          &
  l_fix_ukca_activate_pdf,                                                     &
  l_fix_ukca_activate_vert_rep,                                                &
  l_bug_repro_tke_index,                                                       &
  l_fix_ukca_hygroscopicities,                                                 &
  l_fix_ukca_water_content,                                                    &
  ! New settings in Ticket #25
  l_ukca_drydep_off,                                                           &
  l_ukca_wetdep_off,                                                           &
  nlev_above_trop_o3_env,                                                      &
  nlev_ch4_stratloss,                                                          &
  l_tracer_lumping,                                                            &
  env_log_step,                                                                &
  l_aero_rainout,                                                                   &
  l_impc_scav)
! ----------------------------------------------------------------------
! Description:
!   Write given settings and UKCA equivalents for comparison
! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER, OPTIONAL, INTENT(IN) ::  row_length
INTEGER, OPTIONAL, INTENT(IN) ::  rows
INTEGER, OPTIONAL, INTENT(IN) ::  model_levels
INTEGER, OPTIONAL, INTENT(IN) ::  bl_levels
INTEGER, OPTIONAL, INTENT(IN) ::  nlev_ent_tr_mix
INTEGER, OPTIONAL, INTENT(IN) ::  ntype
INTEGER, OPTIONAL, INTENT(IN) ::  npft
INTEGER, OPTIONAL, INTENT(IN) ::  i_brd_leaf
INTEGER, OPTIONAL, INTENT(IN) ::  i_brd_leaf_dec
INTEGER, OPTIONAL, INTENT(IN) ::  i_brd_leaf_eg_trop
INTEGER, OPTIONAL, INTENT(IN) ::  i_brd_leaf_eg_temp
INTEGER, OPTIONAL, INTENT(IN) ::  i_ndl_leaf
INTEGER, OPTIONAL, INTENT(IN) ::  i_ndl_leaf_dec
INTEGER, OPTIONAL, INTENT(IN) ::  i_ndl_leaf_eg
INTEGER, OPTIONAL, INTENT(IN) ::  i_c3_grass
INTEGER, OPTIONAL, INTENT(IN) ::  i_c3_crop
INTEGER, OPTIONAL, INTENT(IN) ::  i_c3_pasture
INTEGER, OPTIONAL, INTENT(IN) ::  i_c4_grass
INTEGER, OPTIONAL, INTENT(IN) ::  i_c4_crop
INTEGER, OPTIONAL, INTENT(IN) ::  i_c4_pasture
INTEGER, OPTIONAL, INTENT(IN) ::  i_shrub
INTEGER, OPTIONAL, INTENT(IN) ::  i_shrub_dec
INTEGER, OPTIONAL, INTENT(IN) ::  i_shrub_eg
INTEGER, OPTIONAL, INTENT(IN) ::  i_urban
INTEGER, OPTIONAL, INTENT(IN) ::  i_lake
INTEGER, OPTIONAL, INTENT(IN) ::  i_soil
INTEGER, OPTIONAL, INTENT(IN) ::  i_ice
INTEGER, ALLOCATABLE, OPTIONAL, INTENT(IN) ::  i_elev_ice(:)
REAL, OPTIONAL, INTENT(IN) ::  dzsoil_layer1
LOGICAL, OPTIONAL, INTENT(IN) ::  l_cal360
REAL, OPTIONAL, INTENT(IN) ::  timestep
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_chem
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_chem_aero
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_mode
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_ageair
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_emissions_off
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_tropopause_level
INTEGER, OPTIONAL, INTENT(IN) ::  i_ageair_reset_method
INTEGER, OPTIONAL, INTENT(IN) ::  max_ageair_reset_level
REAL, OPTIONAL, INTENT(IN) ::  max_ageair_reset_height
LOGICAL, OPTIONAL, INTENT(IN) ::  l_enable_diag_um
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_persist_off
LOGICAL, OPTIONAL, INTENT(IN) ::  l_timer
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_chem_version
INTEGER, OPTIONAL, INTENT(IN) ::  nrsteps
INTEGER, OPTIONAL, INTENT(IN) ::  chem_timestep
INTEGER, OPTIONAL, INTENT(IN) ::  dts0
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_asad_columns
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_asad_full
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_debug_asad
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_intdd
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_ddepo3_ocean
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_ddep_lev1
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_dry_dep_so2wet
LOGICAL, OPTIONAL, INTENT(IN) ::  l_deposition_jules
INTEGER, OPTIONAL, INTENT(IN) ::  nit
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_quasinewton
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_quasinewton_start
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_quasinewton_end
INTEGER, OPTIONAL, INTENT(IN) ::  ukca_chem_seg_size
REAL, OPTIONAL, INTENT(IN) ::  max_z_for_offline_chem
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_topboundary
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_ro2_ntp
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_ro2_perm
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_het_psc
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_hetconfig
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_limit_nat
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_sa_clim
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_trophet
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_classic_hetchem
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_photol
INTEGER, OPTIONAL, INTENT(IN) ::  fastjx_mode
REAL, OPTIONAL, INTENT(IN) ::  fastjx_prescutoff
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_solcyc
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_solcyc_start_year
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_ibvoc
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_inferno
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_inferno_ch4
INTEGER, OPTIONAL, INTENT(IN) ::  i_inferno_emi
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_so2ems_expvolc
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_so2ems_plumeria
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_qch4inter
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_emsdrvn_ch4
REAL, OPTIONAL, INTENT(IN) ::  mode_parfrac
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_enable_seadms_ems
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_dms_flux
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_scale_seadms_ems
REAL, OPTIONAL, INTENT(IN) ::  seadms_ems_scaling
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_linox_scaling
REAL, OPTIONAL, INTENT(IN) ::  lightnox_scale_fac
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_light_param
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_scale_soa_yield_mt
REAL, OPTIONAL, INTENT(IN) ::  soa_yield_scaling_mt
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_scale_soa_yield_isop
REAL, OPTIONAL, INTENT(IN) ::  soa_yield_scaling_isop
LOGICAL, OPTIONAL, INTENT(IN) ::  l_support_ems_vertprof
LOGICAL, OPTIONAL, INTENT(IN) ::  l_support_ems_gridbox_units
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_h2o_feedback
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_conserve_h
LOGICAL, OPTIONAL, INTENT(IN) ::  l_param_conv
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ctile
LOGICAL, OPTIONAL, INTENT(IN) ::  l_zon_av_ozone
INTEGER, OPTIONAL, INTENT(IN) ::  i_strat_lbc_source
LOGICAL, OPTIONAL, INTENT(IN) ::  l_chem_environ_gas_scalars
LOGICAL, OPTIONAL, INTENT(IN) ::  l_chem_environ_co2_fld
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_prescribech4
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_classic_so4
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_classic_soot
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_classic_ocff
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_classic_biogenic
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_classic_seasalt
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_gridbox_volume
LOGICAL, OPTIONAL, INTENT(IN) ::  l_use_gridbox_mass
LOGICAL, OPTIONAL, INTENT(IN) ::  l_environ_z_top
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_improve_drydep
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_h2dd_x
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_drydep_so2_water
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_offox_h2o_fac
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_h2so4_ystore
LOGICAL, OPTIONAL, INTENT(IN) ::  l_improve_aero_drydep
INTEGER, OPTIONAL, INTENT(IN) ::  i_mode_nzts
INTEGER, OPTIONAL, INTENT(IN) ::  ukca_mode_seg_size
INTEGER, OPTIONAL, INTENT(IN) ::  i_mode_setup
LOGICAL, OPTIONAL, INTENT(IN) ::  l_mode_bhn_on
LOGICAL, OPTIONAL, INTENT(IN) ::  l_mode_bln_on
INTEGER, OPTIONAL, INTENT(IN) ::  i_mode_bln_param_method
INTEGER, OPTIONAL, INTENT(IN) ::  i_mode_nucscav
REAL, OPTIONAL, INTENT(IN) ::  mode_activation_dryr
REAL, OPTIONAL, INTENT(IN) ::  mode_incld_so2_rfrac
LOGICAL, OPTIONAL, INTENT(IN) ::  l_cv_rainout
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ddepaer
LOGICAL, OPTIONAL, INTENT(IN) ::  l_dust_mp_slinn_impc_scav
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_primss
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_primsu
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_primdu
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_primbcoc
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_prim_moc
LOGICAL, OPTIONAL, INTENT(IN) ::  l_bcoc_bf
LOGICAL, OPTIONAL, INTENT(IN) ::  l_bcoc_bm
LOGICAL, OPTIONAL, INTENT(IN) ::  l_bcoc_ff
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_scale_biom_aer_ems
REAL, OPTIONAL, INTENT(IN) ::  biom_aer_ems_scaling
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_fine_no3_prod
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_coarse_no3_prod
LOGICAL, OPTIONAL, INTENT(IN) ::  l_no3_prod_in_aero_step
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_scale_sea_salt_ems
REAL, OPTIONAL, INTENT(IN) ::  sea_salt_ems_scaling
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_scale_marine_pom_ems
REAL, OPTIONAL, INTENT(IN) ::  marine_pom_ems_scaling
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_radaer
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_tune_bc
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_activation_scheme
INTEGER, OPTIONAL, INTENT(IN) ::  i_ukca_nwbins
REAL, OPTIONAL, INTENT(IN) ::  sigwmin
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ntpreq_n_activ_sum
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ntpreq_dryd_nuc_sol
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_sfix
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_neg_pvol_wat
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_impscav
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_nacl_density
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_activate_pdf
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_activate_vert_rep
LOGICAL, OPTIONAL, INTENT(IN) ::  l_bug_repro_tke_index
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_hygroscopicities
LOGICAL, OPTIONAL, INTENT(IN) ::  l_fix_ukca_water_content
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_drydep_off
LOGICAL, OPTIONAL, INTENT(IN) ::  l_ukca_wetdep_off
INTEGER, OPTIONAL, INTENT(IN) ::  nlev_above_trop_o3_env
INTEGER, OPTIONAL, INTENT(IN) ::  nlev_ch4_stratloss
LOGICAL, OPTIONAL, INTENT(IN) ::  l_tracer_lumping
INTEGER, OPTIONAL, INTENT(IN) ::  env_log_step
LOGICAL, OPTIONAL, INTENT(IN) ::  l_aero_rainout
LOGICAL, OPTIONAL, INTENT(IN) ::  l_impc_scav

INTEGER ::  row_length_from_ukca
INTEGER ::  rows_from_ukca
INTEGER ::  model_levels_from_ukca
INTEGER ::  bl_levels_from_ukca
INTEGER ::  nlev_ent_tr_mix_from_ukca
INTEGER ::  ntype_from_ukca
INTEGER ::  npft_from_ukca
INTEGER ::  i_brd_leaf_from_ukca
INTEGER ::  i_brd_leaf_dec_from_ukca
INTEGER ::  i_brd_leaf_eg_trop_from_ukca
INTEGER ::  i_brd_leaf_eg_temp_from_ukca
INTEGER ::  i_ndl_leaf_from_ukca
INTEGER ::  i_ndl_leaf_dec_from_ukca
INTEGER ::  i_ndl_leaf_eg_from_ukca
INTEGER ::  i_c3_grass_from_ukca
INTEGER ::  i_c3_crop_from_ukca
INTEGER ::  i_c3_pasture_from_ukca
INTEGER ::  i_c4_grass_from_ukca
INTEGER ::  i_c4_crop_from_ukca
INTEGER ::  i_c4_pasture_from_ukca
INTEGER ::  i_shrub_from_ukca
INTEGER ::  i_shrub_dec_from_ukca
INTEGER ::  i_shrub_eg_from_ukca
INTEGER ::  i_urban_from_ukca
INTEGER ::  i_lake_from_ukca
INTEGER ::  i_soil_from_ukca
INTEGER ::  i_ice_from_ukca
INTEGER, ALLOCATABLE ::  i_elev_ice_from_ukca(:)
REAL ::  dzsoil_layer1_from_ukca
LOGICAL ::  l_cal360_from_ukca
REAL ::  timestep_from_ukca
INTEGER ::  i_ukca_chem_from_ukca
LOGICAL ::  l_ukca_chem_aero_from_ukca
LOGICAL ::  l_ukca_mode_from_ukca
LOGICAL ::  l_ukca_ageair_from_ukca
LOGICAL ::  l_ukca_emissions_off_from_ukca
LOGICAL ::  l_fix_tropopause_level_from_ukca
INTEGER ::  i_ageair_reset_method_from_ukca
INTEGER ::  max_ageair_reset_level_from_ukca
REAL ::  max_ageair_reset_height_from_ukca
LOGICAL ::  l_enable_diag_um_from_ukca
LOGICAL ::  l_ukca_persist_off_from_ukca
LOGICAL ::  l_timer_from_ukca
INTEGER ::  i_ukca_chem_version_from_ukca
INTEGER ::  nrsteps_from_ukca
INTEGER ::  chem_timestep_from_ukca
INTEGER ::  dts0_from_ukca
LOGICAL ::  l_ukca_asad_columns_from_ukca
LOGICAL ::  l_ukca_asad_full_from_ukca
LOGICAL ::  l_ukca_debug_asad_from_ukca
LOGICAL :: l_ukca_intdd_from_ukca
LOGICAL ::  l_ukca_ddepo3_ocean_from_ukca
LOGICAL ::  l_ukca_ddep_lev1_from_ukca
LOGICAL ::  l_ukca_dry_dep_so2wet_from_ukca
LOGICAL ::  l_deposition_jules_from_ukca
INTEGER ::  nit_from_ukca
LOGICAL ::  l_ukca_quasinewton_from_ukca
INTEGER ::  i_ukca_quasinewton_start_from_ukca
INTEGER ::  i_ukca_quasinewton_end_from_ukca
INTEGER ::  ukca_chem_seg_size_from_ukca
REAL ::  max_z_for_offline_chem_from_ukca
INTEGER ::  i_ukca_topboundary_from_ukca
LOGICAL ::  l_ukca_ro2_ntp_from_ukca
LOGICAL ::  l_ukca_ro2_perm_from_ukca
LOGICAL ::  l_ukca_het_psc_from_ukca
INTEGER ::  i_ukca_hetconfig_from_ukca
LOGICAL ::  l_ukca_limit_nat_from_ukca
LOGICAL ::  l_ukca_sa_clim_from_ukca
LOGICAL ::  l_ukca_trophet_from_ukca
LOGICAL ::  l_ukca_classic_hetchem_from_ukca
INTEGER ::  i_scheme_from_photol
INTEGER ::  fastjx_mode_from_photol
REAL ::  fastjx_prescutoff_from_photol
INTEGER ::  i_solcyc_from_photol
INTEGER ::  i_solcyc_start_year_from_photol
LOGICAL ::  l_ukca_ibvoc_from_ukca
LOGICAL ::  l_ukca_inferno_from_ukca
LOGICAL ::  l_ukca_inferno_ch4_from_ukca
INTEGER ::  i_inferno_emi_from_ukca
LOGICAL ::  l_ukca_so2ems_expvolc_from_ukca
LOGICAL ::  l_ukca_so2ems_plumeria_from_ukca
LOGICAL ::  l_ukca_qch4inter_from_ukca
LOGICAL ::  l_ukca_emsdrvn_ch4_from_ukca
REAL ::  mode_parfrac_from_ukca
LOGICAL ::  l_ukca_enable_seadms_ems_from_ukca
INTEGER ::  i_ukca_dms_flux_from_ukca
LOGICAL ::  l_ukca_scale_seadms_ems_from_ukca
REAL ::  seadms_ems_scaling_from_ukca
LOGICAL ::  l_ukca_linox_scaling_from_ukca
REAL ::  lightnox_scale_fac_from_ukca
INTEGER ::  i_ukca_light_param_from_ukca
LOGICAL :: l_ukca_scale_soa_yield_mt_from_ukca
LOGICAL :: l_ukca_scale_soa_yield_isop_from_ukca
REAL :: soa_yield_scaling_mt_from_ukca
REAL :: soa_yield_scaling_isop_from_ukca
LOGICAL ::  l_support_ems_vertprof_from_ukca
LOGICAL ::  l_support_ems_gridbox_units_from_ukca
LOGICAL ::  l_ukca_h2o_feedback_from_ukca
LOGICAL ::  l_ukca_conserve_h_from_ukca
LOGICAL ::  l_param_conv_from_ukca
LOGICAL ::  l_ctile_from_ukca
LOGICAL ::  l_zon_av_ozone_from_ukca
INTEGER ::  i_strat_lbc_source_from_ukca
LOGICAL ::  l_chem_environ_gas_scalars_from_ukca
LOGICAL ::  l_chem_environ_co2_fld_from_ukca
LOGICAL ::  l_ukca_prescribech4_from_ukca
LOGICAL ::  l_use_classic_so4_from_ukca
LOGICAL ::  l_use_classic_soot_from_ukca
LOGICAL ::  l_use_classic_ocff_from_ukca
LOGICAL ::  l_use_classic_biogenic_from_ukca
LOGICAL ::  l_use_classic_seasalt_from_ukca
LOGICAL ::  l_use_gridbox_volume_from_ukca
LOGICAL ::  l_use_gridbox_mass_from_ukca
LOGICAL ::  l_environ_z_top_from_ukca
LOGICAL ::  l_fix_improve_drydep_from_ukca
LOGICAL ::  l_fix_ukca_h2dd_x_from_ukca
LOGICAL ::  l_fix_drydep_so2_water_from_ukca
LOGICAL ::  l_fix_ukca_offox_h2o_fac_from_ukca
LOGICAL ::  l_fix_ukca_h2so4_ystore_from_ukca
LOGICAL ::  l_improve_aero_drydep_from_ukca
INTEGER ::  i_mode_nzts_from_ukca
INTEGER ::  ukca_mode_seg_size_from_ukca
INTEGER ::  i_mode_setup_from_ukca
LOGICAL ::  l_mode_bhn_on_from_ukca
LOGICAL ::  l_mode_bln_on_from_ukca
INTEGER ::  i_mode_bln_param_method_from_ukca
INTEGER ::  i_mode_nucscav_from_ukca
REAL ::  mode_activation_dryr_from_ukca
REAL ::  mode_incld_so2_rfrac_from_ukca
LOGICAL ::  l_cv_rainout_from_ukca
LOGICAL ::  l_ddepaer_from_ukca
LOGICAL ::  l_dust_mp_slinn_impc_scav_from_ukca
LOGICAL ::  l_ukca_primss_from_ukca
LOGICAL ::  l_ukca_primsu_from_ukca
LOGICAL ::  l_ukca_primdu_from_ukca
LOGICAL ::  l_ukca_primbcoc_from_ukca
LOGICAL ::  l_ukca_prim_moc_from_ukca
LOGICAL ::  l_bcoc_bf_from_ukca
LOGICAL ::  l_bcoc_bm_from_ukca
LOGICAL ::  l_bcoc_ff_from_ukca
LOGICAL ::  l_ukca_scale_biom_aer_ems_from_ukca
REAL ::  biom_aer_ems_scaling_from_ukca
LOGICAL ::  l_ukca_fine_no3_prod_from_ukca
LOGICAL ::  l_ukca_coarse_no3_prod_from_ukca
LOGICAL :: l_no3_prod_in_aero_step_from_ukca
LOGICAL ::  l_ukca_scale_sea_salt_ems_from_ukca
REAL ::  sea_salt_ems_scaling_from_ukca
LOGICAL ::  l_ukca_scale_marine_pom_ems_from_ukca
REAL ::  marine_pom_ems_scaling_from_ukca
LOGICAL ::  l_ukca_radaer_from_ukca
INTEGER ::  i_ukca_tune_bc_from_ukca
INTEGER ::  i_ukca_activation_scheme_from_ukca
INTEGER ::  i_ukca_nwbins_from_ukca
REAL ::  sigwmin_from_ukca
LOGICAL ::  l_ntpreq_n_activ_sum_from_ukca
LOGICAL ::  l_ntpreq_dryd_nuc_sol_from_ukca
LOGICAL ::  l_ukca_sfix_from_ukca
LOGICAL ::  l_fix_neg_pvol_wat_from_ukca
LOGICAL ::  l_fix_ukca_impscav_from_ukca
LOGICAL ::  l_fix_nacl_density_from_ukca
LOGICAL ::  l_fix_ukca_activate_pdf_from_ukca
LOGICAL ::  l_fix_ukca_activate_vert_rep_from_ukca
LOGICAL ::  l_bug_repro_tke_index_from_ukca
LOGICAL ::  l_fix_ukca_hygroscopicities_from_ukca
LOGICAL ::  l_fix_ukca_water_content_from_ukca
LOGICAL ::  l_ukca_drydep_off_from_ukca
LOGICAL ::  l_ukca_wetdep_off_from_ukca
INTEGER ::  nlev_above_trop_o3_env_from_ukca
INTEGER ::  nlev_ch4_stratloss_from_ukca
LOGICAL ::  l_tracer_lumping_from_ukca
INTEGER ::  env_log_step_from_ukca
LOGICAL ::  l_aero_rainout_from_ukca
LOGICAL ::  l_impc_scav_from_ukca

CALL ukca_get_config(                                                          &
  ! Context information
  row_length=row_length_from_ukca,                                             &
  rows=rows_from_ukca,                                                         &
  model_levels=model_levels_from_ukca,                                         &
  bl_levels=bl_levels_from_ukca,                                               &
  nlev_ent_tr_mix=nlev_ent_tr_mix_from_ukca,                                   &
  ntype=ntype_from_ukca,                                                       &
  npft=npft_from_ukca,                                                         &
  i_brd_leaf=i_brd_leaf_from_ukca,                                             &
  i_brd_leaf_dec=i_brd_leaf_dec_from_ukca,                                     &
  i_brd_leaf_eg_trop=i_brd_leaf_eg_trop_from_ukca,                             &
  i_brd_leaf_eg_temp=i_brd_leaf_eg_temp_from_ukca,                             &
  i_ndl_leaf=i_ndl_leaf_from_ukca,                                             &
  i_ndl_leaf_dec=i_ndl_leaf_dec_from_ukca,                                     &
  i_ndl_leaf_eg=i_ndl_leaf_eg_from_ukca,                                       &
  i_c3_grass=i_c3_grass_from_ukca,                                             &
  i_c3_crop=i_c3_crop_from_ukca,                                               &
  i_c3_pasture=i_c3_pasture_from_ukca,                                         &
  i_c4_grass=i_c4_grass_from_ukca,                                             &
  i_c4_crop=i_c4_crop_from_ukca,                                               &
  i_c4_pasture=i_c4_pasture_from_ukca,                                         &
  i_shrub=i_shrub_from_ukca,                                                   &
  i_shrub_dec=i_shrub_dec_from_ukca,                                           &
  i_shrub_eg=i_shrub_eg_from_ukca,                                             &
  i_urban=i_urban_from_ukca,                                                   &
  i_lake=i_lake_from_ukca,                                                     &
  i_soil=i_soil_from_ukca,                                                     &
  i_ice=i_ice_from_ukca,                                                       &
  i_elev_ice=i_elev_ice_from_ukca,                                             &
  dzsoil_layer1=dzsoil_layer1_from_ukca,                                       &
  l_cal360=l_cal360_from_ukca,                                                 &
  timestep=timestep_from_ukca,                                                 &
  ! General UKCA configuration options
  i_ukca_chem=i_ukca_chem_from_ukca,                                           &
  l_ukca_chem_aero=l_ukca_chem_aero_from_ukca,                                 &
  l_ukca_mode= l_ukca_mode_from_ukca,                                          &
  l_ukca_ageair=l_ukca_ageair_from_ukca,                                       &
  l_ukca_emissions_off=l_ukca_emissions_off_from_ukca,                         &
  l_fix_tropopause_level=l_fix_tropopause_level_from_ukca,                     &
  i_ageair_reset_method=i_ageair_reset_method_from_ukca,                       &
  max_ageair_reset_level=max_ageair_reset_level_from_ukca,                     &
  max_ageair_reset_height=max_ageair_reset_height_from_ukca,                   &
  l_enable_diag_um=l_enable_diag_um_from_ukca,                                 &
  l_ukca_persist_off=l_ukca_persist_off_from_ukca,                             &
  l_timer=l_timer_from_ukca,                                                   &
  ! Chemistry configuration options
  i_ukca_chem_version=i_ukca_chem_version_from_ukca,                           &
  nrsteps=nrsteps_from_ukca,                                                   &
  chem_timestep=chem_timestep_from_ukca,                                       &
  dts0=dts0_from_ukca,                                                         &
  l_ukca_asad_columns=l_ukca_asad_columns_from_ukca,                           &
  l_ukca_asad_full=l_ukca_asad_full_from_ukca,                                 &
  l_ukca_debug_asad=l_ukca_debug_asad_from_ukca,                               &
  l_ukca_intdd=l_ukca_intdd_from_ukca,                                         &
  l_ukca_ddepo3_ocean=l_ukca_ddepo3_ocean_from_ukca,                           &
  l_ukca_ddep_lev1=l_ukca_ddep_lev1_from_ukca,                                 &
  l_ukca_dry_dep_so2wet=l_ukca_dry_dep_so2wet_from_ukca,                       &
  l_deposition_jules=l_deposition_jules_from_ukca,                             &
  nit=nit_from_ukca,                                                           &
  l_ukca_quasinewton=l_ukca_quasinewton_from_ukca,                             &
  i_ukca_quasinewton_start=i_ukca_quasinewton_start_from_ukca,                 &
  i_ukca_quasinewton_end=i_ukca_quasinewton_end_from_ukca,                     &
  ukca_chem_seg_size=ukca_chem_seg_size_from_ukca,                             &
  max_z_for_offline_chem=max_z_for_offline_chem_from_ukca,                     &
  i_ukca_topboundary=i_ukca_topboundary_from_ukca,                             &
  l_ukca_ro2_ntp=l_ukca_ro2_ntp_from_ukca,                                     &
  l_ukca_ro2_perm=l_ukca_ro2_perm_from_ukca,                                   &
  ! Chemistry - Heterogeneous chemistry
  l_ukca_het_psc=l_ukca_het_psc_from_ukca,                                     &
  i_ukca_hetconfig=i_ukca_hetconfig_from_ukca,                                 &
  l_ukca_limit_nat=l_ukca_limit_nat_from_ukca,                                 &
  l_ukca_sa_clim=l_ukca_sa_clim_from_ukca,                                     &
  l_ukca_trophet=l_ukca_trophet_from_ukca,                                     &
  l_ukca_classic_hetchem=l_ukca_classic_hetchem_from_ukca,                     &
  ! UKCA emissions configuration options
  l_ukca_ibvoc=l_ukca_ibvoc_from_ukca,                                         &
  l_ukca_inferno=l_ukca_inferno_from_ukca,                                     &
  l_ukca_inferno_ch4=l_ukca_inferno_ch4_from_ukca,                             &
  i_inferno_emi=i_inferno_emi_from_ukca,                                       &
  l_ukca_so2ems_expvolc=l_ukca_so2ems_expvolc_from_ukca,                       &
  l_ukca_so2ems_plumeria=l_ukca_so2ems_plumeria_from_ukca,                     &
  l_ukca_qch4inter=l_ukca_qch4inter_from_ukca,                                 &
  l_ukca_emsdrvn_ch4=l_ukca_emsdrvn_ch4_from_ukca,                             &
  mode_parfrac=mode_parfrac_from_ukca,                                         &
  l_ukca_enable_seadms_ems=l_ukca_enable_seadms_ems_from_ukca,                 &
  i_ukca_dms_flux=i_ukca_dms_flux_from_ukca,                                   &
  l_ukca_scale_seadms_ems=l_ukca_scale_seadms_ems_from_ukca,                   &
  seadms_ems_scaling=seadms_ems_scaling_from_ukca,                             &
  l_ukca_linox_scaling=l_ukca_linox_scaling_from_ukca,                         &
  lightnox_scale_fac=lightnox_scale_fac_from_ukca,                             &
  i_ukca_light_param=i_ukca_light_param_from_ukca,                             &
  l_ukca_scale_soa_yield_mt=l_ukca_scale_soa_yield_mt_from_ukca,               &
  soa_yield_scaling_mt=soa_yield_scaling_mt_from_ukca,                         &
  l_ukca_scale_soa_yield_isop=l_ukca_scale_soa_yield_isop_from_ukca,           &
  soa_yield_scaling_isop=soa_yield_scaling_isop_from_ukca,                     &
  l_support_ems_vertprof=l_support_ems_vertprof_from_ukca,                     &
  l_support_ems_gridbox_units=l_support_ems_gridbox_units_from_ukca,           &
  ! UKCA feedback configuration options
  l_ukca_h2o_feedback=l_ukca_h2o_feedback_from_ukca,                           &
  l_ukca_conserve_h=l_ukca_conserve_h_from_ukca,                               &
  ! UKCA environmental driver configuration options
  l_param_conv=l_param_conv_from_ukca,                                         &
  l_ctile=l_ctile_from_ukca,                                                   &
  l_zon_av_ozone=l_zon_av_ozone_from_ukca,                                     &
  i_strat_lbc_source=i_strat_lbc_source_from_ukca,                             &
  l_chem_environ_gas_scalars=l_chem_environ_gas_scalars_from_ukca,             &
  l_chem_environ_co2_fld=l_chem_environ_co2_fld_from_ukca,                     &
  l_ukca_prescribech4=l_ukca_prescribech4_from_ukca,                           &
  l_use_classic_so4=l_use_classic_so4_from_ukca,                               &
  l_use_classic_soot=l_use_classic_soot_from_ukca,                             &
  l_use_classic_ocff=l_use_classic_ocff_from_ukca,                             &
  l_use_classic_biogenic=l_use_classic_biogenic_from_ukca,                     &
  l_use_classic_seasalt=l_use_classic_seasalt_from_ukca,                       &
  l_use_gridbox_volume=l_use_gridbox_volume_from_ukca,                         &
  l_use_gridbox_mass=l_use_gridbox_mass_from_ukca,                             &
  l_environ_z_top=l_environ_z_top_from_ukca,                                   &
  ! UKCA temporary logicals
  l_fix_improve_drydep=l_fix_improve_drydep_from_ukca,                         &
  l_fix_ukca_h2dd_x=l_fix_ukca_h2dd_x_from_ukca,                               &
  l_fix_drydep_so2_water=l_fix_drydep_so2_water_from_ukca,                     &
  l_fix_ukca_offox_h2o_fac=l_fix_ukca_offox_h2o_fac_from_ukca,                 &
  l_fix_ukca_h2so4_ystore=l_fix_ukca_h2so4_ystore_from_ukca,                   &
  l_improve_aero_drydep=l_improve_aero_drydep_from_ukca,                       &
  ! General GLOMAP configuration options
  i_mode_nzts=i_mode_nzts_from_ukca,                                           &
  ukca_mode_seg_size=ukca_mode_seg_size_from_ukca,                             &
  i_mode_setup=i_mode_setup_from_ukca,                                         &
  l_mode_bhn_on=l_mode_bhn_on_from_ukca,                                       &
  l_mode_bln_on=l_mode_bln_on_from_ukca,                                       &
  i_mode_bln_param_method=i_mode_bln_param_method_from_ukca,                   &
  i_mode_nucscav=i_mode_nucscav_from_ukca,                                     &
  mode_activation_dryr=mode_activation_dryr_from_ukca,                         &
  mode_incld_so2_rfrac=mode_incld_so2_rfrac_from_ukca,                         &
  l_cv_rainout=l_cv_rainout_from_ukca,                                         &
  l_ddepaer=l_ddepaer_from_ukca,                                               &
  l_dust_mp_slinn_impc_scav=l_dust_mp_slinn_impc_scav_from_ukca,               &
  ! GLOMAP emissions configuration options
  l_ukca_primss=l_ukca_primss_from_ukca,                                       &
  l_ukca_primsu=l_ukca_primsu_from_ukca,                                       &
  l_ukca_primdu=l_ukca_primdu_from_ukca,                                       &
  l_ukca_primbcoc=l_ukca_primbcoc_from_ukca,                                   &
  l_ukca_prim_moc=l_ukca_prim_moc_from_ukca,                                   &
  l_bcoc_bf=l_bcoc_bf_from_ukca,                                               &
  l_bcoc_bm=l_bcoc_bm_from_ukca,                                               &
  l_bcoc_ff=l_bcoc_ff_from_ukca,                                               &
  l_ukca_scale_biom_aer_ems=l_ukca_scale_biom_aer_ems_from_ukca,               &
  biom_aer_ems_scaling=biom_aer_ems_scaling_from_ukca,                         &
  l_ukca_fine_no3_prod=l_ukca_fine_no3_prod_from_ukca,                         &
  l_ukca_coarse_no3_prod=l_ukca_coarse_no3_prod_from_ukca,                     &
  l_no3_prod_in_aero_step=l_no3_prod_in_aero_step_from_ukca,                   &
  l_ukca_scale_sea_salt_ems=l_ukca_scale_sea_salt_ems_from_ukca,               &
  sea_salt_ems_scaling=sea_salt_ems_scaling_from_ukca,                         &
  l_ukca_scale_marine_pom_ems=l_ukca_scale_marine_pom_ems_from_ukca,           &
  marine_pom_ems_scaling=marine_pom_ems_scaling_from_ukca,                     &
  ! GLOMAP feedback configuration options
  l_ukca_radaer=l_ukca_radaer_from_ukca,                                       &
  i_ukca_tune_bc=i_ukca_tune_bc_from_ukca,                                     &
  i_ukca_activation_scheme=i_ukca_activation_scheme_from_ukca,                 &
  i_ukca_nwbins=i_ukca_nwbins_from_ukca,                                       &
  sigwmin=sigwmin_from_ukca,                                                   &
  l_ntpreq_n_activ_sum=l_ntpreq_n_activ_sum_from_ukca,                         &
  l_ntpreq_dryd_nuc_sol=l_ntpreq_dryd_nuc_sol_from_ukca,                       &
  l_ukca_sfix=l_ukca_sfix_from_ukca,                                           &
  ! GLOMAP temporary logicals
  l_fix_neg_pvol_wat=l_fix_neg_pvol_wat_from_ukca,                             &
  l_fix_ukca_impscav=l_fix_ukca_impscav_from_ukca,                             &
  l_fix_nacl_density=l_fix_nacl_density_from_ukca,                             &
  l_fix_ukca_activate_pdf=l_fix_ukca_activate_pdf_from_ukca,                   &
  l_fix_ukca_activate_vert_rep=l_fix_ukca_activate_vert_rep_from_ukca,         &
  l_bug_repro_tke_index=l_bug_repro_tke_index_from_ukca,                       &
  l_fix_ukca_hygroscopicities = l_fix_ukca_hygroscopicities_from_ukca,         &
  l_fix_ukca_water_content = l_fix_ukca_water_content_from_ukca,               &
  ! New settings in Ticket #25
  l_ukca_drydep_off=l_ukca_drydep_off_from_ukca,                               &
  l_ukca_wetdep_off=l_ukca_wetdep_off_from_ukca,                               &
  nlev_above_trop_o3_env=nlev_above_trop_o3_env_from_ukca,                     &
  nlev_ch4_stratloss=nlev_ch4_stratloss_from_ukca,                             &
  l_tracer_lumping=l_tracer_lumping_from_ukca,                                 &
  env_log_step=env_log_step_from_ukca,                                         &
  l_aero_rainout=l_aero_rainout_from_ukca,                                               &
  l_impc_scav=l_impc_scav_from_ukca)

  ! -- Photolysis options
  CALL photol_get_config(                                                      &
    i_photol_scheme=i_scheme_from_photol,                                      &
    fastjx_mode=fastjx_mode_from_photol,                                       &
    fastjx_prescutoff=fastjx_prescutoff_from_photol,                           &
    i_solcylc_type=i_solcyc_from_photol,                                       &
    solcylc_start_year=i_solcyc_start_year_from_photol)

  ! Context information
  IF (PRESENT(row_length)) &
    CALL print_setting(row_length, row_length_from_ukca, 'row_length')
  IF (PRESENT(rows)) &
    CALL print_setting(rows, rows_from_ukca, 'rows')
  IF (PRESENT(model_levels)) &
    CALL print_setting(model_levels, model_levels_from_ukca, 'model_levels')
  IF (PRESENT(bl_levels)) &
    CALL print_setting(bl_levels, bl_levels_from_ukca, 'bl_levels')
  IF (PRESENT(nlev_ent_tr_mix)) &
    CALL print_setting(nlev_ent_tr_mix, nlev_ent_tr_mix_from_ukca, 'nlev_ent_tr_mix')
  IF (PRESENT(ntype)) &
    CALL print_setting(ntype, ntype_from_ukca, 'ntype')
  IF (PRESENT(npft)) &
    CALL print_setting(npft, npft_from_ukca, 'npft')
  IF (PRESENT(i_brd_leaf)) &
    CALL print_setting(i_brd_leaf, i_brd_leaf_from_ukca, 'i_brd_leaf')
  IF (PRESENT(i_brd_leaf_dec)) &
    CALL print_setting(i_brd_leaf_dec, i_brd_leaf_dec_from_ukca, 'i_brd_leaf_dec')
  IF (PRESENT(i_brd_leaf_eg_trop)) &
    CALL print_setting(i_brd_leaf_eg_trop, i_brd_leaf_eg_trop_from_ukca, 'i_brd_leaf_eg_trop')
  IF (PRESENT(i_brd_leaf_eg_temp)) &
    CALL print_setting(i_brd_leaf_eg_temp, i_brd_leaf_eg_temp_from_ukca, 'i_brd_leaf_eg_temp')
  IF (PRESENT(i_ndl_leaf)) &
    CALL print_setting(i_ndl_leaf, i_ndl_leaf_from_ukca, 'i_ndl_leaf')
  IF (PRESENT(i_ndl_leaf_dec)) &
    CALL print_setting(i_ndl_leaf_dec, i_ndl_leaf_dec_from_ukca, 'i_ndl_leaf_dec')
  IF (PRESENT(i_ndl_leaf_eg)) &
    CALL print_setting(i_ndl_leaf_eg, i_ndl_leaf_eg_from_ukca, 'i_ndl_leaf_eg')
  IF (PRESENT(i_c3_grass)) &
    CALL print_setting(i_c3_grass, i_c3_grass_from_ukca, 'i_c3_grass')
  IF (PRESENT(i_c3_crop)) &
    CALL print_setting(i_c3_crop, i_c3_crop_from_ukca, 'i_c3_crop')
  IF (PRESENT(i_c3_pasture)) &
    CALL print_setting(i_c3_pasture, i_c3_pasture_from_ukca, 'i_c3_pasture')
  IF (PRESENT(i_c4_grass)) &
    CALL print_setting(i_c4_grass, i_c4_grass_from_ukca, 'i_c4_grass')
  IF (PRESENT(i_c4_crop)) &
    CALL print_setting(i_c4_crop, i_c4_crop_from_ukca, 'i_c4_crop')
  IF (PRESENT(i_c4_pasture)) &
    CALL print_setting(i_c4_pasture, i_c4_pasture_from_ukca, 'i_c4_pasture')
  IF (PRESENT(i_shrub)) &
    CALL print_setting(i_shrub, i_shrub_from_ukca, 'i_shrub')
  IF (PRESENT(i_shrub_dec)) &
    CALL print_setting(i_shrub_dec, i_shrub_dec_from_ukca, 'i_shrub_dec')
  IF (PRESENT(i_shrub_eg)) &
    CALL print_setting(i_shrub_eg, i_shrub_eg_from_ukca, 'i_shrub_eg')
  IF (PRESENT(i_urban)) &
    CALL print_setting(i_urban, i_urban_from_ukca, 'i_urban')
  IF (PRESENT(i_lake)) &
    CALL print_setting(i_lake, i_lake_from_ukca, 'i_lake')
  IF (PRESENT(i_soil)) &
    CALL print_setting(i_soil, i_soil_from_ukca, 'i_soil')
  IF (PRESENT(i_ice)) &
    CALL print_setting(i_ice, i_ice_from_ukca, 'i_ice')
  IF (PRESENT(i_elev_ice)) &
    CALL print_setting(i_elev_ice, i_elev_ice_from_ukca, 'i_elev_ice')
  IF (PRESENT(dzsoil_layer1)) &
    CALL print_setting(dzsoil_layer1, dzsoil_layer1_from_ukca, 'dzsoil_layer1')
  IF (PRESENT(l_cal360)) &
    CALL print_setting(l_cal360, l_cal360_from_ukca, 'l_cal360')
  IF (PRESENT(timestep)) &
    CALL print_setting(timestep, timestep_from_ukca, 'timestep')
  ! General UKCA configuration options
  IF (PRESENT(i_ukca_chem)) &
    CALL print_setting(i_ukca_chem, i_ukca_chem_from_ukca, 'i_ukca_chem')
  IF (PRESENT(l_ukca_chem_aero)) &
    CALL print_setting(l_ukca_chem_aero, l_ukca_chem_aero_from_ukca, 'l_ukca_chem_aero')
  IF (PRESENT(l_ukca_mode)) &
    CALL print_setting(l_ukca_mode, l_ukca_mode_from_ukca, 'l_ukca_mode')
  IF (PRESENT(l_ukca_ageair)) &
    CALL print_setting(l_ukca_ageair, l_ukca_ageair_from_ukca, 'l_ukca_ageair')
  IF (PRESENT(l_ukca_emissions_off)) &
    CALL print_setting(l_ukca_emissions_off, l_ukca_emissions_off_from_ukca, 'l_ukca_emissions_off')
  IF (PRESENT(l_fix_tropopause_level)) &
    CALL print_setting(l_fix_tropopause_level, l_fix_tropopause_level_from_ukca, 'l_fix_tropopause_level')
  IF (PRESENT(i_ageair_reset_method)) &
    CALL print_setting(i_ageair_reset_method, i_ageair_reset_method_from_ukca, 'i_ageair_reset_method')
  IF (PRESENT(max_ageair_reset_level)) &
    CALL print_setting(max_ageair_reset_level, max_ageair_reset_level_from_ukca, 'max_ageair_reset_level')
  IF (PRESENT(max_ageair_reset_height)) &
    CALL print_setting(max_ageair_reset_height, max_ageair_reset_height_from_ukca, 'max_ageair_reset_height')
  IF (PRESENT(l_enable_diag_um)) &
    CALL print_setting(l_enable_diag_um, l_enable_diag_um_from_ukca, 'l_enable_diag_um')
  IF (PRESENT(l_ukca_persist_off)) &
    CALL print_setting(l_ukca_persist_off, l_ukca_persist_off_from_ukca, 'l_ukca_persist_off')
  IF (PRESENT(l_timer)) &
    CALL print_setting(l_timer, l_timer_from_ukca, 'l_timer')
  ! Chemistry configuration options
  IF (PRESENT(i_ukca_chem_version)) &
    CALL print_setting(i_ukca_chem_version, i_ukca_chem_version_from_ukca, 'i_ukca_chem_version')
  IF (PRESENT(nrsteps)) &
    CALL print_setting(nrsteps, nrsteps_from_ukca, 'nrsteps')
  IF (PRESENT(chem_timestep)) &
    CALL print_setting(chem_timestep, chem_timestep_from_ukca, 'chem_timestep')
  IF (PRESENT(dts0)) &
    CALL print_setting(dts0, dts0_from_ukca, 'dts0')
  IF (PRESENT(l_ukca_asad_columns)) &
    CALL print_setting(l_ukca_asad_columns, l_ukca_asad_columns_from_ukca, 'l_ukca_asad_columns')
  IF (PRESENT(l_ukca_asad_full)) &
    CALL print_setting(l_ukca_asad_full, l_ukca_asad_full_from_ukca, 'l_ukca_asad_full')
  IF (PRESENT(l_ukca_debug_asad)) &
    CALL print_setting(l_ukca_debug_asad, l_ukca_debug_asad_from_ukca, 'l_ukca_debug_asad')
  IF (PRESENT(l_ukca_intdd)) &
    CALL print_setting(l_ukca_intdd, l_ukca_intdd_from_ukca, 'l_ukca_intdd')
  IF (PRESENT(l_ukca_ddepo3_ocean)) &
    CALL print_setting(l_ukca_ddepo3_ocean, l_ukca_ddepo3_ocean_from_ukca, 'l_ukca_ddepo3_ocean')
  IF (PRESENT(l_ukca_ddep_lev1)) &
    CALL print_setting(l_ukca_ddep_lev1, l_ukca_ddep_lev1_from_ukca, 'l_ukca_ddep_lev1')
  IF (PRESENT(l_ukca_dry_dep_so2wet)) &
    CALL print_setting(l_ukca_dry_dep_so2wet, l_ukca_dry_dep_so2wet_from_ukca, 'l_ukca_dry_dep_so2wet')
  IF (PRESENT(l_deposition_jules)) &
    CALL print_setting(l_deposition_jules, l_deposition_jules_from_ukca, 'l_deposition_jules')
  IF (PRESENT(nit)) &
    CALL print_setting(nit, nit_from_ukca, 'nit')
  IF (PRESENT(l_ukca_quasinewton)) &
    CALL print_setting(l_ukca_quasinewton, l_ukca_quasinewton_from_ukca, 'l_ukca_quasinewton')
  IF (PRESENT(i_ukca_quasinewton_start)) &
    CALL print_setting(i_ukca_quasinewton_start, i_ukca_quasinewton_start_from_ukca, 'i_ukca_quasinewton_start')
  IF (PRESENT(i_ukca_quasinewton_end)) &
    CALL print_setting(i_ukca_quasinewton_end, i_ukca_quasinewton_end_from_ukca, 'i_ukca_quasinewton_end')
  IF (PRESENT(ukca_chem_seg_size)) &
    CALL print_setting(ukca_chem_seg_size, ukca_chem_seg_size_from_ukca, 'ukca_chem_seg_size')
  IF (PRESENT(max_z_for_offline_chem)) &
    CALL print_setting(max_z_for_offline_chem, max_z_for_offline_chem_from_ukca, 'max_z_for_offline_chem')
  IF (PRESENT(i_ukca_topboundary)) &
    CALL print_setting(i_ukca_topboundary, i_ukca_topboundary_from_ukca, 'i_ukca_topboundary')
  IF (PRESENT(l_ukca_ro2_ntp)) &
    CALL print_setting(l_ukca_ro2_ntp, l_ukca_ro2_ntp_from_ukca, 'l_ukca_ro2_ntp')
  IF (PRESENT(l_ukca_ro2_perm)) &
    CALL print_setting(l_ukca_ro2_perm, l_ukca_ro2_perm_from_ukca, 'l_ukca_ro2_perm')
  ! Chemistry - Heterogeneous chemistry
  IF (PRESENT(l_ukca_het_psc)) &
    CALL print_setting(l_ukca_het_psc, l_ukca_het_psc_from_ukca, 'l_ukca_het_psc')
  IF (PRESENT(i_ukca_hetconfig)) &
    CALL print_setting(i_ukca_hetconfig, i_ukca_hetconfig_from_ukca, 'i_ukca_hetconfig')
  IF (PRESENT(l_ukca_limit_nat)) &
    CALL print_setting(l_ukca_limit_nat, l_ukca_limit_nat_from_ukca, 'l_ukca_limit_nat')
  IF (PRESENT(l_ukca_sa_clim)) &
    CALL print_setting(l_ukca_sa_clim, l_ukca_sa_clim_from_ukca, 'l_ukca_sa_clim')
  IF (PRESENT(l_ukca_trophet)) &
    CALL print_setting(l_ukca_trophet, l_ukca_trophet_from_ukca, 'l_ukca_trophet')
  IF (PRESENT(l_ukca_classic_hetchem)) &
    CALL print_setting(l_ukca_classic_hetchem, l_ukca_classic_hetchem_from_ukca, 'l_ukca_classic_hetchem')
  ! Chemistry - Photolysis
  IF (PRESENT(i_ukca_photol)) &
    CALL print_setting(i_ukca_photol, i_scheme_from_photol, 'i_photol_scheme')
  IF (PRESENT(fastjx_mode)) &
    CALL print_setting(fastjx_mode, fastjx_mode_from_photol, 'fastjx_mode')
  IF (PRESENT(i_ukca_solcyc)) &
    CALL print_setting(i_ukca_solcyc, i_solcyc_from_photol, 'i_solcyc_opt')
  IF (PRESENT(i_ukca_solcyc_start_year)) &
    CALL print_setting(i_ukca_solcyc_start_year, i_solcyc_start_year_from_photol, 'i_solcyc_start_year')
  ! UKCA emissions configuration options
  IF (PRESENT(l_ukca_ibvoc)) &
    CALL print_setting(l_ukca_ibvoc, l_ukca_ibvoc_from_ukca, 'l_ukca_ibvoc')
  IF (PRESENT(l_ukca_inferno)) &
    CALL print_setting(l_ukca_inferno, l_ukca_inferno_from_ukca, 'l_ukca_inferno')
  IF (PRESENT(l_ukca_inferno_ch4)) &
    CALL print_setting(l_ukca_inferno_ch4, l_ukca_inferno_ch4_from_ukca, 'l_ukca_inferno_ch4')
  IF (PRESENT(i_inferno_emi)) &
    CALL print_setting(i_inferno_emi, i_inferno_emi_from_ukca, 'i_inferno_emi')
  IF (PRESENT(l_ukca_so2ems_expvolc)) &
    CALL print_setting(l_ukca_so2ems_expvolc, l_ukca_so2ems_expvolc_from_ukca, 'l_ukca_so2ems_expvolc')
  IF (PRESENT(l_ukca_so2ems_plumeria)) &
    CALL print_setting(l_ukca_so2ems_plumeria, l_ukca_so2ems_plumeria_from_ukca, 'l_ukca_so2ems_plumeria')
  IF (PRESENT(l_ukca_qch4inter)) &
    CALL print_setting(l_ukca_qch4inter, l_ukca_qch4inter_from_ukca, 'l_ukca_qch4inter')
  IF (PRESENT(l_ukca_emsdrvn_ch4)) &
    CALL print_setting(l_ukca_emsdrvn_ch4, l_ukca_emsdrvn_ch4_from_ukca, 'l_ukca_emsdrvn_ch4')
  IF (PRESENT(mode_parfrac)) &
    CALL print_setting(mode_parfrac, mode_parfrac_from_ukca, 'mode_parfrac')
  IF (PRESENT(l_ukca_enable_seadms_ems)) &
    CALL print_setting(l_ukca_enable_seadms_ems, l_ukca_enable_seadms_ems_from_ukca, 'l_ukca_enable_seadms_ems')
  IF (PRESENT(i_ukca_dms_flux)) &
    CALL print_setting(i_ukca_dms_flux, i_ukca_dms_flux_from_ukca, 'i_ukca_dms_flux')
  IF (PRESENT(seadms_ems_scaling)) &
    CALL print_setting(seadms_ems_scaling, seadms_ems_scaling_from_ukca, 'seadms_ems_scaling')
  IF (PRESENT(l_ukca_linox_scaling)) &
    CALL print_setting(l_ukca_linox_scaling, l_ukca_linox_scaling_from_ukca, 'l_ukca_linox_scaling')
  IF (PRESENT(lightnox_scale_fac)) &
    CALL print_setting(lightnox_scale_fac, lightnox_scale_fac_from_ukca, 'lightnox_scale_fac')
  IF (PRESENT(i_ukca_light_param)) &
    CALL print_setting(i_ukca_light_param, i_ukca_light_param_from_ukca, 'i_ukca_light_param')
  IF (PRESENT(l_ukca_scale_soa_yield_mt)) &
    CALL print_setting(l_ukca_scale_soa_yield_mt, l_ukca_scale_soa_yield_mt_from_ukca, 'l_ukca_scale_soa_yield_mt')
  IF (PRESENT(soa_yield_scaling_mt)) &
    CALL print_setting(soa_yield_scaling_mt, soa_yield_scaling_mt_from_ukca, 'soa_yield_scaling_mt')
  IF (PRESENT(l_ukca_scale_soa_yield_isop)) &
    CALL print_setting(l_ukca_scale_soa_yield_isop, l_ukca_scale_soa_yield_isop_from_ukca, 'l_ukca_scale_soa_yield_isop')
  IF (PRESENT(soa_yield_scaling_isop)) &
    CALL print_setting(soa_yield_scaling_isop, soa_yield_scaling_isop_from_ukca, 'soa_yield_scaling_isop')
  IF (PRESENT(l_support_ems_vertprof)) &
    CALL print_setting(l_support_ems_vertprof, l_support_ems_vertprof_from_ukca, 'l_support_ems_vertprof')
  IF (PRESENT(l_support_ems_gridbox_units)) &
    CALL print_setting(l_support_ems_gridbox_units, l_support_ems_gridbox_units_from_ukca, 'l_support_ems_gridbox_units')
  ! UKCA feedback configuration options
  IF (PRESENT(l_ukca_h2o_feedback)) &
    CALL print_setting(l_ukca_h2o_feedback, l_ukca_h2o_feedback_from_ukca, 'l_ukca_h2o_feedback')
  IF (PRESENT(l_ukca_conserve_h)) &
    CALL print_setting(l_ukca_conserve_h, l_ukca_conserve_h_from_ukca, 'l_ukca_conserve_h')
  ! UKCA environmental driver configuration options
  IF (PRESENT(l_param_conv)) &
    CALL print_setting(l_param_conv, l_param_conv_from_ukca, 'l_param_conv')
  IF (PRESENT(l_ctile)) &
    CALL print_setting(l_ctile, l_ctile_from_ukca, 'l_ctile')
  IF (PRESENT(l_zon_av_ozone)) &
    CALL print_setting(l_zon_av_ozone, l_zon_av_ozone_from_ukca, 'l_zon_av_ozone')
  IF (PRESENT(i_strat_lbc_source)) &
    CALL print_setting(i_strat_lbc_source, i_strat_lbc_source_from_ukca, 'i_strat_lbc_source')
  IF (PRESENT(l_chem_environ_gas_scalars)) &
    CALL print_setting(l_chem_environ_gas_scalars, l_chem_environ_gas_scalars_from_ukca, 'l_chem_environ_gas_scalars')
  IF (PRESENT(l_chem_environ_co2_fld)) &
    CALL print_setting(l_chem_environ_co2_fld, l_chem_environ_co2_fld_from_ukca, 'l_chem_environ_co2_fld')
  IF (PRESENT(l_ukca_prescribech4)) &
    CALL print_setting(l_ukca_prescribech4, l_ukca_prescribech4_from_ukca, 'l_ukca_prescribech4')
  IF (PRESENT(l_use_classic_so4)) &
    CALL print_setting(l_use_classic_so4, l_use_classic_so4_from_ukca, 'l_use_classic_so4')
  IF (PRESENT(l_use_classic_soot)) &
    CALL print_setting(l_use_classic_soot, l_use_classic_soot_from_ukca, 'l_use_classic_soot')
  IF (PRESENT(l_use_classic_ocff)) &
    CALL print_setting(l_use_classic_ocff, l_use_classic_ocff_from_ukca, 'l_use_classic_ocff')
  IF (PRESENT(l_use_classic_biogenic)) &
    CALL print_setting(l_use_classic_biogenic, l_use_classic_biogenic_from_ukca, 'l_use_classic_biogenic')
  IF (PRESENT(l_use_classic_seasalt)) &
    CALL print_setting(l_use_classic_seasalt, l_use_classic_seasalt_from_ukca, 'l_use_classic_seasalt')
  IF (PRESENT(l_use_gridbox_volume)) &
    CALL print_setting(l_use_gridbox_volume, l_use_gridbox_volume_from_ukca, 'l_use_gridbox_volume')
  IF (PRESENT(l_use_gridbox_mass)) &
    CALL print_setting(l_use_gridbox_mass, l_use_gridbox_mass_from_ukca, 'l_use_gridbox_mass')
  IF (PRESENT(l_environ_z_top)) &
    CALL print_setting(l_environ_z_top, l_environ_z_top_from_ukca, 'l_environ_z_top')
  ! UKCA temporary logicals
  IF (PRESENT(l_fix_improve_drydep)) &
    CALL print_setting(l_fix_improve_drydep, l_fix_improve_drydep_from_ukca, 'l_fix_improve_drydep')
  IF (PRESENT(l_fix_ukca_h2dd_x)) &
    CALL print_setting(l_fix_ukca_h2dd_x, l_fix_ukca_h2dd_x_from_ukca, 'l_fix_ukca_h2dd_x')
  IF (PRESENT(l_fix_drydep_so2_water)) &
    CALL print_setting(l_fix_drydep_so2_water, l_fix_drydep_so2_water_from_ukca, 'l_fix_drydep_so2_water')
  IF (PRESENT(l_fix_ukca_offox_h2o_fac)) &
    CALL print_setting(l_fix_ukca_offox_h2o_fac, l_fix_ukca_offox_h2o_fac_from_ukca, 'l_fix_ukca_offox_h2o_fac')
  IF (PRESENT(l_fix_ukca_h2so4_ystore)) &
    CALL print_setting(l_fix_ukca_h2so4_ystore, l_fix_ukca_h2so4_ystore_from_ukca, 'l_fix_ukca_h2so4_ystore')
  IF (PRESENT(l_improve_aero_drydep)) &
    CALL print_setting(l_improve_aero_drydep, l_improve_aero_drydep_from_ukca, 'l_improve_aero_drydep')
  ! General GLOMAP configuration options
  IF (PRESENT(i_mode_nzts)) &
    CALL print_setting(i_mode_nzts, i_mode_nzts_from_ukca, 'i_mode_nzts')
  IF (PRESENT(ukca_mode_seg_size)) &
    CALL print_setting(ukca_mode_seg_size, ukca_mode_seg_size_from_ukca, 'ukca_mode_seg_size')
  IF (PRESENT(i_mode_setup)) &
    CALL print_setting(i_mode_setup, i_mode_setup_from_ukca, 'i_mode_setup')
  IF (PRESENT(l_mode_bhn_on)) &
    CALL print_setting(l_mode_bhn_on, l_mode_bhn_on_from_ukca, 'l_mode_bhn_on')
  IF (PRESENT(l_mode_bln_on)) &
    CALL print_setting(l_mode_bln_on, l_mode_bln_on_from_ukca, 'l_mode_bln_on')
  IF (PRESENT(i_mode_bln_param_method)) &
    CALL print_setting(i_mode_bln_param_method, i_mode_bln_param_method_from_ukca, 'i_mode_bln_param_method')
  IF (PRESENT(i_mode_nucscav)) &
    CALL print_setting(i_mode_nucscav, i_mode_nucscav_from_ukca, 'i_mode_nucscav')
  IF (PRESENT(mode_activation_dryr)) &
    CALL print_setting(mode_activation_dryr, mode_activation_dryr_from_ukca, 'mode_activation_dryr')
  IF (PRESENT(mode_incld_so2_rfrac)) &
    CALL print_setting(mode_incld_so2_rfrac, mode_incld_so2_rfrac_from_ukca, 'mode_incld_so2_rfrac')
  IF (PRESENT(l_cv_rainout)) &
    CALL print_setting(l_cv_rainout, l_cv_rainout_from_ukca, 'l_cv_rainout')
  IF (PRESENT(l_ddepaer)) &
    CALL print_setting(l_ddepaer, l_ddepaer_from_ukca, 'l_ddepaer')
  IF (PRESENT(l_dust_mp_slinn_impc_scav)) &
    CALL print_setting(l_dust_mp_slinn_impc_scav, l_dust_mp_slinn_impc_scav_from_ukca, 'l_dust_mp_slinn_impc_scav')
  ! GLOMAP emissions configuration options
  IF (PRESENT(l_ukca_primss)) &
    CALL print_setting(l_ukca_primss, l_ukca_primss_from_ukca, 'l_ukca_primss')
  IF (PRESENT(l_ukca_primsu)) &
    CALL print_setting(l_ukca_primsu, l_ukca_primsu_from_ukca, 'l_ukca_primsu')
  IF (PRESENT(l_ukca_primdu)) &
    CALL print_setting(l_ukca_primdu, l_ukca_primdu_from_ukca, 'l_ukca_primdu')
  IF (PRESENT(l_ukca_primbcoc)) &
    CALL print_setting(l_ukca_primbcoc, l_ukca_primbcoc_from_ukca, 'l_ukca_primbcoc')
  IF (PRESENT(l_ukca_prim_moc)) &
    CALL print_setting(l_ukca_prim_moc, l_ukca_prim_moc_from_ukca, 'l_ukca_prim_moc')
  IF (PRESENT(l_bcoc_bf)) &
    CALL print_setting(l_bcoc_bf, l_bcoc_bf_from_ukca, 'l_bcoc_bf')
  IF (PRESENT(l_bcoc_bm)) &
    CALL print_setting(l_bcoc_bm, l_bcoc_bm_from_ukca, 'l_bcoc_bm')
  IF (PRESENT(l_bcoc_ff)) &
    CALL print_setting(l_bcoc_ff, l_bcoc_ff_from_ukca, 'l_bcoc_ff')
  IF (PRESENT(l_ukca_scale_biom_aer_ems)) &
    CALL print_setting(l_ukca_scale_biom_aer_ems, l_ukca_scale_biom_aer_ems_from_ukca, 'l_ukca_scale_biom_aer_ems')
  IF (PRESENT(biom_aer_ems_scaling)) &
    CALL print_setting(biom_aer_ems_scaling, biom_aer_ems_scaling_from_ukca, 'biom_aer_ems_scaling')
  IF (PRESENT(l_ukca_fine_no3_prod)) &
    CALL print_setting(l_ukca_fine_no3_prod, l_ukca_fine_no3_prod_from_ukca, 'l_ukca_fine_no3_prod')
  IF (PRESENT(l_ukca_coarse_no3_prod)) &
    CALL print_setting(l_ukca_coarse_no3_prod, l_ukca_coarse_no3_prod_from_ukca, 'l_ukca_coarse_no3_prod')
  IF (PRESENT(l_no3_prod_in_aero_step)) &
    CALL print_setting(l_no3_prod_in_aero_step, l_no3_prod_in_aero_step_from_ukca, 'l_no3_prod_in_aero_step')
  IF (PRESENT(l_ukca_scale_sea_salt_ems)) &
    CALL print_setting(l_ukca_scale_sea_salt_ems, l_ukca_scale_sea_salt_ems_from_ukca, 'l_ukca_scale_sea_salt_ems')
  IF (PRESENT(sea_salt_ems_scaling)) &
    CALL print_setting(sea_salt_ems_scaling, sea_salt_ems_scaling_from_ukca, 'sea_salt_ems_scaling')
  IF (PRESENT(l_ukca_scale_marine_pom_ems)) &
    CALL print_setting(l_ukca_scale_marine_pom_ems, l_ukca_scale_marine_pom_ems_from_ukca, 'l_ukca_scale_marine_pom_ems')
  IF (PRESENT(marine_pom_ems_scaling)) &
    CALL print_setting(marine_pom_ems_scaling, marine_pom_ems_scaling_from_ukca, 'marine_pom_ems_scaling')
  ! GLOMAP feedback configuration options
  IF (PRESENT(l_ukca_radaer)) &
    CALL print_setting(l_ukca_radaer, l_ukca_radaer_from_ukca, 'l_ukca_radaer')
  IF (PRESENT(i_ukca_tune_bc)) &
    CALL print_setting(i_ukca_tune_bc, i_ukca_tune_bc_from_ukca, 'i_ukca_tune_bc')
  IF (PRESENT(i_ukca_activation_scheme)) &
    CALL print_setting(i_ukca_activation_scheme, i_ukca_activation_scheme_from_ukca, 'i_ukca_activation_scheme')
  IF (PRESENT(i_ukca_nwbins)) &
    CALL print_setting(i_ukca_nwbins, i_ukca_nwbins_from_ukca, 'i_ukca_nwbins')
  IF (PRESENT(sigwmin)) &
    CALL print_setting(sigwmin, sigwmin_from_ukca, 'sigwmin')
  IF (PRESENT(l_ntpreq_n_activ_sum)) &
    CALL print_setting(l_ntpreq_n_activ_sum, l_ntpreq_n_activ_sum_from_ukca, 'l_ntpreq_n_activ_sum')
  IF (PRESENT(l_ntpreq_dryd_nuc_sol)) &
    CALL print_setting(l_ntpreq_dryd_nuc_sol, l_ntpreq_dryd_nuc_sol_from_ukca, 'l_ntpreq_dryd_nuc_sol')
  IF (PRESENT(l_ukca_sfix)) &
    CALL print_setting(l_ukca_sfix, l_ukca_sfix_from_ukca, 'l_ukca_sfix')
  ! GLOMAP temporary logicals
  IF (PRESENT(l_fix_neg_pvol_wat)) &
    CALL print_setting(l_fix_neg_pvol_wat, l_fix_neg_pvol_wat_from_ukca, 'l_fix_neg_pvol_wat')
  IF (PRESENT(l_fix_ukca_impscav)) &
    CALL print_setting(l_fix_ukca_impscav, l_fix_ukca_impscav_from_ukca, 'l_fix_ukca_impscav')
  IF (PRESENT(l_fix_nacl_density)) &
    CALL print_setting(l_fix_nacl_density, l_fix_nacl_density_from_ukca, 'l_fix_nacl_density')
  IF (PRESENT(l_fix_ukca_activate_pdf)) &
    CALL print_setting(l_fix_ukca_activate_pdf, l_fix_ukca_activate_pdf_from_ukca, 'l_fix_ukca_activate_pdf')
  IF (PRESENT(l_fix_ukca_activate_vert_rep)) &
    CALL print_setting(l_fix_ukca_activate_vert_rep, l_fix_ukca_activate_vert_rep_from_ukca, 'l_fix_ukca_activate_vert_rep')
  IF (PRESENT(l_bug_repro_tke_index)) &
    CALL print_setting(l_bug_repro_tke_index, l_bug_repro_tke_index_from_ukca, 'l_bug_repro_tke_index')
  IF (PRESENT(l_fix_ukca_hygroscopicities)) &
    CALL print_setting(l_fix_ukca_hygroscopicities, l_fix_ukca_hygroscopicities_from_ukca, 'l_fix_ukca_hygroscopicities')
  IF (PRESENT(l_fix_ukca_water_content)) &
    CALL print_setting(l_fix_ukca_water_content, l_fix_ukca_water_content_from_ukca, 'l_fix_ukca_water_content')
  ! New settings in Ticket #25
  IF (PRESENT(l_ukca_drydep_off)) &
    CALL print_setting(l_ukca_drydep_off, l_ukca_drydep_off_from_ukca, 'l_ukca_drydep_off')
  IF (PRESENT(l_ukca_wetdep_off)) &
    CALL print_setting(l_ukca_wetdep_off, l_ukca_wetdep_off_from_ukca, 'l_ukca_wetdep_off')
  IF (PRESENT(nlev_above_trop_o3_env)) &
    CALL print_setting(nlev_above_trop_o3_env, nlev_above_trop_o3_env_from_ukca, 'nlev_above_trop_o3_env')
  IF (PRESENT(nlev_ch4_stratloss)) &
    CALL print_setting(nlev_ch4_stratloss, nlev_ch4_stratloss_from_ukca, 'nlev_ch4_stratloss')
  IF (PRESENT(l_tracer_lumping)) &
    CALL print_setting(l_tracer_lumping, l_tracer_lumping_from_ukca, 'l_tracer_lumping')
  IF (PRESENT(env_log_step)) &
    CALL print_setting(env_log_step, env_log_step_from_ukca, 'env_log_step')
  IF (PRESENT(l_aero_rainout)) &
    CALL print_setting(l_aero_rainout, l_aero_rainout_from_ukca, 'l_aero_rainout')
  IF (PRESENT(l_impc_scav)) &
    CALL print_setting(l_impc_scav, l_impc_scav_from_ukca, 'l_impc_scav')

RETURN

END SUBROUTINE print_config_settings

! ----------------------------------------------------------------------
SUBROUTINE print_setting_integer(parent_value, ukca_value, ukca_varname)
! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: parent_value
INTEGER, INTENT(IN) :: ukca_value
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname

CHARACTER(LEN=30) :: varname_txt
CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_SETTING_INTEGER'

varname_txt = ukca_varname
WRITE(umMessage,'(A30,2(1X,I11),1X,L1)') &
      varname_txt, parent_value, ukca_value, ukca_value == parent_value
CALL umPrint(umMessage, src=RoutineName)

RETURN

END SUBROUTINE print_setting_integer 

! ----------------------------------------------------------------------
SUBROUTINE print_setting_integer_vec(parent_value, ukca_value, ukca_varname)
! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER, ALLOCATABLE, INTENT(IN) :: ukca_value(:)
INTEGER, ALLOCATABLE, INTENT(IN) :: parent_value(:)
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname

INTEGER :: i
LOGICAL :: match
CHARACTER(LEN=30) :: varname_txt
CHARACTER(LEN=11) :: parent_txt
CHARACTER(LEN=11) :: ukca_txt
CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_SETTING_INTEGER_VEC'

varname_txt = ukca_varname
IF (ALLOCATED(parent_value)) THEN
  parent_txt = 'ALLOCATED'
ELSE
  parent_txt = 'UNALLOCATED'
END IF
IF (ALLOCATED(ukca_value)) THEN
  ukca_txt = 'ALLOCATED'
ELSE
  ukca_txt = 'UNALLOCATED'
END IF
match = (ukca_txt == parent_txt)
IF (ALLOCATED(ukca_value)) THEN
  DO i = 1, SIZE(ukca_value)
    IF (ukca_value(i) /= parent_value(i)) match = .FALSE.
  END DO
END IF 
WRITE(umMessage,'(A30,2(1X,A11),1X,L1)') &
      varname_txt, parent_txt, ukca_txt, match 
CALL umPrint(umMessage, src=RoutineName)

RETURN

END SUBROUTINE print_setting_integer_vec

! ----------------------------------------------------------------------
SUBROUTINE print_setting_real(parent_value, ukca_value, ukca_varname)
! ----------------------------------------------------------------------

IMPLICIT NONE

REAL, INTENT(IN) :: parent_value
REAL, INTENT(IN) :: ukca_value
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname

CHARACTER(LEN=30) :: varname_txt
CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_SETTING_REAL'

varname_txt = ukca_varname
WRITE(umMessage,'(A30,2(1X,ES11.4),1X,L1)') &
      varname_txt, parent_value, ukca_value, ukca_value == parent_value
CALL umPrint(umMessage, src=RoutineName)

RETURN

END SUBROUTINE print_setting_real

! ----------------------------------------------------------------------
SUBROUTINE print_setting_logical(parent_value, ukca_value, ukca_varname)
! ----------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT(IN) :: parent_value
LOGICAL, INTENT(IN) :: ukca_value
CHARACTER(LEN=*), INTENT(IN) :: ukca_varname

CHARACTER(LEN=30) :: varname_txt
CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_SETTING_LOGICAL'

varname_txt = ukca_varname
WRITE(umMessage,'(A30,2(1X,L11),1X,L1)') &
      varname_txt, parent_value, ukca_value, ukca_value .EQV. parent_value
CALL umPrint(umMessage, src=RoutineName)

RETURN

END SUBROUTINE print_setting_logical
 
END MODULE box_ukca_setup_mod

