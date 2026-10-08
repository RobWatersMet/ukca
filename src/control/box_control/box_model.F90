! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine: BOX_MODEL                     -------------------------
!
!    Purpose: High level control program for the UKCA box model
!             (master routine).  Calls lower level control routines
!             according to top level switch settings. Called by
!             top level routine UKCA_SHELL which provides dimension sizes
!             for dynamic allocation of data arrays.
!             Based on u_model_4a and atm_step_4A
!
!    -------------------------------------------------------------------
MODULE box_model_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BOX_MODEL_MOD'

!++SAN tracer_ukca now allocated and populated in this routine
REAL, ALLOCATABLE, PUBLIC :: tracer_ukca(:,:,:,:)

!++SAN ntp_data now made in this routine and passed to box_ukca
!      Code is based on ukca_um_interf_mod
REAL, ALLOCATABLE, PUBLIC :: ntp_data(:,:,:,:)

!++SAN Defining here some variables that are defined elsewhere in the 
!      UM, but put here for ease of use for the box model
!++ Environmental variables
! 1D real
REAL, ALLOCATABLE    :: soil_moisture_layer1 (:)

! 2D Logical
LOGICAL, ALLOCATABLE :: land_sea_mask (:,:)
LOGICAL, ALLOCATABLE :: l_tile_active (:,:)

! 2D Real
REAL, ALLOCATABLE    :: conv_cloud_lwp (:,:)
REAL, ALLOCATABLE    :: tstar(:,:)
REAL, ALLOCATABLE    :: zbl(:,:)
REAL, ALLOCATABLE    :: rough_length(:,:)
REAL, ALLOCATABLE    :: seaice_frac(:,:)
REAL, ALLOCATABLE    :: pstar(:,:)
REAL, ALLOCATABLE    :: zhsc(:,:)
REAL, ALLOCATABLE    :: u_scalar_10m(:,:)
REAL, ALLOCATABLE    :: u_s(:,:)
REAL, ALLOCATABLE    :: surf_albedo(:,:)
REAL, ALLOCATABLE    :: frac_types(:,:)
REAL, ALLOCATABLE    :: laift_lp(:,:)
REAL, ALLOCATABLE    :: canhtft_lp(:,:)
REAL, ALLOCATABLE    :: tstar_tile(:,:)
REAL, ALLOCATABLE    :: z0tile_lp(:,:)
REAL, ALLOCATABLE    :: surf_hf(:,:)

REAL, ALLOCATABLE    :: dms_sea_conc(:,:)

! 2D Integer
INTEGER, ALLOCATABLE :: kent(:,:)
INTEGER, ALLOCATABLE :: kent_dsc(:,:)
INTEGER, ALLOCATABLE :: conv_cloud_base(:,:)
INTEGER, ALLOCATABLE :: conv_cloud_top(:,:)

! 3D real
REAL, ALLOCATABLE    :: theta(:,:,:)
REAL, ALLOCATABLE    :: q(:,:,:)
REAL, ALLOCATABLE    :: qcf(:,:,:)
REAL, ALLOCATABLE    :: conv_cloud_amount(:,:,:)
REAL, ALLOCATABLE    :: rho_r2(:,:,:)
REAL, ALLOCATABLE    :: qcl(:,:,:)
REAL, ALLOCATABLE    :: exner_rho_levels(:,:,:)
REAL, ALLOCATABLE    :: area_cloud_fraction(:,:,:)
REAL, ALLOCATABLE    :: cloud_frac(:,:,:)
REAL, ALLOCATABLE    :: cloud_liq_frac(:,:,:)
REAL, ALLOCATABLE    :: exner_theta_levels(:,:,:)
REAL, ALLOCATABLE    :: p_rho_levels(:,:,:)
REAL, ALLOCATABLE    :: p_theta_levels(:,:,:)
REAL, ALLOCATABLE    :: t_theta_levels(:,:,:)
REAL, ALLOCATABLE    :: rhokh_mix(:,:,:)
REAL, ALLOCATABLE    :: dtrdz_charney_grid(:,:,:)
REAL, ALLOCATABLE    :: rhokh_rdz(:,:,:)
REAL, ALLOCATABLE    :: dtrdz(:,:,:)
REAL, ALLOCATABLE    :: we_lim(:,:,:)
REAL, ALLOCATABLE    :: t_frac(:,:,:)
REAL, ALLOCATABLE    :: zrzi(:,:,:)
REAL, ALLOCATABLE    :: we_lim_dsc(:,:,:)
REAL, ALLOCATABLE    :: t_frac_dsc(:,:,:)
REAL, ALLOCATABLE    :: zrzi_dsc(:,:,:)
REAL, ALLOCATABLE    :: ls_rain3d(:,:,:)
REAL, ALLOCATABLE    :: ls_snow3d(:,:,:)
REAL, ALLOCATABLE    :: autoconv(:,:,:)
REAL, ALLOCATABLE    :: accretion(:,:,:)
REAL, ALLOCATABLE    :: pv_on_theta_mlevs(:,:,:)
REAL, ALLOCATABLE    :: stcon(:,:,:)

REAL, ALLOCATABLE    :: conv_rain3d(:,:,:)
REAL, ALLOCATABLE    :: conv_snow3d(:,:,:)

REAL, ALLOCATABLE    :: so4_sa_clim(:,:,:)

REAL, ALLOCATABLE    :: rim_cry(:,:,:)
REAL, ALLOCATABLE    :: rim_agg(:,:,:)
REAL, ALLOCATABLE    :: vertvel(:,:,:)
REAL, ALLOCATABLE    :: bl_tke(:,:,:)


!---------------------------------------------------------------------------
! Define namelist(s)
!---------------------------------------------------------------------------

!++SAN - code that could be added if include this namelist in future
!s NAMELIST/INT_DRY_DEP_PARS/                                 &
!s soil_moisture_layer1_in, l_tile_active_in, frac_types_in,  &
!s laift_lp_in, canhtft_lp_in, tstar_tile_in, z0tile_lp_in,   &
!s surf_hf_in, stcon_in

! PRIVATE :: INT_DRY_DEP_PARS

!s NAMELIST/AEROSOL_PARS/                                 &
! rim_cry_in, rim_agg_in, vertvel_in, bl_tke_in
! PRIVATE :: AEROSOL_PARS

!=============================================================================


CONTAINS

!++SAN - debate whether or not to use a dump field in similar way for input
!++      to run box model. Run subroutine without arguments for now
!s SUBROUTINE box_model(dump_unit)
SUBROUTINE box_model()

!++SAN - deleted unnecesary dependencies
!++SAN - replacing atm_step_4A with box model version of caller to UKCA
!s USE atm_step_4A_mod,        ONLY: atm_step_4A
USE box_ukca_mod, ONLY: box_ukca

!++SAN all requests to UKCA side routines must go through API
USE ukca_api_mod,           ONLY: ukca_get_tracer_varlist,                    &
                                  ukca_get_ntp_varlist,                       &
                                  ukca_maxlen_fieldname,                      &
                                  ukca_get_environment_varlist,               &
                                  ukca_set_environment,                       &
                                  ukca_maxlen_message,                        &
                                  ukca_maxlen_procname,                       &
                                  ukca_get_config

!++SAN - added d1_defs for now, needed for some codes and lengths of arrays
USE ukca_d1_defs, ONLY: n_ntp

!++SAN added code for setting concentrations of long-lived trace gases
USE ukca_trace_gas_mixratio, ONLY: ukca_set_trace_gas_mixratio

!++SAN Add call to allocate empty stash arrays
USE stash_array_mod, ONLY: allocate_stash_arrays, sf, nsects, nitems

!++SAN need to lump halogens in first call
USE ukca_transform_halogen_mod,                                                &
                            ONLY: ukca_transform_halogen

!++SAN Might need to find a way to move this to the API:
USE ukca_config_defs_mod, ONLY: n_use_tracers, n_chem_tracers, n_aero_tracers, &
      n_mode_tracers, n_nonchem_tracers

USE ukca_option_mod, ONLY:                                                     &
     l_ukca, l_ukca_aie1, l_ukca_radaer, l_ukca_radaer_sustrat,                &
     i_ukca_scenario, i_ukca_scenario_um, l_ukca_prescribech4,                 &
     l_ukca_set_trace_gases

!++SAN constants for calculating exner
USE planet_constants_mod,   ONLY: repsilon, g, two_omega,                     &
            kappa, p_zero, sclht, planet_radius, r, c_virtual

USE missing_data_mod,       ONLY: rmdi, imdi

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE exitchek_mod,           ONLY: exitchek

!s USE gas_calc_all_mod,       ONLY: gas_calc_all
USE get_env_var_mod,        ONLY: get_env_var
USE get_wallclock_time_mod, ONLY: get_wallclock_time
USE incrtime_mod,           ONLY: incrtime
!++SAN initial_4A needs replacing with initial_box
!++ UPDATE - initial_box was so simple, decide to remove and move the calls 
!            up a level to box_model
!s USE initial_box_mod,         ONLY: initial_box
USE io_configuration_mod,   ONLY: print_runtime_info
USE io,                     ONLY: io_timestep

USE model_time_mod,         ONLY: secs_per_stepim, stepim, target_end_stepim
!++SAN Adding more global variables that need defining
USE nlsizes_namelist_mod,   ONLY: row_length, rows, global_row_length, model_levels, &
                                  tr_ukca 
USE nlstcall_mod,           ONLY: lpp, lnc, ldump, lmean, lprint, lancillary,  &
                                  lexit, lboundary, ltimer, lstashdumptimer

!++SAN Will replace this file with a box version with just the variables needed
!      to call the box model
!++UPDATE - will just declare inside this routine, as is done for many other variables
!s USE atm_fields_mod,         ONLY: tracer_ukca ! , q

!++SAN Added namelist options for Box model I/O
USE nlstcall_box_namelist_mod, ONLY: ukca_box_nml, tracer_in_filename,         &
                                     tracer_out_filename, tracer_nullval,      &
                                     flux_out_filename, rate_out_filename,     &
                                     ntp_out_filename

USE read_environment_pars_mod, ONLY: read_environment_pars

USE run_info,               ONLY: start_time
USE settsctl_mod,           ONLY: settsctl
USE submodel_mod,           ONLY: submodel_partition_list, atmos_im

USE timer_mod,              ONLY: timer
!++SAN moved inittime up one level
USE inittime_mod,           ONLY: inittime

!++SAN moved trigonometir setting up
USE trignometric_mod,       ONLY: FV_cos_theta_latitude,                       &
                              true_latitude, true_longitude
USE cderived_mod,           ONLY: delta_lambda, delta_phi

!++SAN adding more global variables that need defining
USE level_heights_mod,      ONLY: r_rho_levels, r_theta_levels,                &
                                eta_theta_levels

!++SAN adding helpful conversions
USE conversions_mod,        ONLY: pi_over_180

!s USE turb_diff_mod,          ONLY: l_print_pe
USE um_parcore,             ONLY: mype
USE umPrintMgr,             ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE parkind1,               ONLY: jprb, jpim
USE yomhook,                ONLY: lhook, dr_hook

USE rad_input_mod, ONLY: c11mmr, c12mmr, c113mmr, c114mmr,ch4mmr,co2_mmr,     &
     n2ommr, o2mmr,hcfc22mmr, hfc125mmr, hfc134ammr

USE box_output_chem_diags_mod, ONLY: box_output_chem_diags

IMPLICIT NONE

REAL :: time_end_run

!    Interface and arguments: ------------------------------------------
!++SAN - debate whether to use this same framework for having an input
!        file to initial variables. Comment out for now
!sINTEGER, INTENT(IN) :: dump_unit  ! unit attached to input data file
!s                                  ! (astart or checkpoint_dump_im)

! ----------------------------------------------------------------------
!
!
!  Local variables
!
INTEGER :: internal_model    ! Work - Internal model identifier
INTEGER :: submodel          ! Work - Submodel id for dump partition
INTEGER :: meanlev           ! Work - Mean level indicator
INTEGER :: iabort            ! Work - Internal return code
INTEGER :: g_theta_field_size
                             ! Sizes for MPP dynamic allocation
                          ! in A-O coupling routines
INTEGER :: land_points       ! 1 for land
!
LOGICAL :: lexitNOW          ! Work - Immediate exit indicator

! 3-D fields of species to be passed down to radiation
INTEGER, PARAMETER :: ngrgas = 8
INTEGER, SAVE :: grgas_addr(ngrgas)

! Error reporting
INTEGER :: icode       ! =0 normal exit; >0 error exit
INTEGER :: istatus
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*) :: RoutineName

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
PARAMETER (   RoutineName='BOX_MODEL')

!++SAN Local variables for holding chemistry information read from chem_tracer_pars
REAL, ALLOCATABLE                 :: inp_tracer_mmrs(:)
CHARACTER(LEN=10), ALLOCATABLE    :: inp_tracer_names(:)

! Copies of internal UKCA config variables
LOGICAL :: l_ukca_intdd
INTEGER :: ntype
INTEGER :: npft
LOGICAL :: l_ukca_strat
LOGICAL :: l_ukca_strattrop
LOGICAL :: l_ukca_stratcfc
LOGICAL :: l_ukca_cristrat

!---------------------------------------------------------------------------
!++SAN scalar versions of variables for reading in namelist input
REAL    :: latitude_in             = rmdi
REAL    :: longitude_in            = rmdi
LOGICAL :: land_sea_mask_in        = .false.
REAL    :: conv_cloud_lwp_in       = rmdi
REAL    :: tstar_in                = rmdi
REAL    :: zbl_in                  = rmdi
REAL    :: rough_length_in         = rmdi
REAL    :: seaice_frac_in          = rmdi
REAL    :: pstar_in                = rmdi
REAL    :: zhsc_in                 = rmdi
REAL    :: u_scalar_10m_in         = rmdi
REAL    :: u_s_in                  = rmdi
REAL    :: surf_albedo_in          = rmdi
REAL    :: dms_sea_conc_in         = rmdi
INTEGER :: kent_in                 = imdi
INTEGER :: kent_dsc_in             = imdi
INTEGER :: conv_cloud_base_in      = imdi
INTEGER :: conv_cloud_top_in       = imdi
REAL    :: theta_in                = rmdi
REAL    :: q_in                    = rmdi
REAL    :: qcf_in                  = rmdi
REAL    :: qcl_in                  = rmdi
REAL    :: conv_cloud_amount_in    = rmdi
REAL    :: rho_r2_in               = rmdi
REAL    :: exner_rho_levels_in     = rmdi
REAL    :: area_cloud_fraction_in  = rmdi
REAL    :: cloud_frac_in           = rmdi
REAL    :: cloud_liq_frac_in       = rmdi
REAL    :: exner_theta_levels_in   = rmdi
REAL    :: p_rho_levels_in         = rmdi
REAL    :: p_theta_levels_in       = rmdi
REAL    :: t_theta_levels_in       = rmdi
REAL    :: rhokh_mix_in            = rmdi
REAL    :: dtrdz_charney_grid_in   = rmdi
REAL    :: rhokh_rdz_in            = rmdi
REAL    :: dtrdz_in                = rmdi
REAL    :: we_lim_in               = rmdi
REAL    :: t_frac_in               = rmdi
REAL    :: zrzi_in                 = rmdi
REAL    :: we_lim_dsc_in           = rmdi
REAL    :: t_frac_dsc_in           = rmdi
REAL    :: zrzi_dsc_in             = rmdi
REAL    :: ls_rain3d_in            = rmdi
REAL    :: ls_snow3d_in            = rmdi
REAL    :: autoconv_in             = rmdi
REAL    :: accretion_in            = rmdi
REAL    :: pv_on_theta_mlevs_in    = rmdi
REAL    :: conv_rain3d_in          = rmdi
REAL    :: conv_snow3d_in          = rmdi
REAL    :: so4_sa_clim_in          = rmdi

!++SAN Variables that should be in int dry dep namelist, if used
!s REAL    :: soil_moisture_layer1_in = rmdi
!s LOGICAL :: l_tile_active_in        = .false.
!s REAL    :: frac_types_in           = rmdi
!s REAL    :: laift_lp_in             = rmdi
!s REAL    :: canhtft_lp_in           = rmdi
!s REAL    :: tstar_tile_in           = rmdi
!s REAL    :: z0tile_lp_in            = rmdi
!s REAL    :: surf_hf_in              = rmdi
!s REAL    :: stcon_in                = rmdi

!++SAN variables that should be in aerosol namelist, if used
!s REAL    :: rim_cry_in              = rmdi
!s REAL    :: rim_agg_in              = rmdi
!s REAL    :: vertvel_in              = rmdi
!s REAL    :: bl_tke_in               = rmdi

! List of environmental driver fields required for the UKCA configuration
CHARACTER(LEN=ukca_maxlen_fieldname), SAVE, POINTER :: environ_varnames(:)

! List of tracers required for the UKCA configuration
CHARACTER(LEN=ukca_maxlen_fieldname), POINTER :: tracer_varnames(:)
INTEGER              :: n_tracer_required

! List of non-transported prognostics required for the UKCA configuration
CHARACTER(LEN=ukca_maxlen_fieldname), POINTER :: ntp_varnames(:)
INTEGER              :: n_ntp_required

! Format string for outputting of box model
CHARACTER(LEN=80) :: out_format

! File for outputting hex values
CHARACTER(LEN=80) :: tracer_hex_filename
INTEGER :: ppos
INTEGER, ALLOCATABLE :: int_tracers(:)

!++SAN variables for tracer mapping
INTEGER :: nlines             !SAN to count number of lines in input file
INTEGER :: io_status
INTEGER :: ll                 ! Line loop variable
INTEGER :: i                  ! Loop variable
INTEGER :: j                  ! Loop variable
LOGICAL :: tracer_match

! file units
INTEGER, PARAMETER :: tracer_out_unit = 81
INTEGER, PARAMETER :: flux_out_unit   = 82
INTEGER, PARAMETER :: rate_out_unit   = 83
INTEGER, PARAMETER :: tracer_hex_unit = 84
INTEGER, PARAMETER :: ntp_out_unit    = 85

!++SAN local variables for calculating rho_r2
REAL    :: alt
REAL    :: theta_v
REAL    :: rho

!++SAN for lumping species with ukca_transform_halogen
LOGICAL, PARAMETER :: unlump_species = .TRUE.
LOGICAL, PARAMETER :: lump_species   = .FALSE.

! UKCA error reporting variables
CHARACTER(LEN=ukca_maxlen_message)  :: ukca_errmsg    ! Error return message
CHARACTER(LEN=ukca_maxlen_procname) :: ukca_errproc   ! Routine in which error
                                                      ! was trapped
! UM error reporting variables
INTEGER                             :: errcode=0      ! Error flag (0 = OK)
CHARACTER(LEN=*), PARAMETER :: errproc_suffix = ' in UKCA'

INTEGER, PARAMETER :: npft_tr_mix = 3 ! Number of plant functional types
                                      ! considered by tr_mix

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0
cmessage=''

! ----------------------------------------------------------------------
!  0. Start Timer call for BOX_MODEL
!
!++SAN
CALL umPrint('SAN BOX_MODEL: Startting box model',src='box_model')

IF (ltimer) CALL timer('BOX_MODEL ',5)

! Retrieve internal UKCA configuration variables that are required by box model
CALL ukca_get_config(l_ukca_intdd=l_ukca_intdd,                                &
                     ntype=ntype,                                              &
                     npft=npft,                                                &
                     l_ukca_strat=l_ukca_strat,                                &
                     l_ukca_strattrop=l_ukca_strattrop,                        &
                     l_ukca_stratcfc=l_ukca_stratcfc,                          &
                     l_ukca_cristrat=l_ukca_cristrat)

! ----------------------------------------------------------------------
!  0.1 General initialisation of control and physical data blocks
!

icode=0

CALL timer('INITIAL ',5)

! Hardcode some values for the box, using the grid spacing for 1.875*1.25 (N96)
delta_lambda  = 0.0327249 ! EW (x) grid spacing in radians
delta_phi     = 0.0218166 ! NS (y) grid spacing in radians

!++SAN allocate a dummy stash array
nsects = 52
nitems = 999
CALL allocate_stash_arrays()
! Make sure all of sf is false, so that it is not used
sf(:,:) = .FALSE.

!++SAN - replace with box model version of initialisation
CALL umPrint('SAN BOX_MODEL: Initialising box model',src='box_model')
!++SAN - initialising box model is actually very simple (for now)
!        Decide to only make a couple of calls in this routine rather
!        than seperate module. If this grows to become quite cumbersome,
!        may decide to make a new initial_box routine which does all of
!        these step on first call
submodel=submodel_partition_list(1)
!  Initialise the model time and check that history file data time
CALL umPrint('SAN BOX_MODEL: calling inittime',src='box_model')
CALL inittime(submodel)
!  Set timestep control switches for initial step
CALL umPrint('SAN BOX_MODEL: calling settsctl',src='box_model')
CALL settsctl (internal_model,meanlev,icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,'(A)') 'Failure in call to SETTSCTL'
  CALL umPrint(umMessage,src='initial_box')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

! Get list of environment fields required by the current UKCA configuration
NULLIFY(environ_varnames)
!++SAN
CALL umPrint('SAN; BOX_MODEL: ukca_get_environment_varlist',src='box_model')

CALL ukca_get_environment_varlist(environ_varnames, errcode,                 &
                                  error_message=ukca_errmsg,                 &
                                  error_routine=ukca_errproc)

! ----------------------------------------------------------------------
!  0.2 Allocate all of the memory needed for the environment variables
!      and set defaults for the arrays.
!      These can be checked against the list required by the model

WRITE(umMessage,*) 'SAN; BOX_MODEL: ALLOCATING ALL of the things'
CALL umPrint(umMessage, src='box_model')

!++SAN Allocate tracer and npt arrays, needed for running UKCA
CALL umPrint(umMessage, src='box_model')
DO i = 1, SIZE(environ_varnames)
  WRITE(umMessage,*) 'SAN; BOX_MODEL: ALLOCATING ', TRIM(environ_varnames(i))
  CALL umPrint(umMessage, src='box_model')
  SELECT CASE (TRIM(environ_varnames(i)))

    !++SAN Allocate ALL of the environment variables
  CASE("land_sea_mask")
    IF (.NOT. ALLOCATED(land_sea_mask)) &
        ALLOCATE(land_sea_mask(row_length, rows))
    land_sea_mask(:,:)         = .FALSE.

  CASE("conv_cloud_lwp")
    IF (.NOT. ALLOCATED(conv_cloud_lwp)) &
        ALLOCATE(conv_cloud_lwp(row_length, rows))
    conv_cloud_lwp(:,:)        = 0.0

  CASE("tstar")
    IF (.NOT. ALLOCATED(tstar)) &
        ALLOCATE(tstar(row_length, rows))
    tstar(:,:)                 = 0.0

  CASE("zbl")
    IF (.NOT. ALLOCATED(zbl)) &
        ALLOCATE(zbl(row_length, rows))
    zbl(:,:)                   = 0.0

  CASE("rough_length")
    IF (.NOT. ALLOCATED(rough_length)) &
        ALLOCATE(rough_length(row_length, rows))
    rough_length(:,:)          = 0.0

  CASE("seaice_frac")
    IF (.NOT. ALLOCATED(seaice_frac)) &
        ALLOCATE(seaice_frac(row_length, rows))
    seaice_frac(:,:)           = 0.0

  CASE("pstar")
    IF (.NOT. ALLOCATED(pstar)) &
        ALLOCATE(pstar(row_length, rows))
    pstar(:,:)                 = 0.0

  CASE("zhsc")
    IF (.NOT. ALLOCATED(zhsc)) &
        ALLOCATE(zhsc(row_length, rows))
    zhsc(:,:)                  = 0.0

  CASE("u_scalar_10m")
    IF (.NOT. ALLOCATED(u_scalar_10m)) &
        ALLOCATE(u_scalar_10m(row_length, rows))
    u_scalar_10m(:,:)          = 0.0

  CASE("u_s")
    IF (.NOT. ALLOCATED(u_s)) &
        ALLOCATE(u_s(row_length, rows))
    u_s(:,:)                   = 0.0

  CASE("kent")
    IF (.NOT. ALLOCATED(kent)) &
        ALLOCATE(kent(row_length, rows))
    kent(:,:)                  = 0

  CASE("kent_dsc")
    IF (.NOT. ALLOCATED(kent_dsc)) &
        ALLOCATE(kent_dsc(row_length, rows))
    kent_dsc(:,:)              = 0

  CASE("theta")
    IF (.NOT. ALLOCATED(theta)) &
        ALLOCATE(theta(row_length, rows, model_levels))
    theta(:,:,:)               = 0.0

  CASE("q")
    IF (.NOT. ALLOCATED(q)) &
        ALLOCATE(q(row_length, rows, model_levels))
    q(:,:,:)                   = 0.0

  CASE("qcf")
    IF (.NOT. ALLOCATED(qcf)) &
        ALLOCATE(qcf(row_length, rows, model_levels))
    qcf(:,:,:)                 = 0.0

  CASE("qcl")
    IF (.NOT. ALLOCATED(qcl)) &
        ALLOCATE(qcl(row_length, rows, model_levels))
    qcl(:,:,:)                 = 0.0

  CASE("conv_cloud_amount")
    IF (.NOT. ALLOCATED(conv_cloud_amount)) &
        ALLOCATE(conv_cloud_amount(row_length, rows, model_levels))
    conv_cloud_amount(:,:,:)   = 0.0

  CASE("rho_r2")
    IF (.NOT. ALLOCATED(rho_r2)) &
        ALLOCATE(rho_r2(row_length, rows, model_levels+1))
    rho_r2(:,:,:)              = 0.0

  CASE("exner_rho_levels")
    IF (.NOT. ALLOCATED(exner_rho_levels)) &
        ALLOCATE(exner_rho_levels(row_length, rows, model_levels+1))
    exner_rho_levels(:,:,:)    = 0.0

  CASE("area_cloud_fraction")
    IF (.NOT. ALLOCATED(area_cloud_fraction)) &
        ALLOCATE(area_cloud_fraction(row_length, rows, model_levels))
    area_cloud_fraction(:,:,:) = 0.0

  CASE("cloud_frac")
    IF (.NOT. ALLOCATED(cloud_frac)) &
        ALLOCATE(cloud_frac(row_length, rows, model_levels))
    cloud_frac(:,:,:)          = 0.0

  CASE("cloud_liq_frac")
    IF (.NOT. ALLOCATED(cloud_liq_frac)) &
        ALLOCATE(cloud_liq_frac(row_length, rows, model_levels))
    cloud_liq_frac(:,:,:)      = 0.0

  CASE("exner_theta_levels")
    IF (.NOT. ALLOCATED(exner_theta_levels)) &
        ALLOCATE(exner_theta_levels(row_length, rows, model_levels))
    exner_theta_levels(:,:,:)  = 0.0

  CASE("p_rho_levels")
    IF (.NOT. ALLOCATED(p_rho_levels)) &
        ALLOCATE(p_rho_levels(row_length, rows, model_levels+1))
    p_rho_levels(:,:,:)      = 0.0

  CASE("p_theta_levels")
    IF (.NOT. ALLOCATED(p_theta_levels)) &
        ALLOCATE(p_theta_levels(row_length, rows, model_levels))
    p_theta_levels(:,:,:)      = 0.0

  CASE("rhokh_mix")
    IF (.NOT. ALLOCATED(rhokh_mix)) &
        ALLOCATE(rhokh_mix(row_length, rows, model_levels))
    rhokh_mix(:,:,:)           = 0.0

  CASE("dtrdz_charney_grid")
    IF (.NOT. ALLOCATED(dtrdz_charney_grid)) &
        ALLOCATE(dtrdz_charney_grid(row_length, rows, model_levels))
    dtrdz_charney_grid(:,:,:)  = 0.0

  ! SAN - new variables needed for vn12.1
  CASE("rhokh_rdz")
    IF (.NOT. ALLOCATED(rhokh_rdz)) &
        ALLOCATE(rhokh_rdz(row_length, rows, model_levels))
    rhokh_rdz(:,:,:)  = 0.0

  CASE("dtrdz")
    IF (.NOT. ALLOCATED(dtrdz)) &
        ALLOCATE(dtrdz(row_length, rows, model_levels))
    dtrdz(:,:,:)  = 0.0

    !++SAN variables related to the tracer mixing fields have npft_tr_mix
    !      as third dimension - number of plant functional types
  CASE("we_lim")
    IF (.NOT. ALLOCATED(we_lim)) &
        ALLOCATE(we_lim(row_length, rows, npft_tr_mix))
    we_lim(:,:,:)              = 0.0

  CASE("t_frac")
    IF (.NOT. ALLOCATED(t_frac)) &
        ALLOCATE(t_frac(row_length, rows, npft_tr_mix))
    t_frac(:,:,:)              = 0.0

  CASE("zrzi")
    IF (.NOT. ALLOCATED(zrzi)) &
        ALLOCATE(zrzi(row_length, rows, npft_tr_mix))
    zrzi(:,:,:)                = 0.0

  CASE("we_lim_dsc")
    IF (.NOT. ALLOCATED(we_lim_dsc)) &
        ALLOCATE(we_lim_dsc(row_length, rows, npft_tr_mix))
    we_lim_dsc(:,:,:)          = 0.0

  CASE("t_frac_dsc")
    IF (.NOT. ALLOCATED(t_frac_dsc)) &
        ALLOCATE(t_frac_dsc(row_length, rows, npft_tr_mix))
    t_frac_dsc(:,:,:)          = 0.0

  CASE("zrzi_dsc")
    IF (.NOT. ALLOCATED(zrzi_dsc)) &
        ALLOCATE(zrzi_dsc(row_length, rows, npft_tr_mix))
    zrzi_dsc(:,:,:)            = 0.0

  CASE("ls_rain3d")
    IF (.NOT. ALLOCATED(ls_rain3d)) &
        ALLOCATE(ls_rain3d(row_length, rows, model_levels))
    ls_rain3d(:,:,:)           = 0.0 !SAN - in future, might be possible to set this for wet dep?

  CASE("ls_snow3d")
    IF (.NOT. ALLOCATED(ls_snow3d)) &
        ALLOCATE(ls_snow3d(row_length, rows, model_levels))
    ls_snow3d(:,:,:)           = 0.0

  CASE("autoconv")
    IF (.NOT. ALLOCATED(autoconv)) &
        ALLOCATE(autoconv(row_length, rows, model_levels))
    autoconv(:,:,:)            = 0.0

  CASE("accretion")
    IF (.NOT. ALLOCATED(accretion)) &
        ALLOCATE(accretion(row_length, rows, model_levels))
    accretion(:,:,:)           = 0.0

  CASE("pv_on_theta_mlevs")
    IF (.NOT. ALLOCATED(pv_on_theta_mlevs)) &
        ALLOCATE(pv_on_theta_mlevs(row_length, rows, model_levels))
    pv_on_theta_mlevs(:,:,:)   = 0.0

  CASE("conv_cloud_base")
    IF (.NOT. ALLOCATED(conv_cloud_base)) &
        ALLOCATE(conv_cloud_base(row_length, rows))
    conv_cloud_base(:,:)       = 0

  CASE("conv_cloud_top")
    IF (.NOT. ALLOCATED(conv_cloud_top)) &
        ALLOCATE(conv_cloud_top(row_length, rows))
    conv_cloud_top(:,:)        = 0

  CASE("conv_rain3d")
    IF (.NOT. ALLOCATED(conv_rain3d)) &
        ALLOCATE(conv_rain3d(row_length, rows, model_levels))
    conv_rain3d(:,:,:)         = 0.0

  CASE("conv_snow3d")
    IF (.NOT. ALLOCATED(conv_snow3d)) &
        ALLOCATE(conv_snow3d(row_length, rows, model_levels))
    conv_snow3d(:,:,:)         = 0.0

  CASE("so4_sa_clim")
    IF (.NOT. ALLOCATED(so4_sa_clim)) &
        ALLOCATE(so4_sa_clim(row_length, rows, model_levels))
    so4_sa_clim(:,:,:)         = 0.0 !SAN Might need setting in future...

    !++SAN Allocate environment variables needed for running with aerosols
  CASE("dms_sea_conc")
    IF (.NOT. ALLOCATED(dms_sea_conc)) &
        ALLOCATE(dms_sea_conc(row_length, rows))
    dms_sea_conc(:,:)          = 0.0

        !++ Aerosol environment variables
  CASE("rim_cry")
    IF (.NOT. ALLOCATED(rim_cry)) &
        ALLOCATE(rim_cry(row_length, rows, model_levels))
    rim_cry(:,:,:)             = 0.0

  CASE("rim_agg")
    IF (.NOT. ALLOCATED(rim_agg)) &
        ALLOCATE(rim_agg(row_length, rows, model_levels))
    rim_agg(:,:,:)             = 0.0

  CASE("vertvel")
    IF (.NOT. ALLOCATED(vertvel)) &
        ALLOCATE(vertvel(row_length, rows, model_levels))
    vertvel(:,:,:)             = 0.0

  CASE("bl_tke")
    IF (.NOT. ALLOCATED(bl_tke)) &
        ALLOCATE(bl_tke(row_length, rows, model_levels))
    bl_tke(:,:,:)              = 0.0
  END SELECT
END DO

! Allocate more variables...
IF (.NOT. ALLOCATED(r_rho_levels)) THEN
  ALLOCATE(r_rho_levels(row_length, rows, 1:model_levels))
  r_rho_levels(:,:,:)          = 0.0
END IF
IF (.NOT. ALLOCATED(r_theta_levels)) THEN
  ALLOCATE(r_theta_levels(row_length, rows, 0:model_levels))
  r_theta_levels(:,:,:)        = 0.0
END IF
IF (.NOT. ALLOCATED(eta_theta_levels)) THEN
  ALLOCATE(eta_theta_levels(0:model_levels))
  eta_theta_levels(:)          = 0.0
END IF
IF (.NOT. ALLOCATED(t_theta_levels)) THEN
  ALLOCATE(t_theta_levels(row_length, rows, model_levels))
  t_theta_levels(:,:,:)        = 0.0
END IF
IF (.NOT. ALLOCATED(surf_albedo)) THEN
  ALLOCATE(surf_albedo(row_length, rows))
  surf_albedo(:,:)             = 0.0
END IF

  !++SAN Allocate all of the trig variables used by UKCA
IF (.NOT. ALLOCATED(true_latitude)) THEN
  ALLOCATE(true_latitude(row_length, rows))
  true_latitude(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(true_longitude)) THEN
  ALLOCATE(true_longitude(row_length, rows))
  true_longitude(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(FV_cos_theta_latitude)) THEN
  ALLOCATE(FV_cos_theta_latitude(row_length, rows))
  FV_cos_theta_latitude(:,:) = 0.0
END IF

!++SAN Allocate SF array to max possible size of STASH items
IF (.NOT. ALLOCATED(sf)) ALLOCATE(sf(999,52))



! ----------------------------------------------------------------------
!  1.1 Read in data to define environmental parameters
!
CALL read_environment_pars(icode, iomessage, ukca_box_nml,                    &
  latitude_in, longitude_in, land_sea_mask_in, conv_cloud_lwp_in, tstar_in,   &
  zbl_in, rough_length_in, seaice_frac_in, pstar_in, zhsc_in,                 &
  u_scalar_10m_in, u_s_in, surf_albedo_in, dms_sea_conc_in, kent_in,          &
  kent_dsc_in, conv_cloud_base_in, conv_cloud_top_in, q_in, qcf_in, qcl_in,   &
  conv_cloud_amount_in, rho_r2_in, area_cloud_fraction_in, cloud_frac_in,     &
  cloud_liq_frac_in, p_theta_levels_in, t_theta_levels_in, rhokh_mix_in,      &
  dtrdz_charney_grid_in, rhokh_rdz_in, dtrdz_in, we_lim_in, t_frac_in,        &
  zrzi_in, we_lim_dsc_in, t_frac_dsc_in, zrzi_dsc_in, ls_rain3d_in,           &
  ls_snow3d_in, autoconv_in, accretion_in, pv_on_theta_mlevs_in,              &
  conv_rain3d_in, conv_snow3d_in, so4_sa_clim_in)

IF (icode /= 0) THEN
  CALL Ereport(RoutineName,icode,iomessage)
END IF

!++SAN check these have been read in properly
WRITE(umMessage,*) 'SAN; BOX_MODEL: latitude_in = ', latitude_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: longitude_in = ', longitude_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: t_theta_levels_in = ', t_theta_levels_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: p_theta_levels_in = ', p_theta_levels_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: land_sea_mask_in = ', land_sea_mask_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: tstar_in = ', tstar_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: zbl_in = ', zbl_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: pstar_in = ', pstar_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: q_in = ', q_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: qcf_in = ', qcf_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: qcl_in = ', qcl_in
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: surf_albedo_in = ', surf_albedo_in
CALL umPrint(umMessage, src='box_model')
!--SAN

!++SAN Do some checks on the input values
IF (p_theta_levels_in  >= pstar_in ) THEN
  WRITE(Cmessage, *) 'Surface pressure must be greater than pressure'
  CALL umPrint(Cmessage, src='box_model')
  WRITE(umMessage, *) 'Surface pressure = ', pstar_in 
  CALL umPrint(umMessage, src='box_model')
  WRITE(umMessage, *) 'Model level pressure = ', p_theta_levels_in 
  CALL umPrint(umMessage, src='box_model')
  icode = 1
  CALL ereport(RoutineName, errcode, Cmessage)
END IF

!++Map variables which have been read in to correct internal fields
WRITE(umMessage,*) 'SAN; BOX_MODEL: Mapping variables'
CALL umPrint(umMessage, src='box_model')
!++SAN convert lat and lon into radians
IF (latitude_in /= rmdi) THEN
  true_latitude(:,:) = latitude_in * Pi_over_180
  !++SAN do a dummy calculation of FV_cos_theta_latitude
  FV_cos_theta_latitude(:,:) = COS(latitude_in * Pi_over_180)
END IF
IF (longitude_in /= rmdi) true_longitude(:,:) = longitude_in * Pi_over_180

! Map across all the other variables
IF (land_sea_mask_in .AND. ALLOCATED(land_sea_mask)) THEN
  WRITE(umMessage,*) '   Setting land_sea_mask to ', land_sea_mask_in
  CALL umPrint(umMessage, src='box_model')
  land_sea_mask(:,:)             = land_sea_mask_in
END IF
IF (conv_cloud_lwp_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting conv_cloud_lwp to ', conv_cloud_lwp_in
  CALL umPrint(umMessage, src='box_model')
  conv_cloud_lwp(:,:)   = conv_cloud_lwp_in
END IF
IF (tstar_in /= rmdi .AND. ALLOCATED(tstar)) THEN
  WRITE(umMessage,*) '   Setting tstar to ', tstar_in
  CALL umPrint(umMessage, src='box_model')
  tstar(:,:)                     = tstar_in
END IF
IF (zbl_in /= rmdi .AND. ALLOCATED(zbl)) THEN
  WRITE(umMessage,*) '   Setting zbl to ', zbl_in
  CALL umPrint(umMessage, src='box_model')
  zbl(:,:)                         = zbl_in
END IF
IF (rough_length_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting rough_length to ', rough_length_in
  CALL umPrint(umMessage, src='box_model')
  rough_length(:,:)       = rough_length_in
END IF
IF (seaice_frac_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting seaice_frac to ', seaice_frac_in
  CALL umPrint(umMessage, src='box_model')
  seaice_frac(:,:)         = seaice_frac_in
END IF
IF (pstar_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting pstar to ', pstar_in
  CALL umPrint(umMessage, src='box_model')
  pstar(:,:)                     = pstar_in
END IF
IF (zhsc_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting zhsc to ', zhsc_in
  CALL umPrint(umMessage, src='box_model')
  zhsc(:,:)                       = zhsc_in
END IF
IF (u_scalar_10m_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting u_scalar_10m to ', u_scalar_10m_in
  CALL umPrint(umMessage, src='box_model')
  u_scalar_10m(:,:)       = u_scalar_10m_in
END IF
IF (u_s_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting u_s to ', u_s_in
  CALL umPrint(umMessage, src='box_model')
  u_s(:,:)                         = u_s_in
END IF
IF (surf_albedo_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting surf_albedo to ', surf_albedo_in
  CALL umPrint(umMessage, src='box_model')
  surf_albedo(:,:)         = surf_albedo_in
END IF
IF (kent_in /= imdi) THEN
  WRITE(umMessage,*) '   Setting kent to ', kent_in
  CALL umPrint(umMessage, src='box_model')
  kent(:,:)                       = kent_in
END IF
IF (kent_dsc_in /= imdi) THEN
  WRITE(umMessage,*) '   Setting kent_dsc to ', kent_dsc_in
  CALL umPrint(umMessage, src='box_model')
  kent_dsc(:,:)               = kent_dsc_in
END IF
IF (theta_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting theta to ', theta_in
  CALL umPrint(umMessage, src='box_model')
  theta(:,:,:)                   = theta_in
END IF
IF (q_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting q to ', q_in
  CALL umPrint(umMessage, src='box_model')
  q(:,:,:)                           = q_in
END IF
IF (qcf_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting qcf to ', qcf_in
  CALL umPrint(umMessage, src='box_model')
  qcf(:,:,:)                       = qcf_in
END IF
IF (qcl_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting qcl to ', qcl_in
  CALL umPrint(umMessage, src='box_model')
  qcl(:,:,:)                       = qcl_in
END IF
IF (conv_cloud_amount_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting conv_cloud_amount to ', conv_cloud_amount_in
  CALL umPrint(umMessage, src='box_model')
  conv_cloud_amount(:,:,:)  = conv_cloud_amount_in
END IF
IF (rho_r2_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting rho_r2 to ', rho_r2_in
  CALL umPrint(umMessage, src='box_model')
  rho_r2(:,:,:)                 = rho_r2_in
END IF
IF (exner_rho_levels_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting exner_rho_levels to ', exner_rho_levels_in
  CALL umPrint(umMessage, src='box_model')
  exner_rho_levels(:,:,:)  = exner_rho_levels_in
END IF
IF (area_cloud_fraction_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting area_cloud_fraction to ', area_cloud_fraction_in
  CALL umPrint(umMessage, src='box_model')
  area_cloud_fraction(:,:,:)  = area_cloud_fraction_in
END IF
IF (cloud_frac_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting cloud_frac to ', cloud_frac_in
  CALL umPrint(umMessage, src='box_model')
  cloud_frac(:,:,:)         = cloud_frac_in
END IF
IF (cloud_liq_frac_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting cloud_liq_frac to ', cloud_liq_frac_in
  CALL umPrint(umMessage, src='box_model')
  cloud_liq_frac(:,:,:) = cloud_liq_frac_in
END IF
IF (exner_theta_levels_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting exner_theta_levels to ', exner_theta_levels_in
  CALL umPrint(umMessage, src='box_model')
  exner_theta_levels(:,:,:)  = exner_theta_levels_in
END IF
IF (p_rho_levels_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting p_theta_levels to ', p_rho_levels_in
  CALL umPrint(umMessage, src='box_model')
  p_rho_levels(:,:,:)   = p_rho_levels_in
END IF
IF (p_theta_levels_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting p_theta_levels to ', p_theta_levels_in
  CALL umPrint(umMessage, src='box_model')
  p_theta_levels(:,:,:) = p_theta_levels_in
END IF
IF (t_theta_levels_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting t_theta_levels to ', t_theta_levels_in
  CALL umPrint(umMessage, src='box_model')
  t_theta_levels(:,:,:) = t_theta_levels_in
END IF
IF (rhokh_mix_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting rhokh_mix to ', rhokh_mix_in
  CALL umPrint(umMessage, src='box_model')
  rhokh_mix(:,:,:)           = rhokh_mix_in
END IF
IF (dtrdz_charney_grid_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting dtrdz_charney_grid to ', dtrdz_charney_grid_in
  CALL umPrint(umMessage, src='box_model')
  dtrdz_charney_grid(:,:,:)  = dtrdz_charney_grid_in
END IF
IF (rhokh_rdz_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting rhokh_rdz to ', rhokh_rdz_in
  CALL umPrint(umMessage, src='box_model')
  rhokh_rdz(:,:,:)  = rhokh_rdz_in
END IF
IF (dtrdz_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting dtrdz to ', dtrdz_in
  CALL umPrint(umMessage, src='box_model')
  dtrdz(:,:,:)  = dtrdz_in
END IF

IF (we_lim_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting we_lim to ', we_lim_in
  CALL umPrint(umMessage, src='box_model')
  we_lim(:,:,:)                 = we_lim_in
END IF
IF (t_frac_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting t_frac to ', t_frac_in
  CALL umPrint(umMessage, src='box_model')
  t_frac(:,:,:)                 = t_frac_in
END IF
IF (zrzi_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting zrzi to ', zrzi_in
  CALL umPrint(umMessage, src='box_model')
  zrzi(:,:,:)                     = zrzi_in
END IF
IF (we_lim_dsc_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting we_lim_dsc to ', we_lim_dsc_in
  CALL umPrint(umMessage, src='box_model')
  we_lim_dsc(:,:,:)         = we_lim_dsc_in
END IF
IF (t_frac_dsc_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting t_frac_dsc to ', t_frac_dsc_in
  CALL umPrint(umMessage, src='box_model')
  t_frac_dsc(:,:,:)         = t_frac_dsc_in
END IF
IF (zrzi_dsc_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting zrzi_dsc to ', zrzi_dsc_in
  CALL umPrint(umMessage, src='box_model')
  zrzi_dsc(:,:,:)             = zrzi_dsc_in
END IF
IF (ls_rain3d_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting ls_rain3d to ', ls_rain3d_in
  CALL umPrint(umMessage, src='box_model')
  ls_rain3d(:,:,:)           = ls_rain3d_in
END IF
IF (ls_snow3d_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting ls_snow3d to ', ls_snow3d_in
  CALL umPrint(umMessage, src='box_model')
  ls_snow3d(:,:,:)           = ls_snow3d_in
END IF
IF (autoconv_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting autoconv to ', autoconv_in
  CALL umPrint(umMessage, src='box_model')
  autoconv(:,:,:)             = autoconv_in
END IF
IF (accretion_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting accretion to ', accretion_in
  CALL umPrint(umMessage, src='box_model')
  accretion(:,:,:)           = accretion_in
END IF
IF (pv_on_theta_mlevs_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting pv_on_theta_mlevs to ', pv_on_theta_mlevs_in
  CALL umPrint(umMessage, src='box_model')
  pv_on_theta_mlevs(:,:,:)   = pv_on_theta_mlevs_in
END IF
IF (conv_cloud_base_in /= imdi) THEN
  WRITE(umMessage,*) '   Setting conv_cloud_base to ', conv_cloud_base_in
  CALL umPrint(umMessage, src='box_model')
  conv_cloud_base(:,:) = conv_cloud_base_in
END IF
IF (conv_cloud_top_in /= imdi) THEN
  WRITE(umMessage,*) '   Setting conv_cloud_top to ', conv_cloud_top_in
  CALL umPrint(umMessage, src='box_model')
  conv_cloud_top(:,:)   = conv_cloud_top_in
END IF
IF (conv_rain3d_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting conv_rain3d to ', conv_rain3d_in
  CALL umPrint(umMessage, src='box_model')
  conv_rain3d(:,:,:)    = conv_rain3d_in
END IF
IF (conv_snow3d_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting conv_snow3d to ', conv_snow3d_in
  CALL umPrint(umMessage, src='box_model')
  conv_snow3d(:,:,:)    = conv_snow3d_in
END IF
IF (so4_sa_clim_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting so4_sa_clim to ', so4_sa_clim_in
  CALL umPrint(umMessage, src='box_model')
  so4_sa_clim(:,:,:)    = so4_sa_clim_in
END IF
IF (dms_sea_conc_in /= rmdi) THEN
  WRITE(umMessage,*) '   Setting dms_sea_conc to ', dms_sea_conc_in
  CALL umPrint(umMessage, src='box_model')
  dms_sea_conc(:,:)     = dms_sea_conc_in
END IF

!++SAN - calculating derived atmospheric variables:
WRITE(umMessage,*) 'SAN; BOX_MODEL: Calculating derived variables'
CALL umPrint(umMessage, src='box_model')
! Calculate exner_theta_levels based on pressure
exner_theta_levels(1,1,1) = ( p_theta_levels(1,1,1) / p_zero ) ** kappa
! Calculate theta from T, p and exner:
theta(1,1,1) = t_theta_levels(1,1,1) / exner_theta_levels(1,1,1) 
! Set rho and theta levels to be the same in box model
p_rho_levels(1,1,:) = p_theta_levels(1,1,1)

IF (ALLOCATED(rho_r2) .AND. .NOT. ALLOCATED(exner_rho_levels))                 &
  ALLOCATE(exner_rho_levels(row_length, rows, model_levels))
IF (ALLOCATED(exner_rho_levels))                                               &
  exner_rho_levels(1,1,:) = exner_theta_levels(1,1,1)

!++SAN Additional variables for running with interactive dry dep
!      Set to dummy values for now
IF (l_ukca_intdd) THEN

    IF (land_sea_mask_in) THEN
        land_points = 1
    ELSE
        land_points = 0
    END IF

    IF (.NOT. ALLOCATED(soil_moisture_layer1))                                 &
        ALLOCATE(soil_moisture_layer1(land_points))
    IF (.NOT. ALLOCATED(l_tile_active))                                        &
        ALLOCATE(l_tile_active(land_points, ntype))
    IF (.NOT. ALLOCATED(frac_types))                                           &
        ALLOCATE(frac_types(land_points, ntype))
    IF (.NOT. ALLOCATED(laift_lp))                                             &
        ALLOCATE(laift_lp(land_points, npft))
    IF (.NOT. ALLOCATED(canhtft_lp))                                           &
        ALLOCATE(canhtft_lp(land_points, npft))
    IF (.NOT. ALLOCATED(tstar_tile))                                           &
        ALLOCATE(tstar_tile(land_points, ntype))
    IF (.NOT. ALLOCATED(z0tile_lp))                                            &
        ALLOCATE(z0tile_lp(land_points, ntype))
    IF (.NOT. ALLOCATED(surf_hf))                                              &
        ALLOCATE(surf_hf(row_length, rows))
    IF (.NOT. ALLOCATED(stcon))                                                &
        ALLOCATE(stcon(model_levels, row_length, ntype))

  !++SAN set other environment variables to dummy values, these aren't important for box
  WRITE(umMessage,*) 'SAN; BOX_MODEL: Setting interactive DD variables'

  ! ntype needs to be given a value to work
  soil_moisture_layer1(1)  = 0.0
  l_tile_active(1,:)       = .FALSE.
  frac_types(1,:)          = 0.0
  laift_lp(1,:)            = 0.0
  canhtft_lp(1,:)          = 0.0
  tstar_tile(1,:)          = tstar_in
  z0tile_lp(1,:)           = 0.0
  surf_hf(1,1)             = 0.0
  stcon(1,1,1)             = 0.0

  !++SAN print size of l_tile_active
  WRITE(umMessage,*) 'SAN: l_tile_active 1 =  ', SIZE(l_tile_active, DIM=1)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,*) 'SAN: l_tile_active 2 =  ', SIZE(l_tile_active, DIM=2)
  CALL umPrint(umMessage,src=RoutineName)
  !--SAN

END IF ! l_ukca_intdd

!-----------------------------------------------------------------------

!++SAN Calculate model height from p, using scale height
alt = -sclht * LOG(p_theta_levels(1,1,1)/pstar(1,1))
r_theta_levels(1,1,1) = alt + planet_radius
r_theta_levels(1,1,0) = planet_radius !S Set r_thera levels 0 to sfc?
r_rho_levels(1,1,:) = r_theta_levels(1,1,1)
! Calculate eta, if box level is by definition the top of model
! then eta = 0 at surface and 1 at location of box
eta_theta_levels(0) = 0.
eta_theta_levels(1) = 1.

!++SAN calculate density rho from p and t, using formula from cal_rho
theta_v = theta(1,1,1) * ( 1.0 + c_virtual * (q(1,1,1)+qcl(1,1,1)+qcf(1,1,1)))
IF (ALLOCATED(rho_r2)) THEN
  rho = p_rho_levels(1,1,1) / (r * theta_v * exner_rho_levels(1,1,1))
  rho_r2(1,1,:) = rho * r_rho_levels(1,1,:) * r_rho_levels(1,1,:)
END IF

!++SAN Doing some temporary print statements here to find out if
!      values for key variables have been set or not
WRITE(umMessage,*) 'SAN; BOX_MODEL: alt = ', alt
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: r_rho_levels(1,1,1) = ', r_rho_levels(1,1,1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: r_theta_levels(1,1,1) = ', r_theta_levels(1,1,1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: eta_theta_levels(1) = ', eta_theta_levels(1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: p_theta_levels(1,1,1) = ', p_theta_levels(1,1,1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: exner_theta_levels(1,1,1) = ', exner_theta_levels(1,1,1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: t_theta_levels(1,1,1) = ', t_theta_levels(1,1,1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: theta = ', theta(1,1,1)
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: theta_v = ', theta_v
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; BOX_MODEL: c_virtual = ', c_virtual
CALL umPrint(umMessage, src='box_model')
IF (ALLOCATED(exner_rho_levels)) THEN
  WRITE(umMessage,*) 'SAN; BOX_MODEL: exner_rho_levels(1,1,1) = ', exner_rho_levels(1,1,1)
  CALL umPrint(umMessage, src='box_model')
END IF
IF (ALLOCATED(rho_r2)) THEN
  WRITE(umMessage,*) 'SAN; BOX_MODEL: rho = ', rho
  CALL umPrint(umMessage, src='box_model')
  WRITE(umMessage,*) 'SAN; BOX_MODEL: rho_r2(1,1,1) = ', rho_r2(1,1,1)
  CALL umPrint(umMessage, src='box_model')
END IF


IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)

!++SAN Set all of sf to false
sf(:,:) = .FALSE.


! ----------------------------------------------------------------------
!  1.2 Read in data to define environmental parameters
!
WRITE(umMessage,*) 'SAN; BOX_MODEL: Opening chem tracer input: ', tracer_in_filename
CALL umPrint(umMessage, src='box_model')
OPEN(114,FILE=tracer_in_filename)
! First find out how many lines there are in the file
!=> Might want to edit so it can handle a header
nlines = 0
DO
  READ(114,*,iostat=io_status)
  IF (io_status/=0) EXIT
  nlines=nlines+1
END DO
CLOSE(114)
WRITE(umMessage,*) 'SAN; BOX_MODEL: nlines in chem_tracer_pars = ', nlines
CALL umPrint(umMessage, src='box_model')

! Allocate chemistry arrays based on 
IF (.NOT. ALLOCATED(inp_tracer_mmrs)) ALLOCATE(inp_tracer_mmrs(nlines))
IF (.NOT. ALLOCATED(inp_tracer_names)) ALLOCATE(inp_tracer_names(nlines))

! Loop through file, mapping data
WRITE(umMessage,*) 'SAN; Input chem MMRs:'
CALL umPrint(umMessage, src='box_model')
OPEN(115,FILE=tracer_in_filename)
DO ll = 1, nlines
  READ(115,'(e10.4,1X,A10)') inp_tracer_mmrs(ll), inp_tracer_names(ll)
  WRITE(umMessage,*) ' * ', ll, '. ', inp_tracer_names(ll), ' = ', inp_tracer_mmrs(ll)
  CALL umPrint(umMessage, src='box_model')
END DO
CLOSE(115)

! Load up the tracer varnames needed for the UKCA run
!++SAN
CALL umPrint('SAN; BOX_MODEL: set up tracer mapping',src='box_model')

! Initialise list pointer for safety
NULLIFY(tracer_varnames)

! Get list of tracers required by the current UKCA configuration
CALL ukca_get_tracer_varlist(tracer_varnames, errcode,                       &
                             error_message=ukca_errmsg,                      &
                             error_routine=ukca_errproc)
IF (errcode > 0) THEN
  cmessage = ukca_errmsg
  CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
END IF


n_tracer_required = SIZE(tracer_varnames)

!++SAN print all tracer varnames:
IF (PrintStatus >= PrStatus_Oper) THEN
  WRITE(umMessage,'(A33,I3)') 'SAN; box_model n_tracer_required: ', n_tracer_required
  CALL umPrint(umMessage, src='box_model')
  WRITE(umMessage,*) 'SAN; box_model tracer_varnames: '
  CALL umPrint(umMessage, src='box_model')
  DO i = 1, n_tracer_required
    WRITE(umMessage,'(A3,I3,A2,A10)') ' * ', i, '. ', tracer_varnames(i)
    CALL umPrint(umMessage, src='box_model')
  END DO
END IF

!++SAN Will try to split this mapping between chemical and aerosol species
!      First, double check whether the number of chemical/aerosol species
!      has been set by this point
WRITE(umMessage,*) 'SAN; box_model n_chem_tracers = ', n_chem_tracers
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; box_model n_aero_tracers = ', n_aero_tracers
CALL umPrint(umMessage, src='box_model')
WRITE(umMessage,*) 'SAN; box_model n_mode_tracers = ', n_mode_tracers
CALL umPrint(umMessage, src='box_model')

! Allocate arrays to file with UKCA data
tr_ukca = n_tracer_required
IF (.NOT. ALLOCATED(tracer_ukca)) &
    ALLOCATE(tracer_ukca(row_length, rows, model_levels, tr_ukca))

! Initialise tracer array to small values
tracer_ukca(:,:,:,:) = tracer_nullval

DO i = 1, n_tracer_required
  tracer_match = .FALSE.
  DO ll=1, nlines
    IF (TRIM(tracer_varnames(i)) == TRIM(inp_tracer_names(ll))) THEN
      tracer_match = .TRUE.
      tracer_ukca(:,:,:,i) = inp_tracer_mmrs(ll)
      EXIT
    ELSE
      CONTINUE
    END IF
  END DO

  ! Check whether tracer has been matched
  IF (tracer_match) THEN
    WRITE(umMessage,'(A3,I3,A13,A10,A3,E12.4)') ' * ',                          &
        i, ' MAPPED;     ', tracer_varnames(i), ' = ', tracer_ukca(1,1,1,i) 
    CALL umPrint(umMessage, src='box_model')
  ELSE
    WRITE(umMessage,'(A3,I3,A13,A10,A3,E12.4)') ' * ',                          &
        i, ' NOT MAPPED; ', tracer_varnames(i), ' = ', tracer_ukca(1,1,1,i) 
    CALL umPrint(umMessage, src='box_model')
  END IF
END DO


!++ Do the same thing with NTP array !!!
CALL umPrint('SAN; BOX_MODEL: set up NTP mapping',src='box_model')
! Initialise list pointers for safety
NULLIFY(ntp_varnames)

! Get list of NTP fields required by the current UKCA configuration
CALL ukca_get_ntp_varlist(ntp_varnames, errcode, error_message=ukca_errmsg,  &
                          error_routine=ukca_errproc)
IF (errcode > 0) THEN
  cmessage = ukca_errmsg
  CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
END IF

n_ntp_required = SIZE(ntp_varnames)

!++SAN print all NTP varnames:
IF (PrintStatus >= PrStatus_Oper) THEN
  WRITE(umMessage,'(A33,I3)') 'SAN; BOX_MODEL n_ntp_required: ', n_ntp_required
  CALL umPrint(umMessage, src='box_model')
  WRITE(umMessage,*) 'SAN; BOX_MODEL ntp_varnames: '
  CALL umPrint(umMessage, src='box_model')
  DO i = 1, n_ntp_required
    WRITE(umMessage,'(A3,I3,A2,A10)') ' * ', i, '. ', ntp_varnames(i)
    CALL umPrint(umMessage, src='box_model')
  END DO
END IF

! Allocate arrays to fill with UKCA data
n_ntp = n_ntp_required
IF (.NOT. ALLOCATED(ntp_data)) &
    ALLOCATE(ntp_data(row_length, rows, model_levels, n_ntp))

! Initialise NTP array to small values
ntp_data(:,:,:,:) = tracer_nullval

! Map data for relevant tracers input to ntp array
DO i = 1, n_ntp_required
  ! For NO2, BrO and HCl, these need to be mapped from their initial tracer values
  IF (TRIM(ntp_varnames(i)) == 'NO2' .OR. TRIM(ntp_varnames(i)) == 'BrO' .OR. &
      TRIM(ntp_varnames(i)) == 'HCl') THEN
    DO j=1, n_tracer_required
      IF (TRIM(ntp_varnames(i)) == TRIM(tracer_varnames(j))) THEN
        ntp_data(:,:,:,i) = tracer_ukca(:,:,:,j)
        EXIT
      ELSE 
        CONTINUE
      END IF
    END DO
  END IF
END DO

! ----------------------------------------------------------------------
!  1.3 Set up  outputting of tracer and ntp arrays
!

!      By default file will be written to '~/cylc-run/u-<suiteID>/work/1/ukca/'
!s OPEN(81,FILE='tracer_out.csv')
!++SAN UPDATE - output now defined in namelist
!++ Also, only output on thread 0
IF (mype == 0) THEN
    WRITE(umMessage,*) 'BOX_MODEL: Writing to ',tracer_out_filename
    CALL umPrint(umMessage,src='box_model')
    OPEN(tracer_out_unit,FILE=tracer_out_filename)
    WRITE(out_format,  '(A13, I0, A13)') '("# COL ",I4,', n_tracer_required, '(",",1X,I15))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_out_unit, TRIM(ADJUSTL(out_format))) 1,(i+1, i=1,n_tracer_required)
    WRITE(out_format, '(A9, I0, A13)') '("# ",A8,', n_tracer_required, '(",",1X,A15))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_out_unit, TRIM(ADJUSTL(out_format))) 'Timestep', &
         (TRIM(ADJUSTL(tracer_varnames(i))), i=1,n_tracer_required)
    ! Write initial tracer mass mixing ratios
    CALL umPrint('SAN BOX_MODEL: Writing tracers to output file',src='box_model')
    WRITE(out_format, '(A5, I0, A18)') '(I10,', n_tracer_required, '(",",1X,ES15.6E3))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_out_unit, TRIM(ADJUSTL(out_format))) 0, tracer_ukca

    WRITE(umMessage,*) 'BOX_MODEL: Writing to ',ntp_out_filename
    CALL umPrint(umMessage,src='box_model')
    OPEN(ntp_out_unit,FILE=ntp_out_filename)
    WRITE(out_format,  '(A13, I0, A13)') '("# COL ",I4,', n_ntp_required, '(",",1X,I15))'
    CALL umPrint(out_format, src='box_model')
    WRITE(ntp_out_unit, TRIM(ADJUSTL(out_format))) 1,(i+1, i=1,n_ntp_required)
    WRITE(out_format, '(A9, I0, A13)') '("# ",A8,', n_ntp_required, '(",",1X,A15))'
    CALL umPrint(out_format, src='box_model')
    WRITE(ntp_out_unit, TRIM(ADJUSTL(out_format))) 'Timestep', &
         (TRIM(ADJUSTL(ntp_varnames(i))), i=1,n_ntp_required)
    ! Write initial ntp mass mixing ratios
    CALL umPrint('SAN BOX_MODEL: Writing tracers to output file',src='box_model')
    WRITE(out_format, '(A5, I0, A18)') '(I10,', n_ntp_required, '(",",1X,ES15.6E3))'
    CALL umPrint(out_format, src='box_model')
    WRITE(ntp_out_unit, TRIM(ADJUSTL(out_format))) 0, ntp_data
END IF
! define hex filename
ppos = scan(trim(adjustl(tracer_out_filename)),".", BACK= .true.)
if ( ppos > 0 ) tracer_hex_filename = tracer_out_filename(1:ppos)//"hex"
! make integer array of size of tracer array
IF (.NOT. ALLOCATED(int_tracers)) ALLOCATE(int_tracers(1:n_tracer_required))
int_tracers = 0
int_tracers = TRANSFER(tracer_ukca, int_tracers)
IF (mype == 0) THEN
    WRITE(umMessage,*) 'BOX_MODEL: Writing to ',tracer_hex_filename
    CALL umPrint(umMessage,src='box_model')
    OPEN(tracer_hex_unit,FILE=tracer_hex_filename)
    WRITE(out_format,  '(A13, I0, A13)') '("# COL ",I4,', n_tracer_required, '(",",1X,I16))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_hex_unit, TRIM(ADJUSTL(out_format))) 1,(i+1, i=1,n_tracer_required)
    WRITE(out_format, '(A9, I0, A13)') '("# ",A8,', n_tracer_required, '(",",1X,A16))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_hex_unit, TRIM(ADJUSTL(out_format))) 'Timestep', &
         (TRIM(ADJUSTL(tracer_varnames(i))), i=1,n_tracer_required)
    ! Write initial tracer mass mixing ratios
    CALL umPrint('SAN BOX_MODEL: Writing tracers to output file',src='box_model')
    WRITE(out_format, '(A5, I0, A18)') '(I10,', n_tracer_required, '(",",1X,Z0))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_hex_unit, TRIM(ADJUSTL(out_format))) 0, int_tracers
END IF
! open the chemical diagnostic files - can't write header just yet
WRITE(umMessage,*) 'BOX_MODEL: Opening ',flux_out_filename
CALL umPrint(umMessage,src='box_model')
OPEN(flux_out_unit,FILE=flux_out_filename)
WRITE(umMessage,*) 'BOX_MODEL: Opening ',rate_out_filename
CALL umPrint(umMessage,src='box_model')
OPEN(rate_out_unit,FILE=rate_out_filename)


CALL umPrint('SAN BOX_MODEL: End of initialisation',src='box_model')
CALL timer('INITIAL ',6)


! ----------------------------------------------------------------------
!  2. Check for nothing-to-do
!
WRITE(umMessage,'(A,L3,3(1X,I0),1X,L3)') 'SAN BOX_MODEL: Checking lexit ',     &
     lexitNOW, atmos_im,                                                       &
     stepim(atmos_im), target_end_stepim(atmos_im),                            &
     (stepim(atmos_im) >= target_end_stepim(atmos_im))
CALL umPrint(umMessage,src='box_model')

CALL exitchek(lexitNOW)

! ----------------------------------------------------------------------
!  3. Start group of timesteps

IF (.NOT. lexitnow) THEN

  WRITE(umMessage,'(A,F10.2,A)') 'Model running with timestep ',               &
                                 secs_per_stepim(atmos_im),' seconds'
  CALL umPrint(umMessage,src='box_model')

END IF
! ----------------------------------------------------------------------
!  3. Start group of timesteps
!
CALL umPrint('SAN BOX_MODEL: Starting timestep loop',src='box_model')

DO WHILE (.NOT. lexitnow) ! Keep looping until the exit flag is set

  ! ----------------------------------------------------------------------
  !  3.1. Start main timestep loop
  !
  !  3.1.1 Increment model time ..
  CALL umPrint('SAN BOX_MODEL: Calling incrtime',src='box_model')
  CALL incrtime (internal_model)

  ! Print out all tracers again to double check has worked properly
  DO i = 1, n_tracer_required
    ! Check whether tracer has been matched
    WRITE(umMessage,'(A3,I3,A2,A10,A3,E12.4)') ' * ',                          &
      i, '. ', tracer_varnames(i), ' = ', tracer_ukca(1,1,1,i) 
    CALL umPrint(umMessage, src='box_model')
  END DO

  !++SAN map the non-varying environment fields to what can be read by UKCA
  !++    This is very long and cumbersome code, so may be moved into a separate
  !      routine in time. It needs to be done on every timestep.
  !      At some point, should be possible to add some time-varying fields,
  !      e.g. by reading the location at each time from a separate file
  !      in order to have a trajectory following setup
  IF (ALLOCATED(land_sea_mask)) THEN
    CALL umPrint('SAN: Mapping land_sea_mask', src='box_model')
    CALL ukca_set_environment('land_sea_mask', land_sea_mask,                  &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(tstar)) THEN
    CALL umPrint('SAN: Mapping tstar', src='box_model')
    CALL ukca_set_environment('tstar', tstar,                    &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(zbl)) THEN
    CALL umPrint('SAN: Mapping zbl', src='box_model')
    CALL ukca_set_environment('zbl', zbl,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF


  IF (ALLOCATED(pstar)) THEN
    CALL umPrint('SAN: Mapping pstar', src='box_model')
    CALL ukca_set_environment('pstar', pstar,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(surf_albedo)) THEN
    CALL umPrint('SAN: Mapping surf_albedo', src='box_model')
    CALL ukca_set_environment('surf_albedo', surf_albedo,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(conv_cloud_lwp)) THEN
    CALL umPrint('SAN: Mapping conv_cloud_lwp', src='box_model')
    CALL ukca_set_environment('conv_cloud_lwp', conv_cloud_lwp,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(rough_length)) THEN
    CALL umPrint('SAN: Mapping rough_length', src='box_model')
    CALL ukca_set_environment('rough_length', rough_length,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(seaice_frac)) THEN
    CALL umPrint('SAN: Mapping seaice_frac', src='box_model')
    CALL ukca_set_environment('seaice_frac', seaice_frac,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(zhsc)) THEN
    CALL umPrint('SAN: Mapping zhsc', src='box_model')
    CALL ukca_set_environment('zhsc', zhsc,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(u_scalar_10m)) THEN
    CALL umPrint('SAN: Mapping u_scalar_10m', src='box_model')
    CALL ukca_set_environment('u_scalar_10m', u_scalar_10m,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(u_s)) THEN
    CALL umPrint('SAN: Mapping u_s', src='box_model')
    CALL ukca_set_environment('u_s', u_s,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(kent)) THEN
    CALL umPrint('SAN: Mapping kent', src='box_model')
    CALL ukca_set_environment('kent', kent,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(kent_dsc)) THEN
    CALL umPrint('SAN: Mapping kent_dsc', src='box_model')
    CALL ukca_set_environment('kent_dsc', kent_dsc,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(theta)) THEN
    CALL umPrint('SAN: Mapping theta', src='box_model')
    CALL ukca_set_environment('theta', theta,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  !++SAN map q
  IF (ALLOCATED(q)) THEN
    CALL umPrint('SAN: Mapping q', src='box_model')
    CALL ukca_set_environment('q', q,                                            &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(qcf)) THEN
    CALL umPrint('SAN: Mapping qcf', src='box_model')
    CALL ukca_set_environment('qcf', qcf,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(conv_cloud_amount)) THEN
    CALL umPrint('SAN: Mapping conv_cloud_amount', src='box_model')
    CALL ukca_set_environment('conv_cloud_amount', conv_cloud_amount,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(rho_r2)) THEN
    CALL umPrint('SAN: Mapping rho_r2', src='box_model')
    CALL ukca_set_environment('rho_r2', rho_r2,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(qcl)) THEN
    CALL umPrint('SAN: Mapping qcl', src='box_model')
    CALL ukca_set_environment('qcl', qcl,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(exner_rho_levels)) THEN
    CALL umPrint('SAN: Mapping exner_rho_levels', src='box_model')
    CALL ukca_set_environment('exner_rho_levels', exner_rho_levels,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(area_cloud_fraction)) THEN
    CALL umPrint('SAN: Mapping area_cloud_fraction', src='box_model')
    CALL ukca_set_environment('area_cloud_fraction', area_cloud_fraction,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(cloud_frac)) THEN
    CALL umPrint('SAN: Mapping cloud_frac', src='box_model')
    CALL ukca_set_environment('cloud_frac', cloud_frac,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(cloud_liq_frac)) THEN
    CALL umPrint('SAN: Mapping cloud_liq_frac', src='box_model')
    CALL ukca_set_environment('cloud_liq_frac', cloud_liq_frac,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(exner_theta_levels)) THEN
    CALL umPrint('SAN: Mapping exner_theta_levels', src='box_model')
    CALL ukca_set_environment('exner_theta_levels', exner_theta_levels,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(p_rho_levels)) THEN
    CALL umPrint('SAN: Mapping p_rho_levels', src='box_model')
    CALL ukca_set_environment('p_rho_levels', p_rho_levels,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(p_theta_levels)) THEN
    CALL umPrint('SAN: Mapping p_theta_levels', src='box_model')
    CALL ukca_set_environment('p_theta_levels', p_theta_levels,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(rhokh_mix)) THEN
    CALL umPrint('SAN: Mapping rhokh_mix', src='box_model')
    CALL ukca_set_environment('rhokh_mix', rhokh_mix,                          &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(dtrdz_charney_grid)) THEN
    CALL umPrint('SAN: Mapping dtrdz_charney_grid', src='box_model')
    CALL ukca_set_environment('dtrdz_charney_grid', dtrdz_charney_grid,        &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(rhokh_rdz)) THEN
    CALL umPrint('SAN: Mapping rhokh_rdz', src='box_model')
    CALL ukca_set_environment('rhokh_rdz', rhokh_rdz,                          &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(dtrdz)) THEN
    CALL umPrint('SAN: Mapping dtrdz', src='box_model')
    CALL ukca_set_environment('dtrdz', dtrdz,                                  &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(we_lim)) THEN
    CALL umPrint('SAN: Mapping we_lim', src='box_model')
    CALL ukca_set_environment('we_lim', we_lim,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(t_frac)) THEN
    CALL umPrint('SAN: Mapping t_frac', src='box_model')
    CALL ukca_set_environment('t_frac', t_frac,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(zrzi)) THEN
    CALL umPrint('SAN: Mapping zrzi', src='box_model')
    CALL ukca_set_environment('zrzi', zrzi,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(we_lim_dsc)) THEN
    CALL umPrint('SAN: Mapping we_lim_dsc', src='box_model')
    CALL ukca_set_environment('we_lim_dsc', we_lim_dsc,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(t_frac_dsc)) THEN
    CALL umPrint('SAN: Mapping t_frac_dsc', src='box_model')
    CALL ukca_set_environment('t_frac_dsc', t_frac_dsc,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(zrzi_dsc)) THEN
    CALL umPrint('SAN: Mapping zrzi_dsc', src='box_model')
    CALL ukca_set_environment('zrzi_dsc', zrzi_dsc,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(ls_rain3d)) THEN
    CALL umPrint('SAN: Mapping ls_rain3d', src='box_model')
    CALL ukca_set_environment('ls_rain3d', ls_rain3d,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(ls_snow3d)) THEN
    CALL umPrint('SAN: Mapping ls_snow3d', src='box_model')
    CALL ukca_set_environment('ls_snow3d', ls_snow3d,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(autoconv)) THEN
    CALL umPrint('SAN: Mapping autoconv', src='box_model')
    CALL ukca_set_environment('autoconv', autoconv,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(accretion)) THEN
    CALL umPrint('SAN: Mapping accretion', src='box_model')
    CALL ukca_set_environment('accretion', accretion,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(pv_on_theta_mlevs)) THEN
    CALL umPrint('SAN: Mapping pv_on_theta_mlevs', src='box_model')
    CALL ukca_set_environment('pv_on_theta_mlevs', pv_on_theta_mlevs,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(conv_cloud_base)) THEN
    CALL umPrint('SAN: Mapping conv_cloud_base', src='box_model')
    CALL ukca_set_environment('conv_cloud_base', conv_cloud_base,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(conv_cloud_top)) THEN
    CALL umPrint('SAN: Mapping conv_cloud_top', src='box_model')
    CALL ukca_set_environment('conv_cloud_top', conv_cloud_top,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(conv_rain3d)) THEN
    CALL umPrint('SAN: Mapping conv_rain3d', src='box_model')
    CALL ukca_set_environment('conv_rain3d', conv_rain3d,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(conv_snow3d)) THEN
    CALL umPrint('SAN: Mapping conv_snow3d', src='box_model')
    CALL ukca_set_environment('conv_snow3d', conv_snow3d,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(so4_sa_clim)) THEN
    CALL umPrint('SAN: Mapping so4_sa_clim', src='box_model')
    CALL ukca_set_environment('so4_sa_clim', so4_sa_clim,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(dms_sea_conc)) THEN
    CALL umPrint('SAN: Mapping dms_sea_conc', src='box_model')
    CALL ukca_set_environment('dms_sea_conc', dms_sea_conc,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(rim_cry)) THEN
    CALL umPrint('SAN: Mapping rim_cry', src='box_model')
    CALL ukca_set_environment('rim_cry', rim_cry,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(rim_agg)) THEN
    CALL umPrint('SAN: Mapping rim_agg', src='box_model')
    CALL ukca_set_environment('rim_agg', rim_agg,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(vertvel)) THEN
    CALL umPrint('SAN: Mapping vertvel', src='box_model')
    CALL ukca_set_environment('vertvel', vertvel,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (ALLOCATED(bl_tke)) THEN
    CALL umPrint('SAN: Mapping bl_tke', src='box_model')
    CALL ukca_set_environment('bl_tke', bl_tke,                    &
                              errcode, error_message=ukca_errmsg,                &
                              error_routine=ukca_errproc)
    IF (errcode > 0) THEN
      cmessage = ukca_errmsg
      CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
    END IF
  END IF

  IF (l_ukca_intdd) THEN
    IF (ALLOCATED(soil_moisture_layer1)) THEN
      CALL umPrint('SAN: Mapping soil_moisture_layer1', src='box_model')
      CALL ukca_set_environment('soil_moisture_layer1', soil_moisture_layer1,      &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF ! l_ukca_intdd

    IF (ALLOCATED(l_tile_active)) THEN
      CALL umPrint('SAN: Mapping l_tile_active', src='box_model')
      CALL ukca_set_environment('l_tile_active', l_tile_active,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(frac_types)) THEN
      CALL umPrint('SAN: Mapping frac_types', src='box_model')
      CALL ukca_set_environment('frac_types', frac_types,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(laift_lp)) THEN
      CALL umPrint('SAN: Mapping laift_lp', src='box_model')
      CALL ukca_set_environment('laift_lp', laift_lp,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(canhtft_lp)) THEN
      CALL umPrint('SAN: Mapping canhtft_lp', src='box_model')
      CALL ukca_set_environment('canhtft_lp', canhtft_lp,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(tstar_tile)) THEN
      CALL umPrint('SAN: Mapping tstar_tile', src='box_model')
      CALL ukca_set_environment('tstar_tile', tstar_tile,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(z0tile_lp)) THEN
      CALL umPrint('SAN: Mapping z0tile_lp', src='box_model')
      CALL ukca_set_environment('z0tile_lp', z0tile_lp,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(surf_hf)) THEN
      CALL umPrint('SAN: Mapping surf_hf', src='box_model')
      CALL ukca_set_environment('surf_hf', surf_hf,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

    IF (ALLOCATED(stcon)) THEN
      CALL umPrint('SAN: Mapping stcon', src='box_model')
      CALL ukca_set_environment('stcon', stcon,                    &
                                errcode, error_message=ukca_errmsg,                &
                                error_routine=ukca_errproc)
      IF (errcode > 0) THEN
        cmessage = ukca_errmsg
        CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
      END IF
    END IF

  END IF
  !--SAN

  ! Code to calculate mixing ratio of well-mixed greenhouse gases
  ! if scenarios for their time variation have been prescribed.
  !++SAN - comment out for now, will need to decide how this is done in box model...
  !s CALL gas_calc_all()
  !++SAN set trace gas mixing ratios.
  WRITE(umMessage,*) 'SAN; BOX_MODEL: Setting Trace gas mixing ratios'
  CALL umPrint(umMessage,src='box_model')
  WRITE(umMessage,*) 'ch4mmr: ', ch4mmr, '; co2_mmr: ', co2_mmr
  CALL umPrint(umMessage,src='box_model')
  WRITE(umMessage,*) 'n2ommr: ', n2ommr, '; o2mmr: ', o2mmr
  CALL umPrint(umMessage,src='box_model')
  WRITE(umMessage,*) 'c11mmr: ', c11mmr, '; c12mmr: ', c12mmr
  CALL umPrint(umMessage,src='box_model')
  WRITE(umMessage,*) 'hcfc22mmr: ', hcfc22mmr, '; hfc125mmr: ', hfc125mmr
  CALL umPrint(umMessage,src='box_model')
  WRITE(umMessage,*) 'hfc134ammr: ', hfc134ammr
  CALL umPrint(umMessage,src='box_model')

  CALL umPrint(umMessage, src='box_model')
  !      These use the same calls as done in atmos_physics1.F90 in the UM
  !  Write trace gas mixing ratios to module for use in UKCA. The value can be
  !  either that which was passed to this routine or the value calculated from
  !  the time interpolation, depending on whether L_CLMCHFCG is true or false.
  IF (L_ukca) THEN
    CALL ukca_set_trace_gas_mixratio(                                            &
          ch4mmr, co2_mmr, n2ommr, o2mmr,                                        &
          c11mmr, c12mmr, c113mmr, c114mmr,                                      &
          hcfc22mmr, hfc125mmr, hfc134ammr)
  END IF


  !  3.1.2 .. set timestep control switches
  CALL umPrint('SAN BOX_MODEL: Calling settsctl',src='box_model')
  CALL settsctl (                                                              &
           internal_model,meanlev,icode,cmessage)
  IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)

  ! Send an EOT message to IOS
  ! He may or may not act on this to purge outstanding items
  CALL umPrint('SAN BOX_MODEL: Calling io_timestep',src='box_model')
  CALL io_timestep(stepim(atmos_im))

  !        Integrate atmosphere (box model) or ocean by 1 timestep
  IF (internal_model == atmos_im) THEN

    !++SAN - this will need to be replaced with the box_ukca_mod routine
    !s CALL atm_step_4a (                                                         &
    !++SAN Call box_ukca_mod directly
    !++ ntp_data now made in this routine and passed to box_ukca
    !++ tracer_ukca and ntp_data will need to be made in this routine
    CALL umPrint('SAN BOX_MODEL: Calling box_ukca',src='box_model')
    CALL box_ukca(tracer_ukca, q, ntp_data, tracer_varnames, n_tracer_required,&
                  environ_varnames, p_theta_levels, theta, exner_theta_levels, &
                  p_rho_levels, qcl,qcf)

  END IF        ! internal_model = atmos_im

  ! Print out all tracers again to double check has worked properly
  DO i = 1, n_tracer_required
    ! Check whether tracer has been matched
    WRITE(umMessage,'(A3,I3,A2,A10,A3,E12.4)') ' * ',                          &
      i, '. ', tracer_varnames(i), ' = ', tracer_ukca(1,1,1,i) 
    CALL umPrint(umMessage, src='box_model')
  END DO
  !--SAN

  !++SAN Add some basic output code here
  IF (mype == 0) THEN ! Only on thread 0
    CALL umPrint('SAN BOX_MODEL: Writing tracers to output file',src='box_model')
    WRITE(out_format, '(A5, I0, A18)') '(I10,', n_tracer_required, '(",",1X,ES15.6E3))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_out_unit, TRIM(ADJUSTL(out_format))) stepim(atmos_im), tracer_ukca
    CALL umPrint('SAN BOX_MODEL: Finished writing tracers to output file',src='box_model')

    ! ntp output
    CALL umPrint('SAN BOX_MODEL: Writing ntp to output file',src='box_model')
    WRITE(out_format, '(A5, I0, A18)') '(I10,', n_ntp_required, '(",",1X,ES15.6E3))'
    CALL umPrint(out_format, src='box_model')
    WRITE(ntp_out_unit, TRIM(ADJUSTL(out_format))) stepim(atmos_im), ntp_data
    CALL umPrint('SAN BOX_MODEL: Finished writing ntp to output file',src='box_model')

    ! hex output
    int_tracers = TRANSFER(tracer_ukca, int_tracers)
    WRITE(out_format, '(A5, I0, A18)') '(I10,', n_tracer_required, '(",",1X,Z0))'
    CALL umPrint(out_format, src='box_model')
    WRITE(tracer_hex_unit, TRIM(ADJUSTL(out_format))) stepim(atmos_im), int_tracers
    CALL umPrint('SAN BOX_MODEL: Finished writing tracers to output file',src='box_model')
  END IF

  ! output ASAD chemical diagnostic information here
  CALL box_output_chem_diags(stepim(atmos_im), flux_out_unit, rate_out_unit)
  
  !---------------------------------------------------------------------
  ! Ensure atmos coupling data in prognostic areas of D1 are up to date.
  ! The logic here is that cpl_update_step will be TRUE on the timestep
  ! BEFORE coupling is due to take place - SAN DELETE

  !  3.1.5 If printed output time, call print control routine
  IF (lprint) THEN

    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,*) RoutineName,':Warning, Printing of climate global '   &
                ,'and zonal diagnostics no longer supported'
      CALL umPrint(umMessage,src='box_model')
    END IF  ! PrintStatus test

  END IF

  !  3.1.9 If exit check time, check for immediate exit
  IF (lexit) THEN

    CALL exitchek(lexitNOW)

    IF (lexitNOW) THEN
      IF (.NOT. ldump) THEN
        CYCLE ! Restart the loop now (this will exit since lexitnow is true)
      END IF
    END IF
  END IF
  !
  !       End main timestep loop
  ! ----------------------------------------------------------------------
END DO

!++SAN Close main output files
IF (mype == 0) THEN
  CALL umPrint('SAN BOX_MODEL: Closing main output file',src='box_model')
  CLOSE(tracer_out_unit, IOSTAT=icode)
  IF (icode /= 0) THEN
    WRITE(umMessage,*) 'box_model: error closing output stream ', tracer_out_unit, icode
    CALL umPrint(umMessage,src='box_model')
  ELSE
    WRITE(umMessage,*) 'box_model: successfully closed output stream ', tracer_out_unit, icode
    CALL umPrint(umMessage,src='box_model')
  END IF

  CALL umPrint('SAN BOX_MODEL: Closing ntp output file',src='box_model')
  CLOSE(ntp_out_unit, IOSTAT=icode)
  IF (icode /= 0) THEN
    WRITE(umMessage,*) 'box_model: error closing output stream ', ntp_out_unit, icode
    CALL umPrint(umMessage,src='box_model')
  ELSE
    WRITE(umMessage,*) 'box_model: successfully closed output stream ', ntp_out_unit, icode
    CALL umPrint(umMessage,src='box_model')
  END IF

  CLOSE(tracer_hex_unit, IOSTAT=icode)
  IF (icode /= 0) THEN
    WRITE(umMessage,*) 'box_model: error closing hex output stream ', tracer_hex_unit, icode
    CALL umPrint(umMessage,src='box_model')
  ELSE
    WRITE(umMessage,*) 'box_model: successfully closed hex output stream ', tracer_hex_unit, icode
    CALL umPrint(umMessage,src='box_model')
  END IF
END IF
! flux file
CALL umPrint('SAN BOX_MODEL: Closing flux output file',src='box_model')
CLOSE(flux_out_unit, IOSTAT=icode)
IF (icode /= 0) THEN
   WRITE(umMessage,*) 'box_model: error closing output stream ', flux_out_unit, icode
   CALL umPrint(umMessage,src='box_model')
ELSE
   WRITE(umMessage,*) 'box_model: successfully closed output stream ', flux_out_unit, icode
   CALL umPrint(umMessage,src='box_model')
END IF
! rate file
CALL umPrint('SAN BOX_MODEL: Closing rate output file',src='box_model')
CLOSE(rate_out_unit, IOSTAT=icode)
IF (icode /= 0) THEN
   WRITE(umMessage,*) 'box_model: error closing output stream ', rate_out_unit, icode
   CALL umPrint(umMessage,src='box_model')
ELSE
   WRITE(umMessage,*) 'box_model: successfully closed output stream ', rate_out_unit, icode
   CALL umPrint(umMessage,src='box_model')
END IF
!s CLOSE(82)

IF (print_runtime_info) THEN
  !s IF ( L_print_pe .OR. mype == 0 ) THEN
  IF ( mype == 0 ) THEN
    CALL umPrint( '',src='box_model')
    time_end_run = get_wallclock_time()
    WRITE(umMessage,'(A,A,F10.3,A)')                                           &
      'box_model: Info: Exiting last box model timestep ',                   &
      ' at time=',time_end_run - Start_time,' seconds'
    CALL umPrint(umMessage,src='box_model', level=-PrintStatus)
    CALL umPrint( '',src='box_model')
  END IF ! mype == 0
END IF  ! print_runtime_info

! ----------------------------------------------------------------------
!  4. Exit processing: Output error messages and perform tidy-up
!
CALL umPrint('SAN BOX_MODEL: exit processing',src='box_model')

!  4.1 Exit processing: If abnormal completion, output error message
iabort = icode
IF (icode /= 0) THEN

  CALL Ereport(RoutineName,icode,Cmessage)

END IF

!++SAN deallocate arrays which are now made in this routine
CALL umPrint('SAN BOX_MODEL: deallocating NTP and tracer arrays',src='box_model')
IF (ALLOCATED(ntp_data)) DEALLOCATE(ntp_data)
IF (ALLOCATED(tracer_ukca)) DEALLOCATE(tracer_ukca)

!++SAN Deallocate trigs
CALL umPrint('SAN BOX_MODEL: deallocating trigs',src='box_model')
IF (ALLOCATED(true_latitude)) DEALLOCATE(true_latitude)
IF (ALLOCATED(true_longitude)) DEALLOCATE(true_longitude)
IF (ALLOCATED(FV_cos_theta_latitude)) DEALLOCATE(FV_cos_theta_latitude)

!++SAN deallocate arrays which are now made in this routine
CALL umPrint('SAN BOX_MODEL: deallocating environmental fields',src='box_model')
IF (ALLOCATED(land_sea_mask)) DEALLOCATE(land_sea_mask)
IF (ALLOCATED(conv_cloud_lwp)) DEALLOCATE(conv_cloud_lwp)
IF (ALLOCATED(tstar)) DEALLOCATE(tstar)
IF (ALLOCATED(zbl)) DEALLOCATE(zbl)
IF (ALLOCATED(rough_length)) DEALLOCATE(rough_length)
IF (ALLOCATED(seaice_frac)) DEALLOCATE(seaice_frac)
IF (ALLOCATED(pstar)) DEALLOCATE(pstar)
IF (ALLOCATED(zhsc)) DEALLOCATE(zhsc)
IF (ALLOCATED(u_scalar_10m)) DEALLOCATE(u_scalar_10m)
IF (ALLOCATED(u_s)) DEALLOCATE(u_s)
IF (ALLOCATED(kent)) DEALLOCATE(kent)
IF (ALLOCATED(kent_dsc)) DEALLOCATE(kent_dsc)
IF (ALLOCATED(theta)) DEALLOCATE(theta)
IF (ALLOCATED(q)) DEALLOCATE(q)
IF (ALLOCATED(qcf)) DEALLOCATE(qcf)
IF (ALLOCATED(conv_cloud_amount)) DEALLOCATE(conv_cloud_amount)
IF (ALLOCATED(rho_r2)) DEALLOCATE(rho_r2)
IF (ALLOCATED(qcl)) DEALLOCATE(qcl)
IF (ALLOCATED(exner_rho_levels)) DEALLOCATE(exner_rho_levels)
IF (ALLOCATED(area_cloud_fraction)) DEALLOCATE(area_cloud_fraction)
IF (ALLOCATED(cloud_frac)) DEALLOCATE(cloud_frac)
IF (ALLOCATED(cloud_liq_frac)) DEALLOCATE(cloud_liq_frac)
IF (ALLOCATED(exner_theta_levels)) DEALLOCATE(exner_theta_levels)
IF (ALLOCATED(r_theta_levels)) DEALLOCATE(r_theta_levels)
IF (ALLOCATED(r_rho_levels)) DEALLOCATE(r_rho_levels)
IF (ALLOCATED(eta_theta_levels)) DEALLOCATE(eta_theta_levels)
IF (ALLOCATED(p_rho_levels)) DEALLOCATE(p_rho_levels)
IF (ALLOCATED(p_theta_levels)) DEALLOCATE(p_theta_levels)
IF (ALLOCATED(t_theta_levels)) DEALLOCATE(t_theta_levels)
IF (ALLOCATED(rhokh_mix)) DEALLOCATE(rhokh_mix)
IF (ALLOCATED(dtrdz_charney_grid)) DEALLOCATE(dtrdz_charney_grid)
IF (ALLOCATED(rhokh_rdz)) DEALLOCATE(rhokh_rdz)
IF (ALLOCATED(dtrdz)) DEALLOCATE(dtrdz)
IF (ALLOCATED(we_lim)) DEALLOCATE(we_lim)
IF (ALLOCATED(t_frac)) DEALLOCATE(t_frac)
IF (ALLOCATED(zrzi)) DEALLOCATE(zrzi)
IF (ALLOCATED(we_lim_dsc)) DEALLOCATE(we_lim_dsc)
IF (ALLOCATED(t_frac_dsc)) DEALLOCATE(t_frac_dsc)
IF (ALLOCATED(zrzi_dsc)) DEALLOCATE(zrzi_dsc)
IF (ALLOCATED(ls_rain3d)) DEALLOCATE(ls_rain3d)
IF (ALLOCATED(ls_snow3d)) DEALLOCATE(ls_snow3d)
IF (ALLOCATED(autoconv)) DEALLOCATE(autoconv)
IF (ALLOCATED(accretion)) DEALLOCATE(accretion)
IF (ALLOCATED(pv_on_theta_mlevs)) DEALLOCATE(pv_on_theta_mlevs)
IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)

IF (ALLOCATED(conv_cloud_base)) DEALLOCATE(conv_cloud_base)
IF (ALLOCATED(conv_cloud_top)) DEALLOCATE(conv_cloud_top)
IF (ALLOCATED(conv_rain3d)) DEALLOCATE(conv_rain3d)
IF (ALLOCATED(conv_snow3d)) DEALLOCATE(conv_snow3d)
IF (ALLOCATED(so4_sa_clim)) DEALLOCATE(so4_sa_clim)

IF (l_ukca_intdd) THEN
  IF (ALLOCATED(soil_moisture_layer1)) DEALLOCATE(soil_moisture_layer1)
  IF (ALLOCATED(l_tile_active)) DEALLOCATE(l_tile_active)
  IF (ALLOCATED(frac_types)) DEALLOCATE(frac_types)
  IF (ALLOCATED(laift_lp)) DEALLOCATE(laift_lp)
  IF (ALLOCATED(canhtft_lp)) DEALLOCATE(canhtft_lp)
  IF (ALLOCATED(tstar_tile)) DEALLOCATE(tstar_tile)
  IF (ALLOCATED(z0tile_lp)) DEALLOCATE(z0tile_lp)
  IF (ALLOCATED(surf_hf)) DEALLOCATE(surf_hf)
  IF (ALLOCATED(stcon)) DEALLOCATE(stcon)
END IF

!++SAN Environment arrays needed for aerosol
IF (ALLOCATED(dms_sea_conc)) DEALLOCATE(dms_sea_conc)
IF (ALLOCATED(rim_cry)) DEALLOCATE(rim_cry)
IF (ALLOCATED(rim_agg)) DEALLOCATE(rim_agg)
IF (ALLOCATED(vertvel)) DEALLOCATE(vertvel)
IF (ALLOCATED(bl_tke)) DEALLOCATE(bl_tke)

IF (ALLOCATED(sf)) DEALLOCATE(sf)

! ----------------------------------------------------------------------
!  5. Complete Timer call and return
!
CALL umPrint('SAN BOX_MODEL: Complete timer',src='box_model')
icode=iabort
IF (ltimer) CALL timer('BOX_MODEL ',6)

CALL umPrint('SAN BOX_MODEL: Exiting BOX_MODEL',src='box_model')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE box_model
END MODULE box_model_mod
