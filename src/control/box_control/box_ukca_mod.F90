! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  UM module for coupling with UKCA containing a wrapper subroutine for
!  executing a single UKCA time step in box model mode.
!  Based on code originally in atmos_ukca_mod
!
! Method:
!
!  On the first call only, get the lists of required non-transported
!  prognostics (NTPs) and environmental driver fields from UKCA and set up
!  details of requirements for these fields from the D1 array.
!
!++SAN - won't be using D1 array with the box model. However, this
!        code will be operating the main I/O routines for the box model,
!        writing the tracer concentrations and chemical fluxes at each
!        timestep
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

MODULE box_ukca_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='BOX_UKCA_MOD'

PUBLIC :: box_ukca

! Oxidants Data structure
TYPE, PUBLIC :: ukca_ox_struct
  CHARACTER (LEN=256)  :: file_name    ! Name of source file
  CHARACTER (LEN=80)   :: var_name     ! Name of variable in file
  CHARACTER (LEN=10)   :: tracer_name  ! Emitted species

  INTEGER              :: update_freq  ! Update frequency (hours)
  INTEGER              :: update_type  ! 1 serial, 2 periodic, ...
  INTEGER              :: last_update  ! Num anc update intervals at last update

  LOGICAL              :: l_update     ! True if field updated in a given tstep

  ! Allocatable fields of derived types are not allowed in F95 (they are
  ! an extension of Standard F95). As a consequence we use pointers for
  ! variables that need to be allocated within the emissions structure.
  REAL, POINTER:: values (:,:,:) => NULL()  ! emission data
  REAL, POINTER:: diags  (:,:,:) => NULL()  ! emission diagnostics
END TYPE ukca_ox_struct

CONTAINS

! ----------------------------------------------------------------------

SUBROUTINE box_ukca(tracer_ukca_um, q, ntp_data, tracer_varnames,              &
  n_tracer_required, environ_varnames, p_theta_levels, theta,                  &
  exner_theta_levels, p_rho_levels, qcl, qcf)

!++SAN - deleted any modules no longer needed, especially STASH stuff
! vn12.2 - photolysis variables renamed:
!s                                   ukca_get_photol_varlist,                  &
!s                                  photol_varname_len,                        &
USE ukca_api_mod,           ONLY: ukca_step,                                   &
                                  ukca_get_ntp_varlist,                        &
                                  ukca_get_environment_varlist,                &
                                  ukca_get_photol_reaction_data,               &
                                  ukca_photol_varname_len,                     &
                                  ukca_set_environment,                        &
                                  ukca_maxlen_message,                         &
                                  ukca_maxlen_procname, ukca_maxlen_fieldname, &
                                  ukca_chem_strattrop, ukca_chem_strat,        &
                                  ukca_chem_cristrat
!==NLA
!USE ukca_phot2d,            ONLY: ntphot
USE ukca_option_mod,        ONLY: l_ukca_chem, l_ukca_mode,                    &
                                  l_ukca_dry_dep_so2wet,                       &
                                  l_ukca_offline, l_ukca_offline_be,           &
                                  l_ukca_use_background_aerosol,               &
                                  i_ukca_scenario, i_ukca_scenario_rcp,        &
                                  i_ukca_chem, max_z_for_offline_chem,         &
                                  i_ukca_photol, chem_timestep

USE photol_api_mod,         ONLY: photol_setup

!++SAN ntp_data now made and allocated in the parent box_model
!sUSE ukca_um_interf_mod,     ONLY: ukca_um_d1_initialise, putd1flds, &
!s                                  ukca_um_dealloc,                    &
!s                                  exner_theta_levels
!++SAN here, including enough stuff from ukca_setd1defs to make empty set of
!      D1 codes, in order to be able to run ukca_main1 with segmentation faults
USE ukca_setd1defs_mod,     ONLY: ukca_setd1defs
USE ukca_d1_defs,           ONLY: UkcaD1Codes, Nukca_D1items, n_ntp,           &
                                  n_mode_diags, n_strat_fluxdiags,             &
                                  item1_nitrate_diags, itemN_nitrate_diags
USE missing_data_mod,      ONLY: rmdi, imdi
!--SAN

USE atm_fields_bounds_mod,  ONLY: tdims, tdims_l, tdims_s
USE level_heights_mod,      ONLY: eta_theta_levels
USE dump_headers_mod,       ONLY: a_realhd, rh_z_top_theta
USE ukca_nc_emiss_mod,      ONLY: ukca_set_emissions_from_nc
USE ukca_read_offline_oxidants_mod, ONLY: ukca_set_oxidants_from_nc

!++SAN - for box model, new options to read oxidant values from namelist
USE nlstcall_box_namelist_mod, ONLY: ukca_box_nml
USE read_offline_oxidant_mmr_mod, ONLY: read_offline_oxidant_mmr

!s USE ukca_read_aerosol_mod,  ONLY: ukca_read_aerosol
!++SAN - ukca_um_photol_ctl added in vn12.2
!        Box model version used here
USE ukca_box_photol_ctl_mod, ONLY: photol_ctl
USE ukca_scenario_rcp_mod,  ONLY: ukca_scenario_rcp
USE nlsizes_namelist_mod,   ONLY: tr_ukca, model_levels,                       &
                                  row_length, rows, bl_levels,                 &
                                  global_row_length, global_rows,              &
                                  tr_levels, sm_levels, ntiles
USE dump_headers_mod,       ONLY: a_inthd
USE level_heights_mod,      ONLY: r_rho_levels, r_theta_levels
USE planet_constants_mod,   ONLY: planet_radius
USE cderived_mod,           ONLY: delta_lambda, delta_phi
USE trignometric_mod,       ONLY: true_latitude, true_longitude,               &
                                  FV_cos_theta_latitude
USE conversions_mod,        ONLY: recip_pi_over_180
!s USE land_soil_dimensions_mod, ONLY: land_points
USE submodel_mod,           ONLY: submodel_for_sm, atmos_sm, atmos_im
!s USE ukca_all_tracers_copy_mod, ONLY: ukca_make_tr_map,                         &
!s                                      ukca_tracers_copy_from_um,                &
!s                                      ukca_tracers_copy_to_um,                  &
!s                                      ukca_q_increment_diag
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE model_time_mod,         ONLY: i_year, i_month, i_day, i_hour, i_minute,    &
                                  i_second, i_day_number, stepim,              &
                                  secs_per_stepim, previous_time
USE solpos_mod,             ONLY: solpos
USE umPrintMgr,             ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE bmass_variable_hilem_mod, ONLY: calc_interf_z

USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook



IMPLICIT NONE

! Subroutine arguments

! UKCA Tracers coming from the UM code via subroutine interface
REAL, INTENT(IN OUT) :: tracer_ukca_um                                         &
              (1:row_length, 1:rows, 1:model_levels, 1:tr_ukca)
! Specific humidity coming from the UM code via subroutine interface
REAL, INTENT(IN OUT) :: q                                                   &
               (1:row_length, 1:rows, 1:model_levels)

!++SAN new array in/out for ntp_data
! Non-transported prognostics from the parent model. Dimensions: X,Y,Z,N
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data(:, :, :, :)
!++SAN Also inputting tracer varnames and length
CHARACTER(LEN=*), POINTER, INTENT(IN) :: tracer_varnames(:)
INTEGER, INTENT(IN)                   :: n_tracer_required
CHARACTER(LEN=ukca_maxlen_fieldname), POINTER, INTENT(IN) :: environ_varnames(:)

REAL, INTENT(IN)    :: p_theta_levels (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN)    :: theta (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN)    :: exner_theta_levels (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN)    :: p_rho_levels (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN)    :: qcl (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN)    :: qcf (1:row_length, 1:rows, 1:model_levels)
!--SAN

! Local variables

! List of tracers required for the UKCA configuration
!s CHARACTER(LEN=ukca_maxlen_fieldname), POINTER :: tracer_varnames(:)

! List of non-transported prognostics required for the UKCA configuration
!s CHARACTER(LEN=ukca_maxlen_fieldname), POINTER :: ntp_varnames(:)

! List of environmental driver fields required for the UKCA configuration
!s CHARACTER(LEN=ukca_maxlen_fieldname), SAVE, POINTER :: environ_varnames(:)

! List of photolysis file/variable names required for photolysis calcs
CHARACTER(LEN=10), SAVE, POINTER :: ratj_data(:,:)
CHARACTER(LEN=ukca_photol_varname_len), POINTER :: ratj_varnames(:)
INTEGER, SAVE :: jppj

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER :: current_time(7)

! Model time at previous time step
INTEGER :: i_year_previous
INTEGER :: i_hour_previous
INTEGER :: i_minute_previous
INTEGER :: i_second_previous
INTEGER :: i_day_number_previous
REAL :: secondssincemidnight

INTEGER :: m_atm_modl                   ! Sub model
INTEGER :: first_constant_r_rho_level   ! First rho level at which height
                                        ! is constant

REAL :: z_first_constant_r_rho_level    ! Height of first rho level at which
                                        ! height is constant
REAL, SAVE :: z_top_of_model            ! Top of model SAN - needed?

INTEGER :: n_nitrate_diags            ! Number of mode diagnostics
INTEGER :: n_environ_items            ! No. of environment fields to retrieve

! Environmental driver fields not held in D1

! Fields related to geographic location at points on the horizontal grid.
! Note that 2D fields are required to allow for the use of rotated grids where
! both true latitude and true longitude vary along both rows and columns.
REAL, SAVE, ALLOCATABLE :: latitude(:,:)     ! True latitude (degrees N)
REAL, SAVE, ALLOCATABLE :: longitude(:,:)    ! True longitude (degrees E, 0-360)
REAL, SAVE, ALLOCATABLE :: sin_latitude(:,:) ! SIN(true latitude)
REAL, SAVE, ALLOCATABLE :: cos_latitude(:,:) ! COS(true latitude)
REAL, SAVE, ALLOCATABLE :: tan_latitude(:,:) ! TAN(true latitude)
REAL, SAVE, ALLOCATABLE :: surf_area(:,:)    ! Gridbox surface area (m^2)

! Altitude of model grid-box interfaces
REAL, SAVE, ALLOCATABLE :: interf_z(:,:,:)

! Sulphate aerosol surface area density
!s REAL, SAVE, ALLOCATABLE :: so4_sa_clim(:,:,:)

! UKCA tracer array to be populated from tracer_ukca_um and q.
! This is a working array to hold the tracer fields in the order required by
! UKCA. It is allocated in the subroutine ukca_tracers_copy_from_um before
! the UKCA time step and deallocated in ukca_tracers_copy_to_um afterwards.
REAL, ALLOCATABLE :: ukca_tracer_data(:,:,:,:)

!++SAN loop variables for tracer mapping
INTEGER :: i
INTEGER :: j
INTEGER :: n
INTEGER :: l
INTEGER :: ll

REAL :: sindec    ! Sin solar declination
REAL :: eq_time   ! Equation of time

! Dummy variables for solpos call
REAL :: scs
REAL :: sindec_obs
REAL :: eqt_obs

LOGICAL :: l_planet_obs = .FALSE. ! True to calc. angles towards distant
                                  ! observer in solpos call (always false as
                                  ! this functionality is not used)

!++SAN dimension sizes for emission files, currently hardcoded
!s INTEGER :: glob_emiss_row_length = 96
!s INTEGER :: glob_emiss_rows       = 72
!s INTEGER :: glob_emiss_levels     = 38

!++SAN local variables for handling offline oxidants. Based on code in 
!      ukca_read_offline_oxidants
!++SAN - for now, hardcoding offline oxidant NetCDF option to false
LOGICAL :: l_read_offox_ncdf = .FALSE.
! Oxidant fieldnames
INTEGER, PARAMETER :: max_offl_ox = 5
CHARACTER(LEN=5), PARAMETER:: fldname_h2o2_offl = 'H2O2'
CHARACTER(LEN=5), PARAMETER:: fldname_ho2_offl  = 'HO2'
CHARACTER(LEN=5), PARAMETER:: fldname_no3_offl  = 'NO3'
CHARACTER(LEN=5), PARAMETER:: fldname_o3_offl   = 'O3'
CHARACTER(LEN=5), PARAMETER:: fldname_oh_offl   = 'OH'
CHARACTER(LEN=5), PARAMETER:: fldname_oxid(max_offl_ox) = (/fldname_h2o2_offl, &
                                fldname_ho2_offl, fldname_no3_offl,            &
                                fldname_o3_offl, fldname_oh_offl/)
! Number and names of 'active' oxidants
CHARACTER(LEN=10), SAVE :: oxid_varnames(max_offl_ox)
INTEGER, SAVE      :: n_offline
INTEGER            :: nlines
INTEGER            :: io

! Super array of offline oxidants
TYPE(ukca_ox_struct), SAVE, ALLOCATABLE :: offline_ox(:)

!++SAN scalar versions of offline oxidants for reading namelist input
REAL    :: o3_in         = rmdi
REAL    :: oh_in         = rmdi
REAL    :: ho2_in        = rmdi
REAL    :: h2o2_in       = rmdi
REAL    :: no3_in        = rmdi

REAL, ALLOCATABLE   :: tmpdat1 (:,:,:)         ! Temporal data arrays
!--SAN

LOGICAL, SAVE :: l_first_call = .TRUE.        ! True only on 1st call
LOGICAL, SAVE :: l_calc_plev_diags = .FALSE.  ! True if pressure level
                                              ! diagnostics are required

LOGICAL, SAVE :: l_req_sin_declination = .FALSE.
                                  ! True if UKCA requires sin declination
LOGICAL, SAVE :: l_req_equation_of_time = .FALSE.
                                  ! True if UKCA requires equation of time

LOGICAL :: l_chem_step            ! True if this is a chemistry timestep
                                  ! (required for photolysis)

! No such diagnostics in box model
INTEGER :: nmax_strat_fluxdiags = 0
INTEGER :: nmax_mode_diags = 0

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='BOX_UKCA'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!++SAN
CALL umPrint('SAN; BOX_UKCA: in box_ukca',src='box_ukca')

IF (l_first_call) THEN

  !++SAN
  CALL umPrint('SAN; BOX_UKCA: first call',src='box_ukca')


  ! Get list of NTP fields required by the current UKCA configuration
  !s CALL ukca_get_ntp_varlist(ntp_varnames, errcode, error_message=ukca_errmsg,  &
  !s                           error_routine=ukca_errproc)
!s  IF (errcode > 0) THEN
!s    cmessage = ukca_errmsg
!s    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
!s  END IF

!++SAN - environ_varnames now read in as argument
!s  NULLIFY(environ_varnames)
  !++SAN
  ! Get list of environment fields required by the current UKCA configuration
!s  CALL ukca_get_environment_varlist(environ_varnames, errcode,                 &
!s                                    error_message=ukca_errmsg,                 &
!s                                    error_routine=ukca_errproc)
!s  IF (errcode > 0) THEN
!s    cmessage = ukca_errmsg
!s    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
!s  END IF

  ! Define D1 definitions, emissions and diagnostics required
  ! (Also check compatibility of the UM's UKCA tracer array indices with
  ! radiative feedback if this is active)
  !++SAN no D1
  !s CALL ukca_setd1defs(row_length, rows, bl_levels, tr_levels,                  &
  !s                     land_points, sm_levels, ntiles,                          &
  !s                     ntp_varnames, environ_varnames, l_calc_plev_diags)
  !++SAN Allocate D1 codes and set to null values, needed for ukca_main1 to run
  !      This is dummy code based on code from ukca_setd1defs
  !      D1 is not actually used for output

  n_environ_items = SIZE(environ_varnames)
  n_nitrate_diags = itemN_nitrate_diags -  item1_nitrate_diags + 1
  
  WRITE(umMessage,*) 'SAN; n_environ_items = ', n_environ_items
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; n_ntp = ', n_ntp
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; nmax_mode_diags = ', nmax_mode_diags
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; n_nitrate_diags = ', n_nitrate_diags
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; nmax_strat_fluxdiags = ', nmax_strat_fluxdiags
  CALL umPrint(umMessage, src='box_ukca_mod')

  Nukca_D1items = n_ntp + n_environ_items + nmax_mode_diags +                    &
                  n_nitrate_diags + nmax_strat_fluxdiags
  WRITE(umMessage,*) 'SAN; Nukca_D1items = ', Nukca_D1items
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; Allocating UkcaD1Codes'
  CALL umPrint(umMessage, src='box_ukca_mod')

  IF (.NOT. ALLOCATED(UkcaD1Codes)) ALLOCATE(UkcaD1Codes(Nukca_D1items))

  WRITE(umMessage,*) 'SAN; Looping through = ', Nukca_D1items
  CALL umPrint(umMessage, src='box_ukca_mod')

  DO i=1, Nukca_D1items
    UkcaD1Codes(i)%section=imdi
    UkcaD1Codes(i)%item=imdi
    UkcaD1Codes(i)%n_levels=imdi
    UkcaD1Codes(i)%address=imdi
    UkcaD1Codes(i)%length=imdi
    UkcaD1Codes(i)%halo_type=imdi
    UkcaD1Codes(i)%grid_type=imdi
    UkcaD1Codes(i)%field_type=imdi
    UkcaD1Codes(i)%len_dim1=imdi
    UkcaD1Codes(i)%len_dim2=imdi
    UkcaD1Codes(i)%len_dim3=imdi
    UkcaD1Codes(i)%prognostic=.TRUE.
    UkcaD1Codes(i)%required=.FALSE.
    UkcaD1Codes(i)%fieldname=''
  END DO

  ! Check for environment fields not sourced from D1
  DO n = 1, SIZE(environ_varnames)
    WRITE(umMessage,*) 'SAN; CASE ', TRIM(environ_varnames(n))
    CALL umPrint(umMessage, src='box_ukca_mod')
    SELECT CASE (TRIM(environ_varnames(n)))
    CASE ('sin_declination')
      l_req_sin_declination = .TRUE.
    CASE ('equation_of_time')
      l_req_equation_of_time = .TRUE.
      ! Set location-related fields based on true latitude and longitude
    CASE ('latitude')
      ALLOCATE(latitude(row_length,rows))
      latitude(:,:) = true_latitude(:,:) * recip_pi_over_180
    CASE ('longitude')
      ALLOCATE(longitude(row_length,rows))
      longitude(:,:) = true_longitude(:,:) * recip_pi_over_180
    CASE ('sin_latitude')
      ALLOCATE(sin_latitude(row_length,rows))
      sin_latitude(:,:) = SIN(true_latitude(:,:))
    CASE ('cos_latitude')
      ALLOCATE(cos_latitude(row_length,rows))
      cos_latitude(:,:) = COS(true_latitude(:,:))
    CASE ('tan_latitude')
      ALLOCATE(tan_latitude(row_length,rows))
      tan_latitude(:,:) = SIN(true_latitude(:,:)) / COS(true_latitude(:,:))
      ! Allocate sulphate aerosol surface area density if required
    !++SAN comment out aerosol climatology for now
    !s CASE ('so4_sa_clim')
      !s ALLOCATE(so4_sa_clim(row_length,rows,model_levels))
      ! Allocate array for 2D photolysis data if required
      ! Altitude of grid-cell interfaces
    CASE ('interf_z')
      ALLOCATE(interf_z(row_length,rows,0:model_levels))
      CALL calc_interf_z(row_length, rows, model_levels, interf_z)
      ! Grid box surface area (in m^2).
    CASE ('grid_surf_area')
      ALLOCATE(surf_area(row_length,rows))
      DO j = 1, rows
        DO i = 1, row_length
          surf_area(i,j) = r_theta_levels(i,j,0) * r_theta_levels(i,j,0) *     &
                           delta_lambda * delta_phi * FV_cos_theta_latitude(i,j)
        END DO
      END DO
    END SELECT
  END DO

  !++SAN
  CALL umPrint('SAN; BOX_UKCA: end first_call',src='box_ukca')

END IF  ! l_first_call

! For the explicit Backward-Euler Offline Oxidants scheme, ensure that the
! height to which the chemistry is integrated is above the first constant rho
! level. This is necessary to retain processor bit comparability.
IF (l_ukca_offline_be) THEN
  first_constant_r_rho_level = a_inthd(24)
  z_first_constant_r_rho_level =                                               &
    MAXVAL(r_rho_levels(:,:,first_constant_r_rho_level)) - planet_radius
  IF (max_z_for_offline_chem < z_first_constant_r_rho_level) THEN
    cmessage = ' max_z_for_offline_chem is set too low for processor'//        &
               ' bit comparability'
    errcode = 1
    WRITE(umMessage,'(A74,2E12.4)') cmessage,                                  &
          max_z_for_offline_chem, z_first_constant_r_rho_level
    CALL umPrint(umMessage,src=RoutineName)
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

m_atm_modl = submodel_for_sm(atmos_im)   ! submodel code

! Set solar position parameters in UKCA if required

IF (l_req_sin_declination .OR. l_req_equation_of_time) THEN

  i_year_previous = previous_time(1)
  i_hour_previous = previous_time(4)
  i_minute_previous = previous_time(5)
  i_second_previous = previous_time(6)
  i_day_number_previous = previous_time(7)
  sindec  = 0.0
  eq_time = 0.0

  secondssincemidnight = REAL(i_hour_previous * 3600 +                         &
                              i_minute_previous * 60 +                         &
                              i_second_previous)

  CALL solpos (i_day_number_previous, i_year_previous,                         &
       secondssincemidnight, secs_per_stepim(atmos_im), l_planet_obs,          &
       eq_time, sindec, scs, sindec_obs, eqt_obs)

  errcode = 0

  IF (l_req_sin_declination)                                                   &
    CALL ukca_set_environment('sin_declination', sindec,                       &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)

  IF (errcode <= 0 .AND. l_req_equation_of_time)                               &
    CALL ukca_set_environment('equation_of_time', eq_time,                     &
                              errcode, error_message=ukca_errmsg,              &
                              error_routine=ukca_errproc)

  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF

END IF

! Set additional UKCA input fields if required

IF (ALLOCATED(latitude)) THEN
  CALL ukca_set_environment('latitude', latitude,                              &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF
END IF

IF (ALLOCATED(longitude)) THEN
  CALL ukca_set_environment('longitude', longitude,                            &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF
END IF

IF (ALLOCATED(sin_latitude)) THEN
  CALL ukca_set_environment('sin_latitude', sin_latitude,                      &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF
END IF

IF (ALLOCATED(cos_latitude)) THEN
  CALL ukca_set_environment('cos_latitude', cos_latitude,                      &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF
END IF

IF (ALLOCATED(tan_latitude)) THEN
  CALL ukca_set_environment('tan_latitude', tan_latitude,                      &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF
END IF

IF (ALLOCATED(interf_z)) THEN
  CALL ukca_set_environment('interf_z', interf_z,                              &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF
END IF

IF (ALLOCATED(surf_area)) THEN
  CALL ukca_set_environment('grid_surf_area', surf_area,                       &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF

END IF

!++SAN comment out aerosol climatology for now
!s IF (ALLOCATED(so4_sa_clim)) THEN
  ! Read climatological aerosol data into so4_sa_clim
!s   z_top_of_model = a_realhd(rh_z_top_theta)
!s  CALL ukca_read_aerosol(i_year, i_month, row_length, rows,                    &
!s           model_levels, latitude,                                             &
!s           eta_theta_levels(1:model_levels) * z_top_of_model,                  &
!s           l_ukca_use_background_aerosol, so4_sa_clim)
  ! Pass climatological aerosol data into UKCA
!s  CALL ukca_set_environment('so4_sa_clim', so4_sa_clim,                        &
!s                            errcode, error_message=ukca_errmsg,                &
!s                            error_routine=ukca_errproc)
!s  IF (errcode > 0) THEN
!s    cmessage = ukca_errmsg
!s    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
!s  END IF
!s END IF


!++SAN Hardcode z_top_of_model to be = 85km
! Maybe not needed?
z_top_of_model = 85000.0

! Load emissions from NetCDF files. If this is the first call, the routine below
! will call other routines to check validity of files and setup the data
! structure for receiving NetCDF emissions
!++SAN - in first instance, no emissions in box model. In future version, 
!        should be able to read in emissions for single grid cell into
!        box model, and/or pass in dummy fields which would mean that
!        ukca_emiss_ctl can run without errors inside ukca_main1
!s WRITE(umMessage,*) 'BOX_UKCA_MOD, SAN: Calling ukca_set_emissions_from_nc'
!s CALL umPrint(umMessage, src='box_ukca_mod')
!s CALL ukca_set_emissions_from_nc(i_year, i_month, i_day, i_hour, i_minute,      &
!s                 i_second, row_length, rows, model_levels,                      &
!s                 global_row_length, global_rows)
!++SAN see if I can hack this by adding custom sizes for the emission files
!s CALL ukca_set_emissions_from_nc(i_year, i_month, i_day, i_hour, i_minute,      &
!s                 i_second, row_length, rows, glob_emiss_levels,                      &
!s                 glob_emiss_row_length, glob_emiss_rows)


! If this is an offline_oxidants configuration, read in oxidant values from
! NetCDF files
!++SAN - Not in first working version of box model. Will need to be added when
!        when making offline oxidant version
IF (l_ukca_offline .OR. l_ukca_offline_be) THEN

  ! put reading of offline oxidents on a logical switch - will need to make
  ! this a namelist option. 
  IF (l_read_offox_ncdf) THEN
    WRITE(umMessage,*) 'SAN; Setting offline oxidents'
    CALL umPrint(umMessage, src='box_ukca_mod')
    !++SAN changes to this code needed for it to work with box model -
    !      Currently fails due to dimension mismatch between files and model
    CALL ukca_set_oxidants_from_nc(i_year, i_month, i_day, i_hour, i_minute,     &
                  i_second, row_length, rows, model_levels,                      &
                  global_row_length, global_rows, environ_varnames)
  ELSE
    IF (l_first_call) THEN
      n_offline = 0
      oxid_varnames(:) = ' '       ! Names of oxidant species
      DO l = 1, SIZE(environ_varnames)
        IF ( ANY(fldname_oxid == TRIM(environ_varnames(l))) ) THEN
          n_offline = n_offline + 1
          oxid_varnames(n_offline) = TRIM(environ_varnames(l))
        END IF
      END DO
      WRITE(umMessage,*) 'SAN;  n_offline:', n_offline
      CALL umPrint(umMessage, src='box_ukca_mod')

      WRITE(umMessage,*) 'SAN;  Allocating offline_ox'
      CALL umPrint(umMessage, src='box_ukca_mod')
      !++SAN Allocate offline oxidant array
      IF (.NOT. ALLOCATED(offline_ox)) ALLOCATE(offline_ox(n_offline))
      offline_ox(:)%file_name     = ' '
      offline_ox(:)%var_name      = ' '
      offline_ox(:)%tracer_name   = ' '

      offline_ox(:)%update_freq   = 0
      offline_ox(:)%update_type   = 0
      offline_ox(:)%last_update   = 0
      offline_ox(:)%l_update      = .FALSE.

      WRITE(umMessage,*) 'SAN;  Allocating offline_ox%values & diags'
      CALL umPrint(umMessage, src='box_ukca_mod')

      DO i = 1, n_offline
        ALLOCATE (offline_ox(i)%values(row_length, rows, model_levels))
        ALLOCATE (offline_ox(i)%diags (row_length, rows, model_levels))
      END DO

      ! Read in oxidant values from namelist
      CALL read_offline_oxidant_mmr(errcode, cmessage, ukca_box_nml,          &
        o3_in, oh_in, ho2_in, h2o2_in, no3_in)

      IF (errcode /= 0) THEN
        CALL Ereport(RoutineName,errcode,cmessage)
      END IF

      ! Assign values to chemistry arrays based on inputs
      DO i = 1, n_offline
        offline_ox(i)%values(:,:,:) = 0.0 !SAN - may want to set this to default values
        offline_ox(i)%var_name    = TRIM(oxid_varnames(i))
        offline_ox(i)%tracer_name = TRIM(oxid_varnames(i))
        IF (oxid_varnames(i) == 'O3' .AND. o3_in /= rmdi)                     &
          offline_ox(i)%values(:,:,:) = o3_in
        IF (oxid_varnames(i) == 'OH' .AND. oh_in /= rmdi)                     &
          offline_ox(i)%values(:,:,:) = oh_in
        IF (oxid_varnames(i) == 'HO2' .AND. ho2_in /= rmdi)                   &
          offline_ox(i)%values(:,:,:) = ho2_in
        IF (oxid_varnames(i) == 'H2O2' .AND. h2o2_in /= rmdi)                 &
          offline_ox(i)%values(:,:,:) = h2o2_in
        IF (oxid_varnames(i) == 'NO3' .AND. no3_in /= rmdi)                   &
          offline_ox(i)%values(:,:,:) = no3_in
      END DO

    END IF  ! l_first_call

    !++SAN Map data for offline oxidants to environment fields,
    !      code from ukca_read_offline_oxidants
    ! Populate the oxidant environment variables with values read into the
    ! offline_ox structure from the netCDF input fields. This is done at every
    ! timestep, with the same values being sent until the species is to be updated
    ! from the file.
    ! Need to use a temporary data array to match ukca_set_environment argument type
    IF (ALLOCATED(tmpdat1)) DEALLOCATE(tmpdat1)
    DO l = 1, n_offline
      ALLOCATE(tmpdat1( model_levels, row_length, rows))
      tmpdat1(:,:,:) = offline_ox(l)%values(:,:,:)
      CALL ukca_set_environment(TRIM(offline_ox(l)%tracer_name), tmpdat1, errcode,  &
                             error_message=ukca_errmsg, error_routine=ukca_errproc)
      DEALLOCATE(tmpdat1)
      IF ( errcode /= 0 ) THEN
        WRITE(cmessage,'(A)') TRIM(ukca_errproc)//': '//TRIM(ukca_errmsg)
        CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
      END IF
    END DO    ! loop over cdf_offline_flds
  END IF  ! l_read_offox_ncdf

END IF  ! l_ukca_offline

! If using a Stratospheric scheme, values of GHG and CFC species are required as
! Lower Boundary Conditions. One source of the values, if ukca_scenario = 'RCP',
! is an external file, which is read by the following routine and the values
! interpolated to the current time.
!++SAN - not neccesary for box model
!++SAN   might need to find another way to fill in later...
!sIF ((i_ukca_chem == ukca_chem_strat .OR. i_ukca_chem == ukca_chem_strattrop    &
!s     .OR.  i_ukca_chem == ukca_chem_cristrat)                                  &
!s     .AND. i_ukca_scenario == i_ukca_scenario_rcp ) THEN
!s  CALL ukca_scenario_rcp(i_year, i_day_number)
!sEND IF

! Copy the required tracers from the UM's UKCA tracer array to ukca_tracer_data
! in the order expected by UKCA.
! Include the UM's specific humidity tracer if required
!++SAN - this is replaced with a much simpler one-to-one transfer
!        of data from parent to daughter in box model
!s CALL ukca_tracers_copy_from_um(tracer_ukca_um, q, ukca_tracer_data)
! Allocate arrays here. Need to deallocate later
WRITE(umMessage,*) 'SAN; Allocating tracer fields'
CALL umPrint(umMessage, src='box_ukca_mod')

IF (.NOT. ALLOCATED(ukca_tracer_data))                                         &
     ALLOCATE(ukca_tracer_data(row_length, rows, model_levels, tr_ukca))
WRITE(umMessage,*) 'SAN; Initialising data to zero'
CALL umPrint(umMessage, src='box_ukca_mod')

ukca_tracer_data(:,:,:,:) = 0.0

!++SAN - 1-1 match between parent and daughter tracer array in BM
!        all need to do is feed q to the ukca_tracer_data
DO i = 1, n_tracer_required
  IF (TRIM(tracer_varnames(i)) == 'H2O') THEN
    ukca_tracer_data(:,:,:,i) = q(:,:,:)
  ELSE
    ukca_tracer_data(:,:,:,i) = tracer_ukca_um(:,:,:,i)
  END IF
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A3,I3,A2,A10,A3,E12.4)') ' * ',                          &
        i, '. ', tracer_varnames(i), ' = ', ukca_tracer_data(1,1,1,i) 
    CALL umPrint(umMessage, src='box_model')
  END IF
END DO

!++SAN
CALL umPrint('SAN; Call ukca_step()',src='box_ukca')

! Set the current_time array to pass to photol_ctl and UKCA
current_time(1) = i_year
current_time(2) = i_month
current_time(3) = i_day
current_time(4) = i_hour
current_time(5) = i_minute
current_time(6) = i_second
current_time(7) = i_day_number

! If required, calculate photolysis rates using chosen photolysis scheme
! and pass result to UKCA
IF (l_ukca_chem .AND. i_ukca_photol /= 0) THEN
  ! fetch photolysis reaction data
  CALL ukca_get_photol_reaction_data(ratj_data, ratj_varnames)
  jppj = SIZE(ratj_varnames)

  ! Is this a chemistry timestep?
  l_chem_step = MOD(stepim(atmos_im),                                          &
                     chem_timestep / INT(secs_per_stepim(atmos_im))) == 0

  ! calculate photolysis rates
  !++SAN calls box model version of routine in ukca_box_photol_ctl
  CALL photol_ctl(row_length, rows, model_levels, jppj,                        &
                  ratj_data, ratj_varnames,                                    &
                  longitude, sin_latitude, cos_latitude, tan_latitude,         &
                  p_theta_levels(1:row_length, 1:rows, 1:model_levels),        &
                  theta(1:row_length, 1:rows, 1:model_levels),                 &
                  exner_theta_levels(1:row_length, 1:rows, 1:model_levels),    &
                  p_rho_levels(1:row_length, 1:rows, 1:model_levels),          &
                  q(1:row_length, 1:rows, 1:model_levels),                  &
                  qcl(1:row_length, 1:rows, 1:model_levels),                   &
                  qcf(1:row_length, 1:rows, 1:model_levels),                   &
                  secondssincemidnight, sindec, eq_time, current_time,         &
                  ukca_tracer_data(1:row_length, 1:rows, 1:model_levels, :),   &
                  tracer_varnames, z_top_of_model, l_chem_step)

END IF

!++SAN Checking all inputs before calling ukca_step
IF (PrintStatus >= PrStatus_Oper) THEN
  CALL umPrint('SAN; BOX_UKCA: Checking all inputs to ukca_step',src='box_ukca')
  WRITE(umMessage,*) 'SAN: stepim(atmos_im): ', stepim(atmos_im)
  CALL umPrint(umMessage,src='box_ukca')
  WRITE(umMessage,'(A18,6(I4,1X))') 'SAN: CURRENT TIME: ',                          &
        current_time(1), current_time(2), current_time(3), current_time(4),        &
        current_time(5), current_time(6)
  CALL umPrint(umMessage,src='box_ukca')
  WRITE(umMessage,*) 'SAN; ukca_tracer_data: ', ukca_tracer_data(1,1,1,1)
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; ntp_data: ', ntp_data(1,1,1,1)
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; errcode: ', errcode
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; previous_time: ', previous_time
  CALL umPrint(umMessage, src='box_ukca_mod')
  WRITE(umMessage,*) 'SAN; ukca_errmsg: ', ukca_errmsg
  CALL umPrint(umMessage, src='box_ukca_mod')
END IF
!--SAN

CALL ukca_step(stepim(atmos_im), current_time, ukca_tracer_data, ntp_data,      &
                r_theta_levels, r_rho_levels,                                   &
                errcode, previous_time=previous_time,                           &
                error_message=ukca_errmsg, error_routine=ukca_errproc)
IF (errcode > 0) THEN
  cmessage = ukca_errmsg
  CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
END IF

!---------------------------------------------------------------------------------

! Return updated tracers from ukca_tracer_data to the UM's UKCA tracer array
! Also update the UM's humidity tracer if required
!++SAN - Reverse the previous call
WRITE(umMessage,*) 'BOX_UKCA, SAN: Outputting tracer fields'
CALL umPrint(umMessage, src='box_ukca_mod')
DO i = 1, n_tracer_required
  IF (TRIM(tracer_varnames(i)) == 'H2O') THEN
    q(:,:,:) = ukca_tracer_data(:,:,:,i)
  ELSE
    tracer_ukca_um(:,:,:,i) = ukca_tracer_data(:,:,:,i)
  END IF
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A3,I3,A2,A10,A3,E12.4)') ' * ',                          &
        i, '. ', tracer_varnames(i), ' = ', tracer_ukca_um(1,1,1,i) 
    CALL umPrint(umMessage, src='box_model')
  END IF
END DO

WRITE(umMessage,*) 'BOX_UKCA, SAN: Deallocating ukca_tracer_data'
CALL umPrint(umMessage, src='box_ukca_mod')

IF (ALLOCATED(ukca_tracer_data)) DEALLOCATE (ukca_tracer_data)

!sCALL umPrint('BOX_UKCA, SAN: Deallocating UkcaD1Codes', src='box_ukca_mod')
! Doesn't work here because will be needed again on next timestep...
!s IF (ALLOCATED(UkcaD1Codes)) DEALLOCATE (UkcaD1Codes)

WRITE(umMessage,*) 'BOX_UKCA, SAN; Finished deallocating'
CALL umPrint(umMessage, src='box_ukca_mod')

! Deallocate ukca_um_interf_mod variables
!++SAN - comment out for now?...
!s CALL ukca_um_dealloc()

!++SAN
WRITE(umMessage,*) 'SAN; Exiting BOX_UKCA'
CALL umPrint(umMessage, src='box_ukca_mod')

! Reset first call flag
l_first_call = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE box_ukca

END MODULE box_ukca_mod
