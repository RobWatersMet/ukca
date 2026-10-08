! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module handling the reading in and processing of photolysis data,
!  before this data is passed to UKCA.
!  Modified for use with the box model version of UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_box_photol_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_BOX_PHOTOL_CTL_MOD'

CONTAINS

SUBROUTINE photol_ctl(row_length, rows, model_levels, jppj,                    &
                      ratj_data, ratj_varnames,                                &
                      longitude, sin_latitude, cos_latitude, tan_latitude,     &
                      p_theta_levels, theta,                                   &
                      exner_theta_levels, p_rho_levels,                        &
                      q, qcl, qcf,                                             &
                      seconds_since_midnight, sin_declination,                 &
                      equation_of_time, current_time,                          &
                      tracers, tracer_varnames, z_top_of_model,                &
                      l_in_chem_timestep)
!
! Purpose: Subroutine to handle reading in and processing of photolysis
! data before this data is passed to UKCA.
!
!          Called from BOX_UKCA.
!
! ---------------------------------------------------------------------
!

USE ukca_api_mod,       ONLY: ukca_get_config, ukca_set_environment,           &
!                             ukca_calc_ozonecol,                              &
                              ukca_photol_varname_len,                         &
                              ukca_maxlen_message, ukca_maxlen_procname,       &
                              ukca_maxlen_fieldname

!==NLA
!USE ukca_phot2d,        ONLY: ukca_photin, ntphot
USE ukca_option_mod,    ONLY: i_ukca_photol, l_ukca_strat,                     &
                              l_ukca_stratcfc, l_ukca_strattrop,               &
                              l_ukca_asad_columns,                             &
                              fastjx_prescutoff, fastjx_mode,                  &
                              l_environ_jo2, l_environ_jo2b,                   &
                              chem_timestep

USE photol_config_specification_mod,  ONLY:  photol_config, fjx_mode_fastjx,   &
                              fjx_mode_merged, fjx_mode_2Donly,                &
                              photolysis_off => i_scheme_nophot,               &
                              photolysis_2d => i_scheme_phot2d,                &
                              photolysis_fastjx => i_scheme_fastjx
USE photol_constants_mod,  ONLY: rhour_per_day => const_rhour_per_day,         &
                              pi => const_pi, pi_over_180 => const_pi_over_180,&
                              recip_pi_over_180 => const_recip_pi_over_180

USE ukca_solang_mod,    ONLY: ukca_solang
!s USE ukca_um_interf_mod, ONLY: area_cloud_fraction, pstar, fland,               &
!s                              conv_cloud_lwp, conv_cloud_amount,               &
!s                              conv_cloud_base, conv_cloud_top,                 &
!s                              surf_albedo, land_sea_mask,                      &
!s                              so4_aitken, so4_accum
!++SAN in box model - these are taken from box_model instead
USE ukca_environment_fields_mod, ONLY: area_cloud_fraction, pstar, fland,      &
                                 conv_cloud_lwp, conv_cloud_amount,            &
                                 conv_cloud_base, conv_cloud_top,              &
                                 surf_albedo, land_sea_mask,                   &
                                 so4_aitken, so4_accum

!++SAN box model versions of strat routines
! - for now, removing as no photolysis
!s USE ukca_box_strat_photol_mod,  ONLY: strat_photol
!s USE ukca_box_dissoc_mod,        ONLY: strat_photol_init
!s USE rad_ctl_mod,        ONLY: rad_ctl_jo2, rad_ctl_jo2b
!==NLA USE qsat_mod,           ONLY: qsat_wat_mix
USE ukca_um_legacy_mod, ONLY: rh_z_top_theta
USE model_time_mod,     ONLY: secs_per_stepim, stepim
USE submodel_mod,       ONLY: atmos_im

! Dimensions for tracer array
USE nlsizes_namelist_mod, ONLY: tr_ukca

!++SAN Get ratj_defs for the multiplicative factor jfacta
USE ukca_chem_defs_mod,   ONLY: ratj_defs

!++SAN Added namelist options for Box model I/O
USE nlstcall_box_namelist_mod, ONLY: photol_jrate_in_filename
! l_ukca_simple_photol_input,               &

!++SAN adding diagnostic print statements
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE um_parcore,             ONLY: mype

! UM profiling and error handling
USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim
USE errormessagelength_mod,   ONLY: errormessagelength
USE ereport_mod,              ONLY: ereport

IMPLICIT NONE

! -- Input arguments -- !
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: jppj
CHARACTER(LEN=10), POINTER, INTENT(IN) :: ratj_data(:,:)
CHARACTER(LEN=ukca_photol_varname_len), POINTER, INTENT(IN) :: ratj_varnames(:)
REAL, INTENT(IN)    :: longitude(row_length, rows)
REAL, INTENT(IN)    :: sin_latitude(row_length, rows)
REAL, INTENT(IN)    :: cos_latitude(row_length, rows)
REAL, INTENT(IN)    :: tan_latitude(row_length, rows)
REAL, INTENT(IN)    :: p_theta_levels(row_length, rows, model_levels)
REAL, INTENT(IN)    :: theta(row_length, rows, model_levels)
REAL, INTENT(IN)    :: exner_theta_levels(row_length, rows, model_levels)
REAL, INTENT(IN)    :: p_rho_levels(row_length, rows, model_levels)
REAL, INTENT(IN)    :: q(row_length, rows, model_levels)
REAL, INTENT(IN)    :: qcl(row_length, rows, model_levels)
REAL, INTENT(IN)    :: qcf(row_length, rows, model_levels)
REAL, INTENT(IN)    :: seconds_since_midnight
REAL, INTENT(IN)    :: sin_declination
REAL, INTENT(IN)    :: equation_of_time
INTEGER, INTENT(IN) :: current_time(7)
REAL, INTENT(IN)    :: tracers(row_length, rows, model_levels, tr_ukca + 1)
CHARACTER(LEN=ukca_maxlen_fieldname), POINTER, INTENT(IN) :: tracer_varnames(:)
REAL, INTENT(IN)    :: z_top_of_model
LOGICAL, INTENT(IN) :: l_in_chem_timestep   ! chemical timestep flag

!++SAN - no stash in box model
!s INTEGER, INTENT(IN) :: len_stashwork ! Length of diagnostics array
!s REAL, INTENT(IN OUT) :: stashwork(len_stashwork)  ! STASH workspace

! -- Local variables -- !
REAL, ALLOCATABLE, SAVE :: photol_rates_fastjx(:,:,:,:)
REAL, ALLOCATABLE, SAVE :: photol_rates_2d(:,:,:,:)
REAL, ALLOCATABLE       :: photol_rates(:,:,:,:)
REAL, ALLOCATABLE, SAVE :: photol_rates_save(:,:,:,:) !++SAN
REAL, ALLOCATABLE       :: land_fraction(:,:)
REAL, ALLOCATABLE       :: p_layer_boundaries(:,:,:)
REAL, ALLOCATABLE       :: rel_humid_frac(:,:,:)
REAL, ALLOCATABLE       :: t_theta_levels(:,:,:)
REAL, ALLOCATABLE       :: qsmr(:,:,:)
REAL                    :: photol1d(row_length * rows, jppj)
REAL                    :: photol2d(row_length, rows, jppj)
REAL                    :: photol_strat(row_length, rows, jppj)
! Temporary SO4 aerosol arrays
REAL                    :: photol_so4_aitken(row_length, rows, model_levels)
REAL                    :: photol_so4_accum(row_length, rows, model_levels)
! cos of solar zenith angle, used for strat_photol
REAL                    :: cos_zenith_angle(row_length, rows)
! temporary copy of ozone tracer
REAL                    :: ozone_tracer(row_length, rows, model_levels)
! ozonecol for stratospheric chemistry
REAL                    :: ozonecol(row_length, rows, model_levels)
! variables for calculating land_fraction
INTEGER, ALLOCATABLE    :: land_index(:)
INTEGER                 :: land_points
! top of model
!s REAL                    :: z_top_of_model
! loop variables
INTEGER                 :: i
INTEGER                 :: j
INTEGER                 :: k
INTEGER                 :: l
! size of theta levels
INTEGER                 :: theta_field_size
! index of ozone tracer
INTEGER                 :: n_o3
! required constants (from ukca_main/ukca_constants)
REAL, PARAMETER         :: min_SO4_val = 1.0e-25
REAL, PARAMETER         :: c_o3        = 1.657
REAL, SAVE              :: fxb
REAL, SAVE              :: fxc
! tan(declination)
REAL                    :: tan_declin
! local time variables
REAL                    :: tloc(row_length, rows) ! local time
REAL                    :: daylen(row_length, rows) ! local daylength
REAL                    :: cs_hour_ang(row_length, rows) ! cosine hour angle
REAL                    :: tgmt ! GMT time (decimal representation)
REAL                    :: r_secs_per_step
INTEGER                 :: i_day_number
INTEGER                 :: i_hour
REAL                    :: r_minute
!++SAN variables for reading in j-rates from file
CHARACTER(LEN=ukca_photol_varname_len) :: jlabel_dummy
REAL                    :: photol_rate_dummy
LOGICAL                 :: l_exist
INTEGER                 :: nlines
INTEGER                 :: io_status
LOGICAL                 :: j_match

!++SAN Array containing quantum yields (or multiplication factor)
REAL,  ALLOCATABLE      :: jfacta(:)
LOGICAL, ALLOCATABLE    :: l_jmatch_arr(:)       ! Check j rate has been matched

! if ukca_photin called, set to true (to prevent 2D photolysis data being
! read twice)
LOGICAL, SAVE           :: l_ukca_photin_called = .FALSE.
! first time calling photol_ctl
LOGICAL, SAVE           :: l_first_call = .TRUE.

! error handling
INTEGER                             :: errcode = 0    ! Error code: ereport
CHARACTER(LEN=errormessagelength)   :: cmessage       ! Error message
CHARACTER(LEN=ukca_maxlen_message)  :: ukca_errmsg    ! Error return message
CHARACTER(LEN=ukca_maxlen_procname) :: ukca_errproc   ! Routine in which error
                                                      ! was trapped
CHARACTER(LEN=*), PARAMETER :: errproc_suffix = ' in UKCA'
                                                      ! Text to append to
                                                      ! routine name
! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! -- body of subroutine begins here -- !

WRITE(umMessage,*) 'SAN: In PHOTOL_CTL'
CALL umPrint(umMessage,src=RoutineName)

! get number of seconds in parent model timestep as real
r_secs_per_step = REAL(chem_timestep)
! set other time variables
i_day_number = current_time(7)
i_hour = current_time(4)
r_minute = REAL(current_time(5)) - r_secs_per_step / 60.0
! r_minute shifted for UKCA compatibility - see args of ukca_chemistry_ctl

! allocate array to store Fast-JX rates
IF (.NOT. ALLOCATED(photol_rates_fastjx)) THEN
  WRITE(umMessage,*) 'SAN: Allocating photol_rates_fastjx'
  CALL umPrint(umMessage,src=RoutineName)
  ALLOCATE(photol_rates_fastjx(row_length, rows, model_levels, jppj))
END IF
photol_rates_fastjx = 0.0
! allocate main photolysis rates array
IF (.NOT. ALLOCATED(photol_rates)) THEN
  WRITE(umMessage,*) 'SAN: Allocating photol_rates'
  CALL umPrint(umMessage,src=RoutineName)
  ALLOCATE(photol_rates(row_length, rows, model_levels, jppj))
END IF
photol_rates = 0.0

!++SAN allocate array for saving input photolysis rates
! Set to zero on first call
IF (.NOT. ALLOCATED(photol_rates_save)) THEN
  WRITE(umMessage,*) 'SAN: Allocating photol_rates_save'
  CALL umPrint(umMessage,src=RoutineName)
  ALLOCATE(photol_rates_save(row_length, rows, model_levels, jppj))
END IF
IF (l_first_call) photol_rates_save = 0.0
!--SAN

theta_field_size = row_length * rows

! if we're not in a chemical timestep, we only need to send fields to UKCA
IF (.NOT. l_in_chem_timestep) THEN
  !++SAN
  WRITE(umMessage,*) 'SAN: Not in chem timestep'
  CALL umPrint(umMessage,src=RoutineName)
  !--SAN

  ! Pass dummy photolysis rates array to UKCA
  CALL ukca_set_environment('photol_rates', photol_rates,                      &
                            errcode, error_message=ukca_errmsg,                &
                            error_routine=ukca_errproc)
  IF (errcode > 0) THEN
    cmessage = ukca_errmsg
    CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
  END IF

  !s IF (l_first_call .AND. l_strat_chem) THEN
  !s     ! Set up arrays for later strat_photol calculations
  !s     CALL strat_photol_init(theta_field_size)
  !s   END IF

  l_first_call = .FALSE.
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!++SAN - New photolysis option for boxmodel.
!        Offline photolysis which is either on or off,
!        with rates predefined and read in from a file.
!s IF(i_ukca_photol /= ukca_photolysis_off .AND. l_ukca_simple_photol_input) THEN

!++SAN for now - just hardcode to do on first call
IF (l_first_call) THEN

   fxb         = 23.45/recip_pi_over_180
   fxc         = rhour_per_day/pi

  IF (.NOT. ALLOCATED(jfacta)) ALLOCATE(jfacta(jppj))
  IF (.NOT. ALLOCATED(l_jmatch_arr)) ALLOCATE(l_jmatch_arr(jppj))
  jfacta(:) = 0.0
  l_jmatch_arr(:) = .FALSE.


  WRITE(umMessage,'(A)') 'SAN: Reading photolysis rates from file'
  CALL umPrint(umMessage,src=RoutineName)
  
  ! Check file exists
  IF (mype == 0) THEN
    INQUIRE (FILE=TRIM(photol_jrate_in_filename), EXIST=l_exist)
    IF (.NOT. l_exist) THEN
      errcode  = 1
      cmessage = TRIM(photol_jrate_in_filename)//': photolysis rates file does not exist'
      CALL ereport(RoutineName,errcode,cmessage)
    END IF 
  END IF

  ! Open file to read in photolysis rates
  IF (mype == 0) THEN
    WRITE(umMessage,'(2A)') 'SAN: Opening ', photol_jrate_in_filename
    CALL umPrint(umMessage,src=RoutineName)
    OPEN( 23, FILE=TRIM(photol_jrate_in_filename) )

    ! Find out length of photolysis file
    nlines = 0
    DO
      READ(23,*,iostat=io_status)
      IF (io_status/=0) EXIT
      nlines=nlines+1
    END DO

    CLOSE(23)

    ! Open again and map to photolysis array
    OPEN(24, FILE=TRIM(photol_jrate_in_filename))

    ! Map each input photolysis rate to file data
    DO l = 1, nlines

      READ(24,'(e10.4,1X,A10)') photol_rate_dummy, jlabel_dummy
      j_match = .FALSE.

      ! Loop through all j rates and load to photol_rates
      DO j = 1, jppj
        IF (TRIM(jlabel_dummy) == TRIM(ratj_varnames(j))) THEN
          WRITE(umMessage,'(A,A7,e10.4)') 'Mapping ', ratj_varnames(j), photol_rate_dummy
          CALL umPrint(umMessage,src=RoutineName)
          photol_rates_save(:,:,:,j) = photol_rate_dummy
          j_match = .TRUE.
          l_jmatch_arr(j) = .TRUE.
        END IF
      END DO

      IF (j_match .eqv. .FALSE.) THEN
        WRITE(umMessage,'(2A)') 'WARNING: failed to match input ', jlabel_dummy
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    END DO

    ! Close input jrate file
    CLOSE(24)

    ! Check if any entries iin the photol_rates array have failed to match
    DO j = 1, jppj
      IF (.NOT. l_jmatch_arr(j)) THEN
        WRITE(umMessage,'(A,I3,1X,A7)') 'WARNING: No input jrate data for ',   &
            j, ratj_varnames(j)
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    END DO

    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(A)') 'Applying "quantum yield"'
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Apply "quantum yield" (really just a multiplicative factor)
    ! To each j value. Code based on section from fastjx_specs
    DO i=1,jppj
      jfacta(i)=ratj_defs(i)%jfacta/100.0e0
      photol_rates_save(:,:,:,i) = photol_rates_save(:,:,:,i)*jfacta(i)
      IF (PrintStatus >= PrStatus_Oper) THEN
        WRITE(umMessage,'(I6,E12.3,A12,E12.3)')                                 &
            i,jfacta(i),ratj_varnames(i),photol_rates_save(:,:,:,i)
        CALL umPrint(umMessage,src='fastjx_specs')
      END IF
    END DO

    WRITE(umMessage,'(A)') 'SAN: Finished reading in photolysis J-rates'
    CALL umPrint(umMessage,src=RoutineName)
    !--SAN

  END IF ! mype == 0

  WRITE(umMessage,'(A)') 'SAN: Mapping photolysis J-rates'
  CALL umPrint(umMessage,src=RoutineName)
  ! Map to photol_rates
  photol_rates = photol_rates_save

ELSE ! not l_first_call

  ! On later timesteps, just set photol_rates array to saved array
  ! Map to photol_rates
  WRITE(umMessage,'(A)') 'SAN: Mapping photolysis J-rates'
  CALL umPrint(umMessage,src=RoutineName)
  photol_rates = photol_rates_save

END IF ! l_first_call

!++SAN Apply diurnal cycle to the phototlysis rates
!++SAN
WRITE(umMessage,*) 'SAN: Calculating tloc and day length'
CALL umPrint(umMessage,src=RoutineName)
!--SAN

! Calculate tloc and day length here
tgmt = REAL(i_hour) + r_minute / 60.0 + r_secs_per_step * 0.5 / 3600.0
tloc = tgmt + 24.0 * longitude / 360.0
WHERE (tloc > 24.0) tloc = tloc - 24.0

daylen = 0.0
tan_declin = TAN(fxb * SIN(pi_over_180 * (266.0 + i_day_number)))
IF (tan_declin == 0.0) THEN
  cs_hour_ang = 0.0
ELSE
  cs_hour_ang = MAX(-1.0, MIN(1.0, -tan_latitude * tan_declin))
END IF
daylen = fxc * ACOS(cs_hour_ang)

!++SAN
WRITE(umMessage,'(A,E10.3)') 'tloc =       ', tloc(1,1)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,E10.3)') 'tan_declin = ', tan_declin
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,E10.3)') 'cs_hour_ang = ', cs_hour_ang(1,1)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,E10.3)') 'daylen =      ', daylen(1,1)
CALL umPrint(umMessage,src=RoutineName)
!--SAN

!++SAN - note always calculate cos_zenith_angle in box model
WRITE(umMessage,*) 'SAN: Calculating ukca_solang'
CALL umPrint(umMessage,src=RoutineName)
!--SAN

! Calculate cos_zenith_angle
CALL ukca_solang(sin_declination, seconds_since_midnight,                    &
             r_secs_per_step, equation_of_time,                            &
             sin_latitude, cos_latitude, longitude, theta_field_size,      &
             cos_zenith_angle)

!++SAN
WRITE(umMessage,*) 'cos_zenith_angle = ', cos_zenith_angle(1,1)
CALL umPrint(umMessage,src=RoutineName)
!--SAN

WRITE(umMessage,*) 'SAN: Starting k loop'
CALL umPrint(umMessage,src=RoutineName)

!++SAN - column call redundant in box model
!s IF (l_ukca_asad_columns) THEN
! Fill photol_rates array to pass to UKCA
DO k = 1, model_levels
  WRITE(umMessage,*) 'SAN: reset photol2d'
  CALL umPrint(umMessage,src=RoutineName)

  ! reset photol2d to zero here
  photol2d = 0.0

  WRITE(umMessage,*) 'SAN: Simple input photolysis, applying diurnal cycle'
  CALL umPrint(umMessage,src=RoutineName)

  ! Apply simple diurnal cycle using cos_zenith_angle
  photol2d(:,:,:) = photol_rates(:,:,k,:)
  DO j = 1, rows
    DO i = 1, row_length
      IF (cos_zenith_angle(i,j) > 0) THEN
        DO l = 1, jppj
          photol2d(i,j,l) = photol2d(i,j,l) * cos_zenith_angle(i,j)
        END DO
      ELSE
        photol2d(i,j,:) = 0.0
      END IF
    END DO
  END DO

  WRITE(umMessage,*) 'SAN: Passing back photol2d'
  CALL umPrint(umMessage,src=RoutineName)

  ! construct photol_rates array to pass to UKCA
  photol_rates(:,:,k,:) = photol2d(:,:,:)
END DO

! pass proper diagnostics to STASH array
!++SAN no diags in box model
!++SAN check values of photol rates
DO j=1,jppj
  !++SAN
  WRITE(umMessage,*) 'photol_rates', j, ratj_varnames(j), ' = ',        &
    photol_rates(1,1,1,j)
  CALL umPrint(umMessage,src=RoutineName)
  !--SAN
END DO
!--SAN

! Pass interpolated photolysis rates array to UKCA
WRITE(umMessage,*) 'SAN: ukca_set_environment(photol_rates)'
CALL umPrint(umMessage,src=RoutineName)
CALL ukca_set_environment('photol_rates', photol_rates,                        &
                          errcode, error_message=ukca_errmsg,                  &
                          error_routine=ukca_errproc)
IF (errcode > 0) THEN
  cmessage = ukca_errmsg
  CALL ereport(TRIM(ukca_errproc) // errproc_suffix, errcode, cmessage)
END IF

l_first_call = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE photol_ctl


!++SAN - No need to mixing ozone tracer in UKCA box model
!s SUBROUTINE photol_mix_ozone_tracer(row_length, rows, model_levels,             &
!s                                    ozone_tracer)
!
! Purpose: Wrapper around a call to tr_mix (called for tracers generally in
!          ukca_add_emiss) required to mix ozone tracer before the
!          calculation of ozonecol.
!          This is just here for bit comparability reasons, and should be
!          temporary.
!
!          Called from PHOTOL_CTL.
!
! ---------------------------------------------------------------------
!
!s RETURN
!s END SUBROUTINE photol_mix_ozone_tracer


!++SAN - No photolysis diagnostics in UKCA box model
!s SUBROUTINE photol_diags(row_length, rows, model_levels, jppj,               &
!s                              ratj_varnames, stashwork, len_stashwork,       &
!s                              photol_rates)
!
! Purpose: Subroutine to copy required parts of photolysis rates generated by
!          photol_ctl into output photolysis diagnostics in the
!          section 50 STASHwork array.
!          Modelled after ukca_chem_diags.
!
!          Called from PHOTOL_CTL.
!
! ---------------------------------------------------------------------


END MODULE ukca_box_photol_ctl_mod

