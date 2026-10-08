! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: SETTSCTL -------------------------------------------------
!
!    Purpose: Sets timestep loop control switches and STASHflags.
!             Note STEP on entry is the values at the N+1 (ie. updated)
!             timelevel; control switches are therefore set using
!             A_STEP/O_STEP, whereas physical switches are set using
!             A_STEP-1/O_STEP-1 to ensure correct synchronisation.
!             Note also that step on entry is N (i.e. not updated) when
!             called from INITIAL.
!
!
!++SAN: Cut down to just do simple timestep loop control for UKCA box model,
!       which only uses a single timestep (essentially the chemistry step)
!
!    Programming standard: UM Doc Paper 3, version 8.2
!
!    Project task: C0
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!     ------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Top Level

MODULE settsctl_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETTSCTL_MOD'

CONTAINS
SUBROUTINE settsctl (internal_model,meanlev,icode,cmessage)

USE rad_input_mod,     ONLY: l_radiation
!s USE set_rad_steps_mod, ONLY: set_l_rad_step

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
!++SAN deleted unnecessary dependencies

USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,                  &
     dumpfreqim, dumptimesim, dumptimes_len1,                                  &
     meanfreqim, meanfreq_len1,                                                &
     greg_dump_freq, greg_monthly_dump, end_of_run_dump, l_meaning_sequence
USE conversions_mod, ONLY: isec_per_day
USE submodel_mod, ONLY: n_internal_model, atmos_im

USE nlstcall_mod, ONLY: model_basis_time,                                      &
                         ancil_reftime,                                        &
                         lpp,                                                  &
                         lnc,                                                  &
                         ldump,                                                &
                         lmean,                                                &
                         lprint,                                               &
                         lexit,                                                &
                         lancillary,                                           &
                         lboundary,                                            &
                         lassimilation,                                        &
                         lclimrealyr,                                          &
                         l_fastrun,                                            &
                         lcal360

!s USE file_manager, ONLY: init_file_loop, um_file_type

USE history, ONLY: model_data_time, offset_dumpsim

USE model_time_mod, ONLY:                                                      &
    ancillary_stepsim,  &
    basis_time_days, basis_time_secs, bndary_offsetim, boundary_stepsim,       &
    i_day, i_hour, i_minute, i_month, i_second, i_year, iau_dtresetstep,       &
    stepim, target_end_stepim
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
USE tim2step_mod, ONLY: tim2step
USE time2sec_mod, ONLY: time2sec

IMPLICIT NONE

INTEGER :: internal_model   ! IN  - internal model identifier
INTEGER :: meanlev          ! OUT - Mean level indicator
INTEGER :: icode            ! Out - Return code
CHARACTER(LEN=errormessagelength) :: cmessage  ! Out - Return error message

! ----------------------------------------------------------------------

!  Local variables

INTEGER :: i                ! Loop counters
INTEGER :: time_idx         ! Loop index over selected dump times
INTEGER :: step             ! A_STEP or O_STEP for atmos/ocean
INTEGER :: section,il,im,ntab,it ! STASH variables
INTEGER :: modl             ! Int model no, read from STASH list arra

!s INTEGER :: ancil_ref_days,ancil_ref_secs
!s INTEGER :: ancil_offset_steps ! offset of ref. from basis time
INTEGER :: secs_per_step    ! seconds per timestep
INTEGER :: months_in        ! Number of months into forecast
INTEGER :: secs_per_period
INTEGER :: steps_per_period
INTEGER :: dumpfreq
INTEGER :: offset_dumps
INTEGER :: exitfreq
INTEGER :: target_end_step
INTEGER :: ancillary_steps
INTEGER :: boundary_steps
INTEGER :: bndary_offset
INTEGER :: dumptimes(dumptimes_len1)
INTEGER :: meanfreq (meanfreq_len1)
INTEGER :: elapsed_days
INTEGER :: elapsed_secs
INTEGER :: elapsed_steps
INTEGER :: step_compare

!s TYPE(um_file_type), POINTER :: um_file

LOGICAL :: iau_resetdt ! If .TRUE., reset data time.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETTSCTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0

! Initialise internal_model
internal_model = 1 ! UKCA

!  1. Set timestep loop top-level control switches
!
!  1.0  Initialise control switches which are shared in coupled runs
!

step=                        stepim(internal_model)
secs_per_period=  secs_per_periodim(internal_model)
steps_per_period=steps_per_periodim(internal_model)
secs_per_step=secs_per_period/steps_per_period
dumpfreq=                dumpfreqim(internal_model)
offset_dumps=        offset_dumpsim(internal_model)
target_end_step=  target_end_stepim(internal_model)
ancillary_steps=  ancillary_stepsim(internal_model)
boundary_steps =   boundary_stepsim(internal_model)
bndary_offset  =   bndary_offsetim(internal_model)

exitfreq=dumpfreq

DO i=1,dumptimes_len1
  dumptimes(i)=       dumptimesim(i,internal_model)
END DO ! i
DO i=1,meanfreq_len1
  meanfreq(i)=         meanfreqim(i,internal_model)
END DO ! i
meanlev=0

lassimilation=.FALSE.
ldump=        .FALSE.
lexit=        .FALSE.
lmean=        .FALSE.
lprint=       .FALSE.
lancillary=   .FALSE.
lboundary=    .FALSE.
!
!  1.1  Set up PPfile and NCfile switches for the timestep
!++SAN comment all this out for now - don't think will be using these
!      file systems for box model
!
!s lpp =.FALSE.
!s lnc =.FALSE.

!s DO i = 1, 2
!s   NULLIFY(um_file)
!s   SELECT CASE(i)
!s   CASE (1)
!s     um_file => init_file_loop(handler="portio")
!s   CASE (2)
!s     um_file => init_file_loop(handler="netcdf")
!s   END SELECT
!s   DO WHILE (ASSOCIATED(um_file))
!s     IF (um_file % meta % is_output_file) THEN

      ! Ensure we fallback to switched off
!s       um_file % meta % initialise = .FALSE.

!s       step_compare = step
      ! If init_file_type == 'c' compare init_start_step with previous timestep
!s       IF (ASSOCIATED(um_file % pp_meta)) THEN
!s         IF (um_file % pp_meta % init_file_type == 'c') THEN
!s           step_compare = step - 1
!s         END IF
!s       END IF

      ! Allow mix of real-month & regular reinit. periods on
      ! different units:
!s       IF (um_file % meta % init_steps <  0) THEN ! Real-month reinit.

        ! Select files to be reinitialised on this timestep:
        ! First initialisation may be on any timestep, not just 0
!s         IF (step_compare == um_file % meta % init_start_step) THEN

!s           um_file % meta % initialise =.TRUE.

!s         ELSE IF (.NOT. lcal360                                                 &
!s                  .AND. (TRIM(um_file % id) /= "ppvar")                         &
!s                  .AND. (i_day     ==  1)                                       &
!s                  .AND. (i_hour    ==  secs_per_step/3600           )           &
!s                  .AND. (i_minute  ==  (MOD(secs_per_step,3600))/60 )           &
!s                  .AND. (i_second  ==   MOD(secs_per_step,60)       )           &
!s                  .AND. (step /= 1)) THEN

          ! Gregorian calendar:
          ! Reinitialisation can only be selected when day=1 and
          ! hours, minutes and seconds equals exactly one time-step
          ! This code relies on the rule that the time-step length
          ! divides into a day.

          ! Additionally initialisation only takes place on 1-month,
          ! 3-month (season) or 12-month boundaries. So calculate months
          ! since start of job, and request reinitialisation at the
          ! appropriate time according to which option selected.
!s           months_in = i_month - model_basis_time(2) +                          &
!s              12 * (i_year - model_basis_time(1))

!s           IF (months_in >= um_file % meta % init_start_step .AND.              &
!s              (um_file % meta % init_end_step <  0 .OR.                         &
!s               months_in <= um_file % meta % init_end_step)) THEN

!s             IF (um_file % meta % init_steps == -1) THEN
!s               um_file % meta % initialise = .TRUE. ! Months
!s             ELSE IF (um_file % meta % init_steps == -3 .AND.                   &
!s                      MOD(months_in -                                           &
!s                       um_file % meta % init_start_step,3)  == 0) THEN
!s               um_file % meta % initialise = .TRUE. ! Seasons
!s             ELSE IF (um_file % meta % init_steps == -12 .AND.                  &
!s                      MOD(months_in -                                           &
!s                       um_file % meta % init_start_step,12) == 0) THEN
!s               um_file % meta % initialise = .TRUE. ! Years
!s             END IF ! Of um_file % meta % init_steps = reinit period

!s           END IF ! Of months_in within reinitialisation limits

!s         END IF  ! Of lcal360 and nftunit etc.

!s       ELSE IF (um_file % meta % init_steps >  0) THEN ! Regular reinit.

        ! Select files to be reinitialised on this timestep:
        ! First initialisation may be on any timestep, not just 0
!s         IF (step_compare == um_file % meta % init_start_step) THEN

!s           um_file % meta % initialise = .TRUE.

          ! Deal with subsequent reinitialisation after end of period
!s         ELSE IF ((step-1) > um_file % meta % init_start_step .AND.             &
!s                 (um_file % meta % init_end_step <  0 .OR.                      &
!s                 (step-1) <= um_file % meta % init_end_step)) THEN

!s           IF (MOD((step-1) - um_file % meta % init_start_step,                 &
!s               um_file % meta % init_steps) == 0)                               &
!s               um_file % meta % initialise = .TRUE.

          ! Do not reinitialise files on step 1 if they were
          ! initialised at step 0:
!s           IF (step == 1 .AND. um_file % meta % init_start_step == 0)           &
!s             um_file % meta % initialise = .FALSE.

          ! Sub model id must be atmos
!s           um_file % meta % initialise =                                        &
!s               um_file % meta % initialise                                      &
!s               .AND. internal_model == atmos_im

!s         END IF

!s       ELSE  ! for files not reinitialised, ie. init_steps == 0

        !           Initialise at step 0
!s         um_file % meta % initialise = step == 0 .AND.                          &
!s         (um_file % meta % init_steps == 0       .OR.                           &
!s         um_file % meta % init_start_step == 0)

!s       END IF  ! Of init_steps lt, gt or == 0, ie. reinit. type
!s     ELSE ! Of um_file % meta % is_output_file
!s       um_file % meta % initialise = .FALSE.
!s     END IF

!s     SELECT CASE(i)
!s     CASE (1)
!s       lpp = lpp .OR. um_file % meta % initialise
!s     CASE (2)
!s       lnc = lnc .OR. um_file % meta % initialise
!s     END SELECT

    ! Increment pointer to next file
!s     um_file => um_file % next
!s   END DO
!s END DO

!
!  1.2   Set switches for general internal models.
!        For coupled models dump related switches can only be set when
!        the last internal model in a submodel has completed its group
!        of timesteps. For coupled models the only safe restart point
!        is at the completion of all groups within a model timestep.

IF (n_internal_model == 1 .OR. MOD(step,steps_per_period) == 0) THEN
  ! if not coupled model, or last step in group
  ! ldump   : Write-up dump on this timestep
  gregdmp: IF (lclimrealyr) THEN
    ! Calculate the elapsed time since time zero to ensure dumping consistency
    ! between NRUNs and CRUNs
    CALL time2sec(i_year, i_month, i_day, i_hour, i_minute, i_second,          &
       0, 0, elapsed_days, elapsed_secs, lcal360)
    elapsed_steps = (elapsed_days*isec_per_day+elapsed_secs) / secs_per_step
    IF (greg_dump_freq > 0) THEN
      ! Regular dumps within cycle period
      ldump = (MOD(elapsed_steps,greg_dump_freq) == 0)
    END IF
    IF (greg_monthly_dump .AND. i_day == 1 .AND. i_hour == 0                   &
       .AND. i_minute == 0 .AND. i_second == 0) THEN
      ! Dump at the end of month for Gregorian calendar
      ldump = .TRUE.
    END IF
    IF (greg_dump_freq == 0 .AND. .NOT. greg_monthly_dump) THEN
      ! Dump from the time list
      DO time_idx=1,dumptimes_len1
        ldump=ldump .OR. (step == dumptimes(time_idx))
      END DO
    END IF
    exitdump: IF (stepim(1) == target_end_stepim(1) .AND. end_of_run_dump) THEN
      ! make sure there is a dump at the end of the run to ensure
      ! restartability
      ldump = .TRUE.
    END IF exitdump
  ELSE
    IF (dumpfreq >  0) THEN
      ldump=       (MOD(step,dumpfreq)    == 0)
    ELSE
      ldump=.FALSE.
      DO time_idx=1,dumptimes_len1
        ldump=ldump .OR. (step == dumptimes(time_idx))
      END DO
    END IF
  END IF gregdmp

  !  LMEAN   : Perform climate-meaning from dumps on this timestep
  lmean = .FALSE.
  IF (l_meaning_sequence) THEN
    IF (dumpfreq >  0 .AND. meanfreq(1) >  0) THEN
      lmean=     (MOD(step,dumpfreq)       == 0)
    END IF
  END IF

  !  LEXIT   : Check for exit condition on this timestep
  IF (exitfreq >  0) THEN ! Implies climate meaning

    lexit=  ( (MOD(step,exitfreq)    == 0)  .OR.                               &
              (step  >=  target_end_step) )

  ELSE                    ! No climate meaning

    lexit=    (step  >=  target_end_step)

  END IF

  !  lancillary: Update ancillary fields on this timestep

  !   Convert ancillary reference time to days & secs
!s  CALL time2sec(ancil_reftime(1),ancil_reftime(2),                             &
!s                ancil_reftime(3),ancil_reftime(4),                             &
!s                ancil_reftime(5),ancil_reftime(6),                             &
!s                0,0,ancil_ref_days,ancil_ref_secs,lcal360)

  !   Compute offset in timesteps of basis time from ancillary ref.time
!s  CALL tim2step(basis_time_days-ancil_ref_days,                                &
!s                basis_time_secs-ancil_ref_secs,                                &
!s                steps_per_period,secs_per_period,ancil_offset_steps)

!s  IF (ancillary_steps >  0 .AND. step >  0)                                    &
!s    lancillary=(MOD(step+ancil_offset_steps,ancillary_steps) == 0)

END IF     ! Test for non-coupled or coupled + last step in group

!  lboundary    : Update boundary fields on this timestep
!++SAN no boundaries in box model
!s IF (boundary_steps >  0)                                                       &
!s    lboundary= (MOD(step+bndary_offset,boundary_steps)  == 0)

!  1.2.1  Set switches for atmosphere timestep
!
!s IF (internal_model == atmos_im) THEN
  !  1.2.2 Set switches for all cases

  ! Energy correction switches
  !++SAN DELETE - no energy correction in box model


  !  lassimilation: Perform data assimilation on this timestep
  !++SAN DELETE - no data assimilation in UKCA box model

  ! If using radiation, set whether this is a radiation timestep.
  ! If not using radiation, do nothing; leave the l_rad_step flags in
  ! their initialised states (false, set in set_rad_steps_mod).
  !++SAN comment out for now, no radiation fro first iteration of box model
  !s IF (l_radiation)  CALL set_l_rad_step(step)

!s END IF ! Test for atmos_im
!
! ----------------------------------------------------------------------
!  2. Set STASHflags to activate diagnostics this timestep
!++SAN DELETE - no stash in box model
!
!     section is section number
!     IE is item number within section
!     IL is item number within STASHlist
!     II is counter within given section/item sublist for repeated items
!     IM is cumulative active item number within section

!   Clear all STASHflags
!++SAN DELETE

!   Loop over all items in STASHlist, enabling flags for diagnostics
!   which are active this step --
!     note that atmosphere/ocean diagnostics must be discriminated
!++SAN DELETE - no STASH in box model


! ----------------------------------------------------------------------
!  3. If errors occurred in setting STASHflags, set error return code
!s IF (icode == 1) cmessage='SETTSCTL: STASH code error'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE settsctl
END MODULE settsctl_mod

