! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
MODULE inittime_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INITTIME_MOD'

CONTAINS
SUBROUTINE inittime(                                                           &
  submodel)
!
!  Routine: INITTIME -------------------------------------------------
!
!  Purpose: Initialises the model time relative to the calendar zero
!           time.  The basis time is converted to a time in seconds
!           since T=0 with respect to the calendar. If the model basis
!           time as specified via the history file does not match the
!           data time in the dump(s), an error is flagged and the
!           routine exits.
!           Also sets derived time information from supplied time
!           information:-
!           (a) Assimilation start/end timesteps;
!           (b) Interface field generation start/end steps and length;
!           (c) Real timestep(s) in seconds.
!
!++SAN simplified for UKCA box model
!
!  Programming standard: UM Doc Paper 3, version 8.2 (25/3/2009)
!
!  External documentation: On-line UM document C0 - The top-level
!                          control system
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE dump_headers_mod, ONLY: a_fixhd, a_lookup
USE submodel_mod, ONLY:                                                        &
    internal_model_list, n_internal_model, atmos_sm, atmos_im
USE nlstcall_mod, ONLY: model_basis_time,                                      &
                         model_analysis_hrs,                                   &
                         model_analysis_mins,                                  &
                         run_target_end,                                       &
                         l_fastrun, lcal360

USE conversions_mod, ONLY: isec_per_day

USE history, ONLY: model_data_time, h_stepim, newrun, original_basis_time

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE nlsizes_namelist_mod, ONLY: a_prog_lookup, len_fixhd !Added by SAN

USE model_time_mod, ONLY:                                                      &
    assim_extrastepsim, assim_firststepim, assim_stepsim, basis_time_days,     &
    basis_time_secs, data_minus_basis_hrs, i_day,                &
    i_day_number, i_hour, i_minute, i_month, i_second, i_year,                 &
    iau_dtresetstep, l_c360dy, previous_time, secs_per_stepim, stepim,         &
    target_end_stepim

USE tim2step_mod, ONLY: tim2step
USE time2sec_mod, ONLY: time2sec
USE sec2time_mod, ONLY: sec2time
USE add_period_to_date_mod, ONLY: add_period_to_date

USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

IMPLICIT NONE


INTEGER, INTENT(IN) ::   submodel     ! IN  - submodel dump identifier

!
!
!  Local variables
!
LOGICAL :: l_dt_error         ! Flag for Data Time error
INTEGER :: ical               ! 1: Gregorian, 2: 360 day calendar

INTEGER ::                                                                     &
  elapsed_days                                                                 &
                              ! Whole days from basis time to VT etc.
, elapsed_secs                                                                 &
                              ! Secs-in-day from basis time to VT etc.
, data_minus_basis_days                                                        &
                              ! Whole days from basis time to DT
, data_minus_basis_secs

INTEGER :: i                  ! Loop counter

INTEGER :: basis_year
INTEGER :: basis_month
INTEGER :: basis_day
INTEGER :: basis_hour
INTEGER :: basis_minute
INTEGER :: basis_second

INTEGER :: target_end_year
INTEGER :: target_end_month
INTEGER :: target_end_day
INTEGER :: target_end_hour
INTEGER :: target_end_minute
INTEGER :: target_end_second

! time2sec needs a common reference point of first day of chosen calendar
! ie. day number of 1st of January in first year.
INTEGER,PARAMETER :: day_one_days = 0
INTEGER,PARAMETER :: day_one_secs = 0

INTEGER ::                                                                     &
  im                                                                           &
                              ! internal model id index
, ii                                                                           &
                              ! internal model loop counter
, step                                                                         &
, steps_per_period                                                             &
, secs_per_period                                                              &
, a_steps_per_hr                                                               &
, target_end_step                                                              &
, h_step

INTEGER :: icode
LOGICAL :: zerobasis         ! =T: model_basis_time=0

CHARACTER(LEN=9) :: DumpCalStr
CHARACTER(LEN=9) :: ModelCalStr
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INITTIME'
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=1024) :: cMessageLong

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!
!----------------------------------------------------------------------
! 0. If at start of integration (ie. step numbers both at zero),
!     set model_basis_time from validity time of start dump if run is
!     an assimilation or pseudo-assimilation or if model_basis_time is
!     zero (time checking assumed between model time and observations).
!     Check that timestep definitions are valid.
!     Set MODEL_DATA_TIME from start dump data time in all cases.
!
! Check for new run (NRUN) or continuation run (CRUN)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!++SAN
WRITE(umMessage, *) ' INITTIME; SAN starting inittime'
CALL umPrint(umMessage,src='inittime-inittim1')
!--SAN

icode = 0

! initialse STEP to zero, as is used to set stepim(internal_model) below
step = 0

newrun=.TRUE.
DO ii=1,n_internal_model
  im=internal_model_list(ii)
  IF (h_stepim(im) /= 0) THEN    ! CRUN, continuation run
    newrun=.FALSE.
  END IF
END DO ! ii over internal models

! Check for model_basis_time zero
zerobasis=.TRUE.
DO i=1,6
  IF (model_basis_time(i) /= 0) THEN
    zerobasis=.FALSE.
  END IF
END DO

! Compute model_analysis_hrs from model_analysis_mins
model_analysis_hrs =  REAL(model_analysis_mins)/60.0

!++SAN
WRITE(umMessage, '(A, F10.5)') ' INITTIME; SAN model_analysis_hrs = ', model_analysis_hrs
CALL umPrint(umMessage,src='inittime-inittim1')
!--SAN

IF (newrun) THEN

  ! Store model basis time from newrun in history restart module
  original_basis_time(:)=model_basis_time(:)

  !++SAN a_fixhd not used, as we don't have a start dump for box model
  !      Set model_data_time to be model_basis_time
  DO i=1,6
    ! model_data_time(i)=a_fixhd(20+i)
    model_data_time(i)=model_basis_time(i)
  END DO
  
  WRITE(umMessage,'(A,6(I10,1x))')                                             &
      ' INITTIME; model_data_time= ',model_data_time
  CALL umPrint(umMessage,src='inittime-inittim1')
  WRITE(umMessage,'(A,6(I10,1x))')                                             &
      ' INITTIME; model_basis_time= ',model_basis_time
  CALL umPrint(umMessage,src='inittime-inittim1')

  !----------------------------------------------------------------------
  ! 1. Check model_basis_time against dump validity time in header(s)
  !    (first step only - skipped if assimilation or pseudo-assimilation)
  !++SAN DELETE

END IF

!----------------------------------------------------------------------
! 2. Check that model calendar matches the dump header(s)
!
IF (lcal360) THEN
  ical=2
  ModelCalStr = '360 day'
ELSE
  ical=1
  ModelCalStr = 'Gregorian'
END IF

!++SAN - no dump header, delete read from a_fixhd

!
!  Create copy of logical LCAL360 in model_time_mod
!
l_c360dy = lcal360
!
!----------------------------------------------------------------------
! 3. Initialise model time relative to first day of calendar, using
!    30 day month calendar for climate models, full calendar including
!    leap years for other models
!
basis_year   = model_basis_time(1)
basis_month  = model_basis_time(2)
basis_day    = model_basis_time(3)
basis_hour   = model_basis_time(4)
basis_minute = model_basis_time(5)
basis_second = model_basis_time(6)

CALL time2sec(basis_year,basis_month,basis_day,                                &
   basis_hour,basis_minute,basis_second,                                       &
   day_one_days,day_one_secs,basis_time_days,basis_time_secs,lcal360)


! 3.1 Set initial time day number in dump header(s) only at start
!     The i_year outputs from here are not needed later.
IF (h_stepim(atmos_im) == 0) THEN
  CALL sec2time(day_one_days,day_one_secs,basis_time_days,basis_time_secs,     &
     i_year,i_month,i_day,i_hour,i_minute,i_second,                            &
     i_day_number,lcal360)
  a_fixhd(27) = i_day_number
END IF
!
!----------------------------------------------------------------------
! 4. Initialise incremental step counter(s) to accord with basis time
!    and model restart time, and flag any difference which might occur
!    due to altering the timestep partway through an integration.
!
DO ii=1,n_internal_model
  im=internal_model_list(ii)
  secs_per_period =secs_per_periodim(im)
  steps_per_period=steps_per_periodim(im)
  h_step   = h_stepim(im)
  !++SAN at the moment, these are all set to zero
  !      Expect I should hardcode these to 0
  !      Might be a possibility we could read in the integration if we
  !      Set up the capacity to do restarts with the box model, but I'm
  !      not convinced that would be useful...
  IF (im == atmos_im) THEN
    i_year   = a_fixhd(28)
    i_month  = a_fixhd(29)
    i_day    = a_fixhd(30)
    i_hour   = a_fixhd(31)
    i_minute = a_fixhd(32)
    i_second = a_fixhd(33)
  END IF   ! atmos_im

  ! Calculate elapsed time of run according to dump header.
  CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,                 &
     day_one_days,day_one_secs,                                                &
     elapsed_days,elapsed_secs,lcal360)
  elapsed_days = elapsed_days - basis_time_days
  elapsed_secs = elapsed_secs - basis_time_secs
  IF (elapsed_secs < 0) THEN
    elapsed_secs = elapsed_secs + isec_per_day
    elapsed_days = elapsed_days - 1
  END IF

  CALL time2sec(model_data_time(1),model_data_time(2),                         &
     model_data_time(3),model_data_time(4),                                    &
     model_data_time(5),model_data_time(6),                                    &
     day_one_days,day_one_secs,                                                &
     data_minus_basis_days,data_minus_basis_secs,                              &
     lcal360)

  data_minus_basis_days = data_minus_basis_days - basis_time_days
  data_minus_basis_secs = data_minus_basis_secs - basis_time_secs
  IF (data_minus_basis_secs < 0) THEN
    data_minus_basis_secs = data_minus_basis_secs + isec_per_day
    data_minus_basis_days = data_minus_basis_days - 1
  END IF

  data_minus_basis_hrs= data_minus_basis_secs/3600+                            &
     24*data_minus_basis_days
  !
  !++SAN DELETE forecast code for box model

  stepim(im)    = step
  h_stepim(im)  = step
END DO ! ii over N_INTERNAL_MODEL

!----------------------------------------------------------------------
! 5. Set target end steps from target end time using
!     relative time convention
!

! Target_end date = basis date plus run_target_end date
target_end_year   = basis_year
target_end_month  = basis_month
target_end_day    = basis_day
target_end_hour   = basis_hour
target_end_minute = basis_minute
target_end_second = basis_second

CALL add_period_to_date(                                                       &
   target_end_year,target_end_month,target_end_day,                            &
   target_end_hour,target_end_minute,target_end_second,                        &
   run_target_end(1),run_target_end(2),run_target_end(3),                      &
   run_target_end(4),run_target_end(5),run_target_end(6),                      &
   lcal360)

!++SAN
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN target_end_year= ',target_end_year
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN target_end_month= ',target_end_month
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN target_end_day= ',target_end_day
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN target_end_hour= ',target_end_hour
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN target_end_minute= ',target_end_minute
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN target_end_second= ',target_end_second
CALL umPrint(umMessage,src='inittime')
!--SAN

! Convert the target_end date into the requested end date
! relative to basis date (elapsed date)
CALL time2sec(target_end_year,target_end_month,target_end_day,                 &
   target_end_hour,target_end_minute,target_end_second,                        &
   day_one_days,day_one_secs,elapsed_days,elapsed_secs,lcal360)

elapsed_days = elapsed_days - basis_time_days
elapsed_secs = elapsed_secs - basis_time_secs
IF (elapsed_secs < 0) THEN
  elapsed_secs = elapsed_secs + isec_per_day
  elapsed_days = elapsed_days - 1
END IF

IF (elapsed_days <  0) THEN
  icode=1
  cmessage="INITTIME: Negative run length requested"
  GO TO 9999
END IF
! This prevents an infinite resubmission cycle
IF (MOD(elapsed_secs,secs_per_period/steps_per_period) /= 0) THEN
  icode=1
  cmessage="INITTIME: Run length not integral no. of timesteps. See output"
  CALL umPrint('INITTIME: Run length does not divide into timesteps',          &
      src='inittime-inittim1')
  WRITE(umMessage,'(A,I10,A,I10,A)') 'run length ',elapsed_days,' days ',      &
     elapsed_secs,' seconds '
  CALL umPrint(umMessage,src='inittime-inittim1')
  WRITE(umMessage,'(A,I10)')'TIMESTEP ',secs_per_period/steps_per_period
  CALL umPrint(umMessage,src='inittime-inittim1')
  CALL umPrint('Modify run_target_end and resubmit',src='inittime-inittim1')
  GO TO 9999
END IF
!
DO ii=1,n_internal_model
  im=internal_model_list(ii)
  secs_per_period =secs_per_periodim(im)
  steps_per_period=steps_per_periodim(im)
  CALL tim2step(elapsed_days,elapsed_secs,                                     &
     steps_per_period,secs_per_period,target_end_step)
  target_end_stepim(im)=target_end_step
END DO ! ii over N_INTERNAL_MODEL
!----------------------------------------------------------------------
! 7. Set assimilation start timestep, length in steps, and overlap into
!     forecast from basic control information
!++SAN DELETE
!

!----------------------------------------------------------------------
! 7. If running in IAU mode, calculate step on which data time must
!    be reset. A_STEPS_PER_HR was set in the previous section.
!++SAN DELETE

!----------------------------------------------------------------------
! 9. Calculate length of timestep in seconds in model_time_mod
!
secs_per_stepim(atmos_sm) = REAL(secs_per_periodim(atmos_sm))/                 &
   REAL(steps_per_periodim(atmos_sm))
!----------------------------------------------------------------------
! 10. Set current time in model_time_mod according to submodel
!
!++SAN double check these are still zero'd out. if so, can delete the lines...
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN i_year= ',i_year
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN i_month= ',i_month
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN i_day= ',i_day
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN i_hour= ',i_hour
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN i_minute= ',i_minute
CALL umPrint(umMessage,src='inittime')
WRITE(umMessage,'(A, I10)')                                             &
      ' INITTIME; SAN i_second= ',i_second
CALL umPrint(umMessage,src='inittime')

IF (submodel == atmos_sm) THEN
  i_year   = a_fixhd(28)
  i_month  = a_fixhd(29)
  i_day    = a_fixhd(30)
  i_hour   = a_fixhd(31)
  i_minute = a_fixhd(32)
  i_second = a_fixhd(33)
END IF
!--SAN ...up to here
previous_time(1)=i_year
previous_time(2)=i_month
previous_time(3)=i_day
previous_time(4)=i_hour
previous_time(5)=i_minute
previous_time(6)=i_second
previous_time(7)=0        ! Not yet set

9999 CONTINUE
IF (icode  /=  0) THEN
  CALL ereport(routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!----------------------------------------------------------------------
END SUBROUTINE inittime
END MODULE inittime_mod

