! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: INCRTIME -------------------------------------------------
!
!    Purpose: Increments the model time by one atmosphere timestep.
!             Also updates timestamps in dump LOOKUP headers of
!             PROGNOSTIC fields (diagnostic LOOKUP headers
!             are updated exclusively by STASH).
!++SAN - simplified for box model, no dump LOOKUP headers or diagnostics
!
!    Programming standard: UM Doc Paper 3
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level
MODULE incrtime_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INCRTIME_MOD'

CONTAINS
SUBROUTINE incrtime(internal_model)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE IAU_mod, ONLY:                                                             &
    l_iau

!++SAN - deleted a few dependencies so as to compile easier in box model
USE umPrintMgr, ONLY: PrintStatus, PrStatus_Oper, ummessage, umprint
USE dump_headers_mod, ONLY: a_fixhd, a_lookup
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE nlstcall_mod, ONLY: model_analysis_hrs,                                    &
                         model_analysis_mins,                                  &
                         l_fastrun, lcal360

USE history, ONLY: h_stepim

USE nlsizes_namelist_mod, ONLY: a_prog_lookup

USE model_time_mod, ONLY:                                                      &
    basis_time_days, basis_time_secs, data_minus_basis_hrs, forecast_hrs,      &
    i_day, i_day_number, i_hour, i_minute, i_month, i_second, i_year,          &
    previous_time, stepim

USE errormessagelength_mod, ONLY: errormessagelength

USE stp2time_mod, ONLY: stp2time
USE sec2time_mod, ONLY: sec2time

IMPLICIT NONE

!
!   Arguments
!
INTEGER ::    internal_model    ! IN : internal model identifier

!
! ----------------------------------------------------------------------
!  Local variables
!
INTEGER ::                                                                     &
    elapsed_days,                                                              &
                             ! Elapsed days  since basis time
    elapsed_secs,                                                              &
                             ! Elapsed secs  since basis time
    elapsed_days_prev,                                                         &
                             ! Elapsed days, end of previous step
    elapsed_secs_prev,                                                         &
                             ! Elapsed secs, end of previous step
    i                       ! Loop index
INTEGER ::                                                                     &
                             ! Local scalars of internal model
 step                                                                          &
                             !  arrays.
, steps_per_period                                                             &
, secs_per_period

! model_analysis_hrs replaced by model_analysis_mins in cntlall.h -
! Requires ELAPSED_HRS changed to REAL
REAL :: elapsed_hrs             ! Elapsed hours since basis time

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INCRTIME'

!
! ----------------------------------------------------------------------
!  1. General timestep, increment STEP by one and update atmos
!     elapsed seconds (integer) relative to basis time
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
step            =stepim(internal_model)
steps_per_period=steps_per_periodim(internal_model)
secs_per_period =secs_per_periodim(internal_model)

CALL stp2time(step,steps_per_period,secs_per_period,                           &
             elapsed_days_prev,elapsed_secs_prev)
step = step+1
CALL stp2time(step,steps_per_period,secs_per_period,                           &
             elapsed_days,elapsed_secs)

stepim(internal_model) = step
h_stepim(internal_model) = step

!
!  1.1 If integrating backwards, negate elapsed times.
!++SAN DELETE - no backwards integration in box model
!

!
!  1.2 Set FORECAST_HRS - the number of hours relative to the current
!                         data time.
!
model_analysis_hrs =  REAL(model_analysis_mins)/60.0
elapsed_hrs = REAL(elapsed_secs)/3600 + elapsed_days*24

IF ((l_fastrun .OR. l_iau) .AND.                                               &
    elapsed_hrs >= model_analysis_hrs) THEN
  ! At or beyond data time reset step:
  forecast_hrs = elapsed_hrs - model_analysis_hrs
ELSE
  ! Data time of initial dump still valid:
  forecast_hrs = elapsed_hrs - data_minus_basis_hrs
END IF

IF ( PrintStatus  >  PrStatus_Oper ) THEN
  CALL umPrint( '',src='incrtime' )
  WRITE(umMessage,'(A,I10)')   'incrtime: ELAPSED SECS       ',                &
      elapsed_secs
  CALL umPrint(umMessage,src='incrtime')
  WRITE(umMessage,'(A,F10.5)') 'incrtime: ELAPSED HRS        ',                &
      elapsed_hrs
  CALL umPrint(umMessage,src='incrtime')
  WRITE(umMessage,'(A,F10.5)') 'incrtime: FORECAST HRS       ',                &
      forecast_hrs
  CALL umPrint(umMessage,src='incrtime')
  WRITE(umMessage,'(A,F10.5)') 'incrtime: model_analysis_hrs ',                &
      model_analysis_hrs
  CALL umPrint(umMessage,src='incrtime')
  CALL umPrint('',src='incrtime')
END IF

! ----------------------------------------------------------------------
!  2. Convert elapsed seconds since basis time to calendar time/date
!
CALL umPrint('incrtime: SAN calling sec2time',src='incrtime')
CALL sec2time(elapsed_days_prev,elapsed_secs_prev,                             &
             basis_time_days,basis_time_secs,                                  &
             previous_time(1),previous_time(2),previous_time(3),               &
             previous_time(4),previous_time(5),previous_time(6),               &
             previous_time(7),lcal360)
CALL umPrint('incrtime: SAN calling sec2time',src='incrtime')
CALL sec2time(elapsed_days,elapsed_secs,                                       &
             basis_time_days,basis_time_secs,                                  &
             i_year,i_month,i_day,i_hour,i_minute,i_second,                    &
             i_day_number,lcal360)

!++SAN
WRITE(umMessage,'(A, I10)')                                             &
      ' incrtime; SAN i_year= ',i_year
CALL umPrint(umMessage,src='incrtime')
WRITE(umMessage,'(A, I10)')                                             &
      ' incrtime; SAN i_month= ',i_month
CALL umPrint(umMessage,src='incrtime')
WRITE(umMessage,'(A, I10)')                                             &
      ' incrtime; SAN i_day= ',i_day
CALL umPrint(umMessage,src='incrtime')
WRITE(umMessage,'(A, I10)')                                             &
      ' incrtime; SAN i_hour= ',i_hour
CALL umPrint(umMessage,src='incrtime')
WRITE(umMessage,'(A, I10)')                                             &
      ' incrtime; SAN i_minute= ',i_minute
CALL umPrint(umMessage,src='incrtime')
WRITE(umMessage,'(A, I10)')                                             &
      ' incrtime; SAN i_second= ',i_second
CALL umPrint(umMessage,src='incrtime')
!--SAN

! ----------------------------------------------------------------------
!  3. Copy date/time information into the dump header and update
!     VALIDITY TIME and FORECAST PERIOD in prognostic field LOOKUP
!     headers.
!
a_fixhd(28) = i_year
a_fixhd(29) = i_month
a_fixhd(30) = i_day
a_fixhd(31) = i_hour
a_fixhd(32) = i_minute
a_fixhd(33) = i_second
a_fixhd(34) = i_day_number

CALL umPrint('incrtime: SAN Leaving incrtime',src='incrtime')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
! ----------------------------------------------------------------------
END SUBROUTINE incrtime
END MODULE incrtime_mod
