! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: Provide top level interface to model/application initialisation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

!++SAN - Now providing top_level interface for UKCA
MODULE ukca_box_config

! DEPENDS ON: exceptions

! for dependencies
USE f_shum_thread_utils_mod           ! needed by FCM build system in this form
USE um_types, ONLY: integer_64
USE io, ONLY: ioInit, ioShutDown
USE application_description_runtypes, ONLY:                                    &
!s    exe_UM, exe_RCF, exe_scm,                                                  &
    exe_UM, exe_RCF, exe_scm, exe_UKCA,                                        &
    exe_combine, exe_merge,                                                    &
    exe_hreset,                                                                &
    exe_setup, exe_hprint,                                                     &
    exe_pptoanc, exe_pickup,                                                   &
    exe_crmstyle_coarse_grid
USE application_description, ONLY:                                             &
    setApplicationDesc,                                                        &
    getExeType,                                                                &
    isParallel,                                                                &
    isSmallExec
USE app_banner, ONLY:                                                          &
    reportApplication
USE umPrintMgr, ONLY: PrintStatus, PrDiag, umPrintExceptionHandler,            &
                      umPrintSetTarget, umPrintSetLevel,   &
                      umPrintFinalise, & !++SAN - adding umprint
                      umprint
!++SAN UM vn12.0 uses this, likely not needed for box model:
!s #if !defined(PPTOANC)
USE umPrintMgr_nml_mod, ONLY: umPrintLoadOptions
!s #endif

USE UM_ParCore, ONLY:                                                          &
    mype,                                                                      &
    nproc_max
USE file_manager, ONLY: init_file_manager
USE fort2c_exceptions_interfaces, ONLY:                                        &
    signal_set_verbose, signal_rank, signal_command,                           &
    signal_trap, signal_add_handler, signal_traceback, signal_off
IMPLICIT NONE

CONTAINS

SUBROUTINE appInit(exe)

USE file_manager, ONLY: assign_file_unit

IMPLICIT NONE
INTEGER, INTENT(IN) :: exe

INTEGER :: f_unit

! Note that the UM gc initialisation is complex due to
! (a) threading
! (b) oasis
! (c) flume
!
! ... So we will only init gcom for small execs, the UM must do
! its own thing before calling appInit.

!++SAN - 24/02/2021 Crashes at some point in this routine
! Adding diagnostic print statements
! These will be deleted as and when issues are fixed
CALL umPrint( '',src='ukca_box_config')
CALL umPrint( 'SAN - in appInit',src='ukca_box_config')
!--SAN

!++SAN UKCA-Box is a simple exe, so this should be fine
IF (exe/=exe_UM) THEN
  CALL umPrint( 'SAN Calling GC_init',src='ukca_box_config')

  CALL gc_init(' ',mype, nproc_max)
END IF

CALL umPrint( '',src='ukca_box_config')
CALL umPrint( 'SAN calling init_file_manager',src='ukca_box_config')
! Init file manager
CALL init_file_manager()

! Init output management, parallel executables (except SCM) will redirect
! output, other execs will stick with stdout
!++SAN Does not apply to UKCA-Box - seriel job
!++SAN adding the case for UKCA-Box as well, so we can use the same output streams as UM...
IF (exe==exe_UM  .OR.                                                          &
    exe==exe_RCF .OR.                                                          &
    exe==exe_crmstyle_coarse_grid .OR.                                         &
    exe==exe_UKCA) THEN !++SAN
  CALL umPrint( 'SAN calling assign_file_unit',src='ukca_box_config')
  CALL assign_file_unit("STDOUT Stream",                                       &
        f_unit, handler="fortran", id="stdout_reserved_unit")
  CALL umPrint( 'SAN calling umPrintSetTarget',src='ukca_box_config')
  CALL umPrintSetTarget(out_unit=f_unit)
END IF

! Init app registry
CALL umPrint( 'SAN setApplicationDesc',src='ukca_box_config')
CALL setApplicationDesc(exe)
!s CALL umPrint( 'SAN umPrintLoadOptions',src='ukca_box_config')
!s #if !defined(PPTOANC)
CALL umPrintLoadOptions()
!s #endif
CALL umPrint( 'SAN umPrintSetLevel',src='ukca_box_config')
CALL umPrintSetLevel()

! Init exception handling
CALL umPrint( 'SAN umSetApplicationExceptions',src='ukca_box_config')
CALL umSetApplicationExceptions(exe)
CALL umPrint( 'SAN reportApplication',src='ukca_box_config')
CALL reportApplication()

! Initialise IO
CALL umPrint( 'SAN ioInit',src='ukca_box_config')
CALL ioInit()

CALL umPrint( 'SAN Successful completion of appInit!',src='ukca_box_config')

END SUBROUTINE appInit

SUBROUTINE appTerminate()
USE ereport_mod, ONLY: ereport_finalise
USE fort2c_exceptions_interfaces, ONLY: signal_unregister_callbacks
IMPLICIT NONE

CALL ioShutdown()

CALL ereport_finalise()

CALL umPrintFinalise()

IF (getExeType()/=exe_um) THEN
  CALL gc_exit()
END IF

! free memory structures for callbacks added with signal_add_handler()
CALL signal_unregister_callbacks()

END SUBROUTINE appTerminate

! Abort handler for um_abort_mod to call
SUBROUTINE gcom_signal_abort(errcode)

USE um_types, ONLY: integer_64
USE UM_ParCore, ONLY: mype, nproc_max
USE fort2c_exceptions_interfaces, ONLY: signal_controlled_exit

IMPLICIT NONE

INTEGER, INTENT(IN) :: errcode

! DEPENDS ON: exceptions
CALL signal_controlled_exit(INT(errcode, integer_64))

CALL GC_Abort(mype, nproc_max, "um_abort called")

END SUBROUTINE gcom_signal_abort

! turn on signal handling accordingly
SUBROUTINE umSetApplicationExceptions(exe)

USE um_abort_mod, ONLY: set_abort_handler

IMPLICIT NONE
INTEGER, INTENT(IN)     :: exe
INTEGER(KIND=integer_64) :: trapping_option
CHARACTER(LEN=256)      :: command

! All of the signal_* C interfaces use int64_t for integer arguments.

IF ( exe==exe_UM .OR. exe==exe_RCF ) THEN
  trapping_option=signal_traceback ! no core, but traceback
ELSE
  trapping_option=signal_off       ! off, completely
END IF

IF (PrintStatus>=PrDiag) CALL signal_set_verbose()

! Tell the signal handler the mpi rank, so that it can tag output.
CALL signal_rank(INT(mype,integer_64))

! Tell the signal handler the program name. Not all fortran can do this
! however it is not essential.
CALL GET_COMMAND_ARGUMENT(0,command)
CALL signal_command(TRIM(command),INT(LEN_TRIM(command),integer_64))

! Set up OS exception handling
! Must be in a serial region because some c calls used are not reentrant.
CALL signal_trap(trapping_option)

    ! Note that registered handlers will be called both on applciation
    ! exceptions (if specified by signal_trap), and on 'clean' shutdowns
    ! via ereport (or umPrintError if from umPrintMgr) - all exit should go
    ! through one of these.
CALL signal_add_handler(umPrintExceptionHandler)

! Set up abort function
CALL set_abort_handler(gcom_signal_abort)

END SUBROUTINE umSetApplicationExceptions
END MODULE ukca_box_config

