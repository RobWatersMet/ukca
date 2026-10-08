! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: UKCA_SHELL -------------------------------------------------
!
!  Purpose: Control routine for the Box Model.
!           Acquires size information needed for dynamic allocation of
!           configuration-dependent arrays and calls U_MODEL (the
!           master control routine) to allocate the arrays and perform
!           the top-level control functions and timestepping.
!
!  SAN - Stripped down to only call UKCA code
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level
MODULE ukca_shell_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SHELL_MOD'

CONTAINS
SUBROUTINE ukca_shell

!++SAN - add version of UM_config
USE ukca_box_config, ONLY:                                                         &
    appInit,                                                                   &
    appTerminate,                                                              &
    exe_UKCA
USE ukca_config_specification_mod, ONLY: ukca_config

USE mpl, ONLY:                                                                 &
    mpl_max_processor_name,                                                    &
    mpl_thread_single
!s    mpl_thread_multiple,                                                       &
!s    mpl_thread_serialized,                                                     &
!s    mpl_thread_funneled,                                                       &

USE run_info, ONLY: set_start_time,start_time

!$ USE omp_lib, ONLY: omp_get_num_threads, openmp_version
           ! Note OpenMP sentinel

USE filenamelength_mod, ONLY:                                                  &
    filenamelength

USE file_manager, ONLY:                                                        &
    assign_file_unit, release_file_unit, um_file_type, init_file_loop

USE yomhook, ONLY: lhook, dr_hook
USE drhook_control_mod, ONLY: drhook_control_enable, drhook_control_disable

USE IOS_Constants,      ONLY:                                                  &
    IOS_maxServers
USE ereport_mod, ONLY: ereport,ereport_finalise
USE UM_ParCore, ONLY: mype, nproc_max, nproc_um_npes=>nproc
USE umPrintMgr, ONLY: Printstatus, ummessage, umprint, PrStatus_Diag,          &
prstatus_normal, str, umprinterror

USE gcom_mod, ONLY: gc_alltoall_multi, gc_alltoall_version

USE nlsizes_namelist_mod, ONLY:                                                &
    global_row_length, global_rows, max_intf_model_levels,                     &
    max_lbcrow_length, max_lbcrows, model_levels, n_intf_a,                    &
    river_row_length, river_rows, row_length, rows,                            &
!++SAN - added more key variables to be given hard coded values
    theta_field_size, n_rows, model_levels, bl_levels, tr_levels,              &
    sm_levels, ozone_levels, ntiles

USE errormessagelength_mod, ONLY: errormessagelength
USE get_env_var_mod,   ONLY: get_env_var
USE io_configuration_mod, ONLY: print_runtime_info

USE um_submodel_init_mod, ONLY: um_submodel_init

USE get_wallclock_time_mod, ONLY: get_wallclock_time

USE timer_mod, ONLY: timer

USE readcntl_mod, ONLY: readcntl
USE readlsta_mod, ONLY: readlsta

!++SAN - only call the UKCA routines
USE ukca_option_mod, ONLY: l_ukca
!s USE atmos_ukca_setup_mod, ONLY: atmos_ukca_setup
!++SAN setup replaced with version for the box model
USE box_ukca_setup_mod, ONLY: box_ukca_setup

!++SAN CALL for box_model (replaces u_model_4A)
USE box_model_mod, ONLY: box_model

!++SAN temporary include for print statements:
USE nlstcall_mod, ONLY: model_basis_time

IMPLICIT NONE

!
!  Local parameters
!
CHARACTER(LEN=*) :: RoutineName
PARAMETER (RoutineName = 'UKCA_SHELL')
!
!  Local variables
!
INTEGER :: icode       ! Work - Internal return code
INTEGER :: istatus     ! RETURN STATUS FROM OPEN

CHARACTER(LEN=filenamelength) :: shared_filename    = "dummy filename" !Namelist
CHARACTER(LEN=filenamelength) :: atmoscntl_filename = "dummy filename" !Namelist
CHARACTER(LEN=errormessagelength) :: cmessage ! Work - Internal error message
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER :: atm_nprocx          ! number of procs EW for atmosphere
INTEGER :: atm_nprocy          ! number of procs NS for atmosphere
INTEGER :: length              ! length of env var contents
INTEGER :: err                 ! error return from subroutine calls
INTEGER :: atmoscntl_unit      ! Unit to use for ATMOSCNTL file
INTEGER :: shared_unit         ! Unit to use for SHARED file

CHARACTER(LEN=10) :: c_thread           ! to get nproc_x and nproc_y from

CHARACTER(LEN=8) :: ch_date2   !  Date returned from date_and_time
CHARACTER(LEN=10) :: ch_time2  !  Time returned from date_and_time

! Variables for IO Server setup
LOGICAL                      :: isIOServer
INTEGER                      :: numIOServers
CHARACTER(LEN=32)            :: c_io_pes
CHARACTER(LEN=10)            :: thread_level_setc
INTEGER                      :: thread_level_set

CHARACTER(LEN=mpl_max_processor_name) :: env_myhost
INTEGER                               :: env_myhost_len

REAL :: time_end_run


!-----------------------------------------------------------------------

cmessage = ' '
numIOServers = 0  ! Initialise IO server variable

! ----------------------------------------------------------------------
!----------------------------------------------------------------------
! 1.0 Initialise Message Passing Libraries
!
!++SAN checking icode going in
WRITE(umMessage,*) 'SAN UKCA_SHELL: SAN icode = ',icode
CALL umPrint(umMessage, src='readcntl')
!--SAN

CALL umPrint( 'SAN Initialising IO',src='ukca_shell')

! Get the atmosphere decomposition
! Cannot call get_env_var before GCOM initialisation.
CALL GET_ENVIRONMENT_VARIABLE('UM_THREAD_LEVEL',c_thread,length,err)
IF (err  /=  0 .OR. length == 0) THEN
  CALL umPrint('Warning: Environment variable UM_THREAD_LEVEL has ' //         &
      'not been set.',src='um_shell')
  CALL umPrint('Setting thread_level to multiple',src='um_shell')
  thread_level_setc = 'MULTIPLE'
ELSE
  READ(UNIT=c_thread,FMT='(A10)') thread_level_setc
END IF

!++SAN - For now, hardcode thead level to single,
! because only one thread on the box model
!s SELECT CASE (thread_level_setc)
!sCASE ('SINGLE')
CALL umPrint( 'SAN Setting to MPL thread single',src='ukca_shell')
  thread_level_set = mpl_thread_single
!s CASE DEFAULT
!s  WRITE(umMessage,'(A,A,A)') 'Warning: Thread level ', thread_level_setc,      &
!s       ' not recognised, setting to MULTIPLE.'
!s  CALL umPrint(umMessage,src='um_shell')
!s  thread_level_set = mpl_thread_multiple
!sEND SELECT

! The total number of processors required (nproc_max) is determined
!  by a call to gc_init/gc_init_thread:

!++SAN don't need to worry about calling OASIS
!   Standard UM GCOM initialisation when OASIS is not used
CALL gc_init_thread(mype,nproc_max, thread_level_set)

  !   Permit calls to DrHook, then call the first top-level caliper.
!s  CALL drhook_control_enable()

!s IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! initialise timers after MPI initialisation
CALL set_start_time()

CALL umPrint('SAN I am PE '//TRIM(str(mype)),src='ukca_shell')
CALL umPrint( 'SAN Initialising APP exe_UKCA',src='ukca_shell')
!s CALL appInit(exe_UM)
CALL appInit(exe_UKCA)
CALL umPrint( 'SAN Out of appInit',src='ukca_shell')

IF (print_runtime_info) THEN
!s  IF ( L_print_pe .OR. mype == 0 ) THEN
  IF (mype == 0 ) THEN
    CALL umPrint( '',src='ukca_shell')
    CALL umPrint('ukca_shell: Info: Start time set ',src='ukca_shell')
    CALL umPrint( '',src='um_shell')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

! Determine the number of processors the model has been configured to run on:
CALL umPrint( 'SAN Setting the number of processors',src='ukca_shell')
!++SAN - Hardcode to uce a 1x1 grid
atm_nprocx=1
atm_nprocy=1


! Calculate total number of atmos processors required:
nproc_um_npes = atm_nprocx * atm_nprocy

! Get number of I/O PEs from the environment
CALL get_env_var('FLUME_IOS_NPROC',c_io_pes,allow_missing=.TRUE.,length=length)
IF (length < 0) THEN
      ! If not specified, try to work out a valid number
  numIOServers=nproc_max-nproc_um_npes
  icode=-10
  WRITE(cmessage,'(A,A,I4)')                                                   &
      'FLUME_IOS_NPROC environment variable not set',                          &
      ', I/O PE count set to ',numIOServers
  CALL ereport(routinename,icode,cmessage)
ELSE
  READ (UNIT=c_io_pes,FMT='(I5)') numIOServers
END IF

IF ( numIOServers < 0 .OR. numIOServers > IOS_maxServers ) THEN
  icode=-10
  WRITE(cmessage,'(A,I4)')                                                     &
      'I/O PE count is outside allowed range: ',numIOServers
  CALL ereport(routinename,icode,cmessage)

      ! try to work out a valid number
  numIOServers=nproc_max-nproc_um_npes
  IF ( numIOServers < 0 .OR. numIOServers > IOS_maxServers ) THEN
    icode=10
    WRITE(cmessage,'(A)')                                                      &
        'A valid I/O PE count could not be set'
    CALL ereport(routinename,icode,cmessage)
  END IF
ELSE
  CALL umPrint('Enabling '//TRIM(str(numIOServers))//' I/O PEs',               &
      pe=0,src='ukca_shell', level=-PrintStatus)
END IF

IF (nproc_max  <   0) THEN
  CALL umPrint( 'Parallel initialisation failed',src='um_shell')
  GO TO 9999
ELSE

  ! Set GCOM to use the alternative version of RALLTOALLE
  ! throughout the run
  CALL gc_setopt(gc_alltoall_version, gc_alltoall_multi, err)

  CALL umPrint('',src='ukca_shell')
  CALL umPrint('**************************** PROCESSOR '//                     &
      'INFORMATION ****************************',src='ukca_shell')
  CALL umPrint('',src='um_shell')
  CALL umPrint(TRIM(str(nproc_max))//' Processors initialised.',               &
      src='um_shell')
  CALL MPL_Get_processor_name(env_myhost, env_myhost_len, err)
  IF (err /= 0) THEN
    CALL umPrint('I am PE '//TRIM(str(mype)),src='ukca_shell')
  ELSE
    CALL umPrint('I am PE '//TRIM(str(mype))//' on '//TRIM(env_myhost),        &
        src='ukca_shell')
  END IF
  ! Only want OpenMP section executing if OpenMP is compiled in,
  ! so protect by sentinal
!$OMP PARALLEL DEFAULT(NONE) SHARED(PrintStatus)
!$OMP MASTER
!$  WRITE(umMessage,'(A,I2,A)') 'I am running with ',                          &
!$      omp_get_num_threads(),' thread(s).'
!$  CALL umPrint(umMessage,src='ukca_shell', level=-PrintStatus)
!$  WRITE(umMessage,'(A,I6)') 'OpenMP Specification: ',openmp_version
!$  CALL umPrint(umMessage,src='ukca_shell', level=-PrintStatus)
!$OMP END MASTER
!$OMP END PARALLEL
END IF

#if defined(IBM_XL_FORTRAN)
! On IBM force buffering of Fortran I/O for initialisation
CALL setrteopts('buffering=enable')
#endif
!
CALL timer(RoutineName,1,starttime=Start_time)

!----------------------------------------------------------------------
!
!    Open the two main namelist files on PE 0 only.
!    All runtime control variables are subsequently read in from here.
!
CALL get_env_var("ATMOSCNTL",atmoscntl_filename)
CALL assign_file_unit(atmoscntl_filename, atmoscntl_unit, handler="fortran",   &
                      id="atmoscntl")

IF (mype == 0) THEN
  OPEN(UNIT=atmoscntl_unit,FILE=atmoscntl_filename, ACTION='READ',             &
       IOSTAT=istatus, IOMSG=iomessage)

  IF (istatus /= 0) THEN
    icode=500
    CALL umPrint( ' ERROR OPENING ATMOSPHERE-ONLY NAMELIST FILE',src='ukca_shell')
    WRITE(umMessage,'(A,A)') ' FILENAME =',TRIM(atmoscntl_filename)
    CALL umPrint(umMessage,src='ukca_shell')
    WRITE(umMessage,'(A,I0)') ' IOSTAT =',istatus
    CALL umPrint(umMessage,src='ukca_shell')
    WRITE(umMessage,'(A,A)') ' IOMSG =',TRIM(iomessage)
    CALL umPrint(umMessage,src='ukca_shell')
    GO TO 9999
  END IF
END IF

CALL get_env_var("SHARED_NLIST",shared_filename)
CALL assign_file_unit(shared_filename, shared_unit, handler="fortran",         &
                      id="shared")

IF (mype == 0) THEN
  OPEN(UNIT=shared_unit,FILE=shared_filename, ACTION='READ', IOSTAT=istatus,   &
       IOMSG=iomessage)

  IF (istatus /= 0) THEN
    icode=510
    CALL umPrint( ' ERROR OPENING SHARED NAMELIST FILE',src='um_shell')
    WRITE(umMessage,'(A,A)') ' FILENAME =',TRIM(shared_filename)
    CALL umPrint(umMessage,src='ukca_shell')
    WRITE(umMessage,'(A,I0)') ' IOSTAT =',istatus
    CALL umPrint(umMessage,src='ukca_shell')
    WRITE(umMessage,'(A,A)') ' IOMSG =',TRIM(iomessage)
    CALL umPrint(umMessage,src='ukca_shell')
    GO TO 9999
  END IF
END IF
! ------------------------------------------------------------------
!  0.1 Get submodel/internal model components of model run.
!

icode=0
CALL UM_Submodel_Init(icode)
IF (icode  /=  0) THEN
  cmessage = 'Error calling UM_Submodel_init'
  GO TO 9999
END IF

CALL ukca_Shell_banner('Start')

!----------------------------------------------------------------------
!++SAN no readhist for now in box model
  !  Read history files for NRUN or CRUN.
!s  CALL readhist ( icode,cmessage )
!s  IF (icode > 0) GO TO 9999

!++SAN Need to run readcntl in order to read in namelist settings
CALL umPrint('SAN UKCA_SHELL: calling readcntl',src='ukca_shell')
!--SAN

!  Read Control file on standard input.
!
CALL readcntl ( icode,cmessage )

IF (icode  >   0) GO TO 9999

!----------------------------------------------------------------------

!  Call READLSTA to read namelists to control atmosphere integration
!  and diagnostic point print.
!++SAN Think we need this to read UKCA namelists
CALL umPrint('SAN UKCA_SHELL: calling readlsta',src='ukca_shell')

CALL readlsta()

CALL umPrint('********************************************'//                &
      '***********************************',src='ukca_shell')

!----------------------------------------------------------------------
!  1.1 Get configuration-dependent sizes needed for dynamic allocation.
!

! Decompose atmosphere data and find new local data size
!++SAN for now, hardcoding these settings
CALL umPrint('ukca_shell SAN: Hardcoding model array sizes',src='ukca_shell')
global_row_length = 1
global_rows       = 1
model_levels      = 1
river_rows        = 1
river_row_length  = 1
row_length        = 1
rows              = 1
!++SAN additional global variables added to also be defined in this routine
theta_field_size  = row_length*rows
n_rows            = 1
model_levels      = 1
bl_levels         = 1
tr_levels         = 1
sm_levels         = 1
ozone_levels      = 1
ntiles            = 0

n_intf_a          = 1
max_intf_model_levels = 1
max_lbcrow_length = 1
max_lbcrows       = 1



!-----------------------------------------------------------------------
! 1.1.1 Set up the UKCA model if required.
! This step is dependent on the domain decomposition and is required
! prior to STASH request processing.
! Processing in ukca_setup includes checking that UKCA logicals are
! consistent and setting of internal UKCA values from UKCA and other
! namelists.
CALL umPrint('SAN UKCA_SHELL Section 1.1.1.',src='ukca_shell')
!++SAN atmos_ukca_setup replaced with call to box_ukca_setup for box model
IF (l_ukca) THEN
  CALL umPrint('SAN Calling box_ukca_setup', src='ukca_shell', level=-PrintStatus)
  CALL box_ukca_setup()
END IF

!++SAN check config options
CALL umPrint(umMessage, src='ukca_shell')
! Hardcode options that are needed to be set ot a particular value to run
! the box model version of UKCA
IF (ukca_config%l_environ_z_top) THEN
  WRITE(umMessage,*) 'SAN UKCA_SHELL: l_environ_z_top must be FALSE in box model'
  CALL umPrint(umMessage, src='ukca_shell')
  ukca_config%l_environ_z_top = .FALSE.
  WRITE(umMessage,*) 'SAN UKCA_SHELL: l_environ_z_top = ', ukca_config%l_environ_z_top
  CALL umPrint(umMessage, src='ukca_shell')
END IF
!--SAN


!-----------------------------------------------------------------------
! 1.2 Call STASH_PROC: top level control routine for processing of
!                      STASH requests and STASH addressing.
!++SAN - avoid STASH requests for now...

! ----------------------------------------------------------------------
!  1.3 Calculate addresses of super arrays passed down for dynamic 
!      allocation.
!++SAN - avoid d1 allocation...

! ----------------------------------------------------------------------
!  2. Call U_MODEL_4A master routine to allocate the main data arrays
!     and do the calculations.
!
!++SAN - this is important next stage - replaced previous call to u_model_4a
!++ with call to box_model, which in turn calls box_ukca_mod (replacement
!++ of atmos_ukca_mod). These are designed to do the minimum setup needed
!++ for one timestep of the UKCA box model to run. box_ukca_mod calls ukca_step,
!++ which is itself linked to ukca_main1, and from there the UKCA model runs
!++ as it does for the full 3D model.
CALL umPrint('SAN UKCA_SHELL Calling box_model()',src='ukca_shell')
CALL box_model ()

9999  CONTINUE

! Namelist files were only open on PE 0, so should only be closed on such.
!++SAN Need to run readcntl in order to read in namelist settings
WRITE(umMessage,*) 'SAN UKCA_SHELL: close atmoscntl_unit and shared unit'
CALL umPrint(umMessage,src='ukca_shell')
!--SAN
IF (mype == 0) THEN
  CLOSE(atmoscntl_unit)
  CLOSE(shared_unit)
END IF
CALL release_file_unit(atmoscntl_unit, handler="fortran")
CALL release_file_unit(shared_unit, handler="fortran")

CALL ukca_Shell_banner('End')

IF (icode /= 0) THEN
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

IF (print_runtime_info) THEN
  IF ( mype == 0 ) THEN
    CALL umPrint( '',src='ukca_shell', level=-PrintStatus)
    time_end_run = get_wallclock_time()
    WRITE(umMessage,'(A,A,F10.3,A)')                                           &
      'ukca_shell: Info: End model run',                                         &
      ' at time=',time_end_run - Start_time,' seconds'
    CALL umPrint(umMessage,src='ukca_shell', level=-PrintStatus)
    CALL umPrint( '',src='ukca_shell', level=-PrintStatus)
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

CALL timer(RoutineName,2)

!++SAN - work out which of these needed for box model
!s CALL Halo_Exchange_Finalise()

! End the memory usage collection and print the final report.
!s CALL memory_usage_collect(sync_arg = .TRUE.)
!s CALL memory_usage_report()
!s CALL memory_usage_fini()

!s CALL appTerminate()

! Final top-level DrHook caliper, then prevent any further calls to DrHook.
!s IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!s CALL drhook_control_disable()

!++SAN
CALL umPrint('SAN - END UKCA_SHELL', src='ukca_shell', level=-PrintStatus)
!--SAN

RETURN

CONTAINS

SUBROUTINE ukca_shell_banner(stampname)

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) ::  stampname
IF (mype == 0) THEN
  CALL DATE_AND_TIME(ch_date2, ch_time2)
  CALL umPrint('',src='ukca_shell', level=-PrintStatus)
  CALL umPrint('********************************************'//                &
                '***********************************',src='um_shell',          &
         level=-PrintStatus)
  WRITE(umMessage,'(23A)')                                                     &
         '**************** ',stampname,' of UM RUN Job : ',                    &
         ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),                    &
          ' on ',                                                              &
          ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),                   &
          ' *****************'
  CALL umPrint(umMessage,src='ukca_shell', level=-PrintStatus)
!s  WRITE(umMessage,'(3A)')                                                      &
!s        '**************** Based upon UM release vn', um_version_char,          &
!s        '             *****************'
!s  CALL umPrint(umMessage,src='ukca_shell', level=-PrintStatus)
  CALL umPrint('*****************************************'//                   &
               '**************************************',src='ukca_shell',        &
               level=-PrintStatus)
  CALL umPrint('',src='ukca_shell', level=-PrintStatus)
END IF

END SUBROUTINE ukca_Shell_banner

END SUBROUTINE ukca_shell
END MODULE ukca_shell_mod

