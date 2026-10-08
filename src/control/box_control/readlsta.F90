! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Read run-time control namelist information configurable
!   in the gui, as required to set up parametrization constants and
!   logical switches needed by physics and dynamics schemes for the
!   Atmosphere model.
!
! SAN: Updated to just read UKCA relevant namelists
!
! Method:  Sequential read of namelists.
!          Optional print out of namelists
!          Check input namelists
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! System component: Control Atmos
!

! Subroutine Interface:
MODULE readlsta_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READLSTA_MOD'

CONTAINS
SUBROUTINE Readlsta()

!++SAN - Deleted out all unnecesary modules for UKCa-box

!S USE atmos_max_sizes,                  ONLY: model_levels_max
!++SAN - trying readding convection options - Needed to pass some of the UKCA checks

USE cv_run_mod,                       ONLY: l_dcpl_cld4pc2,                    &
read_nml_run_convection, print_nlist_run_convection, check_run_convection


USE ereport_mod,                      ONLY: ereport

USE errormessagelength_mod,           ONLY: errormessagelength

USE file_manager,                     ONLY: get_file_unit_by_id


!++SAN suspect I will need the general physics inputs?
! Module for general physics options
USE gen_phys_inputs_mod,              ONLY:                                    &
    read_nml_gen_phys_inputs,                                                  &
    gen_phys_inputs,                                                           &
    print_nlist_gen_phys_inputs

! module for GLOMAP_CLIM options
USE glomap_clim_option_mod,           ONLY:                                    &
    print_nlist_run_glomap_clim,                                               &
    read_nml_run_glomap_clim,                                                  &
    check_run_glomap_clim

!++SAN Delete?
!sUSE lw_rad_input_mod, ONLY: lw_input

USE model_domain_mod,                 ONLY:                                    &
    check_nml_model_domain,                                                    &
    read_nml_model_domain,                                                     &
    print_nlist_model_domain

!++SAN starting chemistry namelists?
USE ozone_inputs_mod,                 ONLY:                                    &
    print_nlist_run_ozone,                                                     &
    check_run_ozone,                                                           &
    read_nml_run_ozone

USE parkind1,                         ONLY:                                    &
    jprb,                                                                      &
    jpim

USE pc2_constants_mod,                ONLY: i_cld_pc2

USE planet_constants_mod,             ONLY:                                    &
    planet_constants,                                                          &
    set_planet_constants,                                                      &
    read_nml_planet_constants,                                                 &
    print_nlist_planet_constants

USE rad_input_mod,                    ONLY:                                    &
    check_run_radiation,                                                       &
    print_nlist_run_radiation,                                                 &
    read_nml_run_radiation,                                                    &
    l_radiation

! module for aerosol emissions options
USE run_aerosol_mod,                  ONLY:                                    &
    print_nlist_run_aerosol,                                                   &
    read_nml_run_aerosol

USE science_fixes_mod,   ONLY: read_nml_temp_fixes, print_nlist_temp_fixes,    &
                                   warn_temp_fixes

!s USE set_rad_steps_mod,                ONLY: set_a_radstep

USE sl_input_mod,                     ONLY:                                    &
    print_nlist_run_sl,                                                        &
    check_run_sl,                                                              &
    read_nml_run_sl

!sUSE sw_rad_input_mod, ONLY: sw_input

! module for Segments options (tuning_segments namelist)
USE tuning_segments_mod,              ONLY:                                    &
    check_tuning_segments,                                                     &
    print_nlist_tuning_segments,                                               &
    read_nml_tuning_segments

!++SAN UKCA namelist options!!!
! module for UKCA options
USE ukca_option_mod,                  ONLY:                                    &
    check_run_ukca,                                                            &
    print_nlist_RUN_ukca,                                                      &
    read_nml_run_ukca,                                                         &
    l_ukca

USE UM_ParCore,                       ONLY: mype

USE umPrintMgr,                       ONLY:                                    &
    PrintStatus,                                                               &
    PrStatus_Normal,                                                           &
    umPrint

USE yomhook,                          ONLY:                                    &
    lhook,                                                                     &
    dr_hook

! module for EasyAerosol options
!sUSE easyaerosol_option_mod,           ONLY:                                    &
!s    print_nlist_easyaerosol,                                                   &
!s    read_nml_easyaerosol

!s USE cv_set_dependent_switches_mod, ONLY: cv_set_dependent_switches

IMPLICIT NONE


! Local parameters:
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='READLSTA')

! Local scalars:
INTEGER ::                                                                     &
 ErrorStatus      ! Return code : 0 Normal Exit : >0 Error
CHARACTER(LEN=errormessagelength) ::                                           &
 CMessage         ! Error message if Errorstatus >0

INTEGER :: atmoscntl_unit  ! unit no. for ATMOSCNTL file
INTEGER :: shared_unit     ! unit no. for SHARED file

LOGICAL :: l_print_namelist = .FALSE.  ! print out namelist entries.

CHARACTER(LEN=50000) :: lineBuffer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! -------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0

! Retrieve the unit number of the main namelist files
atmoscntl_unit = get_file_unit_by_id("atmoscntl", handler="fortran")
shared_unit    = get_file_unit_by_id("shared",    handler="fortran")

! height values to split model levels into l/m/h cloud
! not an input at present.
!++SAN comment out for now...
!S CALL h_split_defaults()

IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(lineBuffer,'(A)') '******************** ' // RoutineName //            &
                          ': Atmosphere run-time constants *******************'
  CALL umPrint(lineBuffer,src=RoutineName)
  l_print_namelist = .TRUE.
END IF

! -----------------------------------------
! read in scientific fixes namelist
! controls what fixes are applied
! each of these fixes are anticipated
! to become the default code in the future
! -----------------------------------------

CALL read_nml_temp_fixes(shared_unit)

! JULES switches have been moved to a shared UM-JULES namelist, but have kept
! the calls to the routines here to keep all the temporary fix information
! together.
!++SAN delete JULES options

! Carbon options - SAN delete

! Coupling options - SAN deleted

! Model domain - SAN deleted

! Planet constants - SAN keep for now...
CALL read_nml_planet_constants(shared_unit)
CALL set_planet_constants()
! Print needs to be after set otherwise g etc look unset for Earth simulations.
IF (l_print_namelist) CALL print_nlist_planet_constants()

! Mineral dust modelling - SAN comment out for now...
! Initialise sizes in case not set in namelist
!sCALL dust_size_dist_initialise
!sCALL read_nml_run_dust(shared_unit)
!sIF (l_print_namelist) CALL print_nlist_RUN_Dust()
!sCALL dust_parameters_load
!sCALL dust_parameters_check

! FV-track for cyclone tracking diagnostics - SAN deleted

! GLOMAP_CLIM Sub-model - SAN comment out for now...
!s CALL read_nml_run_glomap_clim(shared_unit)
!s IF (l_print_namelist) CALL print_nlist_RUN_Glomap_Clim()

!++SAN Read in UKCA namelist !!!
! UKCA Sub-model
CALL read_nml_run_ukca(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_UKCA()

! Gravity wave drag physics - SAN delete

! Vera visibility scheme - SAN delete

! Murk aerosol physics - SAN delete

! Convection physics - SAN delete
!++SAN readd... Need convection options to run UKCA
CALL read_nml_run_convection(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Convection()
CALL check_run_convection()


! Boundary layer physics - SAN delete

! Large scale precipitation physics - SAN delete

! Radiation physics
!++SAN - some of this needed to set long lived gas tracer values
CALL read_nml_run_radiation(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Radiation()
CALL check_run_radiation()

! Large scale cloud physics - SAN delete

! Aerosol Modelling - SAN Keep
CALL read_nml_run_aerosol(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Aerosol()
! run_aerosol_check called in readsize because it need model_levels

! LAM configuration - SAN delete

! Ozone - SAN keep?
CALL read_nml_run_ozone(shared_unit)
IF (l_print_namelist) CALL print_nlist_run_ozone()
CALL check_run_ozone()

! Free tracers - SAN delete

! Energy correction physics - SAN delete

! Calc pmsl - SAN delete

! General physics - SAN keep?
CALL read_nml_gen_phys_inputs(shared_unit)
IF (l_print_namelist) CALL print_nlist_gen_phys_inputs()

! LBC options - SAN delete

! UM nudging - SAN delete

! Generalised integration and GCR dynamics - SAN delete

! Generalised integration and GCR dynamics - SAN delete

! Idealised Model if required - SAN delete

! Semi-Lagrangian advection dynamics - SAN delete

! Diffusion, divergence damping and filtering dynamics - SAN delete

! Call to COSP - SAN delete

! Diagnostic double call to radiation - SAN comment out for now...?
!sCALL read_nml_radfcdia(atmoscntl_unit)
!sIF (l_print_namelist) CALL print_nlist_RADFCDIA()

! EasyAerosol - SAN comment out for now???
!sCALL read_nml_easyaerosol(shared_unit)
!sIF (l_print_namelist) CALL print_nlist_easyaerosol()

!-------------------------------------------------------
! Set other dependent convective switches valid for whole run - SAN delete

! Set model timesteps per radiation timestep (valid for whole run), and
! the radiation timestep length - SAN delete
!s IF (l_radiation)  CALL set_a_radstep()

!-------------------------------------------------------
! Set up switches based on CASIM, valid for whole run.
! Only call this routine if CASIM is switched on.
! SAN delete

!-------------------------------------------------------
! Check that if L_GLOMAP_CLIM_RADAER is selected certain
! other switches are not set so that they conflict.
!-------------------------------------------------------
CALL check_run_glomap_clim()

!++SAN delete cld_vn stuff...

!++SAN comment out rad input stuff for now...

! Options for the shortwave radiation
!sCALL sw_input

! Options for the longwave radiation
!sCALL lw_input

!++SAN no idea what this stuff is, comment out for now...
!sCALL read_nml_clmchfcg(atmoscntl_unit)
!sCALL clmchfcg_rates()
!sIF (l_print_namelist) CALL print_nlist_clmchfcg()

!sCALL coradoca_defaults

! ACP, ACDIAG - SAN delete

! Read the JULES namelists - SAN delete

! Stochastic physics - SAN delete

! Electric physics - SAN delete

!-----------------------------------------------------
! Below is where we provide a location for inter namelist checks
!-----------------------------------------------------

! If UKCA is active, check UKCA logicals are consistent with physics inputs
! Needs to be after UKCA namelist is read and after gen_phys_inputs
IF (l_ukca) THEN
  CALL check_run_ukca()
END IF

! In the coupled model case, set alpham/c to ssalpham/c and dtice to ssdtice
! if l_ssice_albedo == T. Note, ssalpham/c, ssdtice accessed via rad_input_mod
! Do this after JULES namelist reads, as l_ssice_albedo is in jules_sea_seaice
!++SAN - delete

!--------------------------------------------------------

! Read IAU namelist: - SAN delete

! Segments
CALL read_nml_tuning_segments(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_tuning_segments()
CALL check_tuning_segments()

! Check error condition
IF (ErrorStatus >  0) THEN

  CALL Ereport(ModuleName//':'//RoutineName,ErrorStatus,Cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Readlsta
END MODULE readlsta_mod

