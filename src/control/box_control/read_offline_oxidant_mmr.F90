!    -------------------------------------------------------------------
!    Subroutine: READ_OFFLINE_OXIDANT_MMR_MOD                --------------
!
!    Purpose: Programme to read environmental parameters namelist
!             for UKCA box model.
!
!    Programming standard: @£($@*£$)
!
!    External documentation: @£*£@£%
!
!    -------------------------------------------------------------------
!
! Module for UKCA Box model namelist : ENVIRONMENT_PARS

MODULE read_offline_oxidant_mmr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_OFFLINE_OXIDANT_MMR'

!---------------------------------------------------------------------------
! Define namelist
!---------------------------------------------------------------------------

REAL    :: o3
REAL    :: oh
REAL    :: ho2
REAL    :: h2o2
REAL    :: no3

! Public procedures
PUBLIC :: read_offline_oxidant_mmr

NAMELIST/OFFLINE_OXIDANT_MMR/                                                      &
  o3, oh, ho2, h2o2, no3

PRIVATE :: OFFLINE_OXIDANT_MMR

! ----------------------------------------------------------------------

CONTAINS

SUBROUTINE read_offline_oxidant_mmr(icode, iomessage, ukca_box_nml,              &
  o3_in, oh_in, ho2_in, h2o2_in, no3_in)

USE filenamelength_mod, ONLY: filenamelength

USE missing_data_mod,       ONLY: rmdi, imdi

USE umPrintMgr,             ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

!=============================================================================

IMPLICIT NONE

INTEGER, INTENT(OUT)                           :: icode
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: iomessage
CHARACTER(LEN=filenamelength), INTENT(IN)      :: ukca_box_nml

!++SAN scalar versions of variables for reading in namelist input
REAL, INTENT(OUT)    :: o3_in
REAL, INTENT(OUT)    :: oh_in
REAL, INTENT(OUT)    :: ho2_in
REAL, INTENT(OUT)    :: h2o2_in
REAL, INTENT(OUT)    :: no3_in

! Local variables
INTEGER :: istatus
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_OFFLINE_OXIDANT_MMR'

!=============================================================================

! Initialise values
o3            = rmdi
oh            = rmdi
ho2           = rmdi
h2o2          = rmdi
no3           = rmdi

!++SAN Setting up to read from a namelist
IF (PrintStatus >= PrStatus_Oper) THEN
  WRITE(umMessage,*) 'SAN; READ_OFFLINE_OXIDANT_MMR: Opening box model namelist'
  CALL umPrint(umMessage, src='read_offline_oxidant_mmr')
  WRITE(umMessage,*) 'Namelist file called ', ukca_box_nml
  CALL umPrint(umMessage, src='read_offline_oxidant_mmr')
END IF

OPEN(67,FILE=TRIM(ADJUSTL(ukca_box_nml)), ACTION='READ', IOSTAT=istatus)

IF (istatus /= 0) THEN
  icode = 500
  iomessage = ' Error opening file on unit 67 from '//routinename
  WRITE(umMessage,*) iomessage
  CALL umPrint(umMessage,src='read_offline_oxidant_mmr')
  WRITE(umMessage,*) ' Filename = '//TRIM(ADJUSTL(ukca_box_nml))
  CALL umPrint(umMessage,src='read_offline_oxidant_mmr')
  WRITE(umMessage,*) ' IOstat =', istatus
  CALL umPrint(umMessage,src='read_offline_oxidant_mmr')
ELSE 
  READ (UNIT=67, NML=OFFLINE_OXIDANT_MMR)
END IF

CLOSE (67)

! Map variables
o3_in            = o3
oh_in            = oh
ho2_in           = ho2
h2o2_in          = h2o2
no3_in           = no3

RETURN
END SUBROUTINE read_offline_oxidant_mmr

!=============================================================================
END MODULE read_offline_oxidant_mmr_mod
