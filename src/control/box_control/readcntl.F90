! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!++SAN deleted unnecesary dependencies for box model
!
! Subroutine Interface:
MODULE readcntl_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READCNTL_MOD'

CONTAINS
SUBROUTINE readcntl(icode,cmessage)

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE filenamelength_mod, ONLY:                                                  &
          filenamelength

! MPP - SAN DELETED

USE io, ONLY: setpos, buffin
USE io_constants, ONLY: ioOpenReadOnly, ioNoDelete
!s USE ancilcta_namelist_mod, ONLY: read_nml_ancilcta, check_nml_ancilcta,        &
!s                                  print_nlist_ancilcta
USE model_id_mod, ONLY: itab, configid, read_nml_configid, check_configid
USE nlstgen_mod, ONLY: ppxm, read_nml_nlstcgen
USE check_nlstcgen_mod, ONLY: check_nlstcgen
USE missing_data_mod, ONLY: imdi
USE nlstcall_mod, ONLY: read_nml_nlstcall, print_nlist_nlstcall,               &
                        lcal360,  model_basis_time
                        
!s USE nlstcall_pp_namelist_mod, ONLY: read_nml_nlstcall_pp
!s USE nlstcall_nc_namelist_mod, ONLY: read_nml_nlstcall_nc,                      &
!s                                     read_nml_nlstcall_nc_options,              &
!s                                     l_netcdf
!++SAN new namelist for box model I/O options
USE nlstcall_box_namelist_mod, ONLY: read_nml_nlstcall_box_options

USE file_manager, ONLY: get_file_unit_by_id, assign_file_unit,                 &
                        release_file_unit, um_file_type,                       &
                        get_file_by_unit

USE nlsizes_namelist_mod, ONLY:                                                &
    len_fixhd

USE um_parcore, ONLY: mype, nproc_max

!++SAN adding print statements
USE umPrintMgr, ONLY: Printstatus, ummessage, umprint, PrStatus_Diag,          &
prstatus_normal, str, umprinterror

IMPLICIT NONE
!
! Description:
!  Reads overall, generic and model-specific control variables.
!
! Method:
!  Reads namelist containing all overall, generic and model-specific
!  control variables not in History file.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!
! Declarations:
!

! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
INTEGER :: icode ! return code
CHARACTER (LEN=*), PARAMETER  :: RoutineName='READCNTL'
CHARACTER(LEN=errormessagelength) :: cmessage
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays:

INTEGER   :: fixhdr(len_fixhd)
INTEGER   :: len_io
REAL      :: a
INTEGER, ALLOCATABLE  :: lookup(:,:)
INTEGER   :: atmoscntl_unit ! unit no. for ATMOSCNTL file
INTEGER   :: shared_unit    ! unit no. for SHARED file
INTEGER   :: pp_unit
INTEGER   :: lbc_unit

!++SAN loop variable
INTEGER :: i

! Reading in of the various NLSTCALL namelists
TYPE(um_file_type), POINTER :: um_file

!- End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!++SAN checking icode going in
WRITE(umMessage,*) 'readcntl: SAN icode = ',icode
CALL umPrint(umMessage, src='readcntl')
!--SAN


!  To ensure the namelists work correctly, owing to dependence on variables
!  the namelists must be read in the following order:
!     NLSTCGEN
!     NLSTCALL

      ! Retrieve the unit number of the main namelist file
atmoscntl_unit = get_file_unit_by_id("atmoscntl", handler="fortran")
shared_unit   = get_file_unit_by_id("shared", handler="fortran")

!++SAN Check the contents of these units
WRITE(umMessage,*) 'readcntl: SAN atmoscntl_unit',atmoscntl_unit
CALL umPrint(umMessage, src='readcntl')
WRITE(umMessage,*) 'readcntl: SAN shared_unit',shared_unit
CALL umPrint(umMessage, src='readcntl')
WRITE(umMessage,*) 'readcntl: SAN calling read_nml_nlstcgen'
CALL umPrint(umMessage, src='readcntl')
!--SAN

CALL read_nml_nlstcgen(atmoscntl_unit)

!++SAN comment out pp file options

!++SAN Added read for box model I/O options
WRITE(umMessage,*) 'readcntl: SAN calling read_nml_nlstcall_box_options'
CALL umPrint(umMessage, src='readcntl')
CALL read_nml_nlstcall_box_options(atmoscntl_unit)
IF ( mype == 0 ) REWIND(atmoscntl_unit)

!++SAN
WRITE(umMessage,*) 'readcntl: SAN calling read_nml_nlstcall'
CALL umPrint(umMessage, src='readcntl')

CALL read_nml_nlstcall(shared_unit)
CALL print_nlist_nlstcall()

!++SAN double check the model_basis_time at this level:
WRITE(umMessage,*) 'readcntl: SAN double checking model_basis_time'
CALL umPrint(umMessage, src='readcntl')
DO i=1,6
  WRITE(umMessage,'(A,I0,A,I0)')'  model_basis_time(',i,') = ',               &
       model_basis_time(i)
  CALL umPrint(umMessage,src='nlstcall_mod')
END DO
!--SAN

!++SAN
WRITE(umMessage,*) 'readcntl: SAN calling check_nlstcgen'
CALL umPrint(umMessage, src='readcntl')
!--SAN

! Check nlstcgen has to be called after NLSTCALL is read, as it is dependent
! on a variable in this namelist
CALL check_nlstcgen()

!  Read configuration id, is itab set explictly as an input?
!++SAN DELETE

!++SAN DELETED Read MPP configuration data

!++SAN
WRITE(umMessage,*) 'readcntl: SAN exiting readcntl'
CALL umPrint(umMessage, src='readcntl')
!--SAN

!++SAN checking icode going out
WRITE(umMessage,*) 'readcntl: SAN icode = ',icode
CALL umPrint(umMessage, src='readcntl')
!--SAN


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE readcntl
END MODULE readcntl_mod
