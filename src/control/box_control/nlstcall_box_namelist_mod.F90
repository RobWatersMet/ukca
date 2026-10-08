! *****************************COPYRIGHT*******************************
!
! Copyright 2017-2019 University of Reading
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors
! may be used to endorse or promote products derived from this software without
! specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
! Description: Module based interface to the nlstcall_box_option namelist,
!              which configures box model input/output streams and files
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: NetCDF output
!
! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_box_namelist_mod

USE filenamelength_mod, ONLY: filenamelength
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook,dr_hook
USE parkind1, ONLY: jprb,jpim

USE profilename_length_mod, ONLY: fileid_length


IMPLICIT NONE

! Read variables for nlstcall_box_options
CHARACTER(LEN=filenamelength) :: ukca_box_nml
CHARACTER(LEN=filenamelength) :: tracer_in_filename
CHARACTER(LEN=filenamelength) :: tracer_out_filename
CHARACTER(LEN=filenamelength) :: ntp_out_filename
CHARACTER(LEN=filenamelength) :: flux_out_filename
CHARACTER(LEN=filenamelength) :: rate_out_filename
CHARACTER(LEN=filenamelength) :: photol_jrate_in_filename
REAL                          :: tracer_nullval

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

NAMELIST / nlstcall_box_options / ukca_box_nml,                                &
                                  tracer_in_filename,                          &
                                  tracer_out_filename,                         &
                                  ntp_out_filename,                            &
                                  flux_out_filename,                           &
                                  rate_out_filename,                           &
                                  photol_jrate_in_filename,                    &
                                  tracer_nullval

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLSTCALL_BOX_NAMELIST_MOD'

CONTAINS

! This routine will read as many nlstcall_box_options namelists as it finds
! in the file attached to the provided unit
!-------------------------------------------------------------------------------
SUBROUTINE read_nml_nlstcall_box_options(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type
USE umprintmgr, ONLY: printstatus, prstatus_oper

IMPLICIT NONE

INTEGER, INTENT(IN)  :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL_BOX_OPTIONS'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_chars = 6*filenamelength
INTEGER, PARAMETER :: n_real  = 1

TYPE :: my_namelist
  SEQUENCE
  CHARACTER (LEN=filenamelength) :: ukca_box_nml
  CHARACTER (LEN=filenamelength) :: tracer_in_filename
  CHARACTER (LEN=filenamelength) :: tracer_out_filename
  CHARACTER (LEN=filenamelength) :: ntp_out_filename
  CHARACTER (LEN=filenamelength) :: flux_out_filename
  CHARACTER (LEN=filenamelength) :: rate_out_filename
  CHARACTER (LEN=filenamelength) :: photol_jrate_in_filename
  REAL                           :: tracer_nullval
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_chars_in=n_chars,             &
  n_real_in=n_real)

! Setup defaults
ukca_box_nml               = 'TEST.box'
tracer_in_filename         = 'chem_tracer_pars'
tracer_out_filename        = 'tracer_out.csv'
ntp_out_filename           = 'ntp_out.csv'
flux_out_filename          = 'flux_out.csv'
rate_out_filename          = 'rate_out.csv'
photol_jrate_in_filename   = ''
tracer_nullval             = 1e-15

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=nlstcall_box_options,                                  &
       IOSTAT=ErrorStatus, IOMSG=iomessage)

  CALL check_iostat(errorstatus, "namelist NLSTCALL_BOX_OPTIONS", iomessage)

  my_nml % ukca_box_nml               = ukca_box_nml
  my_nml % tracer_in_filename         = tracer_in_filename
  my_nml % tracer_out_filename        = tracer_out_filename
  my_nml % ntp_out_filename           = ntp_out_filename
  my_nml % flux_out_filename          = flux_out_filename
  my_nml % rate_out_filename          = rate_out_filename
  my_nml % photol_jrate_in_filename   = photol_jrate_in_filename
  my_nml % tracer_nullval             = tracer_nullval

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

    ukca_box_nml               = my_nml % ukca_box_nml
    tracer_in_filename         = my_nml % tracer_in_filename
    tracer_out_filename        = my_nml % tracer_out_filename
    ntp_out_filename           = my_nml % ntp_out_filename
    flux_out_filename          = my_nml % flux_out_filename
    rate_out_filename          = my_nml % rate_out_filename
    photol_jrate_in_filename   = my_nml % photol_jrate_in_filename
    tracer_nullval             = my_nml % tracer_nullval

END IF

IF (printstatus >= prstatus_oper .AND. mype == 0) THEN
  CALL print_nlist_nlstcall_box_options()
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall_box_options

!-------------------------------------------------------------------------------
SUBROUTINE print_nlist_nlstcall_box_options()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLSTCALL_BOX_OPTIONS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlstcall_box_options', src=modulename)

WRITE(linebuffer,"(A,A)") ' ukca_box_nml             = ', TRIM(ukca_box_nml)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' tracer_in_filename       = ', TRIM(tracer_in_filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' tracer_out_filename      = ', TRIM(tracer_out_filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' ntp_out_filename         = ', TRIM(ntp_out_filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' flux_out_filename        = ', TRIM(flux_out_filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' rate_out_filename        = ', TRIM(rate_out_filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' photol_jrate_in_filename = ', TRIM(photol_jrate_in_filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,E12.3)") ' tracer_nullval = ', tracer_nullval
CALL umprint(linebuffer,src=modulename)


CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlstcall_box_options

END MODULE nlstcall_box_namelist_mod
