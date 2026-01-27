! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsno2_mod

   USE parkind1, ONLY: jprb, jpim
   USE yomhook, ONLY: lhook, dr_hook
   IMPLICIT NONE

! Description:
!     Calculate NO2 temperature dependent cross sections,

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ACSNO2_MOD'

CONTAINS
   SUBROUTINE acsno2(temp, ano2)
      USE ukca_parpho_mod, ONLY: jpwav
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface

      INTEGER :: i

      REAL, INTENT(IN)    :: temp
! Contains the cross sections at model
! wavelengths.Contains the result on exit
      REAL, INTENT(IN OUT) :: ano2(jpwav)

! Local variables
! Contains the cross sections at 298K
      REAL, PARAMETER ::  ano2t(jpwav) = [ &
                         (0.0, i=1, 48), &
                         0.000E+00, 0.000E+00, 0.000E+00, 2.670E-19, 2.780E-19, 2.900E-19, &
                         2.790E-19, 2.600E-19, 2.420E-19, 2.450E-19, 2.480E-19, 2.750E-19, &
                         4.145E-19, 4.478E-19, 4.454E-19, 4.641E-19, 4.866E-19, 4.818E-19, &
                         5.022E-19, 4.441E-19, 4.713E-19, 3.772E-19, 3.929E-19, 2.740E-19, &
                         2.778E-19, 1.689E-19, 1.618E-19, 8.812E-20, 7.472E-20, 3.909E-20, &
                         2.753E-20, 2.007E-20, 1.973E-20, 2.111E-20, 2.357E-20, 2.698E-20, &
                         3.247E-20, 3.785E-20, 5.030E-20, 5.880E-20, 7.000E-20, 8.150E-20, &
                         9.720E-20, 1.154E-19, 1.344E-19, 1.589E-19, 1.867E-19, 2.153E-19, &
                         2.477E-19, 2.807E-19, 3.133E-19, 3.425E-19, 3.798E-19, 4.065E-19, &
                         4.313E-19, 4.717E-19, 4.833E-19, 5.166E-19, 5.315E-19, 5.508E-19, &
                         5.644E-19, 5.757E-19, 5.927E-19, 5.845E-19, 6.021E-19, 5.781E-19, &
                         5.999E-19, 5.651E-19, 5.812E-19, 0.000E+00, 0.000E+00, 0.000E+00, &
                         (0.0, i=1, 83)]
!     slope A
      REAL, PARAMETER ::  a(jpwav) = [ &
                         (0.0, i=1, 86), 0.075, 0.082, -0.053, -0.043, -0.031, -0.162, -0.284, &
                         -0.357, -0.536, -0.686, -0.786, -1.105, -1.355, -1.277, -1.612, &
                         -1.890, -1.219, -1.921, -1.095, -1.322, -1.102, -0.806, -0.867, &
                         -0.945, -0.923, -0.738, -0.599, -0.545, -1.129, 0.001, -1.208, &
                         (0.0, i=1, 86)]

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACSNO2'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

      ano2(87:117) = ano2t(87:117) + 1.0E-22*a(87:117)*(temp - 273.0)

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN

   END SUBROUTINE acsno2
END MODULE acsno2_mod
