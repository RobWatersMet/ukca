! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsf11_mod

   USE parkind1, ONLY: jprb, jpim
   USE yomhook, ONLY: lhook, dr_hook
   IMPLICIT NONE

! Description:
!     Calculate T-dependent F11 cross sections
!     The expression is valid for the wavelength range; 174-230 nm
!                            and the temperature range; 210-300 K.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ACSF11_MOD'

CONTAINS
   SUBROUTINE acsf11(t, wavenm, af11)
      USE ukca_parpho_mod, ONLY: jpwav
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface
      REAL, INTENT(IN) :: t                ! Temperature in kelvin
      REAL, INTENT(IN) :: wavenm(jpwav)    ! Wavelength of each interval in nm
! Absorption cross-section of F11 in cm^2 for temperature T for each
! interval wavelength.
      REAL, INTENT(IN OUT) :: af11(jpwav)

! Local variables

      INTEGER :: i
      INTEGER :: jw

      REAL :: tc
      REAL :: arg
! Absorption cross-section of F11 at 298K
! CFCl3: JPL 1992 T=298K values.
      REAL, PARAMETER :: af11t(jpwav) = [ &
                         (0.0, i=1, 36), &
                         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.139E-19, 2.478E-18, &
                         3.182E-18, 3.168E-18, 3.141E-18, 3.098E-18, 3.042E-18, 3.096E-18, &
                         2.968E-18, 2.765E-18, 2.558E-18, 2.319E-18, 2.107E-18, 1.839E-18, &
                         1.574E-18, 1.332E-18, 1.092E-18, 8.911E-19, 7.221E-19, 5.751E-19, &
                         4.389E-19, 3.340E-19, 2.377E-19, 1.700E-19, 1.171E-19, 7.662E-20, &
                         5.082E-20, 3.184E-20, 1.970E-20, 1.206E-20, 8.000E-21, 4.834E-21, &
                         2.831E-21, 1.629E-21, 9.327E-22, 5.209E-22, 3.013E-22, 1.617E-22, &
                         9.035E-23, 5.427E-23, 3.474E-23, 2.141E-23, 9.102E-24, 1.499E-25, &
                         (0.0, i=1, 119)]

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACSF11'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!     Check that temperature is in range.
      tc = MIN(t, 300.0)
      tc = MAX(tc, 210.0)

!     Wavelength loop.
      DO jw = 45, 72
         arg = 1.0E-4*(wavenm(jw) - 184.9)*(tc - 298.0)
         af11(jw) = af11t(jw)*EXP(arg)
      END DO
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN

   END SUBROUTINE acsf11
END MODULE acsf11_mod
