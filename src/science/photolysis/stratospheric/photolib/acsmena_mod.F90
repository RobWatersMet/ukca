! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsmena_mod

   USE parkind1, ONLY: jprb, jpim
   USE yomhook, ONLY: lhook, dr_hook
   IMPLICIT NONE

! Description:
!     Calculate T-dependent MeONO2 cross sections
!     The expression is valid for the wavelength range; 186.1-296.3 nm
!                            and the temperature range; 225-295 K.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ACSMENA_MOD'

CONTAINS
   SUBROUTINE acsmena(t, amena)
      USE ukca_parpho_mod, ONLY: jpwav
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

      REAL, INTENT(IN)    :: t              ! Temperature in kelvin
! Absorption cross-section of COS in cm^2 for temperature T for each
! interval wavelength.
      REAL, INTENT(IN OUT) :: amena(jpwav)
! Local variables

      INTEGER :: i

      REAL :: tc
!     JPL 2002 T=295K values.
      REAL, PARAMETER :: amena298(jpwav) = [ &
                         (0.0, i=1, 75), &
                         6.049E-20, 5.0688E-20, 4.142E-20, 3.77E-20, &
                         3.4972E-20, 3.3116E-20, 3.1512E-20, 2.9788E-20, &
                         2.7758E-20, 2.504E-20, 2.2262E-20, 1.9244E-20, &
                         1.6052E-20, 1.2914E-20, 9.99601E-21, 7.372E-21, &
                         5.13921E-21, 3.3664E-21, 2.092E-21, 1.34E-21, &
                         6.33E-22, 3.16E-22, 1.44E-22, 6.61E-23, &
                         2.74E-23, 1.22E-23, (0.0, i=1, 102)]
!     JPL 2002 T=225K values.
      REAL, PARAMETER :: b(jpwav) = [ &
                         (0.0, i=1, 75), &
                         0.003499, 0.0033888, 0.0032636, 0.003059, &
                         0.0029152, 0.0028256, 0.0028262, 0.0028552, &
                         0.0029182, 0.003032, 0.003164, 0.0033214, &
                         0.0034962, 0.0037098, 0.0039256, 0.004212, &
                         0.0045922, 0.0050392, 0.0056062, 0.00633, &
                         0.00734, 0.00874, 0.00997, 0.0136, &
                         0.0136, 0.0136, (0.0, i=1, 102)]

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACSMENA'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
!     Check that temperature is in range.
      tc = MAX(MIN(t, 330.0), 240.0)

! ln sigma = ln sigma(298K) + b * (T - 298K)
      amena(76:101) = EXP(LOG(amena298(76:101)) + &
                          b(76:101)*(tc - 298.0))

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
   END SUBROUTINE acsmena
END MODULE acsmena_mod
