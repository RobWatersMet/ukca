! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acshno3_mod

   USE parkind1, ONLY: jprb, jpim
   USE yomhook, ONLY: lhook, dr_hook
   IMPLICIT NONE

! Description:
!     Calculate HNO3 temperature dependent cross sections

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ACSHNO3_MOD'

CONTAINS
   SUBROUTINE acshno3(temp, ahno3)
      USE ukca_parpho_mod, ONLY: jpwav
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface
      REAL, INTENT(IN)    :: temp
! Contains the cross sections at model wavelengths.Contains the result on exit.
      REAL, INTENT(IN OUT) :: ahno3(jpwav)

! Local variables

      INTEGER :: i

      REAL :: tc            ! Contains the wavelength intervals.
! Contains the cross sections at 298K
      REAL, PARAMETER :: ahno3t(jpwav) = [(0.0, i=1, 48), &
                                          0.000E+00, 8.305E-18, 1.336E-17, 1.575E-17, 1.491E-17, 1.385E-17, &
                                          1.265E-17, 1.150E-17, 1.012E-17, 8.565E-18, 6.739E-18, 5.147E-18, &
                                          3.788E-18, 2.719E-18, 1.796E-18, 1.180E-18, 7.377E-19, 4.487E-19, &
                                          2.810E-19, 1.826E-19, 1.324E-19, 1.010E-19, 8.020E-20, 6.479E-20, &
                                          5.204E-20, 4.178E-20, 3.200E-20, 2.657E-20, 2.298E-20, 2.086E-20, &
                                          1.991E-20, 1.962E-20, 1.952E-20, 1.929E-20, 1.882E-20, 1.804E-20, &
                                          1.681E-20, 1.526E-20, 1.335E-20, 1.136E-20, 9.242E-21, 7.186E-21, &
                                          5.320E-21, 3.705E-21, 2.393E-21, 1.442E-21, 8.140E-22, 4.131E-22, &
                                          1.970E-22, 9.434E-23, 4.310E-23, 2.204E-23, 1.030E-23, 5.841E-24, &
                                          4.170E-24, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
                                          (0.0, i=1, 95)]
! Contains B coefficient. Intercept B from Burkholder
      REAL, PARAMETER :: b(jpwav) = [(0.0, i=1, 54), &
                                     0.000E+00, 0.000E+00, 1.152E-03, 1.668E-03, 1.653E-03, 1.673E-03, &
                                     1.720E-03, 1.750E-03, 1.817E-03, 1.935E-03, 2.060E-03, 2.168E-03, &
                                     2.178E-03, 2.195E-03, 2.106E-03, 1.987E-03, 1.840E-03, 1.782E-03, &
                                     1.838E-03, 1.897E-03, 1.970E-03, 1.978E-03, 1.855E-03, 1.655E-03, &
                                     1.416E-03, 1.247E-03, 1.162E-03, 1.121E-03, 1.136E-03, 1.199E-03, &
                                     1.315E-03, 1.493E-03, 1.637E-03, 1.767E-03, 1.928E-03, 2.139E-03, &
                                     2.380E-03, 2.736E-03, 3.139E-03, 3.695E-03, 4.230E-03, 5.151E-03, &
                                     6.450E-03, 7.327E-03, 9.750E-03, 1.013E-02, 1.180E-02, 1.108E-02, &
                                     9.300E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
                                     (0.0, i=1, 95)]

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACSHNO3'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
      tc = temp
      tc = MAX(tc, 200.0)
      tc = MIN(tc, 360.0)

      ahno3(50:103) = ahno3t(50:103)*EXP(b(50:103)*(tc - 298.0))
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN

   END SUBROUTINE acshno3
END MODULE acshno3_mod
