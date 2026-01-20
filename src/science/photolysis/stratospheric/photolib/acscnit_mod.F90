! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acscnit_mod

   USE parkind1, ONLY: jprb, jpim
   USE yomhook, ONLY: lhook, dr_hook
   IMPLICIT NONE

! Description:
!     Calculate ClONO2 temperature dependent cross sections,

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ACSCNIT_MOD'

CONTAINS
   SUBROUTINE acscnit(temp, wavenm, acnita, acnitb)
      USE ukca_parpho_mod, ONLY: jpwav
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface
      REAL, INTENT(IN)    :: temp
      REAL, INTENT(IN)    :: wavenm(jpwav)  ! Contains the wavelength intervals
! ACNITA, ACNITB Contain the cross sections times quantum yields
! at model wavelengths.Contains the result on exit.
      REAL, INTENT(IN OUT) :: acnita(jpwav)
      REAL, INTENT(OUT) :: acnitb(jpwav)

! Local variables

      INTEGER :: i

      REAL :: tc
      LOGICAL, SAVE :: first = .TRUE.
!     ACNIT at 296K
      REAL, PARAMETER :: acnitt(jpwav) = [(0.0, i=1, 54), &
                                          4.320E-19, 1.980E-18, 2.911E-18, 3.015E-18, 2.871E-18, 2.785E-18, &
                                          2.780E-18, 2.839E-18, 2.956E-18, 3.097E-18, 3.264E-18, 3.386E-18, &
                                          3.448E-18, 3.392E-18, 3.236E-18, 2.971E-18, 2.640E-18, 2.268E-18, &
                                          1.922E-18, 1.591E-18, 1.314E-18, 1.086E-18, 8.967E-19, 7.444E-19, &
                                          6.074E-19, 5.129E-19, 4.352E-19, 3.703E-19, 3.152E-19, 2.662E-19, &
                                          2.213E-19, 1.840E-19, 1.498E-19, 1.211E-19, 9.519E-20, 7.333E-20, &
                                          5.500E-20, 4.007E-20, 2.969E-20, 2.190E-20, 1.600E-20, 1.142E-20, &
                                          8.310E-21, 6.114E-21, 4.660E-21, 3.657E-21, 3.020E-21, 2.576E-21, &
                                          2.290E-21, 2.079E-21, 2.000E-21, 1.795E-21, 1.590E-21, 1.414E-21, &
                                          1.210E-21, 1.056E-21, 9.090E-22, 7.588E-22, 6.380E-22, 5.376E-22, &
                                          4.440E-22, 3.672E-22, 3.160E-22, 2.314E-22, 1.890E-22, 5.264E-23, &
                                          (0.0, i=1, 83)]
!     Coeffs A1, A2 from Burkholder et al GRL 1994
      REAL, PARAMETER :: a1(jpwav) = [(0.0, i=1, 54), &
                                      1.73E-05, 7.50E-05, 1.00E-04, 8.82E-05, 3.61E-05, -5.88E-05, &
                                      -1.95E-04, -3.44E-04, -5.11E-04, -6.59E-04, -7.85E-04, -8.71E-04, &
                                      -9.03E-04, -8.73E-04, -7.83E-04, -6.38E-04, -4.53E-04, -2.35E-04, &
                                      -6.97E-06, 2.19E-04, 4.16E-04, 5.64E-04, 6.78E-04, 7.81E-04, &
                                      9.08E-04, 1.08E-03, 1.26E-03, 1.44E-03, 1.59E-03, 1.74E-03, &
                                      1.88E-03, 2.03E-03, 2.21E-03, 2.37E-03, 2.55E-03, 2.74E-03, &
                                      2.95E-03, 3.28E-03, 3.72E-03, 4.12E-03, 4.53E-03, 4.98E-03, &
                                      5.40E-03, 5.75E-03, 5.92E-03, 5.85E-03, 5.51E-03, 4.92E-03, &
                                      4.02E-03, 3.27E-03, 2.70E-03, 2.08E-03, 1.33E-03, 7.65E-04, &
                                      3.53E-04, 2.39E-04, 4.10E-04, 7.77E-04, 1.38E-03, 2.15E-03, &
                                      3.38E-03, 4.88E-03, 6.70E-03, 8.10E-03, 9.72E-03, 9.96E-03, &
                                      (0.0, i=1, 83)]
      REAL, PARAMETER :: a2(jpwav) = [(0.0, i=1, 54), &
                                      -1.16E-06, -5.32E-06, -7.85E-06, -8.21E-06, -7.81E-06, -7.52E-06, &
                                      -7.46E-06, -7.62E-06, -7.93E-06, -8.30E-06, -8.69E-06, -9.03E-06, &
                                      -9.25E-06, -9.37E-06, -9.37E-06, -9.27E-06, -9.06E-06, -8.65E-06, &
                                      -7.98E-06, -7.13E-06, -6.36E-06, -6.00E-06, -6.09E-06, -6.45E-06, &
                                      -6.57E-06, -6.03E-06, -5.08E-06, -4.20E-06, -3.50E-06, -2.74E-06, &
                                      -1.87E-06, -9.85E-07, 6.15E-08, 1.13E-06, 2.14E-06, 3.05E-06, &
                                      3.74E-06, 5.29E-06, 7.78E-06, 9.88E-06, 1.20E-05, 1.49E-05, &
                                      1.84E-05, 2.27E-05, 2.70E-05, 3.01E-05, 3.11E-05, 2.86E-05, &
                                      2.07E-05, 1.38E-05, 8.59E-06, 2.01E-06, -7.40E-06, -1.44E-05, &
                                      -1.91E-05, -2.11E-05, -2.05E-05, -1.87E-05, -1.42E-05, -7.14E-06, &
                                      4.47E-06, 1.93E-05, 3.87E-05, 5.57E-05, 7.52E-05, 7.81E-05, &
                                      (0.0, i=1, 83)]
      REAL, SAVE :: quanty(jpwav) ! Quantum yield for Cl + NO3 channel

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACSCNIT'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
      IF (first) THEN
         quanty = MIN(MAX(7.143E-3*wavenm - 1.6, 0.6), 1.0)
         first = .FALSE.
      END IF

      tc = temp
      tc = MAX(tc, 220.0)
      tc = MIN(tc, 298.0)

! Initialise acnitb here as acnita already initialised in fill_spectra
      acnitb(:) = 0.0

! ClONO2 cross section, following JPL (2002)
      acnita(55:120) = acnitt(55:120)*(1.0 + a1(55:120)*(tc - 296.0) &
                                       + a2(55:120)*((tc - 296.0)**2))
! ClONO2 cross section, times quantum yield for ClO + NO2 channel
      acnitb(55:120) = (1.0 - quanty(55:120))*acnita(55:120)
! ClONO2 cross section, times quantum yield for Cl  + NO3 channel
      acnita(55:120) = quanty(55:120)*acnita(55:120)
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN

   END SUBROUTINE acscnit
END MODULE acscnit_mod
