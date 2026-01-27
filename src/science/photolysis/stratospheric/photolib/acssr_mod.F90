! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acssr_mod

   USE parkind1, ONLY: jprb, jpim
   USE yomhook, ONLY: lhook, dr_hook
   IMPLICIT NONE

! Description:
!     O2 cross section in Schumann-Runge region.

! Method:
!     Origin of data: Table
!     7.7 in WMO (1985) gives transmission in the Schumann-Runge window
!     (175-205 nm). d ln T/ d (column oxygen) is displayed here, which
!     is the O2 cross section.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ACSSR_MOD'

CONTAINS
   SUBROUTINE acssr(tc, jpwav, ao2sr)
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface
! Total slant path O2 column above the point in molecules per cm^2.
      REAL, INTENT(IN)    :: tc
      INTEGER, INTENT(IN) :: jpwav !     Number of wavelength intervals
! Absorption cross-section of O2 in cm^2 for each wavelength interval
! in the Schumann-Runge bands.
      REAL, INTENT(IN OUT) :: ao2sr(jpwav)

! Local parameters and variables
!     Number of Schumann-Runge wavelength intervals.
      INTEGER, PARAMETER :: jpwavesr = 17
! Number of O2 column intervals
      INTEGER, PARAMETER :: no2col = 20

!     Oxygen column intervol
      INTEGER :: jo2
      REAL :: frac
      LOGICAL, SAVE :: first = .TRUE.

! O2 cross section in Schumann-Runge bands
      REAL, SAVE :: sr(jpwavesr, no2col)
      REAL, SAVE :: logsr(jpwavesr, no2col)
      REAL, SAVE :: o2col(no2col)
      REAL, SAVE :: logo2col(no2col)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACSSR'

! Main block: Initialize data upon first entry
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
      IF (first) THEN
         sr(:, 1) = [ &
                    2.0596E-19, 1.3273E-19, 6.6113E-20, 5.1048E-20, 4.5025E-20, 1.4983E-20, &
                    8.9858E-21, 2.9953E-21, 2.9932E-21, 1.0000E-40, 1.0000E-40, 1.0000E-40, &
                    1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 2) = [ &
                    1.9882E-19, 1.1754E-19, 5.8335E-20, 5.0977E-20, 4.0015E-20, 1.4504E-20, &
                    6.0382E-21, 3.6212E-21, 1.2071E-21, 1.0000E-40, 1.0000E-40, 1.0000E-40, &
                    1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 3) = [ &
                    1.8704E-19, 9.5822E-20, 5.1962E-20, 4.8400E-20, 3.6629E-20, 1.2836E-20, &
                    6.4058E-21, 2.9866E-21, 1.2793E-21, 8.5257E-22, 4.2614E-22, 1.0000E-40, &
                    1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 4) = [ &
                    1.6133E-19, 6.4581E-20, 4.2558E-20, 4.2452E-20, 2.9852E-20, 1.1903E-20, &
                    6.1489E-21, 2.8427E-21, 1.1954E-21, 4.4801E-22, 5.9735E-22, 1.4923E-22, &
                    1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 5) = [ &
                    1.2039E-19, 3.7382E-20, 2.9637E-20, 3.1606E-20, 2.0773E-20, 1.0706E-20, &
                    5.7048E-21, 2.7256E-21, 1.1954E-21, 5.4263E-22, 5.4263E-22, 1.0845E-22, &
                    1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 6) = [ &
                    8.1089E-20, 2.2320E-20, 1.7193E-20, 2.0043E-20, 1.2370E-20, 8.3474E-21, &
                    4.9413E-21, 2.4020E-21, 1.1291E-21, 4.5832E-22, 5.2091E-22, 1.4555E-22, &
                    4.1557E-23, 1.0000E-40, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 7) = [ &
                    5.2569E-20, 1.4185E-20, 9.5606E-21, 1.1560E-20, 7.1082E-21, 5.4135E-21, &
                    3.7658E-21, 1.8750E-21, 1.0425E-21, 4.5744E-22, 5.1708E-22, 1.2648E-22, &
                    3.3682E-23, 8.4137E-24, 1.0000E-40, 1.0000E-40, 1.0000E-40]
         sr(:, 8) = [ &
                    3.3180E-20, 9.4479E-21, 5.6342E-21, 6.7482E-21, 4.2844E-21, 3.1560E-21, &
                    2.4997E-21, 1.3218E-21, 8.5520E-22, 4.1620E-22, 4.9829E-22, 1.3120E-22, &
                    3.6328E-23, 1.0000E-40, 3.6300E-24, 3.6300E-24, 3.6274E-24]
         sr(:, 9) = [ &
                    2.1783E-20, 6.5782E-21, 3.6155E-21, 4.2548E-21, 2.7253E-21, 1.9454E-21, &
                    1.5785E-21, 9.1690E-22, 6.2044E-22, 3.6428E-22, 4.5727E-22, 1.2660E-22, &
                    3.4726E-23, 3.2997E-24, 1.6499E-24, 1.0000E-40, 1.0000E-40]
         sr(:, 10) = [ &
                     1.6102E-20, 4.7863E-21, 2.5785E-21, 2.8755E-21, 1.8324E-21, 1.2969E-21, &
                     1.0247E-21, 6.3919E-22, 4.0578E-22, 2.7130E-22, 3.8879E-22, 1.1894E-22, &
                     3.2220E-23, 4.6906E-24, 1.5631E-24, 1.0000E-40, 1.0000E-40]
         sr(:, 11) = [ &
                     1.4089E-20, 3.6881E-21, 2.0761E-21, 3.7937E-21, 1.3054E-21, 9.1906E-22, &
                     7.2141E-22, 4.5887E-22, 2.7029E-22, 2.1300E-22, 3.0038E-22, 1.0554E-22, &
                     2.9401E-23, 4.2685E-24, 3.8789E-24, 3.8765E-25, 7.7525E-25]
         sr(:, 12) = [ &
                     1.4608E-20, 3.1322E-21, 1.8465E-21, 8.0188E-22, 8.9807E-22, 6.8966E-22, &
                     5.5373E-22, 3.4768E-22, 1.9326E-22, 1.6128E-22, 2.1394E-22, 8.7938E-23, &
                     2.4211E-23, 5.0370E-24, 1.4084E-24, 6.0272E-25, 2.0093E-25]
         sr(:, 13) = [ &
                     1.4608E-20, 2.9988E-21, 1.7814E-21, 1.3742E-21, 7.9193E-22, 5.3877E-22, &
                     4.6246E-22, 2.8948E-22, 1.4895E-22, 1.2973E-22, 1.5169E-22, 6.8160E-23, &
                     2.0371E-23, 5.7585E-24, 2.4894E-24, 6.4781E-25, 5.3968E-25]
         sr(:, 14) = [ &
                     1.4608E-20, 3.1296E-21, 1.7853E-21, 1.1892E-21, 6.3662E-22, 4.2117E-22, &
                     3.8823E-22, 1.8296E-22, 1.1938E-22, 1.0509E-22, 1.0961E-22, 5.1871E-23, &
                     1.7125E-23, 6.7396E-24, 2.6245E-24, 7.7167E-25, 4.1525E-25]
         sr(:, 15) = [ &
                     1.4608E-20, 2.9777E-21, 1.5883E-21, 9.6037E-22, 5.0665E-22, 3.0382E-22, &
                     2.7925E-22, 1.3249E-22, 8.6811E-23, 7.6305E-23, 7.3698E-23, 3.7378E-23, &
                     1.3751E-23, 5.5964E-24, 2.6709E-24, 7.9697E-25, 4.1381E-25]
         sr(:, 16) = [ &
                     1.4608E-20, 2.9777E-21, 1.2210E-21, 7.0899E-22, 3.6453E-22, 2.0328E-22, &
                     1.7512E-22, 8.4507E-23, 5.6116E-23, 4.8982E-23, 4.3889E-23, 2.4642E-23, &
                     9.6018E-24, 4.3275E-24, 2.4343E-24, 6.6242E-25, 4.1886E-25]
         sr(:, 17) = [ &
                     1.4608E-20, 2.9777E-21, 1.2210E-21, 5.5897E-22, 2.7411E-22, 1.4980E-22, &
                     1.1815E-22, 5.8790E-23, 3.6523E-23, 3.2044E-23, 2.5235E-23, 1.5842E-23, &
                     6.5220E-24, 3.2660E-24, 2.2322E-24, 6.0985E-25, 4.1236E-25]
         sr(:, 18) = [ &
                     1.4608E-20, 2.9777E-21, 1.2210E-21, 5.5897E-22, 2.3540E-22, 1.2884E-22, &
                     1.7942E-22, 4.6834E-23, 2.6335E-23, 2.2983E-23, 1.5504E-23, 1.0513E-23, &
                     4.5950E-24, 2.6892E-24, 2.1383E-24, 5.7770E-25, 3.9394E-25]
         sr(:, 19) = [ &
                     1.4608E-20, 2.9777E-21, 1.2210E-21, 5.5897E-22, 2.3540E-22, 1.1905E-22, &
                     5.0663E-23, 4.0295E-23, 2.1386E-23, 1.7647E-23, 1.0693E-23, 7.3588E-24, &
                     3.3738E-24, 2.3488E-24, 2.0769E-24, 5.6322E-25, 3.9207E-25]
         sr(:, 20) = [ &
                     1.4608E-20, 2.9777E-21, 1.2210E-21, 5.5897E-22, 2.3540E-22, 1.1905E-22, &
                     5.0663E-23, 3.5399E-23, 1.8521E-23, 1.4160E-23, 8.0985E-24, 5.4163E-24, &
                     2.5977E-24, 2.1190E-24, 2.0413E-24, 5.4415E-25, 3.8368E-25]

         logsr = LOG(sr)

         o2col = [ &
                 5.3368E+16, 1.0627E+17, 2.4629E+17, 6.4304E+17, 1.7548E+18, 4.7351E+18, &
                 1.2299E+19, 3.0403E+19, 7.1301E+19, 1.5943E+20, 3.4126E+20, 6.9993E+20, &
                 1.3797E+21, 2.6309E+21, 4.9365E+21, 9.3681E+21, 1.8360E+22, 3.7372E+22, &
                 7.8501E+22, 1.6851E+23]

         logo2col = LOG(o2col)
         first = .FALSE.
      END IF

! Perform linear interpolation in log-log space of cross section in
! Schumann-Runge window. In case the O2 column is outside the window
! covered, take

      IF (tc < o2col(1)) THEN
         jo2 = 1
         frac = 0.0
      ELSE IF (tc >= o2col(no2col)) THEN
         jo2 = no2col - 1
         frac = 1.0
      ELSE
         jo2 = no2col - 1
         DO WHILE (o2col(jo2) > tc)
            jo2 = jo2 - 1
         END DO
         frac = (LOG(tc) - logo2col(jo2))/ &
                (logo2col(jo2 + 1) - logo2col(jo2))
      END IF

      ao2sr(46:45 + jpwavesr) = &
         EXP((1.0 - frac)*logsr(:, jo2) + frac*logsr(:, jo2 + 1))
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN

   END SUBROUTINE acssr
END MODULE acssr_mod
