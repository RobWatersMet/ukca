! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates water content of each mode given component
!    concentrations (in air) using ZSR and binary molalities
!    evaluated using water activity data from Jacobson,
!    "Fundamentals of Atmospheric Modelling", page 610 Table B
!    ("Water Activity Data" Table). Equations are from Chapter
!    18 ("Chemical Equilibrium and Dissolution Processes"
!    Chapter. Be aware that in some editions the chapter numbers
!    may be different.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_water_content_v_mod

   IMPLICIT NONE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_WATER_CONTENT_V_MOD'

CONTAINS

   SUBROUTINE ukca_water_content_v(nv, mask, cl, rh, ions, wc)

!-----------------------------------------------------------------------
!
! Purpose:
! -------
! Calculates water content of each mode given component
! concentrations (in air) using ZSR and binary molalities
! evaluated using water activity data from Jacobson,
! "Fundamentals of Atmospheric Modelling", page 610 Table B
! ("Water Activity Data" Table).
!
! Inputs
! ------
! NV     : Total number of gridboxes in domain
! IONS   : Logical indicating presence of each ion
! CL     : Concentration of each ion (moles/cc air)
! RH     : Relative humidity (fraction)
! MASK   : Logical array for where in domain to do calculation.
!
! Outputs
! -------
! WC     : Water content for aerosol (moles/cm3 of air)
!
! Local Variables
! ---------------
! IC     : Loop variable for cations
! IA     : Loop variable for anions
! AW     : Water activity (local copy of RH fraction)
! CLI    : Internal copy of CL
! CLP    : Ion pair concentrations (moles/cc air)
! MB     : Ion pair solution molalities (moles/cc water)
! N      : Ion stoiciometries for each electrolyte
! Z      : Charge for each ion
! Y      : Coefficients in expressions for binary molalities from
!        : Jacobson page 610 (Table B.10) for each electrolyte
!        : ("Water Activity Data" Table)
! RH_MIN : Lowest rh for which expression is valid
! MOLAL_MAX : Highest molality for which expression is valid.
!
!----------------------------------------------------------------------
      USE ukca_mode_setup, ONLY: ncation, nanion
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ukca_types_mod, ONLY: log_small
      USE ukca_config_specification_mod, ONLY: glomap_config

      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: nv                       ! No of points
      LOGICAL(KIND=log_small), INTENT(IN) :: mask(nv)                 ! Domain mask
      LOGICAL(KIND=log_small), INTENT(IN) :: ions(nv, -nanion:ncation) ! Ion presence
      ! switches
      REAL, INTENT(IN)    :: rh(nv)                   ! Relative humidity
      REAL, INTENT(IN)    :: cl(nv, -nanion:ncation)   ! Ion conc. (mol/cm3 of air)
      REAL, INTENT(OUT)   :: wc(nv)                   ! water content (mol/cm3 of air)

! Local variables
      INTEGER :: i                   ! Loop counter
      INTEGER :: j                   ! Loop counter
      INTEGER :: ic                  ! Cation loop variable
      INTEGER :: ia                  ! Anion  loop variable
      INTEGER :: m                   ! Counter
      INTEGER :: idx(nv)             ! Index

!WATER ACTIVITY (%RH EXPRESSED AS A FRACTION)
      REAL    :: aw(nv)
      REAL    :: dum(nv)
!INTERNAL COPY OF CL
      REAL    :: cli(nv, -nanion:ncation)
!ION PAIR CONCENTRATIONS
      REAL    :: clp(nv, ncation, -nanion:-1)
!ION PAIR BINARY SOLUTION MOLALITIES AT R
      REAL    :: mb(nv, ncation, -nanion:-1)
!ION STOICIOMETRIES FOR EACH ELECTROLYTE
      REAL    :: n(-nanion:ncation)
!ION CHARGES
      REAL    :: z(-nanion:ncation)
      REAL    :: y(3, -4:-1, 0:7)
      REAL    :: rh_min(3, -4:-1)
      REAL    :: molal_max(3, -4:-1)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_WATER_CONTENT_V'

      DATA(z(i), i=-nanion, ncation)/1.0, 1.0, 2.0, 1.0, 0.0, 1.0, 1.0, 1.0/
!                                    Cl NO3 SO4 HSO4 H2O  H NH4  Na

! As this subroutine is in an OMP parallel section,
! the use of the DATA statements cause y to be shared between threads,
! which creates issues when modifying the array outside of a DATA statement.
! Setting it as threadprivate avoids that.

!$OMP THREADPRIVATE(y)

! Set y coefficients for the calculations of solute molalities from
! Jacobson page 610 (Table B.10) for each electrolyte.
! ("Water Activity Data" Table)
! Also set the min rh and max molality for which the expressions are valid.
!     H+ HSO4- (1,-1)
      DATA(y(1, -1, j), j=0, 7)/ &
         3.0391387536E1, -1.8995058929E2, 9.7428231047E2, &
         -3.1680155761E3, 6.1400925314E3, -6.9116348199E3, &
         4.1631475226E3, -1.0383424491E3/
      DATA rh_min(1, -1), molal_max(1, -1)/0.0E0, 30.4E0/

!     2H+ SO42- (1,-2)
      DATA(y(1, -2, j), j=0, 7)/ &
         3.0391387536E1, -1.8995058929E2, 9.7428231047E2, &
         -3.1680155761E3, 6.1400925314E3, -6.9116348199E3, &
         4.1631475226E3, -1.0383424491E3/
      DATA rh_min(1, -2), molal_max(1, -2)/0.0E0, 30.4E0/

!     H+ NO3- (1,-3)
! One of the y co-efficients is incorrect here (j=6)
! Fixed with glomap_config%l_fix_ukca_water_content below.
      DATA(y(1, -3, j), j=0, 7)/ &
         2.306844303E1, -3.563608869E1, -6.210577919E1, &
         5.510176187E2, -1.460055286E3, 1.894467542E3, &
         -1.220611402E2, 3.098597737E2/
      DATA rh_min(1, -3), molal_max(1, -3)/0.0E0, 22.6E0/

!     H+ Cl- (1,-4)
      DATA(y(1, -4, j), j=0, 7)/ &
         1.874637647E1, -2.052465972E1, -9.485082073E1, &
         5.362930715E2, -1.223331346E3, 1.427089861E3, &
         -8.344219112E2, 1.90992437E2/
      DATA rh_min(1, -4), molal_max(1, -4)/0.0E0, 18.5E0/

!     Na+ HSO4- (2,-1)
      DATA(y(2, -1, j), j=0, 7)/ &
         1.8457001681E2, -1.6147765817E3, 8.444076586E3, &
         -2.6813441936E4, 5.0821277356E4, -5.5964847603E4, &
         3.2945298603E4, -8.002609678E3/
      DATA rh_min(2, -1), molal_max(2, -1)/1.9E0, 158.0E0/

!     2Na+ SO42- (2,-2)
      DATA(y(2, -2, j), j=0, 7)/ &
         5.5983158E2, -2.56942664E3, 4.47450201E3, &
         -3.45021842E3, 9.8527913E2, 0.0E0, &
         0.0E0, 0.0E0/
      DATA rh_min(2, -2), molal_max(2, -2)/58.0E0, 13.1E0/

!     Na+ NO3- (2,-3)
      DATA(y(2, -3, j), j=0, 7)/ &
         3.10221762E2, -1.82975944E3, 5.13445395E3, &
         -8.01200018E3, 7.07630664E3, -3.33365806E3, &
         6.5442029E2, 0.0E0/
      DATA rh_min(2, -3), molal_max(2, -3)/30.0E0, 56.8E0/

!     Na+ Cl- (2,-4)
      DATA(y(2, -4, j), j=0, 7)/ &
         5.875248E1, -1.8781997E2, 2.7211377E2, &
         -1.8458287E2, 4.153689E1, 0.0E0, &
         0.0E0, 0.0E0/
      DATA rh_min(2, -4), molal_max(2, -4)/47.0E0, 13.5E0/

!     NH4+ HSO4- (3,-1)
      DATA(y(3, -1, j), j=0, 7)/ &
         2.9997156464E2, -2.8936374637E3, 1.4959985537E4, &
         -4.5185935292E4, 8.110895603E4, -8.4994863218E4, &
         4.7928255412E4, -1.1223105556E4/
      DATA rh_min(3, -1), molal_max(3, -1)/6.5E0, 165.0E0/

!     2NH4+ SO42- (3,-2)
      DATA(y(3, -2, j), j=0, 7)/ &
         1.1065495E2, -3.6759197E2, 5.0462934E2, &
         -3.1543839E2, 6.770824E1, 0.0E0, &
         0.0E0, 0.0E0/
      DATA rh_min(3, -2), molal_max(3, -2)/37.0E0, 29.0E0/

!     NH4+ NO3- (3,-3)
      DATA(y(3, -3, j), j=0, 7)/ &
         3.983916445E3, 1.153123266E4, -2.13956707E5, &
         7.926990533E5, -1.407853405E6, 1.351250086E6, &
         -6.770046795E5, 1.393507324E5/
      DATA rh_min(3, -3), molal_max(3, -3)/62.0E0, 28.0E0/

!     NH4+ Cl- (3,-4)
      DATA(y(3, -4, j), j=0, 7)/ &
         -7.110541604E3, 7.217772665E4, -3.071054075E5, &
         7.144764216E5, -9.840230371E5, 8.03407288E5, &
         -3.603924022E5, 6.856992393E4/
      DATA rh_min(3, -4), molal_max(3, -4)/47.0E0, 23.2E0/

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!Put the correct value into the y array (N.B., can't do with DATA statment).
      IF (glomap_config%l_fix_ukca_water_content) y(1, -3, 6) = -1.220611402E3

      m = 0
      DO i = 1, nv
         IF (mask(i)) THEN
            m = m + 1
            idx(m) = i
         END IF
      END DO

! Write cl to internal variable to be adjusted in this subroutine only
      DO i = -nanion, ncation
         cli(:m, i) = cl(idx(:m), i)
      END DO
! Calculate mole concentrations of hypothetical ion pairs
      DO ic = 1, ncation
         DO ia = -nanion, -1
            ! ..Calculate stoichiometries for each ion pair
            n(ic) = z(ia)
            n(ia) = z(ic)
            IF (ABS(z(ic) - z(ia)) < EPSILON(0.0) .AND. &
                ABS(z(ia) - 1.0) > EPSILON(0.0)) THEN
               n(ic) = n(ic)/z(ic)
               n(ia) = n(ia)/z(ia)
            END IF
            WHERE (ions(idx(:m), ia) .AND. ions(idx(:m), ic))
               ! .. Calculate minimum ion pair concentration and subtract from
               ! .. Ion concentration. Eqn. 18.72 of Jacobson.
               clp(:m, ic, ia) = MIN(cli(:m, ic)/n(ic), cli(:m, ia)/n(ia))
               cli(:m, ic) = cli(:m, ic) - n(ic)*clp(:m, ic, ia)
               cli(:m, ia) = cli(:m, ia) - n(ia)*clp(:m, ic, ia)
            END WHERE
         END DO
      END DO
! Calculate binary electrolyte molalities at given aw Eqn. 18.66 of
! Jacobson.
      IF (glomap_config%l_fix_ukca_water_content) THEN
         DO ic = 1, ncation
            DO ia = -nanion, -1
               !If the water content bug is being fixed then copy the
               !original rh values into aw each loop since they may get erroneously
               !overwritten each time around the loop.
               !We want to apply a different min aw for each anion/cation pair:-
               aw(:m) = rh(idx(:m))
               !Prevent aw from going below rh_min (from Table B10 of Jacobson)
               !for the anion/cation pair
               WHERE (aw(:m) < rh_min(ic, ia)/1.0E2)
                  aw(:m) = rh_min(ic, ia)/1.0E2
               END WHERE
               mb(:m, ic, ia) = 0.0
               WHERE (ions(idx(:m), ia) .AND. ions(idx(:m), ic))
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 0)*aw(:m)**0
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 1)*aw(:m)**1
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 2)*aw(:m)**2
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 3)*aw(:m)**3
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 4)*aw(:m)**4
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 5)*aw(:m)**5
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 6)*aw(:m)**6
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 7)*aw(:m)**7
                  mb(:m, ic, ia) = MIN(mb(:m, ic, ia), molal_max(ic, ia))
               END WHERE
            END DO
         END DO
      ELSE
         ! Copy fractional relative humidity to water activity (local)
         aw(:m) = rh(idx(:m))
         DO ic = 1, ncation
            DO ia = -nanion, -1
               !Prevent aw from going below rh_min (from Table B10 of Jacobson)
               !for the anion/cation pair
               WHERE (aw(:m) < rh_min(ic, ia)/1.0E2)
                  aw(:m) = rh_min(ic, ia)/1.0E2
               END WHERE
               mb(:m, ic, ia) = 0.0
               WHERE (ions(idx(:m), ia) .AND. ions(idx(:m), ic))
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 0)*aw(:m)**0
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 1)*aw(:m)**1
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 2)*aw(:m)**2
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 3)*aw(:m)**3
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 4)*aw(:m)**4
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 5)*aw(:m)**5
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 6)*aw(:m)**6
                  mb(:m, ic, ia) = mb(:m, ic, ia) + y(ic, ia, 7)*aw(:m)**7
                  mb(:m, ic, ia) = MIN(mb(:m, ic, ia), molal_max(ic, ia))
               END WHERE
            END DO
         END DO
      END IF
! Calculate water content (mol/cm3 air); Eqn. 18.71 of Jacobson.
      dum(:m) = 0.0
      DO ic = 1, ncation
         DO ia = -nanion, -1
            WHERE (ions(idx(:m), ia) .AND. ions(idx(:m), ic))
               dum(:m) = dum(:m) + clp(:m, ic, ia)/mb(:m, ic, ia)
            END WHERE
         END DO
      END DO
      wc(idx(:m)) = (1.0/18.0E-3)*dum(:m)

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_water_content_v
END MODULE ukca_water_content_v_mod
