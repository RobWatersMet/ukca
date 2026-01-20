! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    To give calculate heterogeneous rates for UKCA.
!    Contains the following routines:
!     UKCA_HETERO
!     UKCA_SOLIDPHASE
!     UKCA_CALCKPSC
!     UKCA_EQCOMP
!     UKCA_POSITION
!     UKCA_PSCPRES
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!          Originally used in SLIMCAT CTM.
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_hetero_mod

   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim

   IMPLICIT NONE

   REAL, PARAMETER, PRIVATE  :: rho1 = 1.35     ! Density of type 1 PSCs (g/cm3)
   REAL, PARAMETER, PRIVATE  :: rho2 = 0.928    ! Density of type 2 PSCs (g/cm3)
   REAL, PARAMETER, PRIVATE  :: rad1 = 1.0E-4   ! Radius of type 1 PSCs (cm)
   REAL, PARAMETER, PRIVATE  :: rad2 = 10.0E-4  ! Radius of type 2 PSCs (cm)
   REAL, PARAMETER, PRIVATE  :: radsa = 1.0E-5  ! aerosol radius(cm)

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_HETERO_MOD'

CONTAINS

   SUBROUTINE ukca_hetero(n_points, have_nat, stratflag)
! Description:
!
! Changed from version by Peter Breasicke to allow for dynamical limitation of
! reaction rates as functions of the abundances of the reaction partners.
! Also remove a bug with the ordering of reactions, and allow for changing
! of the order of reactions without recoding. Finally, introduce limits
! on pressure, temperature and latitude where heterogeneous chemistry is
! performed.
!
!
! Method: Heterogeneous reaction rates are specified as pseudo-bimolecular
!         reactions.
!
!           Reactions:
!           1) ClONO2 + HCl -> 2Cl + HONO2
!           2) ClONO2 + H2O -> HOCl + HONO2
!           3) HOCl + HCl   -> H2O + 2Cl
!           4) N2O5 + H2O   -> 2HONO2
!           5) N2O5 + HCl   -> ClNO2 + HONO2
!           6) HOBr + HCl   -> H2O + BrCl
!           7) HCl + BrONO2 -> BrCl + HONO2
!           8) BrONO2 + H2O -> HOBr + HONO2
!           9) HOBr+HBr     -> 2Br + H2O
!           10) HOCl+HBr    -> BrCl + H2O
!           11) ClONO2 + HBr-> BrCl + HONO2
!           12) BrONO2 + HBr-> 2Br + HONO2
!           13) N2O5 + HBr  -> BrNO2 + HONO2
!
! Code Description:
! Language: FORTRAN 90 + common extensions.
!
! Declarations:
! These are of the form:-

      USE asad_mod, ONLY: specf, cdt_diag, f, nhrkx, p, peps, rk, &
                          shno3, sph, sph2o, sphno3, t, tnd, wp, za, &
                          fpsc1, sh2o, jpcspf, jphk
      USE ukca_config_constants_mod, ONLY: avogadro, boltzmann
      USE ukca_constants, ONLY: pi, m_clono2, m_hocl, m_brono2, m_hobr, m_n2o5, &
                                m_h2o, m_hno3
      USE ukca_config_specification_mod, ONLY: ukca_config
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points
! logical to indicate whether natpsc formation is
! allowed at this point (based on height above surface)
      LOGICAL, INTENT(IN) :: have_nat(n_points)
! logical to indicate whether point is in stratosphere
      LOGICAL, INTENT(IN) :: stratflag(n_points)

! Local variables

      LOGICAL, SAVE :: gpsa
      LOGICAL, SAVE :: gphocl
      LOGICAL, SAVE :: gppsc
      LOGICAL, SAVE :: gpsimp
      LOGICAL :: L_ukca_sulphur
      LOGICAL :: L_ukca_presaer

      INTEGER :: js
      INTEGER :: jh
      INTEGER :: jl

      INTEGER :: n_hk

! Tracer names
      INTEGER, SAVE :: ih2o = 0
      INTEGER, SAVE :: ihno3 = 0
      INTEGER, SAVE :: ihcl = 0
      INTEGER, SAVE :: iclono2 = 0
      INTEGER, SAVE :: ihocl = 0
      INTEGER, SAVE :: in2o5 = 0
      INTEGER, SAVE :: ihbr = 0
      INTEGER, SAVE :: ibrono2 = 0
      INTEGER, SAVE :: ihobr = 0

! reaction names
      INTEGER, SAVE :: n_clono2_hcl = 0
      INTEGER, SAVE :: n_clono2_h2o = 0
      INTEGER, SAVE :: n_n2o5_h2o = 0
      INTEGER, SAVE :: n_n2o5_hcl = 0
      INTEGER, SAVE :: n_hocl_hcl = 0
      INTEGER, SAVE :: n_brono2_hcl = 0
      INTEGER, SAVE :: n_brono2_h2o = 0
      INTEGER, SAVE :: n_hobr_hcl = 0
      INTEGER, SAVE :: n_hobr_hbr = 0
      INTEGER, SAVE :: n_hocl_hbr = 0
      INTEGER, SAVE :: n_clono2_hbr = 0
      INTEGER, SAVE :: n_brono2_hbr = 0
      INTEGER, SAVE :: n_n2o5_hbr = 0

      LOGICAL, SAVE :: first = .TRUE.
      LOGICAL, SAVE :: first_pass = .TRUE.

      REAL :: zp(n_points)  ! pressure (hPa)
      REAL :: zt(n_points)  ! temperature (K)

! number density (cm-3)
      REAL :: zhno3(n_points)
      REAL :: zh2o(n_points)
      REAL :: zhcl(n_points)
      REAL :: zclono2(n_points)
      REAL :: zn2o5(n_points)
      REAL :: zhocl(n_points)
      REAL :: zbrono2(n_points)
      REAL :: zhobr(n_points)
      REAL :: zhbr(n_points)

!====================
      REAL :: psc1(n_points)
      REAL :: psc2(n_points)
      REAL :: psc3(n_points)
      REAL :: psc4(n_points)
      REAL :: psc5(n_points)
      REAL :: hk(n_points, 5)
!====================

      REAL :: zrate
      REAL :: zdhcl, zdhbr
      REAL :: zfact

! collision frequency term
      REAL :: c_cf(n_points)

! PSC surface area
      REAL :: psc1sa(n_points)
      REAL :: psc2sa(n_points)

! heterogeneous reaction rate
      REAL :: kpsc(n_points, 13)

      REAL :: gam1(13) ! uptake coefficients on ice and NAT
      REAL :: gam2(13) ! uptake coefficients on ice
      REAL :: gam3(13) ! uptake coefficients on sulphate aerosol

! spatially varying uptake coefficients
      REAL :: gam_calc(n_points, 7) ! for non constant uptake coefficients
      REAL :: gam3calc(n_points, 8) !

      REAL :: mm_arr(13) ! relevent molar mass for each reaction

      REAL :: amu        ! atomic mass unit (kg)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_HETERO'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!! The uptake coefficients for reactions 1-5 on NAT and ice are those
!! recommended in "Chemical Kinetics and Photochemical Data for use in
!! Atmospheric Studies: Evaluation Number 18" from the Jet Propulsion
!! Laboratory, October 2015 (abbrev. here as JPL 15-10). In many cases
!! these values differ from those used prior to UM vn11.2
!! The uptake coefficients for reactions 1-3 on sulfate aerosol are calculated
!! in the routine ukca_shi_liqiud_aerosol based on Shi et al. 2001 (Kinetic
!! model for reaction of ClONO2 with H2O and HCl and HOCl with HCl in sulfuric
!! acid solutions, JGR, 106, 24259-24272, 2001)
!! For reactions 6-13, added at UM vn11.2, uptake coefficients are those
!! recommended by IUPAC (see Crowley et al. 2010, Evaluated kinetic and
!! photochemical data for atmospheric chemistry: Volume V - heterogeneous
!! reactions on solid substrates, Atmos. Chem. Phys., 10 9059-9223,
!! https://doi.org/10.5194/acp-10-9059-2010, 2010). Cases where the
!! coefficient requires calculation are indentified by values greater than one
!! in the following arrays. The value refers to the index in the array
!! gam_calc returned by the routine ukca_br_gamma

!! uptake coefficient on PSC type 1 (NAT)
      gam1 = [0.2, & ! ClONO2 + HCl     (JPL 15-10)
              0.004, & ! ClONO2 + H2O     (JPL 15-10)
              0.1, & ! HOCl + HCl       (JPL 15-10)
              0.0004, & ! N2O5 + H2O       (JPL 15-10)
              0.003, & ! N2O5 + HCl       (JPL 15-10)
              0.0, & ! HOBr + HCl
              0.0, & ! BrONO2 + HCl
              0.0, & ! BrONO2 + H2O
              0.0, & ! HOBr + HBr
              0.0, & ! HOCl + HBr
              6.0, & ! ClONO2 + HBr     (to be calculated, IUPAC)
              0.0, & ! BrONO2 + HBr
              7.0]       ! N2O5 + HBr       (tbc IUPAC) (JPL has 0.005)

! uptake coefficient on PSC type 2 (ice)
      gam2 = [0.3, & ! ClONO2 + HCl     (JPL 15-10)
              0.3, & ! ClONO2 + H2O     (JPL 15-10)
              0.2, & ! HOCl + HCl       (JPL 15-10)
              0.02, & ! N2O5 + H2O       (JPL 15-10)
              0.03, & ! N2O5 + HCl       (JPL 15-10)
              0.25, & ! HOBr + HCl       (IUPAC) (JPL has 0.3)
              0.3, & ! BrONO2 + HCl     (IUPAC)
              1.0, & ! BrONO2 + H2O     (tbc IUPAC)
              2.0, & ! HOBr + HBr       (tbc IUPAC)
              3.0, & ! HOCl + HBr       (tbc IUPAC)
              4.0, & ! ClONO2 + HBr     (tbc IUPAC)
              5.0, & ! BrONO2 + HBr     (tbc IUPAC)
              0.0]      ! N2O5 + HBr

! uptake coefficient on sulfate aerosols
      gam3 = [0.0, & ! ClONO2 + HCl     (tbc Shi et al. 2001)
              0.0, & ! ClONO2 + H2O     (tbc Shi et al. 2001)
              0.0, & ! HOCl + HCl       (tbc Shi et al. 2001)
              0.1, & ! N2O5 + H2O
              0.0, & ! N2O5 + HCl
              0.0, & ! HOBr + HCl
              0.9, & ! BrONO2 + HCl     (JPL 15-10)
              0.0, & ! BrONO2 + H2O     (tbc IUPAC)
              0.0, & ! HOBr + HBr
              0.0, & ! HOCl + HBr
              0.0, & ! ClONO2 + HBr
              0.0, & ! BrONO2 + HBr
              0.0]      ! N2O5 + HBr

! molar mass
      mm_arr = [m_clono2, &
                m_clono2, &
                m_hocl, &
                m_n2o5, &
                m_n2o5, &
                m_hobr, &
                m_brono2, &
                m_brono2, &
                m_hobr, &
                m_hocl, &
                m_clono2, &
                m_brono2, &
                m_n2o5]

! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (ukca_hetero_init)
      IF (first_pass) THEN
         IF (first) THEN
            DO js = 1, jpcspf
               SELECT CASE (specf(js))
               CASE ('H2O       ', 'H2OS      ')
                  ih2o = js
               CASE ('HONO2     ')
                  ihno3 = js
               CASE ('HCl       ')
                  ihcl = js
               CASE ('ClONO2    ')
                  iclono2 = js
               CASE ('N2O5      ')
                  in2o5 = js
               CASE ('HOCl      ')
                  ihocl = js
               CASE ('HBr       ')
                  ihbr = js
               CASE ('BrONO2    ')
                  ibrono2 = js
               CASE ('HOBr      ')
                  ihobr = js
               END SELECT
            END DO

            DO jh = 1, jphk
               SELECT CASE (sph(jh, 1))
               CASE ('H2O       ', 'H2OS      ')
                  IF (sph(jh, 2) == 'ClONO2    ') n_clono2_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'N2O5      ') n_n2o5_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'BrONO2    ') n_brono2_h2o = nhrkx(jh)
               CASE ('ClONO2    ')
                  IF (sph(jh, 2) == 'HCl       ') n_clono2_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'H2O       ') n_clono2_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'H2OS      ') n_clono2_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'HBr       ') n_clono2_hbr = nhrkx(jh)
               CASE ('HCl       ')
                  IF (sph(jh, 2) == 'ClONO2    ') n_clono2_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'N2O5      ') n_n2o5_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'HOCl      ') n_hocl_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'HOBr      ') n_hobr_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'BrONO2    ') n_brono2_hcl = nhrkx(jh)
               CASE ('BrONO2    ')
                  IF (sph(jh, 2) == 'HCl       ') n_brono2_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'H2O       ') n_brono2_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'H2OS      ') n_brono2_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'HBr       ') n_brono2_hbr = nhrkx(jh)
               CASE ('HOBr      ')
                  IF (sph(jh, 2) == 'HCl       ') n_hobr_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'HBr       ') n_hobr_hbr = nhrkx(jh)
               CASE ('N2O5      ')
                  IF (sph(jh, 2) == 'H2O       ') n_n2o5_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'H2OS      ') n_n2o5_h2o = nhrkx(jh)
                  IF (sph(jh, 2) == 'HCl       ') n_n2o5_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'HBr       ') n_n2o5_hbr = nhrkx(jh)
               CASE ('HOCl      ')
                  IF (sph(jh, 2) == 'HCl       ') n_hocl_hcl = nhrkx(jh)
                  IF (sph(jh, 2) == 'HBr       ') n_hocl_hbr = nhrkx(jh)
               CASE ('HBr       ')
                  IF (sph(jh, 2) == 'ClONO2    ') n_clono2_hbr = nhrkx(jh)
                  IF (sph(jh, 2) == 'BrONO2    ') n_brono2_hbr = nhrkx(jh)
                  IF (sph(jh, 2) == 'HOBr      ') n_hobr_hbr = nhrkx(jh)
                  IF (sph(jh, 2) == 'HOCl      ') n_hocl_hbr = nhrkx(jh)
                  IF (sph(jh, 2) == 'N2O5      ') n_n2o5_hbr = nhrkx(jh)
               END SELECT
            END DO

            first = .FALSE.
            l_ukca_sulphur = .FALSE.
            L_ukca_presaer = .TRUE.

            ! activate sulphur chemistry
            gpsa = L_ukca_sulphur .OR. L_ukca_presaer
            ! DO include HCl + HOCl reaction on SA aerosols
            gphocl = .TRUE.
            ! Use heterogeneous chemistry on NAT and ice PSCs
            gppsc = .TRUE.
            ! Use full not simplified scheme for PSCs.
            gpsimp = .FALSE.
         END IF
         first_pass = .FALSE.
      END IF
!$OMP END CRITICAL (ukca_hetero_init)

! pressure in hPa here.
      zp(1:n_points) = p(1:n_points)/100.0
      zt(1:n_points) = t(1:n_points)

! copy tracers. Make sure they have been found correctly.
      IF (ihno3 > 0) THEN
         zhno3 = f(:, ihno3)
      ELSE
         zhno3 = 0.0
      END IF
! if water vapour tracer is not present, use special water
! vapour field.
      IF (ih2o > 0) THEN
         zh2o = f(:, ih2o)
      ELSE
         zh2o = wp
      END IF
      IF (ihcl > 0) THEN
         zhcl = f(:, ihcl)
      ELSE
         zhcl = 0.0
      END IF
      IF (iclono2 > 0) THEN
         zclono2 = f(:, iclono2)
      ELSE
         zclono2 = 0.0
      END IF
      IF (in2o5 > 0) THEN
         zn2o5 = f(:, in2o5)
      ELSE
         zn2o5 = 0.0
      END IF
      IF (ihocl > 0) THEN
         zhocl = f(:, ihocl)
      ELSE
         zhocl = 0.0
      END IF
      IF (ihobr > 0) THEN
         zhobr = f(:, ihobr)
      ELSE
         zhobr = 0.0
      END IF
      IF (ihbr > 0) THEN
         zhbr = f(:, ihbr)
      ELSE
         zhbr = 0.0
      END IF
      IF (ibrono2 > 0) THEN
         zbrono2 = f(:, ibrono2)
      ELSE
         zbrono2 = 0.0
      END IF

! Remove tropospheric ice clouds. They would cause model instability!
      WHERE (.NOT. (stratflag)) sph2o = 0.0
!
! calculate the amount of hno3 and h2o in the solid phase and return
! the residual gas phase concentration
!
      CALL ukca_pscpres(zt(1:n_points), zp(1:n_points), tnd(1:n_points), &
                        zh2o(1:n_points), zhno3(1:n_points), 1, n_points, &
                        n_points, have_nat(1:n_points), sph2o(1:n_points))

      IF (ihno3 > 0) f(:, ihno3) = zhno3
      IF (ih2o > 0) f(:, ih2o) = zh2o

! =====================================================================
! =====================================================================
      IF (ukca_config%i_ukca_hetconfig == 0) THEN
         CALL ukca_calckpsc(za(1:n_points), zt(1:n_points), &
                            zh2o(1:n_points), zhcl(1:n_points), &
                            zclono2(1:n_points), zn2o5(1:n_points), &
                            zhocl(1:n_points), &
                            psc1(1:n_points), psc2(1:n_points), &
                            psc3(1:n_points), psc4(1:n_points), &
                            psc5(1:n_points), gpsa, gphocl, &
                            gppsc, gpsimp, n_points, 1, n_points, cdt_diag)
         !
         ! divide rates by h2o or hcl as asad treats psc reactions as bimolecular
         !
         WHERE (zhcl > peps)
            hk(:, 2) = psc1/zhcl
            hk(:, 3) = psc5/zhcl
            hk(:, 5) = psc4/zhcl
         ELSE WHERE
            hk(:, 2) = 0.0
            hk(:, 3) = 0.0
            hk(:, 5) = 0.0
         END WHERE
         WHERE (zh2o > peps)
            hk(:, 1) = psc2/zh2o
            hk(:, 4) = psc3/zh2o
         ELSE WHERE
            hk(:, 1) = 0.0
            hk(:, 4) = 0.0
         END WHERE
         !
         ! copy the relevant hk's to rk's
         ! Introduce dynamical upper limit. Consider A + B -> C. Throughput through
         ! reaction rk*[A]*[B]*dt should be less than 0.5*min([A],[B])
         ! Also introduce flexible numbering (allow for reordering of reactions
         ! in rath.d
         ! Olaf Morgenstern  18/10/2004
         ! Do not do limiting in the case of non-families chemistry
         !
         ! 1. ClONO2 + H2O --> HOCl + HNO3
         IF (n_clono2_h2o > 0) THEN
            rk(:, n_clono2_h2o) = hk(:, 1)
         END IF

         IF (n_clono2_hcl > 0) THEN
            ! 2. ClONO2 + HCl --> Cl2 + HNO3
            rk(:, n_clono2_hcl) = hk(:, 2)
         END IF

         IF (n_hocl_hcl > 0) THEN
            ! 3. HOCl + HCl --> Cl2 + H2O
            rk(:, n_hocl_hcl) = hk(:, 3)
         END IF

         IF (n_n2o5_h2o > 0) THEN
            ! 4. N2O5 + H2O -> 2 HNO3
            rk(:, n_n2o5_h2o) = hk(:, 4)
         END IF

         IF (n_n2o5_hcl > 0) THEN
            ! 5. N2O5 + HCl -> ClNO2 + HNO3
            rk(:, n_n2o5_hcl) = hk(:, 5)
         END IF
         ! =====================================================================
         ! =====================================================================

      ELSE ! New config

         ! Calculate non-constant gammas
         IF (ukca_config%i_ukca_hetconfig == 2) THEN
            CALL ukca_br_gamma(1, n_points, n_points, zt(1:n_points), zhbr(1:n_points), &
                               gam_calc)
         END IF

         ! Calculate uptake coefficients on SA using formulation of Shi et al. 2001
         IF (ukca_config%i_ukca_hetconfig > 0) THEN
            CALL ukca_shi_liquid_aerosol(1, n_points, n_points, &
                                         zt(1:n_points), zh2o(1:n_points), &
                                         zhcl(1:n_points), zclono2(1:n_points), &
                                         gam3calc(1:n_points, :))
         END IF
         IF (gppsc) THEN     !  If heterogeneous processes are on.

            ! Set conversion factor (atomic mass unit in kg)
            amu = 1.0/(avogadro*1000.0)

            ! calculate reaction rate on psc
            IF (ukca_config%i_ukca_hetconfig == 2) THEN
               n_hk = 13
            ELSE
               n_hk = 5        ! leave out Br reactions
            END IF

            !  1.  INITIALISE RATES TO ZERO
            !      ---------- ----- -- ----
            kpsc(1:n_points, :) = 0.0

            DO jh = 1, n_hk ! loop over reactions
               !---------------------------------------------------------------------
               !          2.  GENERAL TERMS IN COLLISION FREQUENCY EXPRESSION
               !              ------- ----- -- --------- --------- ----------
               c_cf(1:n_points) = 0.25*SQRT((8.0*boltzmann)/ &
                                            (pi*mm_arr(jh)*amu)*zt(1:n_points))
               !
               !---------------------------------------------------------------------
               !          3.  PSC RATES
               !              --- -----
               IF (gppsc) THEN
                  !              3.1  SIMPLE PSC SCHEME
                  IF (gpsimp) THEN
                     !         Zero order PSC rates.
                     kpsc(1:n_points, jh) = 4.6E-5*fpsc1(1:n_points)
                     !
                  ELSE
                     !          3.2  CALCULATE SURFACE AREA OF PSC'S
                     ! sa[m2/m3] = 3*M*u[kg]*1000[g/kg]*nd[cm-3]/
                     !             (rho[g/cm3]*rad[cm]/100[cm/m])
                     !           TYPE 1
                     psc1sa(1:n_points) = (m_hno3 + 3.0*m_h2o)*amu*3.0E5* &
                                          shno3(1:n_points)/(rho1*rad1)
                     !           TYPE 2
                     psc2sa(1:n_points) = m_h2o*amu*3.0E5*sh2o(1:n_points)/(rho2*rad2)
                     !
                     !           Rate on type 1 and 2
                     IF (gam1(jh) >= 1.0) THEN  ! reaction uses calculated gamma
                        kpsc(1:n_points, jh) = c_cf(1:n_points)*psc1sa(1:n_points)* &
                                               gam_calc(1:n_points, INT(gam1(jh)))
                     ELSE
                        kpsc(1:n_points, jh) = c_cf(1:n_points)*psc1sa(1:n_points)*gam1(jh)
                     END IF

                     IF (gam2(jh) >= 1.0) THEN  ! reaction uses calculated gamma
                        kpsc(1:n_points, jh) = kpsc(1:n_points, jh) + &
                                               c_cf(1:n_points)*psc2sa(1:n_points)* &
                                               gam_calc(1:n_points, INT(gam2(jh)))
                     ELSE
                        kpsc(1:n_points, jh) = kpsc(1:n_points, jh) + c_cf(1:n_points)* &
                                               psc2sa(1:n_points)*gam2(jh)
                     END IF
                  END IF
               END IF
               !--------------------------------------------------------------------
               !          4.  AEROSOL REACTIONS
               !              ------- ---------
               !        If required, include the reactions on the sulphate aerosols.
               IF (ukca_config%i_ukca_hetconfig > 0) THEN
                  IF (jh <= 3 .OR. jh == 8) THEN  ! reaction uses calculated gamma
                     kpsc(1:n_points, jh) = kpsc(1:n_points, jh) + c_cf(1:n_points)* &
                                            100.0*za(1:n_points)*gam3calc(1:n_points, jh)
                  ELSE                ! reaction uses constant gamma
                     kpsc(1:n_points, jh) = kpsc(1:n_points, jh) + c_cf(1:n_points)* &
                                            100.0*za(1:n_points)*gam3(jh)
                  END IF
               END IF

            END DO ! loop over reactions

            ! NOTE: The following checks for low HCl/HBr are based on the similar
            ! procedure in the ukca_calckpsc routine. The step
            ! zdhcl = MIN(zrate,zhcl)
            ! is not performing it's intended role here as it compares the fractional
            ! abundance "zrate" to the actual concentration "zhcl". Practically, zhcl
            ! will always be greater so zfact will always equal 1.
            ! This section and the equivalent in ukca_calckpsc should be removed in
            ! the future
            !        5.  CHECK FOR LOW HCl
            !            ----- --- --- ---
            DO jl = 1, n_points
               zrate = MAX(1.0, &
                           cdt_diag*(kpsc(jl, 1)*zclono2(jl) + &
                                     kpsc(jl, 5)*zn2o5(jl) + &
                                     kpsc(jl, 3)*zhocl(jl) + &
                                     kpsc(jl, 6)*zhobr(jl) + &
                                     kpsc(jl, 7)*zbrono2(jl)))
               zdhcl = MIN(zrate, zhcl(jl))
               zfact = zdhcl/zrate

               kpsc(jl, 1) = zfact*kpsc(jl, 1)
               kpsc(jl, 5) = zfact*kpsc(jl, 5)
               kpsc(jl, 3) = zfact*kpsc(jl, 3)
               kpsc(jl, 6) = zfact*kpsc(jl, 6)
               kpsc(jl, 7) = zfact*kpsc(jl, 7)
            END DO

            !       6. CHECK FOR LOW HBr
            !          ----- --- --- ---
            IF (ukca_config%i_ukca_hetconfig == 2) THEN
               DO jl = 1, n_points
                  zrate = MAX(1.0, &
                              cdt_diag*(kpsc(jl, 9)*zhobr(jl) + &
                                        kpsc(jl, 10)*zhocl(jl) + &
                                        kpsc(jl, 11)*zclono2(jl) + &
                                        kpsc(jl, 12)*zbrono2(jl) + &
                                        kpsc(jl, 13)*zn2o5(jl)))
                  zdhbr = MIN(zrate, zhbr(jl))
                  zfact = zdhbr/zrate

                  kpsc(jl, 9) = zfact*kpsc(jl, 9)
                  kpsc(jl, 10) = zfact*kpsc(jl, 10)
                  kpsc(jl, 11) = zfact*kpsc(jl, 11)
                  kpsc(jl, 12) = zfact*kpsc(jl, 12)
                  kpsc(jl, 13) = zfact*kpsc(jl, 13)
               END DO
            END IF
         END IF  !  If heterogeneous processes are on

         ! divide rates by h2o/hcl/hbr as asad treats psc reactions as bimolecular
         WHERE (zhcl > peps)
            rk(:, n_clono2_hcl) = kpsc(:, 1)/zhcl
            rk(:, n_hocl_hcl) = kpsc(:, 3)/zhcl
            rk(:, n_n2o5_hcl) = kpsc(:, 5)/zhcl
            rk(:, n_hobr_hcl) = kpsc(:, 6)/zhcl
            rk(:, n_brono2_hcl) = kpsc(:, 7)/zhcl
         ELSE WHERE
            rk(:, n_clono2_hcl) = 0.0
            rk(:, n_hocl_hcl) = 0.0
            rk(:, n_n2o5_hcl) = 0.0
            rk(:, n_hobr_hcl) = 0.0
            rk(:, n_brono2_hcl) = 0.0
         END WHERE

         WHERE (zh2o > peps)
            rk(:, n_clono2_h2o) = kpsc(:, 2)/zh2o
            rk(:, n_n2o5_h2o) = kpsc(:, 4)/zh2o
            rk(:, n_brono2_h2o) = kpsc(:, 8)/zh2o
         ELSE WHERE
            rk(:, n_clono2_h2o) = 0.0
            rk(:, n_n2o5_h2o) = 0.0
            rk(:, n_brono2_h2o) = 0.0
         END WHERE

         WHERE (zhbr > peps)
            rk(:, n_hobr_hbr) = kpsc(:, 9)/zhbr
            rk(:, n_hocl_hbr) = kpsc(:, 10)/zhbr
            rk(:, n_clono2_hbr) = kpsc(:, 11)/zhbr
            rk(:, n_brono2_hbr) = kpsc(:, 12)/zhbr
            rk(:, n_n2o5_hbr) = kpsc(:, 13)/zhbr
         ELSE WHERE
            rk(:, n_hobr_hbr) = 0.0
            rk(:, n_hocl_hbr) = 0.0
            rk(:, n_clono2_hbr) = 0.0
            rk(:, n_brono2_hbr) = 0.0
            rk(:, n_n2o5_hbr) = 0.0
         END WHERE
      END IF

! save the solid phase hno3 to add back after end of the chemistry timestep
      sphno3 = shno3

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_hetero

!==========================================================================
   SUBROUTINE ukca_br_gamma(kstart, kend, kchmlev, t, thbr, gam_br)
! Calculate gamma on ice and NAT for bromine reactions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: kchmlev            ! no. of chem levels
      INTEGER, INTENT(IN) :: kstart             ! chem level to begin loop
      INTEGER, INTENT(IN) :: kend               ! chem level to end loop

      REAL, INTENT(IN) :: t(kchmlev)            ! temperature (K)
      REAL, INTENT(IN) :: thbr(kchmlev)         ! number density HBr (cm-3)

      REAL, INTENT(OUT) :: gam_br(kchmlev, 7)   ! uptake coefficients

      REAL :: tt                                ! local temperature (k)
      REAL :: theta_hbr                         ! surface coverage HBr

      INTEGER :: jl

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_BR_GAMMA'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Equations are IUPAC recommendations and can be found in
! Crowley et al. 2010, Atmos. Chem. Physics. (ref. to page numbers)
! or at IUPAC web link

      DO jl = kstart, kend
         ! parameterizations defined for temperatures above 180K
         tt = MAX(t(jl), 180.0)

         theta_hbr = 4.14E-10*thbr(jl)**0.88

         ! BrONO2 + H2O on ice (pg 9160)
         gam_br(jl, 1) = 5.3E-4*EXP(1100.0/tt)

         ! HOBr + HBr on ice (pg 9158)
         gam_br(jl, 2) = 4.8E-4*EXP(1240.0/tt)

         ! HOCl + HBr on ice
         ! http://iupac.pole-ether.fr/htdocs/datasheets/pdf/HOCl+HBr_V.A1.42.pdf
         ! (3/10/2018)
         ! For numerical stability we have refactored the equation from
         !  gamma = 1.0 / (alpha + 2.7 / theta_hbr )
         ! to the following:
         gam_br(jl, 3) = theta_hbr/((0.3*theta_hbr) + 2.7)

         ! ClONO2 + HBr on ice (pg 9153)
         gam_br(jl, 4) = 0.56*theta_hbr

         ! BrONO2 + HBr on ice (pg 9163)
         gam_br(jl, 5) = 6.6E-3*EXP(700.0/tt)

         ! ClONO2 + HBr on NAT (pg 9219)
         gam_br(jl, 6) = 0.56*theta_hbr

         ! N2O5 + HBr on NAT (pg 9221)
         gam_br(jl, 7) = 0.02*theta_hbr
      END DO

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_br_gamma

!===========================================================================
   SUBROUTINE ukca_shi_liquid_aerosol(kstart, kend, kchmlev, t, th2o, &
                                      thcl, tcnit, gam3arr)

!*********************************************************************
!subroutine description
!
!new method for calculating uptake coeffiecients of heterogeneous
!reactions on liquid sulphate aerosol.  produces differnt gamma
!values to the pre-existing ukca scheme for hocl+hcl and clono2+h2o,
!and also includes calculations for the clono2+hcl reaction, which
!was not previously included in the pre-existing ukca scheme
!
!written by jmk 09/07/15
!
!*********************************************************************
!technical notes and references
!
!calculations of uptake coeffiecients based on the parameterisation
!of Shi et al. (2001) "Kinetic model for reaction of ClONO2 with H2O
!and HCl and HOCl with HCl in sulfuric acid solutions" J. Geophys. Res.
!parameterisation breaks down at temperatures below 185k, so minimum
!temperature is limited to 185k
!

      USE ukca_config_constants_mod, ONLY: boltzmann
      USE asad_mod, ONLY: peps
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: kchmlev             ! no. of chem. levels
      INTEGER, INTENT(IN) :: kstart              ! chem level to begin loop
      INTEGER, INTENT(IN) :: kend                ! chem level to end loop

      REAL, INTENT(IN) :: t(kchmlev)             ! temperature (K)
      REAL, INTENT(IN) :: th2o(kchmlev)          ! number density H2O (cm-3)
      REAL, INTENT(IN) :: thcl(kchmlev)          !                HCl
      REAL, INTENT(IN) :: tcnit(kchmlev)         !                ClONO2

      REAL, INTENT(OUT) :: gam3arr(kchmlev, 8)   ! uptake coefficients

      REAL, PARAMETER :: r_vol = 0.082             ! gas constant/vol (atm K-1 M-1)

      REAL :: tt               ! local temperature             (K)
      REAL :: p0h2o            ! saturation water vapour       (mbar)
      REAL :: ph2o             ! partial pressure water vapour (mbar)
      REAL :: phcl             ! partial pressure hcl          (atm)
      REAL :: pclono2          ! partial pressure clono2       (atm)
      REAL :: aw               ! water activity                (unitless)
      REAL :: ml_h2so4         ! h2so4 molality                (mol kg-1)
      REAL :: wt_h2so4         ! h2so4 weight percentage       (%)
      REAL :: mr_h2so4         ! h2so4 molarity                (mol l-1)
      REAL :: rho_h2so4        ! h2so4 solution density        (g cm-3)
      REAL :: vis_h2so4        ! viscosity of h2so4 solution   (cp)
      REAL :: a1, b1, c1, d1
      REAL :: a2, b2, c2, d2
      REAL :: z1, z2, z3
      REAL :: mf_h2so4         ! h2so4 mole fraction
      REAL :: a
      REAL :: t0
      REAL :: ah               ! acid activity                 (mol l-1)
      REAL :: y1, y2
      REAL :: k_h2o, k_hydr, k_hcl, k_hocl ! rate constant     (s-1)
      REAL :: k_h              ! rate constant                 (M-1 s-1)
      REAL :: c_clono2, c_hocl ! molecular velocity    (cm s-1)
      REAL :: h_clono2, h_hocl, h_hcl ! Henry coeff     (M atm-1)
      REAL :: s_clono2, s_hocl ! Setchenow coeff               (M-1)
      REAL :: d_clono2, d_hocl ! diffusion coeff (cm2 s-1)
      REAL :: g_h2o, g_hcl, g_s, g_s2, g_hcl2, g_clono2, g_b  ! react. prob.
      REAL :: g_s2_b           ! intermediate value for checking
      REAL :: gam_hocl_f_hcl   ! intermediate value for checking
      REAL :: m_hcl            ! hcl conc.                     (M)
      REAL :: l_clono2, l_hocl ! reacto-diffusive length       (cm)
      REAL :: f_clono2, f_hocl, f_hcl
      REAL :: gam_hocl, gam_clono2 ! reaction prob.
      REAL :: small_value
      REAL, PARAMETER :: min_value = 1.0E-15
      INTEGER :: jl

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SHI_LIQUID_AEROSOL'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! initialise to zero
      gam3arr(:, :) = 0.0

! calculate a small value for comparison in IF blocks. Use same as for rafeps:
      small_value = SQRT(peps)

      DO jl = kstart, kend
         !set temperature lower limit to 185k - below this the
         !parameterisation breaks down
         tt = MAX(t(jl), 185.0)

         !*********************************************************************
         !calculations of h2so4 wt% from temp and ph2o, after tabazadeh et al. (1997)
         !*********************************************************************
         !calculate sturation water vapour pressure (mbar)
         p0h2o = EXP(18.452406985 - (3505.1578807/tt) - &
                     (330918.55082/tt**2) + (12725068.262/tt**3))

         !calculate partial pressure of h2o (mbar)
         ! may have small/negative values of th2o - this is trapped below
         ! -ve values will NOT be passed back to the UM
         ph2o = 1.0E4*th2o(jl)*boltzmann*tt

         !partial pressure hcl and clono2 (atm)
         ! trap for v. low number densities (molec cm-3). THIS VALUE MAY BE UNPHYSICAL
         phcl = 1.0E4*thcl(jl)*boltzmann*tt/1013.0

         pclono2 = 1.0E4*tcnit(jl)*boltzmann*tt/1013.0

         !calculate water activity - trap for small/negative values
         IF (ABS(p0h2o) < small_value) THEN
            aw = 0.0
         ELSE
            aw = ph2o/p0h2o
            IF (aw < min_value) aw = min_value
         END IF

         ! determine constants for subsequent reactions, from tabazadeh et al. (1997)
         IF (aw < 0.85) THEN
            IF (aw <= 0.05) THEN
               a1 = 12.37208932
               b1 = -0.16125516114
               c1 = -30.490657554
               d1 = -2.1133114241
               a2 = 13.455394705
               b2 = -0.1921312255
               c2 = -34.285174607
               d2 = -1.7620073078
            ELSE
               a1 = 11.820654354
               b1 = -0.20786404244
               c1 = -4.807306373
               d1 = -5.1727540348
               a2 = 12.891938068
               b2 = -0.23233847708
               c2 = -6.4261237757
               d2 = -4.9005471319
            END IF
         ELSE
            a1 = -180.06541028
            b1 = -0.38601102592
            c1 = -93.317846778
            d1 = 273.88132245
            a2 = -176.95814097
            b2 = -0.36257048154
            c2 = -90.469744201
            d2 = 267.45509988
         END IF

         y1 = a1*(aw**b1) + (c1*aw) + d1

         y2 = a2*(aw**b2) + (c2*aw) + d2

         !calculate molality h2so4 (mol kg-1)
         ml_h2so4 = y1 + (((tt - 190.0)*(y2 - y1))/70.0)

         ! MOLALITY CANNOT BE -VE
         ml_h2so4 = MAX(ml_h2so4, small_value)

         !calculate weight % h2so4
         wt_h2so4 = (9800.0*ml_h2so4)/((98.0*ml_h2so4) + 1000.0)

         !*********************************************************************
         !parameters for the h2so4 solution
         !*********************************************************************

         !set up terms for h2so4 solution density calculation
         z1 = 0.12364 - (5.6E-7*(tt**2.0))

         z2 = -0.02954 + (1.814E-7*(tt**2.0))

         z3 = 2.343E-3 - (1.487E-6*tt) - (1.324E-8*(tt**2.0))

         !calculate h2so4 solution density (g cm-3)
         rho_h2so4 = 1.0 + (z1*ml_h2so4) + (z2*(ml_h2so4**1.5)) + &
                     (z3*(ml_h2so4**2.0))

         !calculate molarity h2so4 (mol l-1)
         mr_h2so4 = (rho_h2so4*wt_h2so4)/9.8

         !calculate h2so4 mole fraction
         mf_h2so4 = wt_h2so4/(wt_h2so4 + ((100.0 - wt_h2so4)*(98.0/18.0)))

         !set up terms for h2so4 viscosity calculation
         IF (wt_h2so4 <= 80.0) THEN ! parameterisation valid for wt% < 80
            a = 169.5 + (5.18*wt_h2so4) - (0.0825*(wt_h2so4**2.0)) + &
                (3.27E-3*(wt_h2so4**3.0))

            t0 = 144.11 + (0.166*wt_h2so4) - (0.015*(wt_h2so4**2.0)) + &
                 (2.18E-4*(wt_h2so4**3.0))

            !calculate viscosity of h2so4 solution
            vis_h2so4 = a*(tt**(-1.43))*EXP(448.0/(tt - t0))

         ELSE ! assume wt% = 80
            vis_h2so4 = 1730.14*(tt**(-1.43))*EXP(448.0/(tt - 173.006))

         END IF

         !calculate acid activity
         ah = EXP(60.51 - (0.095*wt_h2so4) + (0.0077*(wt_h2so4**2.0)) &
                  - (1.61E-5*(wt_h2so4**3.0)) &
                  - ((1.76 + (2.52E-4*(wt_h2so4**2.0)))*(tt**0.5)) &
                  + ((-805.89 + (253.05*(wt_h2so4**0.076)))/(tt**0.5)))

         !*********************************************************************
         !calculate uptake parameters for clono2+h2o and clono2+hcl
         !*********************************************************************

         k_h2o = 1.95E10*EXP(-2800.0/tt)

         k_h = 1.22E12*EXP(-6200.0/tt)

         k_hydr = (k_h2o*aw) + (k_h*ah*aw)

         c_clono2 = 1474.0*(tt**0.5)

         s_clono2 = 0.306 + (24.0/tt)

         h_clono2 = 1.6E-6*EXP((4710.0/tt) - (s_clono2*mr_h2so4))

         d_clono2 = 5E-8*tt/vis_h2so4

         g_h2o = (4*h_clono2*r_vol*tt/c_clono2)*(d_clono2*k_hydr)**0.5

         h_hcl = (0.094 - 0.61*mf_h2so4 + 1.2*mf_h2so4**2)* &
                 EXP(-8.68 + (8515.0 - 10718.0*mf_h2so4**0.7)/tt)

         m_hcl = h_hcl*phcl

         k_hcl = 7.9E11*ah*d_clono2*m_hcl

         l_clono2 = (d_clono2/(k_hydr + k_hcl))**0.5

         f_clono2 = (1/(TANH(radsa/l_clono2))) - l_clono2/radsa

         g_clono2 = f_clono2*g_h2o*((1 + (k_hcl/k_hydr))**0.5)

         g_hcl = (g_clono2*k_hcl)/(k_hcl + k_hydr)

         g_s = 66.12*EXP(-1374.0/tt)*h_clono2*m_hcl

         f_hcl = 1/(1 + (0.612*(g_s + g_hcl)*(pclono2/phcl)))

         g_s2 = f_hcl*g_s

         g_hcl2 = f_hcl*g_hcl

         g_b = g_hcl2 + ((g_clono2*k_hydr)/(k_hcl + k_hydr))

         ! traps for small value division
         IF (ABS(g_s2 + g_b) < small_value) THEN
            ! We've got a divide by zero problem
            gam3arr(jl, 1) = 0.0
            gam3arr(jl, 2) = 0.0
         ELSE
            ! calculate this first to test against
            g_s2_b = (1.0 + (1.0/(g_s2 + g_b)))
            IF (ABS(g_s2_b) < small_value) THEN

               gam3arr(jl, 1) = 0.0
               gam3arr(jl, 2) = 0.0
            ELSE
               gam_clono2 = 1.0/g_s2_b

               gam3arr(jl, 1) = gam_clono2*((g_s2 + g_hcl2)/(g_s2 + g_b))

               gam3arr(jl, 2) = gam_clono2 - gam3arr(jl, 1)

               ! make sure at a minimum value
               IF (gam3arr(jl, 1) < min_value) gam3arr(jl, 1) = min_value

               IF (gam3arr(jl, 2) < min_value) gam3arr(jl, 2) = min_value
            END IF
         END IF

         !*********************************************************************
         !calculate uptake parameters for hocl+hcl
         !*********************************************************************

         s_hocl = 0.0776 + 59.18/tt

         h_hocl = 1.91E-6*EXP((5862.4/tt) - (s_hocl*mr_h2so4))

         d_hocl = 6.4E-8*tt/vis_h2so4

         k_hocl = 1.25E9*ah*d_hocl*m_hcl

         c_hocl = 2009.0*(tt**0.5)

         l_hocl = (d_hocl/k_hocl)**0.5

         f_hocl = (1/(TANH(radsa/l_hocl))) - l_hocl/radsa

         ! check for small value division
         IF ((ABS(f_hocl) < small_value) .OR. (ABS(f_hcl) < small_value)) THEN
            gam3arr(jl, 3) = 0.0
         ELSE
            gam_hocl = (f_hocl*4.0*h_hocl*r_vol*tt/c_hocl)*(d_hocl*k_hocl)**0.5

            ! trap for small value division
            IF (gam_hocl*f_hcl < small_value) THEN
               gam3arr(jl, 3) = 0.0
            ELSE
               gam_hocl_f_hcl = (1.0 + 1.0/(gam_hocl*f_hcl))
               gam3arr(jl, 3) = 1.0/gam_hocl_f_hcl

               ! make sure at a minimum value
               IF (gam3arr(jl, 3) < min_value) gam3arr(jl, 3) = min_value
            END IF
         END IF

         ! Additional reactions on sulfate aerosol
         ! BrONO2 + H2O
         gam3arr(jl, 8) = 1.0/(1.0/0.8 + 1.0/(EXP(29.2 - 0.4*wt_h2so4) + 0.11))

      END DO !(jl=kstart, kend)

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_shi_liquid_aerosol

! ######################################################################
   SUBROUTINE ukca_solidphase(n_points)

! Description:
!
! *** *solidphase* - adds HONO2 and H2O in solid state back to main
!                    arrays.
!
!     Purpose.
!     --------
!         The amount of HONO2 and H2O in the solid phase when PSCs are
!     on is calculated and stored during a chemical step. Only the
!     amount still in the gas phase is used to calc. species tendencies.
!
!     Interface
!     ---------
!         This routine *MUST* be called at the end of each chemical
!     substep to add the solid phase HONO2 and H2O back to the main
!     ASAD arrays.
!
!     Method
!     ------
!          See comments in routine 'hetero'.
!
! The return of ice to the gasphase is disactivated if UM_ICE is set
! because the UM has an explicit ice tracer which we don't want to
! affect. Also, returning HNO3S is disactivated in case NAT PSC
! sedimentation is selected because that is done outside of ASAD.
!
!
! Code Description:
! Language: FORTRAN 90 + common extensions.
!
! Declarations:
! These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
!---------------------------------------------------------------------
!

      USE asad_mod, ONLY: f, sphno3, specf, jpcspf
      USE ereport_mod, ONLY: ereport

      USE errormessagelength_mod, ONLY: errormessagelength
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points

! local variables
      CHARACTER(LEN=errormessagelength) :: cmessage
      INTEGER :: errcode                ! Variable passed to ereport
      INTEGER :: js
      INTEGER, SAVE :: ihno3 = 0
      LOGICAL, SAVE :: firstcall = .TRUE.
      LOGICAL, SAVE :: first_pass = .TRUE.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SOLIDPHASE'

!
!     ----------------------------------------------------------------
!           1.  Add amount of HONO2 back to main ASAD arrays.
!               --- ------ -- ----- ---- -- ---- ---- -------
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (ukca_solidphase_init)
      IF (first_pass) THEN
         IF (firstcall) THEN
            DO js = 1, jpcspf
               IF (specf(js) == 'HONO2     ') ihno3 = js
            END DO
            IF (ihno3 == 0) THEN
               errcode = 1
               cmessage = 'Select HONO2 as advected tracer.'
               CALL ereport('SOLIDPHASE', errcode, cmessage)
            END IF
            firstcall = .FALSE.
         END IF
         first_pass = .FALSE.
      END IF
!$OMP END CRITICAL (ukca_solidphase_init)

      f(1:n_points, ihno3) = f(1:n_points, ihno3) + &
                             sphno3(1:n_points)

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_solidphase

! ######################################################################
   SUBROUTINE ukca_calckpsc(sasa, t, th2o, thcl, tcnit, tn2o5, thocl, &
                            akpsc1, akpsc2, akpsc3, akpsc4, akpsc5, &
                            lpsa, lphocl, lppsc, lpsimp, &
                            kchmlev, kstart, kend, dt)
!
!     CALCKPSC - CALCULATION OF HETEROGENEOUS REACTION RATES
!
!         PETER GOOD, UNIVERSITY OF CAMBRIDGE, 9/11/93
!         BASED ON CODE BY LEE GRENFELL AND MARTYN CHIPPERFIELD
!
!     PURPOSE
!     -------
!         THIS ROUTINE CALCULATES FIRST ORDER RATES FOR REACTIONS
!     OCCURING ON TYPE 1 AND 2 PSC'S, AND AEROSOL.
!
!     INTERFACE
!     ---------
!         ARGUMENTS IN :  SASA   - (SULPHATE) AEROSOL SURFACE AREA
!                                  PER UNIT VOLUME (cm2 cm-3)
!                         T      - TEMPERATURE
!                         TH2O   - H2O NUMBER DENSITY (cm-3)
!                         THCL   - HCl NUMBER DENSITY (cm-3)
!                         TCNIT  - ClONO2 NUMBER DENSITY (cm-3)
!                         TN2O5  - N2O5 NUMBER DENSITY (cm-3)
!                         THOCL  - HOCl NUMBER DENSITY (cm-3)
!                         LPSA   - IF .TRUE., AEROSOL REACTIONS ARE ON
!                         LPHOCL - IF .TRUE., HOCL+HCL OCCURS ON AEROSOL
!                         LPPSC  - IF .TRUE., PSC REACTIONS ARE ON
!                         LPSIMP - IF .TRUE., USE SIMPLE PSC SCHEM
!                         KCHMLEV- LEVEL DIMENSION OF ARRAYS
!                         KSTART - TOP LEVEL OF CHEMISTRY
!                         KEND   - BOTTOM LEVEL OF CHEMISTRY
!                         DT     - MODEL TIMESTEP (s)
!
!     RESULTS
!     -------
!         FIRST ORDER RATE COEFFICIENTS AKPSC1 ... AKPSC5:
!
!           AKPSC1: HCl + ClONO2 -> Cl2 + HNO3 -> 2ClOX + HNO3
!           AKPSC2: ClONO2 + H2O -> HOCl + HNO3
!           AKPSC3: N2O5 + H2O   -> 2HNO3
!           AKPSC4: N2O5 + HCl   -> ClNO2 + HNO3
!           AKPSC5: HOCL + HCL   -> H2O + CL2    -> H2O + 2CLOX
!
!-----------------------------------------------------------------------
!
!     METHOD NOTES
!     ------ -----
!
!         AEROSOL REACTIONS 2 AND 5:
!
!              AKPSC2: ClONO2 + H2O -> HOCl + HNO3
!              The rate of this reaction is a strong function of the
!              sulphate acid wt.percent, and hence of temperature.
!
!              AKPSC5: HOCL + HCL -> CL2 + H2O
!              The rate of this reaction is proportional to the Henry's
!              law coefficients for HCl(aq) + Cl-(aq) and HOCl(aq)
!              That for HCl(aq) + Cl-(aq) is estimated from the sulphuri
!              acid wt.percent and temperature; the other is only known
!              for 60% w/w H2SO4.
!
!        See Cox, MacKenzie, Muller, Peter, Crutzen 1993
!
!----------------------------------------------------------------------
!
!     The collision frequency, v, of gas molecules with the
!     reacting surface is calculated using:
!
!     v=(A/4)*SQRT(8kT/pi*M)
!
!     A, Reacting surface concentration per unit volume.
!     k, Boltzmann's constant.
!     T, Temperature in kelvin.
!     M, Molecular mass of air molecules (=RMM*U).
!
!-----------------------------------------------------------------------
!
      USE asad_mod, ONLY: fpsc1, shno3, sh2o
      USE ukca_config_constants_mod, ONLY: avogadro, boltzmann
      USE ukca_constants, ONLY: pi

      IMPLICIT NONE

! Subroutine interface
      LOGICAL, INTENT(IN) :: lpsa
      LOGICAL, INTENT(IN) :: lphocl
      LOGICAL, INTENT(IN) :: lppsc
      LOGICAL, INTENT(IN) :: lpsimp

      INTEGER, INTENT(IN) :: kchmlev
      INTEGER, INTENT(IN) :: kstart
      INTEGER, INTENT(IN) :: kend

      REAL, INTENT(IN)  :: dt
      REAL, INTENT(IN)  :: t(kchmlev)
      REAL, INTENT(IN)  :: sasa(kchmlev)
      REAL, INTENT(IN)  :: thcl(kchmlev)
      REAL, INTENT(IN)  :: th2o(kchmlev)
      REAL, INTENT(IN)  :: tcnit(kchmlev)
      REAL, INTENT(IN)  :: tn2o5(kchmlev)
      REAL, INTENT(IN)  :: thocl(kchmlev)

      REAL, INTENT(OUT) :: akpsc1(kchmlev)
      REAL, INTENT(OUT) :: akpsc2(kchmlev)
      REAL, INTENT(OUT) :: akpsc3(kchmlev)
      REAL, INTENT(OUT) :: akpsc4(kchmlev)
      REAL, INTENT(OUT) :: akpsc5(kchmlev)

! Local variables
      REAL, PARAMETER :: u = 1.66056E-27
!
!     Bimolecular rate coefficient (HOCl+HCl*)(aq) (mol-1 m3 s-1)
      REAL, PARAMETER :: cpk1 = 1.0E+2
!
      REAL, PARAMETER :: gam1a = 0.3
      REAL, PARAMETER :: gam1b = 0.006
      REAL, PARAMETER :: gam1c = 0.0006
      REAL, PARAMETER :: gam1d = 0.003
      REAL, PARAMETER :: gam1e = 0.3
      REAL, PARAMETER :: gam2a = 0.3
      REAL, PARAMETER :: gam2b = 0.3
      REAL, PARAMETER :: gam2c = 0.03
      REAL, PARAMETER :: gam2d = 0.03
      REAL, PARAMETER :: gam2e = 0.3
      REAL, PARAMETER :: gam3c = 0.1
!
      INTEGER :: jl

      REAL :: zfcnit
      REAL :: zfn2o5
      REAL :: zfhocl
      REAL :: vol
      REAL :: order2
      REAL :: zrate
      REAL :: zdhcl
      REAL :: zfact

      REAL :: gam3b(kchmlev)
      REAL :: hshcl(kchmlev)
      REAL :: hhocl(kchmlev)
      REAL :: ccnit(kchmlev)
      REAL :: cn2o5(kchmlev)
      REAL :: chocl(kchmlev)
      REAL :: psc1sa(kchmlev)
      REAL :: psc2sa(kchmlev)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_CALCKPSC'

!
!
!-----------------------------------------------------------------------
!          1.  INITIALISE RATES TO ZERO
!              ---------- ----- -- ----
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
      akpsc1(kstart:kend) = 0.0
      akpsc2(kstart:kend) = 0.0
      akpsc3(kstart:kend) = 0.0
      akpsc4(kstart:kend) = 0.0
      akpsc5(kstart:kend) = 0.0
!
!     If heterogeneous processes are on.
      IF (lppsc .OR. lpsa) THEN
         !
         !-----------------------------------------------------------------------
         !          2.  GENERAL TERMS IN COLLISION FREQUENCY EXPRESSION
         !              ------- ----- -- --------- --------- ----------
         !
         zfcnit = (8.0*boltzmann)/(pi*97.5*u)
         zfn2o5 = (8.0*boltzmann)/(pi*108.0*u)
         zfhocl = (8.0*boltzmann)/(pi*52.5*u)
         !
         ccnit(kstart:kend) = 0.25*SQRT(zfcnit*t(kstart:kend))
         cn2o5(kstart:kend) = 0.25*SQRT(zfn2o5*t(kstart:kend))
         chocl(kstart:kend) = 0.25*SQRT(zfhocl*t(kstart:kend))
         !
         !-----------------------------------------------------------------------
         !          3.  PSC RATES
         !              --- -----
         !
         IF (lppsc) THEN
            !
            !              3.1  SIMPLE PSC SCHEME
            !
            IF (lpsimp) THEN
               !
               !         Zero order PSC rates.
               akpsc1(kstart:kend) = 4.6E-5*fpsc1(kstart:kend)
               akpsc2(kstart:kend) = 4.6E-5*fpsc1(kstart:kend)
               akpsc3(kstart:kend) = 4.6E-5*fpsc1(kstart:kend)
               !
            ELSE
               !
               !              3.2  CALCULATE SURFACE AREA OF PSC'S
               !
               !           TYPE 1
               psc1sa(kstart:kend) = &
                  1.85*63.0*u*3.0E5*shno3(kstart:kend)/(rho1*rad1)
               !           TYPE 2
               psc2sa(kstart:kend) = &
                  18.0*u*3.0E5*sh2o(kstart:kend)/(rho2*rad2)
               !
               akpsc1(kstart:kend) = &
                  ccnit(kstart:kend)*(psc1sa(kstart:kend)*gam1a + &
                                      psc2sa(kstart:kend)*gam2a)
               akpsc2(kstart:kend) = &
                  ccnit(kstart:kend)*(psc1sa(kstart:kend)*gam1b + &
                                      psc2sa(kstart:kend)*gam2b)
               akpsc3(kstart:kend) = &
                  cn2o5(kstart:kend)*(psc1sa(kstart:kend)*gam1c + &
                                      psc2sa(kstart:kend)*gam2c)
               akpsc4(kstart:kend) = &
                  cn2o5(kstart:kend)*(psc1sa(kstart:kend)*gam1d + &
                                      psc2sa(kstart:kend)*gam2d)
               akpsc5(kstart:kend) = &
                  chocl(kstart:kend)*(psc1sa(kstart:kend)*gam1e + &
                                      psc2sa(kstart:kend)*gam2e)
               !
            END IF
            !
         END IF
         !
         !----------------------------------------------------------------------
         !          4.  AEROSOL REACTIONS
         !              ------- ---------
         !
         !        If required, include the reactions on the sulphate aerosols.
         IF (lpsa) THEN
            !
            CALL ukca_eqcomp(t, th2o, kstart, kend, kchmlev, lphocl, &
                             gam3b, hshcl, hhocl)
            !
            akpsc2(kstart:kend) = akpsc2(kstart:kend) + &
                                  ccnit(kstart:kend)*100.0*sasa(kstart:kend)*gam3b(kstart:kend)
            akpsc3(kstart:kend) = akpsc3(kstart:kend) + &
                                  cn2o5(kstart:kend)*100.0*sasa(kstart:kend)*gam3c
            !
            !           If required, include HOCl + HCl -> Cl2 + H2O
            IF (lphocl) THEN
               DO jl = kstart, kend
                  !
                  !             Estimate specific volume from surface area
                  vol = sasa(jl)*radsa/3.0
                  !
                  !             Add aerosol rate, converted to pseudo first order
                  order2 = 1.0E6*cpk1*vol*avogadro*hshcl(jl)*hhocl(jl)* &
                           (boltzmann*t(jl))**2
                  akpsc5(jl) = akpsc5(jl) + order2*thcl(jl)
               END DO
            END IF
            !
            !----------------------------------------------------------------------
            !
         END IF
         !
         !        5.  CHECK FOR LOW HCL
         !            ----- --- --- ---
         !
         DO jl = kstart, kend
            zrate = MAX(1.0, &
                        dt*(akpsc1(jl)*tcnit(jl) + &
                            akpsc4(jl)*tn2o5(jl) + &
                            akpsc5(jl)*thocl(jl)))
            zdhcl = MIN(zrate, thcl(jl))
            zfact = zdhcl/zrate
            akpsc1(jl) = zfact*akpsc1(jl)
            akpsc4(jl) = zfact*akpsc4(jl)
            akpsc5(jl) = zfact*akpsc5(jl)
         END DO
         !
         !----------------------------------------------------------------------
         !
      END IF

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_calckpsc

! ######################################################################
   SUBROUTINE ukca_eqcomp(t, th2o, kstart, kend, kchmlev, lphocl, &
                          rpcncl, hshcl, hhocl)
!
!-----------------------------------------------------------------------
!
!     EQCOMP    - CALCULATIONS INVOLVING EQUILIBRIUM AEROSOL COMPOSITION
!
!            PETER GOOD, UNIVERSITY OF CAMBRIDGE, 10/11/93
!            BASED ON CODE BY SLIMANE BEKKI
!
!     PURPOSE
!     -------
!         CALCULATE THOSE TERMS IN THE AEROSOL REACTION RATE EXPRESSIONS
!     WHICH DEPEND ON THE AEROSOL COMPOSITION.
!
!     INTERFACE
!     ---------
!         ARGUMENTS IN :
!              T   -    TEMPERATURE (K)
!              TH2O   - H2O NUMBER DENSITY (cm-3)
!              KSTART - INDEX OF FIRST CHEMSITRY LEVEL
!              KEND   - LAST CHEMISTRY LEVEL
!              KCHMLEV- LEVEL DIMENSION OF ARRAYS
!              LPHOCL - IF .TRUE. THEN HOCL+HCL IN AEROSOL IS TURNED ON
!
!     RESULTS
!     -------
!              RPCNCL - REACTION PROBABILITY FOR ClONO2 + H2O
!              HSHCL  - MODIFIED HENRY'S LAW COEFFICIENT (mol/Nm)
!              HHOCL  - HENRY'S LAW COEFFICIENT
!
!     REFERENCES
!     ----------
!         HAMILL & STEELE  (TABLE FOR RPCNCL)
!         ZHANG ET AL. 1993   (DATA FOR EVALUATING HSHCL AND HHOCL)
!
!-----------------------------------------------------------------------
!
      USE ukca_config_constants_mod, ONLY: boltzmann

      IMPLICIT NONE

! Subroutine interface

      INTEGER, INTENT(IN)  :: kchmlev
      INTEGER, INTENT(IN)  :: kstart
      INTEGER, INTENT(IN)  :: kend
      LOGICAL, INTENT(IN)  :: lphocl
      REAL, INTENT(IN)  :: t(kchmlev)
      REAL, INTENT(IN)  :: th2o(kchmlev)
      REAL, INTENT(OUT) :: hshcl(kchmlev)
      REAL, INTENT(OUT) :: rpcncl(kchmlev)
      REAL, INTENT(OUT) :: hhocl(kchmlev)

! Local variables
      INTEGER, PARAMETER :: noph2o = 16
      INTEGER, PARAMETER :: notemp = 28
      INTEGER, PARAMETER :: nocomp = 4
      INTEGER :: jxp1
      INTEGER :: jyp1
      INTEGER :: jx
      INTEGER :: jy
      INTEGER :: ierx
      INTEGER :: iery
      INTEGER :: i
      INTEGER :: jl
      REAL :: sxy
      REAL :: sx1y
      REAL :: sx1y1
      REAL :: sxy1
      REAL :: ta
      REAL :: tb
      REAL :: tt
      REAL :: ua
      REAL :: ub
      REAL :: u
      REAL :: a
      REAL :: b
      REAL :: ajx
      REAL :: ajxp1
      REAL :: bjx
      REAL :: bjxp1
      REAL :: ph2o
!
      REAL :: comp(kchmlev)
      REAL :: compindx(nocomp)
      REAL :: aindex(nocomp)
      REAL :: bindex(nocomp)
      REAL :: pindx(noph2o)
      REAL :: tindx(notemp)
      REAL :: c(noph2o, notemp)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_EQCOMP'

!
!     Composition index for HSHCL
      DATA compindx/35.0, 40.0, 50.0, 60.0/
!
!     Intercept and gradient look-up table for HSHCL
      DATA aindex/-11.35, -11.95, -13.80, -14.68/
      DATA bindex/6.10E+03, 5.88E+03, 5.53E+03, 4.92E+03/
!
!     P:Ambient partial pressure of water (mb)
      DATA pindx/0.1000E-05, 0.1000E-04, 0.5000E-04, 0.1000E-03, &
         0.1500E-03, 0.2000E-03, 0.3000E-03, 0.5000E-03, &
         0.6000E-03, 0.8000E-03, 0.1000E-02, 0.1200E-02, &
         0.1500E-02, 0.2000E-02, 0.3000E-02, 0.1000E-01/
!
!     T:Ambient temperature (K)
      DATA tindx/0.1750E+03, 0.1800E+03, 0.1850E+03, 0.1900E+03, &
         0.1950E+03, 0.2000E+03, 0.2050E+03, 0.2100E+03, &
         0.2150E+03, 0.2200E+03, 0.2250E+03, 0.2300E+03, &
         0.2350E+03, 0.2400E+03, 0.2450E+03, 0.2500E+03, &
         0.2550E+03, 0.2600E+03, 0.2650E+03, 0.2700E+03, &
         0.2750E+03, 0.2800E+03, 0.2850E+03, 0.2900E+03, &
         0.2950E+03, 0.3000E+03, 0.3050E+03, 0.3100E+03/
!
!     C:Composition (wt)
      DATA(c(i, 1), i=1, 16)/ &
         0.9400E+02, 0.7800E+02, 0.5957E+02, 0.4363E+02, &
         0.3430E+02, 0.2768E+02, 0.1836E+02, 0.4000E+02, &
         0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, &
         0.4000E+02, 0.4000E+02, 0.4000E+02, 0.9000E+01/
      DATA(c(i, 2), i=1, 16)/ &
         0.9400E+02, 0.7900E+02, 0.6152E+02, 0.4723E+02, &
         0.3887E+02, 0.3294E+02, 0.2458E+02, 0.1272E+02, &
         0.4000E+02, 0.4000E+02, 0.1142E+02, 0.4000E+02, &
         0.1341E+02, 0.4000E+02, 0.4000E+02, 0.9000E+01/
      DATA(c(i, 3), i=1, 16)/ &
         0.9400E+02, 0.8000E+02, 0.6348E+02, 0.5084E+02, &
         0.4344E+02, 0.3820E+02, 0.3080E+02, 0.1929E+02, &
         0.1196E+02, 0.1681E+02, 0.1728E+02, 0.1300E+02, &
         0.1928E+02, 0.1442E+02, 0.4000E+02, 0.9000E+01/
      DATA(c(i, 4), i=1, 16)/ &
         0.9400E+02, 0.8100E+02, 0.6543E+02, 0.5444E+02, &
         0.4801E+02, 0.4345E+02, 0.3702E+02, 0.2585E+02, &
         0.1538E+02, 0.2541E+02, 0.2315E+02, 0.1806E+02, &
         0.2515E+02, 0.1988E+02, 0.1615E+02, 0.9000E+01/
      DATA(c(i, 5), i=1, 16)/ &
         0.9400E+02, 0.8200E+02, 0.6935E+02, 0.6165E+02, &
         0.5715E+02, 0.5396E+02, 0.4946E+02, 0.4226E+02, &
         0.3935E+02, 0.3402E+02, 0.2902E+02, 0.2313E+02, &
         0.3102E+02, 0.2535E+02, 0.2334E+02, 0.9000E+01/
      DATA(c(i, 6), i=1, 16)/ &
         0.9400E+02, 0.8300E+02, 0.7125E+02, 0.6594E+02, &
         0.6283E+02, 0.6062E+02, 0.5751E+02, 0.5278E+02, &
         0.5073E+02, 0.4693E+02, 0.4369E+02, 0.4086E+02, &
         0.3689E+02, 0.3082E+02, 0.3052E+02, 0.9000E+01/
      DATA(c(i, 7), i=1, 16)/ &
         0.9400E+02, 0.8367E+02, 0.7395E+02, 0.6976E+02, &
         0.6731E+02, 0.6557E+02, 0.6312E+02, 0.5955E+02, &
         0.5811E+02, 0.5561E+02, 0.5344E+02, 0.5144E+02, &
         0.4863E+02, 0.4449E+02, 0.3771E+02, 0.9000E+01/
      DATA(c(i, 8), i=1, 16)/ &
         0.9400E+02, 0.8420E+02, 0.7626E+02, 0.7284E+02, &
         0.7084E+02, 0.6942E+02, 0.6742E+02, 0.6455E+02, &
         0.6341E+02, 0.6147E+02, 0.5983E+02, 0.5838E+02, &
         0.5646E+02, 0.5369E+02, 0.4849E+02, 0.1209E+02/
      DATA(c(i, 9), i=1, 16)/ &
         0.9400E+02, 0.8519E+02, 0.7841E+02, 0.7548E+02, &
         0.7377E+02, 0.7256E+02, 0.7085E+02, 0.6845E+02, &
         0.6752E+02, 0.6594E+02, 0.6462E+02, 0.6347E+02, &
         0.6196E+02, 0.5983E+02, 0.5640E+02, 0.3239E+02/
      DATA(c(i, 10), i=1, 16)/ &
         0.9400E+02, 0.8603E+02, 0.8020E+02, 0.7768E+02, &
         0.7621E+02, 0.7517E+02, 0.7370E+02, 0.7163E+02, &
         0.7083E+02, 0.6949E+02, 0.6839E+02, 0.6743E+02, &
         0.6619E+02, 0.6447E+02, 0.6175E+02, 0.4271E+02/
      DATA(c(i, 11), i=1, 16)/ &
         0.9424E+02, 0.8691E+02, 0.8179E+02, 0.7959E+02, &
         0.7830E+02, 0.7738E+02, 0.7609E+02, 0.7429E+02, &
         0.7360E+02, 0.7244E+02, 0.7148E+02, 0.7066E+02, &
         0.6960E+02, 0.6815E+02, 0.6589E+02, 0.5007E+02/
      DATA(c(i, 12), i=1, 16)/ &
         0.9433E+02, 0.8780E+02, 0.8323E+02, 0.8127E+02, &
         0.8012E+02, 0.7930E+02, 0.7815E+02, 0.7656E+02, &
         0.7595E+02, 0.7493E+02, 0.7410E+02, 0.7338E+02, &
         0.7245E+02, 0.7119E+02, 0.6925E+02, 0.5567E+02/
      DATA(c(i, 13), i=1, 16)/ &
         0.9445E+02, 0.8860E+02, 0.8451E+02, 0.8275E+02, &
         0.8172E+02, 0.8099E+02, 0.7996E+02, 0.7853E+02, &
         0.7798E+02, 0.7708E+02, 0.7633E+02, 0.7570E+02, &
         0.7489E+02, 0.7377E+02, 0.7207E+02, 0.6017E+02/
      DATA(c(i, 14), i=1, 16)/ &
         0.9478E+02, 0.8945E+02, 0.8571E+02, 0.8411E+02, &
         0.8317E+02, 0.8250E+02, 0.8156E+02, 0.8027E+02, &
         0.7977E+02, 0.7896E+02, 0.7829E+02, 0.7772E+02, &
         0.7699E+02, 0.7600E+02, 0.7449E+02, 0.6392E+02/
      DATA(c(i, 15), i=1, 16)/ &
         0.9568E+02, 0.9057E+02, 0.8700E+02, 0.8546E+02, &
         0.8456E+02, 0.8392E+02, 0.8302E+02, 0.8183E+02, &
         0.8138E+02, 0.8063E+02, 0.8002E+02, 0.7951E+02, &
         0.7885E+02, 0.7795E+02, 0.7659E+02, 0.6707E+02/
      DATA(c(i, 16), i=1, 16)/ &
         0.9695E+02, 0.9190E+02, 0.8836E+02, 0.8684E+02, &
         0.8595E+02, 0.8532E+02, 0.8443E+02, 0.8327E+02, &
         0.8284E+02, 0.8215E+02, 0.8158E+02, 0.8111E+02, &
         0.8050E+02, 0.7969E+02, 0.7845E+02, 0.6977E+02/
      DATA(c(i, 17), i=1, 16)/ &
         0.9900E+02, 0.9374E+02, 0.9000E+02, 0.8840E+02, &
         0.8746E+02, 0.8679E+02, 0.8585E+02, 0.8467E+02, &
         0.8425E+02, 0.8357E+02, 0.8303E+02, 0.8258E+02, &
         0.8202E+02, 0.8126E+02, 0.8012E+02, 0.7214E+02/
      DATA(c(i, 18), i=1, 16)/ &
         0.9900E+02, 0.9563E+02, 0.9170E+02, 0.9001E+02, &
         0.8902E+02, 0.8832E+02, 0.8733E+02, 0.8610E+02, &
         0.8566E+02, 0.8497E+02, 0.8444E+02, 0.8399E+02, &
         0.8344E+02, 0.8272E+02, 0.8164E+02, 0.7408E+02/
      DATA(c(i, 19), i=1, 16)/ &
         0.9900E+02, 0.9753E+02, 0.9341E+02, 0.9163E+02, &
         0.9059E+02, 0.8985E+02, 0.8881E+02, 0.8753E+02, &
         0.8707E+02, 0.8637E+02, 0.8585E+02, 0.8540E+02, &
         0.8486E+02, 0.8418E+02, 0.8316E+02, 0.7602E+02/
      DATA(c(i, 20), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9511E+02, 0.9324E+02, &
         0.9215E+02, 0.9138E+02, 0.9029E+02, 0.8896E+02, &
         0.8848E+02, 0.8777E+02, 0.8726E+02, 0.8681E+02, &
         0.8628E+02, 0.8564E+02, 0.8468E+02, 0.7796E+02/
      DATA(c(i, 21), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9681E+02, 0.9486E+02, &
         0.9372E+02, 0.9291E+02, 0.9177E+02, 0.9039E+02, &
         0.8989E+02, 0.8917E+02, 0.8867E+02, 0.8822E+02, &
         0.8770E+02, 0.8710E+02, 0.8620E+02, 0.7990E+02/
      DATA(c(i, 22), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9851E+02, 0.9647E+02, &
         0.9528E+02, 0.9444E+02, 0.9325E+02, 0.9182E+02, &
         0.9130E+02, 0.9057E+02, 0.9008E+02, 0.8963E+02, &
         0.8912E+02, 0.8856E+02, 0.8772E+02, 0.8184E+02/
      DATA(c(i, 23), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9809E+02, &
         0.9685E+02, 0.9597E+02, 0.9473E+02, 0.9325E+02, &
         0.9271E+02, 0.9197E+02, 0.9149E+02, 0.9104E+02, &
         0.9054E+02, 0.9002E+02, 0.8924E+02, 0.8378E+02/
      DATA(c(i, 24), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, &
         0.9842E+02, 0.9750E+02, 0.9621E+02, 0.9468E+02, &
         0.9412E+02, 0.9337E+02, 0.9290E+02, 0.9245E+02, &
         0.9196E+02, 0.9148E+02, 0.9076E+02, 0.8572E+02/
      DATA(c(i, 25), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, &
         0.9900E+02, 0.9900E+02, 0.9769E+02, 0.9611E+02, &
         0.9553E+02, 0.9477E+02, 0.9431E+02, 0.9386E+02, &
         0.9338E+02, 0.9294E+02, 0.9228E+02, 0.8766E+02/
      DATA(c(i, 26), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9754E+02, &
         0.9694E+02, 0.9617E+02, 0.9572E+02, 0.9527E+02, &
         0.9480E+02, 0.9440E+02, 0.9380E+02, 0.8960E+02/
      DATA(c(i, 27), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9897E+02, &
         0.9835E+02, 0.9757E+02, 0.9713E+02, 0.9668E+02, &
         0.9622E+02, 0.9586E+02, 0.9532E+02, 0.9154E+02/
      DATA(c(i, 28), i=1, 16)/ &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, &
         0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, &
         0.9900E+02, 0.9897E+02, 0.9854E+02, 0.9809E+02, &
         0.9764E+02, 0.9732E+02, 0.9684E+02, 0.9348E+02/
!
!-----------------------------------------------------------------------
!          1.  REACTION PROBABILITY FOR CLONO2 + H2O
!              -------- ----------- --- ------ - ---
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

      hshcl(:) = 0.0
      rpcncl(:) = 0.0
      hhocl(:) = 0.0

      DO jl = kstart, kend

         ph2o = 1.0E4*th2o(jl)*boltzmann*t(jl)
         CALL ukca_position(pindx, noph2o, ph2o, jx, ierx, 0)
         CALL ukca_position(tindx, notemp, t(jl), jy, iery, 0)
         !
         !     temperature first criteria
         IF (jy == 0) THEN
            comp(jl) = 9.0
         ELSE IF (iery == 1) THEN
            comp(jl) = 95.0
         ELSE IF (jx == 0) THEN
            comp(jl) = 95.0
         ELSE IF (ierx == 1) THEN
            comp(jl) = 9.0
         ELSE
            jxp1 = jx + 1
            jyp1 = jy + 1
            sxy = c(jx, jy)
            sx1y = c(jxp1, jy)
            sx1y1 = c(jxp1, jyp1)
            sxy1 = c(jx, jyp1)
            ta = ph2o - pindx(jx)
            tb = pindx(jxp1) - pindx(jx)
            tt = ta/tb
            ua = t(jl) - tindx(jy)
            ub = tindx(jyp1) - tindx(jy)
            u = ua/ub
            comp(jl) = (1.0 - tt)*(1.0 - u)*sxy + tt*(1.0 - u)*sx1y + tt*u*sx1y1 + &
                       (1.0 - tt)*u*sxy1

            IF (comp(jl) < 9.0) comp(jl) = 9.0
            IF (comp(jl) > 95.0) comp(jl) = 95.0
            !
         END IF
         !
         !       WMO 1991
         rpcncl(jl) = 10.0**(1.87 - (0.074*comp(jl)))
         IF (rpcncl(jl) > 1.0) rpcncl(jl) = 1.0
      END DO
!
!-----------------------------------------------------------------------
!          2.  HENRY'S LAW COEFFICIENTS HSHCL AND HHOCL
!              ------- --- ------------ ----- --- -----
!
      IF (lphocl) THEN
         DO jl = kstart, kend
            !
            CALL ukca_position(compindx, nocomp, comp(jl), jx, ierx, 0)
            !
            IF (jx == 0) THEN
               a = -11.35
               b = 6.10E+03
            ELSE IF (ierx == 1) THEN
               a = -14.68
               b = 4.92E+03
            ELSE
               jxp1 = jx + 1
               ajx = aindex(jx)
               ajxp1 = aindex(jxp1)
               bjx = bindex(jx)
               bjxp1 = bindex(jxp1)
               ta = comp(jl) - compindx(jx)
               tb = compindx(jxp1) - compindx(jx)
               tt = ta/tb
               a = (1.0 - tt)*ajx + tt*ajxp1
               b = (1.0 - tt)*bjx + tt*bjxp1
            END IF
            !
            hshcl(jl) = EXP(a + b/t(jl))/101.3
            hhocl(jl) = EXP(0.71 + 1633.0/t(jl))/101.3
         END DO
      END IF
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_eqcomp

! ######################################################################
   SUBROUTINE ukca_position(xc, n, x, jx, ier, iorder)
!
!     Auxiliary subroutine for Eqcomp, Slimane Bekki, Nov. 1991
!
!-------------------------------------------------------------------------
!
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iorder
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xc(n)

      INTEGER, INTENT(OUT) :: ier
      INTEGER, INTENT(OUT) :: jx

! Local variables
      INTEGER :: i

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_POSITION'

!
!     --------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
      ier = 1
      IF (iorder == 0) THEN
         IF (x < xc(1)) THEN
            jx = 0
            ier = 0
         ELSE
            DO i = 1, n
               IF (x < xc(i)) THEN
                  ier = 0
                  EXIT
               END IF
            END DO
            jx = i - 1
         END IF
      ELSE IF (x > xc(1)) THEN
         jx = 0
         ier = 0
      ELSE
         DO i = 1, n
            IF (x > xc(i)) THEN
               ier = 0
               EXIT
            END IF
         END DO
         jx = i - 1
      END IF
!
!     --------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_position

! ######################################################################
   SUBROUTINE ukca_pscpres(t, p, tnd, th2o, thno3, &
                           kstart, kend, kchmlev, have_nat, sph2o)

      USE asad_mod, ONLY: shno3, sh2o, fpsc1, fpsc2

      IMPLICIT NONE

! Subroutine interface

      INTEGER, INTENT(IN) :: kstart
      INTEGER, INTENT(IN) :: kend
      INTEGER, INTENT(IN) :: kchmlev
      REAL, INTENT(IN) :: t(kchmlev)
      REAL, INTENT(IN) :: p(kchmlev)
      REAL, INTENT(IN) :: tnd(kchmlev)
      LOGICAL, INTENT(IN) :: have_nat(kchmlev)
      REAL, INTENT(IN OUT) :: th2o(kchmlev)
      REAL, INTENT(IN OUT) :: thno3(kchmlev)
      REAL, INTENT(IN OUT) :: sph2o(kchmlev)

! Local variables
      INTEGER :: jl
      REAL :: zh2ot
      REAL :: ztpsc
      REAL :: zmt
      REAL :: zbt
      REAL :: zhno3eq
      REAL :: zh2oeq

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_PSCPRES'

!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
      shno3(kstart:kend) = 0.0
      sh2o(kstart:kend) = 0.0
      fpsc1(kstart:kend) = 0.0
      fpsc2(kstart:kend) = 0.0
!
      DO jl = kstart, kend
         !
         !     Type 1 PSCs
         !
         zh2ot = 0.75*p(jl)*th2o(jl)/tnd(jl)
         zh2ot = MAX(zh2ot, 1.0E-12)
         ztpsc = t(jl)
         zmt = -2.7836 - 0.00088*ztpsc
         zbt = 38.9855 - 11397.0/ztpsc + 0.009179*ztpsc
         zhno3eq = 10.0**(zmt*LOG10(zh2ot) + zbt)
         !
         !     Convert to number density.
         !
         zhno3eq = zhno3eq*133.3*tnd(jl)/(100.0*p(jl))
         !
         ! only perform calculation if considering NAT in this region
         IF (thno3(jl) > zhno3eq .AND. have_nat(jl)) THEN
            fpsc1(jl) = 1.0
            shno3(jl) = thno3(jl) - zhno3eq
            thno3(jl) = zhno3eq
         ELSE
            fpsc1(jl) = 0.0
         END IF
         !
         !     Type 2 PSCs
         !
         IF (.TRUE.) THEN
            ! just sh2o from volume mixing ratio to number density and set
            ! FPSC2 flag
            IF (sph2o(jl) > 0.0) THEN
               sh2o(jl) = sph2o(jl)*tnd(jl)
               fpsc2(jl) = 1.0
            ELSE
               fpsc2(jl) = 0.0
               sh2o(jl) = 0.0
            END IF
         ELSE
            ! calculate water ice number density locally
            zh2oeq = 610.78*EXP(21.875*(t(jl) - 273.16)/(t(jl) - 7.66))
            !
            !     Convert to number density.
            !
            zh2oeq = zh2oeq*tnd(jl)/(100.0*p(jl))
            !
            IF (th2o(jl) > zh2oeq) THEN
               fpsc2(jl) = 1.0
               sh2o(jl) = th2o(jl) - zh2oeq
               th2o(jl) = zh2oeq
            ELSE
               fpsc2(jl) = 0.0
            END IF
         END IF
      END DO
!
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_pscpres
! =======================================================================
END MODULE ukca_hetero_mod
