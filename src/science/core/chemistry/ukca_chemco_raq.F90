! PSyclone CodeBlock (unsupported code) reason:
!  - PSyclone SymbolTable error: Unable to find the scope of node 'Reference[name:'stratflag']' as none of its ancestors are Container or Schedule nodes.
MODULE ukca_chemco_raq_mod
  IMPLICIT NONE
  INTEGER, SAVE :: in2o5 = 0
  INTEGER, SAVE :: in2o5_ss = 0
  CHARACTER(LEN = *), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CHEMCO_RAQ_MOD'
  CONTAINS
  SUBROUTINE ukca_chemco_raq(nr, n_pnts, tc, m, h2o, o2, clw, fcloud, fdiss, rho_air, rh_frac, so4_aitken, so4_accum, soot_fresh, soot_aged, ocff_fresh, ocff_aged, biogenic, sea_salt_film, sea_salt_jet, stratflag, rc, hrc, H_plus_1d_arr)
    USE ukca_chemco_raq_init_mod, ONLY: ukca_chemco_raq_init, rk, pdep
    USE asad_mod, ONLY: peps, jpeq, jpspec, jphk
    USE ukca_constants, ONLY: m_air, m_n2o5
    USE ereport_mod, ONLY: ereport
    USE errormessagelength_mod, ONLY: errormessagelength
    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE ukca_config_constants_mod, ONLY: rho_water, avogadro
    USE ukca_config_specification_mod, ONLY: ukca_config
    USE ukca_mode_setup, ONLY: cp_su, cp_bc, cp_oc, cp_cl, cp_so
    USE ukca_um_legacy_mod, ONLY: calc_surf_area, l_um_calc_surf_area
    USE ukca_trop_hetchem_mod, ONLY: cal_n2o5, cal_loss1k
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nr
    INTEGER, INTENT(IN) :: n_pnts
    REAL, INTENT(IN) :: tc(n_pnts)
    REAL, INTENT(IN) :: m(n_pnts)
    REAL, INTENT(IN) :: h2o(n_pnts)
    REAL, INTENT(IN) :: o2(n_pnts)
    REAL, INTENT(IN) :: clw(n_pnts)
    REAL, INTENT(IN) :: fcloud(n_pnts)
    REAL, INTENT(IN) :: fdiss(n_pnts, jpspec, jpeq + 1)
    REAL, INTENT(IN) :: rho_air(n_pnts)
    REAL, INTENT(IN) :: rh_frac(n_pnts)
    REAL, INTENT(IN) :: so4_aitken(n_pnts)
    REAL, INTENT(IN) :: so4_accum(n_pnts)
    REAL, INTENT(IN) :: soot_fresh(n_pnts)
    REAL, INTENT(IN) :: soot_aged(n_pnts)
    REAL, INTENT(IN) :: ocff_fresh(n_pnts)
    REAL, INTENT(IN) :: ocff_aged(n_pnts)
    REAL, INTENT(IN) :: biogenic(n_pnts)
    REAL, INTENT(IN) :: sea_salt_film(n_pnts)
    REAL, INTENT(IN) :: sea_salt_jet(n_pnts)
    LOGICAL, INTENT(IN) :: stratflag(n_pnts)
    REAL, INTENT(IN) :: H_plus_1d_arr(n_pnts)
    REAL, INTENT(OUT) :: rc(n_pnts, nr)
    REAL, INTENT(OUT) :: hrc(n_pnts, jphk)
    INTEGER, PARAMETER :: npdep = 10
    REAL, PARAMETER :: o2frac = 0.20948
    REAL, PARAMETER :: n2frac = 0.78084
    REAL, PARAMETER :: qcl_min = 1.0E-12
    INTEGER, PARAMETER :: ih_o3 = 4
    INTEGER, PARAMETER :: ih_hno3 = 10
    INTEGER, PARAMETER :: ih_h2o2 = 11
    INTEGER, PARAMETER :: ih_so2 = 45
    INTEGER :: i, j, k
    INTEGER :: r
    INTEGER :: nklo, nkhi, nres
    REAL :: tfc = 1.0 / 300.0
    REAL :: ar
    REAL :: n
    REAL :: ea
    REAL :: fac1
    REAL :: fac2
    REAL :: fac3
    REAL :: brn
    REAL :: f
    REAL :: fc
    REAL :: fac1v(n_pnts)
    REAL :: faq(n_pnts, jpspec)
    REAL :: vr
    REAL :: sa_so4_ait(n_pnts)
    REAL :: sa_so4_acc(n_pnts)
    REAL :: sa_bc_fresh(n_pnts)
    REAL :: sa_bc_aged(n_pnts)
    REAL :: sa_ocff_fresh(n_pnts)
    REAL :: sa_ocff_aged(n_pnts)
    REAL :: sa_soa(n_pnts)
    REAL :: sa_ss_film(n_pnts)
    REAL :: sa_ss_jet(n_pnts)
    REAL :: wr_so4_ait(n_pnts)
    REAL :: wr_so4_acc(n_pnts)
    REAL :: wr_bc_fresh(n_pnts)
    REAL :: wr_bc_aged(n_pnts)
    REAL :: wr_ocff_fresh(n_pnts)
    REAL :: wr_ocff_aged(n_pnts)
    REAL :: wr_soa(n_pnts)
    REAL :: wr_ss_film(n_pnts)
    REAL :: wr_ss_jet(n_pnts)
    REAL :: gamma_n2o5_so4(n_pnts)
    REAL :: gamma_n2o5_bc(n_pnts)
    REAL :: gamma_n2o5_ocff(n_pnts)
    REAL :: gamma_n2o5_soa(n_pnts)
    REAL :: gamma_n2o5_ss(n_pnts)
    REAL :: arg2
    REAL :: zfc
    REAL :: zr
    REAL :: zi(n_pnts)
    REAL :: zo(n_pnts)
    REAL :: alpha(n_pnts)
    REAL :: t300(n_pnts)
    REAL :: at(7)
    LOGICAL, SAVE :: first = .TRUE.
    LOGICAL, SAVE :: first_pass = .TRUE.
    INTEGER :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    INTEGER(KIND = jpim), PARAMETER :: zhook_in = 0
    INTEGER(KIND = jpim), PARAMETER :: zhook_out = 1
    REAL(KIND = jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMCO_RAQ'
    IF (lhook) CALL dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    IF (first_pass) THEN
      IF (first) THEN
        CALL ukca_chemco_raq_init(nr, npdep)
        IF (ukca_config % l_ukca_classic_hetchem) THEN
          in2o5 = 1
          in2o5_ss = 2
        END IF
        first = .FALSE.
      END IF
      first_pass = .FALSE.
    END IF
    rc(:, :) = 0.0
    hrc(:, :) = 0.0
    DO i = 1, nr
      r = rk(i) % r
      ar = rk(i) % a
      n = rk(i) % n
      ea = rk(i) % e
      IF (r == 1) THEN
        rc(:, i) = ar
      ELSE IF (r == 2) THEN
        rc(:, i) = ar * ((tc(:) * tfc) ** n)
      ELSE IF (r == 3) THEN
        rc(:, i) = ar * EXP(ea / tc(:))
      ELSE IF (r == 4) THEN
        rc(:, i) = ar * ((tc(:) * tfc) ** n) * EXP(ea / tc(:))
      END IF
    END DO
    DO k = 1, npdep
      nklo = pdep(k) % klo
      nkhi = pdep(k) % khi
      nres = pdep(k) % res
      f = pdep(k) % fc
      IF (f > 0.01) THEN
        fc = f
        DO j = 1, n_pnts
          rc(j, nklo) = rc(j, nklo) * m(j)
          brn = 0.75 - 1.27 * LOG10(fc)
          fac1 = rc(j, nklo) / rc(j, nkhi)
          fac2 = rc(j, nklo) / (1.0 + fac1)
          fac3 = 1.0 + ((LOG10(fac1) / brn) ** 2)
          rc(j, nres) = fac2 * (fc ** (1.0 / fac3))
        END DO
      ELSE
        DO j = 1, n_pnts
          fc = EXP(- tc(j) / 250.0) + EXP(- 1050.0 / tc(j))
          rc(j, nklo) = rc(j, nklo) * m(j)
          brn = 0.75 - 1.27 * LOG10(fc)
          fac1 = rc(j, nklo) / rc(j, nkhi)
          fac2 = rc(j, nklo) / (1.0 + fac1)
          fac3 = 1.0 + ((LOG10(fac1) / brn) ** 2)
          rc(j, nres) = fac2 * (fc ** (1.0 / fac3))
        END DO
      END IF
    END DO
    rc(:, 7) = (o2frac * rc(:, 6) + n2frac * rc(:, 7)) * m(:)
    fac1v(:) = rc(:, 51) * m(:)
    rc(:, 35) = rc(:, 35) + fac1v(:) / (1.0 + fac1v(:) / rc(:, 50))
    rc(:, 36) = (rc(:, 36) + rc(:, 37) * m(:)) * (1.0 + rc(:, 38) * h2o(:))
    rc(:, 62) = rc(:, 62) - rc(:, 61)
    rc(:, 70) = rc(:, 70) + (rc(:, 69) * m(:))
    fac1v(:) = rc(:, 74) / (1.0 + rc(:, 74))
    rc(:, 74) = rc(:, 80) * (1.0 - fac1v(:))
    rc(:, 80) = rc(:, 80) * fac1v(:)
    rc(:, 94) = rc(:, 89) + rc(:, 94)
    rc(:, 1) = rc(:, 1) * m(:) * o2(:)
    IF (ukca_config % l_ukca_classic_hetchem) THEN
      IF (.NOT. (l_um_calc_surf_area)) THEN
        errcode = 1
        cmessage = 'UM support unavailable for heterogenous chemistry on CLASSIC particles'
        CALL ereport(RoutineName, errcode, cmessage)
      END IF
      CALL calc_surf_area(n_pnts, rho_air(1 : n_pnts), rh_frac(1 : n_pnts), so4_aitken(1 : n_pnts), so4_accum(1 : n_pnts), soot_fresh(1 : n_pnts), soot_aged(1 : n_pnts), ocff_fresh(1 : n_pnts), ocff_aged(1 : n_pnts), biogenic(1 : n_pnts), sea_salt_film(1 : n_pnts), sea_salt_jet(1 : n_pnts), sa_so4_ait(1 : n_pnts), sa_so4_acc(1 : n_pnts), sa_bc_fresh(1 : n_pnts), sa_bc_aged(1 : n_pnts), sa_ocff_fresh(1 : n_pnts), sa_ocff_aged(1 : n_pnts), sa_soa(1 : n_pnts), sa_ss_film(1 : n_pnts), sa_ss_jet(1 : n_pnts), wr_so4_ait(1 : n_pnts), wr_so4_acc(1 : n_pnts), wr_bc_fresh(1 : n_pnts), wr_bc_aged(1 : n_pnts), wr_ocff_fresh(1 : n_pnts), wr_ocff_aged(1 : n_pnts), wr_soa(1 : n_pnts), wr_ss_film(1 : n_pnts), wr_ss_jet(1 : n_pnts))
      IF (in2o5 > 0) THEN
        gamma_n2o5_so4(:) = cal_n2o5(cp_su, n_pnts, tc(:), rh_frac(:))
        gamma_n2o5_bc(:) = cal_n2o5(cp_bc, n_pnts, tc(:), rh_frac(:))
        gamma_n2o5_ocff(:) = cal_n2o5(cp_oc, n_pnts, tc(:), rh_frac(:))
        gamma_n2o5_soa(:) = cal_n2o5(cp_so, n_pnts, tc(:), rh_frac(:))
        hrc(:, in2o5) = cal_loss1k(n_pnts, sa_so4_ait(:), wr_so4_ait(:), m(:), gamma_n2o5_so4(:), tc(:), m_n2o5)
        hrc(:, in2o5) = hrc(:, in2o5) + cal_loss1k(n_pnts, sa_so4_acc(:), wr_so4_acc(:), m(:), gamma_n2o5_so4(:), tc(:), m_n2o5)
        hrc(:, in2o5) = hrc(:, in2o5) + cal_loss1k(n_pnts, sa_bc_fresh(:), wr_bc_fresh(:), m(:), gamma_n2o5_bc(:), tc(:), m_n2o5)
        hrc(:, in2o5) = hrc(:, in2o5) + cal_loss1k(n_pnts, sa_bc_aged(:), wr_bc_aged(:), m(:), gamma_n2o5_bc(:), tc(:), m_n2o5)
        hrc(:, in2o5) = hrc(:, in2o5) + cal_loss1k(n_pnts, sa_ocff_fresh(:), wr_ocff_fresh(:), m(:), gamma_n2o5_ocff(:), tc(:), m_n2o5)
        hrc(:, in2o5) = hrc(:, in2o5) + cal_loss1k(n_pnts, sa_ocff_aged(:), wr_ocff_aged(:), m(:), gamma_n2o5_ocff(:), tc(:), m_n2o5)
        hrc(:, in2o5) = hrc(:, in2o5) + cal_loss1k(n_pnts, sa_soa(:), wr_soa(:), m(:), gamma_n2o5_soa(:), tc(:), m_n2o5)
        WHERE (stratflag) hrc(:, in2o5) = 0.0
      END IF
      IF (in2o5_ss > 0) THEN
        gamma_n2o5_ss(:) = cal_n2o5(cp_cl, n_pnts, tc(:), rh_frac(:))
        hrc(:, in2o5_ss) = cal_loss1k(n_pnts, sa_ss_film(:), wr_ss_film(:), m(:), gamma_n2o5_ss(:), tc(:), m_n2o5)
        hrc(:, in2o5_ss) = hrc(:, in2o5_ss) + cal_loss1k(n_pnts, sa_ss_jet(:), wr_ss_jet(:), m(:), gamma_n2o5_ss(:), tc(:), m_n2o5)
        WHERE (stratflag) hrc(:, in2o5_ss) = 0.0
      END IF
    END IF
    IF (ukca_config % l_ukca_raqaero) THEN
      faq(:, :) = fdiss(:, :, 1) + fdiss(:, :, 2)
      rc(:, 184) = 9.60E-12 * EXP(- 240.0 / tc(:))
      alpha(:) = 1.106E-31 * EXP(7460.0 / tc(:)) * m(:)
      rc(:, 185) = 3.04E-12 * EXP(350.00 / tc(:)) * alpha(:) / (1 + alpha(:))
      rc(:, 186) = 1.90E-13 * EXP(500.00 / tc(:))
      rc(:, 187) = 5.80E-11
      at(1) = 0.60
      at(2) = 3.00E-31
      at(3) = - 3.30
      at(4) = 0.0
      at(5) = 1.50E-12
      at(6) = 0.00
      at(7) = 0.0
      t300(:) = tc(:) / 300.0
      zo(:) = at(2) * EXP(at(3) * LOG(t300(:))) * EXP(- at(4) / tc(:)) * m(:)
      zi(:) = at(5) * EXP(at(6) * LOG(t300(:))) * EXP(- at(7) / tc(:))
      DO i = 1, n_pnts
        IF (zo(i) < peps) THEN
          rc(i, 188) = zi(i)
        ELSE IF (zi(i) < peps) THEN
          rc(i, 188) = zo(i)
        ELSE
          IF (at(1) <= 1.0) THEN
            zfc = at(1)
          ELSE
            zfc = EXP(- tc(i) / at(1))
          END IF
          zr = zo(i) / zi(i)
          arg2 = 1.0 / (1.0 + ((LOG10(zr)) ** 2))
          rc(i, 188) = (zo(i) / (1.0 + zr)) * EXP(arg2 * LOG(zfc))
          rc(i, 188) = rc(i, 188) * (1.0 - (faq(i, ih_so2) * fcloud(i)))
        END IF
      END DO
      DO i = 1, n_pnts
        IF (clw(i) > qcl_min) THEN
          vr = clw(i) * m(i) * m_air * (1E6 / avogadro) / rho_water
          rc(i, 31) = rc(i, 31) * (1.0 - (fdiss(i, ih_h2o2, 1) + fdiss(i, ih_h2o2, 2)) * fcloud(i))
          rc(i, 35) = rc(i, 35) * (1.0 - (fdiss(i, ih_hno3, 1) + fdiss(i, ih_hno3, 2)) * fcloud(i))
          rc(i, 190) = 7.43E+16 * EXP(- 5280.0 / tc(i)) * fcloud(i) * fdiss(i, ih_so2, 3) * fdiss(i, ih_o3, 1) * 1000.0 / (avogadro * vr)
          rc(i, 191) = 2.1295E+14 * EXP(- 4430.0 / tc(i)) * (H_plus_1d_arr(i) / (1.0 + 13.0 * H_plus_1d_arr(i))) * fcloud(i) * fdiss(i, ih_so2, 2) * fdiss(i, ih_h2o2, 1) * 1000.0 / (avogadro * vr)
          rc(i, 197) = 4.0113E+13 * EXP(- 5530.0 / tc(i)) * fcloud(i) * fdiss(i, ih_so2, 2) * fdiss(i, ih_o3, 1) * 1000.0 / (avogadro * vr)
        ELSE
          rc(i, 190) = 0.0
          rc(i, 191) = 0.0
          rc(i, 197) = 0.0
        END IF
      END DO
      rc(:, 192) = 1.2E-11 * EXP(444.0 / tc(:))
      rc(:, 193) = 1.01E-15 * EXP(- 732.0 / tc(:))
      rc(:, 194) = 1.19E-12 * EXP(490.0 / tc(:))
      rc(:, 195) = 0.0
      rc(:, 196) = 0.0
    END IF
    IF (lhook) CALL dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    RETURN
  END SUBROUTINE ukca_chemco_raq
END MODULE ukca_chemco_raq_mod
