module asad_bimol_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_BIMOL_MOD'
  public

  contains
  subroutine asad_bimol(n_points)
    use asad_mod, only : ab, f, jpbk, jpcspf, jpspb, nbrkx, p, peps, rk, spb, specf, t, t300, tnd, wp
    use asad_findreaction_mod, only : asad_findreaction
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ukca_um_legacy_mod, only : exp_v, oneover_v, powr_v
    use errormessagelength_mod, only : errormessagelength
    use ereport_mod, only : ereport
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, parameter :: i_chem_version_117 = 117
    integer, parameter :: i_chem_version_121 = 121
    integer, intent(in) :: n_points
    integer, save :: iohco = 0
    integer, save :: iohhno3 = 0
    integer, save :: iho2 = 0
    integer, save :: iohc3h8a = 0
    integer, save :: iohc3h8b = 0
    integer, save :: imeoomeooa = 0
    integer, save :: imeoomeoob = 0
    integer, save :: imeooro2a = 0
    integer, save :: imeooro2b = 0
    integer, save :: imeoo2a = 0
    integer, save :: imeoo2b = 0
    integer, save :: imeooho2a = 0
    integer, save :: imeooho2b = 0
    integer, save :: idmsoh = 0
    integer, save :: iso3h2o = 0
    integer, save :: ics2oh = 0
    integer, save :: ih2o = 0
    integer, save :: iho2no = 0
    integer :: jtr
    integer :: j
    integer :: jr
    real, allocatable, dimension(:) :: z1
    real, allocatable, dimension(:) :: z3
    real, allocatable, dimension(:) :: z4
    real, allocatable, dimension(:) :: z5
    real, allocatable, dimension(:) :: z6
    real, allocatable, dimension(:) :: alpha
    real, allocatable, dimension(:) :: ratioa2b
    real, allocatable, dimension(:) :: ratiob2total
    CHARACTER(LEN = 10) :: r1
    CHARACTER(LEN = 10) :: r2
    CHARACTER(LEN = 10) :: prods(jpspb)
    logical, save :: first = .true.
    logical, save :: first_pass = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_BIMOL'
    real, dimension(n_points) :: tmp
    real, dimension(n_points) :: tmp_out
    real, dimension(n_points) :: tmp1
    real, dimension(n_points) :: inv_t
    integer, save :: errcode = 0
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19
    integer :: idx_20
    integer :: idx_21
    integer :: idx_22
    integer :: idx_23
    integer :: idx_24
    integer :: idx_25
    integer :: idx_26
    integer :: idx_27
    integer :: idx_28
    integer :: idx_29
    integer :: idx_30
    integer :: idx_31
    integer :: idx_32
    integer :: idx_33
    integer :: idx_34
    integer :: idx_35
    integer :: idx_36
    integer :: idx_37
    integer :: idx_38
    integer :: idx_39
    integer :: idx_40

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    ALLOCATE(z1(1:n_points))
    ALLOCATE(z3(1:n_points))
    ALLOCATE(z4(1:n_points))
    ALLOCATE(z5(1:n_points))
    ALLOCATE(z6(1:n_points))
    ALLOCATE(ratioa2b(1:n_points))
    ALLOCATE(ratiob2total(1:n_points))
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(z1, dim=1), UBOUND(z1, dim=1), 1
      z1(idx) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_1 = LBOUND(z3, dim=1), UBOUND(z3, dim=1), 1
      z3(idx_1) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_2 = LBOUND(z4, dim=1), UBOUND(z4, dim=1), 1
      z4(idx_2) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_3 = LBOUND(z5, dim=1), UBOUND(z5, dim=1), 1
      z5(idx_3) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_4 = LBOUND(z6, dim=1), UBOUND(z6, dim=1), 1
      z6(idx_4) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_5 = LBOUND(ratioa2b, dim=1), UBOUND(ratioa2b, dim=1), 1
      ratioa2b(idx_5) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_6 = LBOUND(ratiob2total, dim=1), UBOUND(ratiob2total, dim=1), 1
      ratiob2total(idx_6) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_7 = 1, n_points, 1
      t300(idx_7) = t(idx_7) / 300.0
    enddo
    !$acc end kernels
    call oneover_v(n_points, t, inv_t)
    if (first_pass) then
      if (first) then
        first = .false.
        !$acc kernels
        !$acc loop seq
        do jtr = 1, jpcspf, 1
          if (specf(jtr) == 'H2O       ') then
            ih2o = jtr
          end if
        enddo
        !$acc end kernels
        r1(:) = '          '
        r2(:) = r1(:)
        prods(:) = r1(:)
        r1(:) = 'OH        '
        r2(:) = 'C3H8      '
        prods(1) = 'n-PrOO    '
        prods(2) = 'H2O       '
        iohc3h8a = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
        prods(1) = 'i-PrOO    '
        prods(2) = 'H2O       '
        iohc3h8b = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
        r1(:) = 'MeOO      '
        r2(:) = 'MeOO      '
        prods(1) = 'MeOH      '
        prods(2) = 'HCHO      '
        imeoomeooa = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
        prods(1) = 'HO2       '
        prods(2) = 'HCHO      '
        prods(3) = 'HO2       '
        prods(4) = 'HCHO      '
        imeoomeoob = asad_findreaction(r1,r2,prods,4,spb,nbrkx,jpbk + 1,jpspb)
        r1(:) = 'MeOO      '
        r2(:) = 'RO2       '
        prods(1) = 'MeOH      '
        prods(2) = 'HCHO      '
        imeooro2a = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
        prods(1) = 'HO2       '
        prods(2) = 'HCHO      '
        imeooro2b = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
        if ((imeoomeooa /= 0 .AND. imeoomeoob /= 0) .AND. (imeooro2a /= 0 .AND. imeooro2b /= 0)) then
          errcode = 2
          cmessage(:) = 'ASAD_BIMOL: Cannot have both MeOO+MeOO & MeOO+RO2 reactions'
          call ereport('ASAD_BIMOL', errcode, cmessage)
        end if
        imeoo2a = MAX(imeoomeooa, imeooro2a)
        imeoo2b = MAX(imeoomeoob, imeooro2b)
        r1(:) = 'HO2       '
        r2(:) = 'MeOO      '
        prods(1) = 'MeOOH     '
        imeooho2a = asad_findreaction(r1,r2,prods,1,spb,nbrkx,jpbk + 1,jpspb)
        prods(1) = 'HCHO      '
        imeooho2b = asad_findreaction(r1,r2,prods,1,spb,nbrkx,jpbk + 1,jpspb)
        r1(:) = 'HO2       '
        r2(:) = 'NO        '
        prods(1) = 'HONO2     '
        iho2no = asad_findreaction(r1,r2,prods,1,spb,nbrkx,jpbk + 1,jpspb)
        if (ukca_config%l_ukca_chem_aero) then
          r1(:) = 'DMS       '
          r2(:) = 'OH        '
          prods(1) = 'SO2       '
          if (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) then
            prods(2) = 'DMSO      '
            idmsoh = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
          else
            if (ukca_config%l_ukca_strat) then
              prods(2) = 'MSA       '
              idmsoh = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
            else
              if (ukca_config%l_ukca_strattrop .AND. ukca_config%i_ukca_chem_version < i_chem_version_121) then
                if (ukca_config%i_ukca_chem_version >= i_chem_version_117) then
                  prods(2) = 'DMSO      '
                  prods(3) = 'MeOO      '
                  idmsoh = asad_findreaction(r1,r2,prods,3,spb,nbrkx,jpbk + 1,jpspb)
                else
                  prods(2) = 'MSA       '
                  idmsoh = asad_findreaction(r1,r2,prods,2,spb,nbrkx,jpbk + 1,jpspb)
                end if
              else
                if (ukca_config%l_ukca_tropisop) then
                  prods(2) = 'DMSO      '
                  prods(3) = 'MeOO      '
                  idmsoh = asad_findreaction(r1,r2,prods,3,spb,nbrkx,jpbk + 1,jpspb)
                else
                  idmsoh = 0
                end if
              end if
            end if
          end if
        else
          idmsoh = 0
        end if
        if (ukca_config%l_ukca_chem_aero .AND. (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be .OR. ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_strattrop .AND. ukca_config%i_ukca_chem_version < i_chem_version_121 .OR. ukca_config%l_ukca_tropisop) .AND. idmsoh == 0) then
          errcode = 1
          cmessage(:) = ' Rate coefficient for DMS + OH lacks correction'
          call ereport('ASAD_BIMOL', errcode, cmessage)
        end if
        !$acc kernels
        !$acc loop seq
        do j = 1, jpbk, 1
          jr = nbrkx(j)
          if (spb(j,1) == 'OH     ' .AND. spb(j,2) == 'CO     ' .OR. spb(j,1) == 'CO     ' .AND. spb(j,2) == 'OH     ') then
            iohco = jr
          end if
          if (spb(j,1) == 'OH     ' .AND. spb(j,2) == 'HONO2  ' .OR. spb(j,1) == 'HONO2  ' .AND. spb(j,2) == 'OH     ') then
            iohhno3 = jr
          end if
          if (spb(j,1) == 'HO2    ' .AND. spb(j,2) == 'HO2    ') then
            iho2 = jr
          end if
          if (spb(j,1) == 'SO3      ' .AND. spb(j,2) == 'H2O   ') then
            iso3h2o = jr
          end if
          if (spb(j,1) == 'CS2      ' .AND. spb(j,2) == 'OH    ') then
            ics2oh = jr
          end if
        enddo
        !$acc end kernels
      end if
      first_pass = .false.
    end if
    do j = 1, jpbk, 1
      jr = nbrkx(j)
      if (ABS(ab(j,2)) < peps .AND. ABS(ab(j,3)) < peps) then
        rk(1:n_points,jr) = ab(j,1)
      else
        if (ABS(ab(j,2)) < peps) then
          tmp(:) = -ab(j,3) * inv_t(:)
          call exp_v(n_points, tmp, tmp_out)
          rk(1:n_points,jr) = ab(j,1) * tmp_out(:)
        else
          if (ABS(ab(j,3)) < peps) then
            call powr_v(n_points, t300, ab(j,2), tmp1)
            rk(1:n_points,jr) = ab(j,1) * tmp1(:)
          else
            tmp(:) = -ab(j,3) * inv_t(:)
            call exp_v(n_points, tmp, tmp_out)
            call powr_v(n_points, t300, ab(j,2), tmp1)
            rk(1:n_points,jr) = ab(j,1) * tmp1(:) * tmp_out(:)
          end if
        end if
      end if
    enddo
    if (iohco /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_8 = 1, n_points, 1
        rk(idx_8,iohco) = rk(idx_8,iohco) * (1.0 + tnd(idx_8) / 4.2e19)
      enddo
      !$acc end kernels
    end if
    if (iohhno3 /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_9 = 1, n_points, 1
        tmp(idx_9) = 460.0 * inv_t(idx_9)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_10 = LBOUND(z1, dim=1), UBOUND(z1, dim=1), 1
        z1(idx_10) = 2.4e-14 * tmp_out(idx_10 + (1 - LBOUND(z1, dim=1)))
      enddo
      !$acc loop gang vector independent
      do idx_11 = 1, n_points, 1
        tmp(idx_11) = 1335.0 * inv_t(idx_11)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_12 = LBOUND(z3, dim=1), UBOUND(z3, dim=1), 1
        z3(idx_12) = (6.5e-34 * tmp_out(idx_12 + (1 - LBOUND(z3, dim=1)))) * tnd(idx_12 + (1 - LBOUND(z3, dim=1)))
      enddo
      !$acc loop gang vector independent
      do idx_13 = 1, n_points, 1
        tmp(idx_13) = 2199.0 * inv_t(idx_13)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_14 = LBOUND(z4, dim=1), UBOUND(z4, dim=1), 1
        z4(idx_14) = 2.7e-17 * tmp_out(idx_14 + (1 - LBOUND(z4, dim=1)))
      enddo
      !$acc loop gang vector independent
      do idx_15 = 1, n_points, 1
        rk(idx_15,iohhno3) = z1(idx_15 + (LBOUND(z1, dim=1) - 1)) + z3(idx_15 + (LBOUND(z3, dim=1) - 1)) / (1.0 + z3(idx_15 + (LBOUND(z3, dim=1) - 1)) / z4(idx_15 + (LBOUND(z4, dim=1) - 1)))
      enddo
      !$acc end kernels
    end if
    if (ih2o /= 0 .AND. iho2 /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_16 = 1, n_points, 1
        tmp(idx_16) = 2200.0 * inv_t(idx_16)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_17 = 1, n_points, 1
        rk(idx_17,iho2) = rk(idx_17,iho2) * (1.0 + 1.4e-21 * f(idx_17,ih2o) * tmp_out(idx_17))
      enddo
      !$acc end kernels
    else
      if (ih2o == 0 .AND. iho2 /= 0) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_18 = 1, n_points, 1
          tmp(idx_18) = 2200.0 * inv_t(idx_18)
        enddo
        !$acc end kernels
        call exp_v(n_points, tmp, tmp_out)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_19 = 1, n_points, 1
          rk(idx_19,iho2) = rk(idx_19,iho2) * (1.0 + 1.4e-21 * wp(idx_19) * tnd(idx_19) * tmp_out(idx_19))
        enddo
        !$acc end kernels
      end if
    end if
    if (iho2no /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_20 = 1, n_points, 1
        rk(idx_20,iho2no) = rk(idx_20,iho2no) * (530.0 * inv_t(idx_20) + 8.53e-4 * (1e-2 * p(idx_20)) - 1.73) / 100.0
      enddo
      !$acc end kernels
    end if
    if (ih2o /= 0 .AND. iso3h2o /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_21 = 1, n_points, 1
        rk(idx_21,iso3h2o) = rk(idx_21,iso3h2o) * f(idx_21,ih2o)
      enddo
      !$acc end kernels
    else
      if (ih2o == 0 .AND. iso3h2o /= 0) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_22 = 1, n_points, 1
          rk(idx_22,iso3h2o) = rk(idx_22,iso3h2o) * wp(idx_22) * tnd(idx_22)
        enddo
        !$acc end kernels
      end if
    end if
    if (ics2oh /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_23 = 1, n_points, 1
        tmp(idx_23) = 3400.0 * inv_t(idx_23)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_24 = LBOUND(z6, dim=1), UBOUND(z6, dim=1), 1
        z6(idx_24) = 1.81e-3 * tmp_out(idx_24 + (1 - LBOUND(z6, dim=1)))
      enddo
      !$acc loop gang vector independent
      do idx_25 = 1, n_points, 1
        rk(idx_25,ics2oh) = rk(idx_25,ics2oh) / (t(idx_25) + z6(idx_25 + (LBOUND(z6, dim=1) - 1)))
      enddo
      !$acc end kernels
    end if
    if (idmsoh /= 0) then
      ALLOCATE(alpha(1:n_points))
      !$acc kernels
      !$acc loop gang vector independent
      do idx_26 = 1, n_points, 1
        tmp(idx_26) = 7460.0 * inv_t(idx_26)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_27 = LBOUND(alpha, dim=1), UBOUND(alpha, dim=1), 1
        alpha(idx_27) = 1.106e-31 * tmp_out(idx_27 + (1 - LBOUND(alpha, dim=1))) * tnd(idx_27 + (1 - LBOUND(alpha, dim=1)))
      enddo
      !$acc loop gang vector independent
      do idx_28 = 1, n_points, 1
        rk(idx_28,idmsoh) = rk(idx_28,idmsoh) * alpha(idx_28 + (LBOUND(alpha, dim=1) - 1)) / (1.0 + alpha(idx_28 + (LBOUND(alpha, dim=1) - 1)))
      enddo
      !$acc end kernels
    end if
    if (iohc3h8a /= 0 .AND. iohc3h8b /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_29 = 1, n_points, 1
        tmp(idx_29) = -816.0 * inv_t(idx_29)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      call powr_v(n_points, t, -0.64, tmp1)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_30 = 1, n_points, 1
        ratioa2b(idx_30) = 226.0 * tmp1(idx_30) * tmp_out(idx_30)
      enddo
      !$acc loop gang vector independent
      do idx_31 = 1, n_points, 1
        rk(idx_31,iohc3h8a) = rk(idx_31,iohc3h8a) * (ratioa2b(idx_31) / (ratioa2b(idx_31) + 1.0))
      enddo
      !$acc loop gang vector independent
      do idx_32 = 1, n_points, 1
        rk(idx_32,iohc3h8b) = rk(idx_32,iohc3h8b) / (ratioa2b(idx_32) + 1.0)
      enddo
      !$acc end kernels
    end if
    if ((imeoo2a /= 0 .AND. imeoo2b /= 0) .AND. (.NOT.ukca_config%l_ukca_cristrat)) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_33 = 1, n_points, 1
        tmp(idx_33) = 1300.0 * inv_t(idx_33)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_34 = 1, n_points, 1
        ratiob2total(idx_34) = 1.0 / (1.0 + tmp_out(idx_34) / 33.0)
      enddo
      !$acc loop gang vector independent
      do idx_35 = 1, n_points, 1
        rk(idx_35,imeoo2a) = rk(idx_35,imeoo2a) * (1.0 - ratiob2total(idx_35))
      enddo
      !$acc loop gang vector independent
      do idx_36 = 1, n_points, 1
        rk(idx_36,imeoo2b) = rk(idx_36,imeoo2b) * ratiob2total(idx_36)
      enddo
      !$acc end kernels
    end if
    if (imeooho2a /= 0 .AND. imeooho2b /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_37 = 1, n_points, 1
        tmp(idx_37) = -1160.0 * inv_t(idx_37)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do idx_38 = 1, n_points, 1
        ratiob2total(idx_38) = 1.0 / (1.0 + 498.0 * tmp_out(idx_38))
      enddo
      !$acc loop gang vector independent
      do idx_39 = 1, n_points, 1
        rk(idx_39,imeooho2a) = rk(idx_39,imeooho2a) * (1.0 - ratiob2total(idx_39))
      enddo
      !$acc loop gang vector independent
      do idx_40 = 1, n_points, 1
        rk(idx_40,imeooho2b) = rk(idx_40,imeooho2b) * ratiob2total(idx_40)
      enddo
      !$acc end kernels
    end if
    if (ALLOCATED(z1)) then
      DEALLOCATE(z1)
    end if
    if (ALLOCATED(z3)) then
      DEALLOCATE(z3)
    end if
    if (ALLOCATED(z4)) then
      DEALLOCATE(z4)
    end if
    if (ALLOCATED(z5)) then
      DEALLOCATE(z5)
    end if
    if (ALLOCATED(z6)) then
      DEALLOCATE(z6)
    end if
    if (ALLOCATED(ratioa2b)) then
      DEALLOCATE(ratioa2b)
    end if
    if (ALLOCATED(ratiob2total)) then
      DEALLOCATE(ratiob2total)
    end if
    if (ALLOCATED(alpha)) then
      DEALLOCATE(alpha)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_bimol

end module asad_bimol_mod
