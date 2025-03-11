module ukca_chemistry_ctl_be_mod
  implicit none
  real, allocatable, dimension(:,:), private :: rc
  integer, allocatable, dimension(:), save, private :: ibimol
  integer, allocatable, dimension(:), save, private :: itrimol
  integer, allocatable, dimension(:), save, private :: ihetero
  real, save, private :: frac_dms_so2
  real, save, private :: frac_dms_dmso
  real, save, private :: frac_dmso_so2
  real, dimension(3), save, private :: frac_monoterp_sec_org
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_CTL_BE_MOD'
  private

  public :: ukca_chemistry_ctl_be

  contains
  subroutine ukca_chemistry_ctl_be(row_length, rows, model_levels, theta_field_size, tot_n_pnts, ntracers, chem_timestep, k_be_top, p_theta_levels, temp, q, qcl, tracer, cloud_frac, volume, uph2so4inaer, delso2_wet_h2o2, delso2_wet_o3, delh2so4_chem, delso2_drydep, delso2_wetdep, h_plus, h2o2_offline, dryrt, wetrt, nlev_with_ddep, firstcall)
    use ukca_um_legacy_mod, only : deposition_from_ukca_chemistry
    use asad_mod, only : advt, f, jpbk, jpctr, jpdd, jpdw, jphk, jpnr, jpspec, jptk, lvmr, ndepd, ndepw, nitfg, nldepd, prk, speci, tnd, y
    use asad_chem_flux_diags, only : asad_chemical_diagnostics, l_asad_use_chem_diags, l_asad_use_drydep, l_asad_use_flux_rxns, l_asad_use_wetdep
    use asad_ftoy_mod, only : asad_ftoy
    use ukca_cspecies, only : c_species, n_h2o2, nn_h2so4, nn_so2
    use ukca_constants, only : c_h2o
    use ukca_environment_fields_mod, only : surf_wetness
    use ukca_config_specification_mod, only : ukca_config
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use ukca_be_drydep_mod, only : ukca_be_drydep
    use ukca_be_wetdep_mod, only : ukca_be_wetdep
    use ukca_ddepctl_mod, only : ukca_ddepctl
    use ukca_ddeprt_mod, only : ukca_ddeprt
    use ukca_wdeprt_mod, only : ukca_wdeprt
    use ukca_drydep_mod, only : ukca_drydep
    use ukca_wetdep_mod, only : ukca_wetdep
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: theta_field_size
    integer, intent(in) :: tot_n_pnts
    integer, intent(in) :: ntracers
    integer, intent(in) :: uph2so4inaer
    integer, intent(in) :: k_be_top
    integer, dimension(row_length,rows), intent(in) :: nlev_with_ddep
    integer, intent(in) :: chem_timestep
    real, dimension(tot_n_pnts), intent(in) :: p_theta_levels
    real, dimension(tot_n_pnts), intent(in) :: temp
    real, dimension(row_length,rows,model_levels), intent(in) :: volume
    real, dimension(tot_n_pnts), intent(in) :: qcl
    real, dimension(tot_n_pnts), intent(in) :: cloud_frac
    real, dimension(theta_field_size,jpdd), intent(in) :: dryrt
    real, dimension(theta_field_size,model_levels,jpdw), intent(in) :: wetrt
    real, dimension(tot_n_pnts), intent(inout) :: q
    real, dimension(tot_n_pnts,ntracers), intent(inout) :: tracer
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wet_h2o2
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wet_o3
    real, dimension(tot_n_pnts), intent(inout) :: delh2so4_chem
    real, dimension(tot_n_pnts), intent(inout) :: delso2_drydep
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wetdep
    real, dimension(tot_n_pnts), intent(in) :: h_plus
    real, dimension(tot_n_pnts), intent(in) :: h2o2_offline
    logical, intent(in) :: firstcall
    logical :: lflux
    real, allocatable, dimension(:,:) :: dflux
    integer, dimension(theta_field_size) :: nlev_with_ddep2
    integer, save :: nr
    integer, save :: n_be_calls
    integer :: ix
    integer :: jy
    integer :: i
    integer :: j
    integer :: js
    integer :: k
    integer :: l
    integer :: jit
    integer :: kcs
    integer :: kce
    integer :: ierr
    integer :: errcode
    CHARACTER(LEN = 72) :: cmessage
    logical, dimension(theta_field_size) :: blmask
    logical, save :: ofirst = .true.
    real, save :: dts
    real :: secs_per_step
    real, dimension(model_levels,jpspec) :: dpd_dummy
    real, dimension(model_levels,jpspec) :: dpw_dummy
    real, dimension(model_levels,jpnr) :: prk_dummy
    real, dimension(model_levels,jpspec) :: y_dummy
    real, dimension(theta_field_size) :: so2_wetox_h2o2
    real, dimension(theta_field_size) :: so2_wetox_o3
    real, dimension(theta_field_size) :: so2_dryox_oh
    real, allocatable, dimension(:) :: ystore
    real, dimension(theta_field_size,jpctr) :: zftr
    real, dimension(theta_field_size) :: zp
    real, dimension(theta_field_size) :: zt
    real, dimension(theta_field_size) :: zq
    real, dimension(theta_field_size) :: zclw
    real, dimension(theta_field_size) :: zfcloud
    real, dimension(theta_field_size,jpdd) :: zdryrt2
    real, dimension(theta_field_size,model_levels,jpdw) :: zwetrt2
    real, dimension(theta_field_size,jpspec) :: zwetrt
    real, dimension(theta_field_size,jpspec) :: zdryrt
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_CTL_BE'
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
    integer :: widx1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    nr = jpbk + jptk + jphk
    secs_per_step = REAL(chem_timestep)
    if (.NOT.ALLOCATED(dflux)) then
      ALLOCATE(dflux(1:theta_field_size,1:nr))
    end if
    if (.NOT.ALLOCATED(ystore) .AND. uph2so4inaer == 1) then
      ALLOCATE(ystore(1:theta_field_size))
    end if
    !NO!$acc kernels
    ix = 0
    jy = 0
    !NO!$acc loop gang vector independent collapse(2)
    do idx = 1, jpspec, 1
      do idx_1 = 1, model_levels, 1
        dpd_dummy(idx_1,idx) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_2 = 1, jpspec, 1
      do idx_3 = 1, model_levels, 1
        dpw_dummy(idx_3,idx_2) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_4 = 1, jpnr, 1
      do idx_5 = 1, model_levels, 1
        prk_dummy(idx_5,idx_4) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_6 = 1, jpspec, 1
      do idx_7 = 1, model_levels, 1
        y_dummy(idx_7,idx_6) = 0.0
      enddo
    enddo
    !NO!$acc end kernels
    if (firstcall) then
      n_be_calls = INT(secs_per_step / ukca_config%dts0)
      if (n_be_calls == 0) then
        n_be_calls = 1
      end if
      dts = secs_per_step / n_be_calls
      if (printstatus >= prstatus_oper) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A12,I7,A6,E10.1)') 'n_be_calls: ', n_be_calls, ' dts: ', dts
        call umprint(ummessage, src='ukca_chemistry_ctl_be')
      end if
      if (.NOT.lvmr) then
        errcode = 1
        cmessage(:) = ' lvmr must be set to true for BE offline chemistry'
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    end if
    nlev_with_ddep2(:) = RESHAPE(nlev_with_ddep(:,:), [theta_field_size])
    if (.NOT.ALLOCATED(rc)) then
      ALLOCATE(rc(1:theta_field_size,1:nr))
    end if
    do k = 1, k_be_top, 1
      !NO!$acc kernels
      kcs = (k - 1) * theta_field_size + 1
      kce = k * theta_field_size
      !NO!$acc loop gang vector independent collapse(2)
      do idx_8 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
        do idx_9 = LBOUND(zdryrt2, dim=1), UBOUND(zdryrt2, dim=1), 1
          zdryrt2(idx_9,idx_8) = 0.0e0
        enddo
      enddo
      !NO!$acc end kernels
      if (ukca_config%l_ukca_intdd) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_10 = LBOUND(blmask, dim=1), UBOUND(blmask, dim=1), 1
          blmask(idx_10) = k <= nlev_with_ddep2(idx_10)
        enddo
        !NO!$acc loop gang vector independent collapse(2)
        do l = 1, jpdd, 1
          do widx1 = 1, theta_field_size, 1
            if (blmask(widx1)) then
              zdryrt2(widx1,l) = dryrt(widx1,l)
            end if
          enddo
        enddo
        !NO!$acc end kernels
      else
        if (k == 1) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do idx_11 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
            do idx_12 = LBOUND(zdryrt2, dim=1), UBOUND(zdryrt2, dim=1), 1
              zdryrt2(idx_12,idx_11) = dryrt(idx_12,idx_11)
            enddo
          enddo
          !NO!$acc end kernels
        end if
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_13 = LBOUND(zp, dim=1), UBOUND(zp, dim=1), 1
        zp(idx_13) = p_theta_levels(idx_13 + (kcs - LBOUND(zp, dim=1)))
      enddo
      !NO!$acc loop gang vector independent
      do idx_14 = LBOUND(zt, dim=1), UBOUND(zt, dim=1), 1
        zt(idx_14) = temp(idx_14 + (kcs - LBOUND(zt, dim=1)))
      enddo
      !NO!$acc end kernels
      if (ukca_config%l_fix_ukca_offox_h2o_fac) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_15 = LBOUND(zq, dim=1), UBOUND(zq, dim=1), 1
          zq(idx_15) = q(idx_15 + (kcs - LBOUND(zq, dim=1))) / c_h2o
        enddo
        !NO!$acc end kernels
      else
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_16 = LBOUND(zq, dim=1), UBOUND(zq, dim=1), 1
          zq(idx_16) = q(idx_16 + (kcs - LBOUND(zq, dim=1)))
        enddo
        !NO!$acc end kernels
      end if
      if (k <= model_levels) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_17 = LBOUND(zclw, dim=1), UBOUND(zclw, dim=1), 1
          zclw(idx_17) = qcl(idx_17 + (kcs - LBOUND(zclw, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_18 = LBOUND(zfcloud, dim=1), UBOUND(zfcloud, dim=1), 1
          zfcloud(idx_18) = cloud_frac(idx_18 + (kcs - LBOUND(zfcloud, dim=1)))
        enddo
        !NO!$acc end kernels
      else
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_19 = LBOUND(zclw, dim=1), UBOUND(zclw, dim=1), 1
          zclw(idx_19) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_20 = LBOUND(zfcloud, dim=1), UBOUND(zfcloud, dim=1), 1
          zfcloud(idx_20) = 0.0
        enddo
        !NO!$acc end kernels
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do js = 1, jpctr, 1
        do idx_21 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
          zftr(idx_21,js) = tracer(idx_21 + (kcs - LBOUND(zftr, dim=1)),js) / c_species(js)
        enddo
      enddo
      !NO!$acc end kernels
      call ukca_be_wetdep(theta_field_size, wetrt(:,k,:), zwetrt)
      call ukca_be_drydep(k, theta_field_size, nlev_with_ddep2, zdryrt2, zdryrt)
      if (ndepw /= 0) then
        call ukca_wetdep(wetrt(:,k,:), theta_field_size)
      end if
      if (ndepd /= 0) then
        call ukca_drydep(k, zdryrt2, theta_field_size)
      end if
      call ukca_chemco_be_offline(theta_field_size, zt, zp, zq, zfcloud, zclw, h_plus(kcs:kce))
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do js = 1, jpctr, 1
        do idx_22 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
          f(idx_22,js) = zftr(idx_22 + (LBOUND(zftr, dim=1) - LBOUND(f, dim=1)),js) * tnd(idx_22 + (LBOUND(tnd, dim=1) - LBOUND(f, dim=1)))
        enddo
      enddo
      !NO!$acc end kernels
      jit = 0
      call asad_ftoy(ofirst, nitfg, jit, theta_field_size, ix, jy, k)
      lflux = l_asad_use_chem_diags .AND. (l_asad_use_flux_rxns .OR. l_asad_use_wetdep .OR. l_asad_use_drydep)
      if (uph2so4inaer == 1) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_23 = LBOUND(ystore, dim=1), UBOUND(ystore, dim=1), 1
          ystore(idx_23) = y(idx_23 + (LBOUND(y, dim=1) - LBOUND(ystore, dim=1)),nn_h2so4)
        enddo
        !NO!$acc end kernels
      end if
      call ukca_deriv_offline(nr, n_be_calls, theta_field_size, dts, y, zdryrt, zwetrt, so2_wetox_h2o2, so2_wetox_o3, so2_dryox_oh, lflux, dflux)
      if (uph2so4inaer == 1) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_24 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_24,nn_h2so4) = ystore(idx_24 + (LBOUND(ystore, dim=1) - LBOUND(y, dim=1)))
        enddo
        !NO!$acc end kernels
      end if
      do j = 1, jpctr, 1
        do i = 1, jpspec, 1
          if (advt(j) == speci(i)) then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_25 = kcs, kce, 1
              tracer(idx_25,j) = y(idx_25 + (LBOUND(y, dim=1) - kcs),i) / tnd(idx_25 + (LBOUND(tnd, dim=1) - kcs)) * c_species(j)
            enddo
            !NO!$acc end kernels
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
      enddo
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_26 = kcs, kce, 1
        delso2_wet_h2o2(idx_26) = so2_wetox_h2o2(idx_26 + (LBOUND(so2_wetox_h2o2, dim=1) - kcs))
      enddo
      !NO!$acc loop gang vector independent
      do idx_27 = kcs, kce, 1
        delso2_wet_o3(idx_27) = so2_wetox_o3(idx_27 + (LBOUND(so2_wetox_o3, dim=1) - kcs))
      enddo
      !NO!$acc loop gang vector independent
      do idx_28 = kcs, kce, 1
        delh2so4_chem(idx_28) = so2_dryox_oh(idx_28 + (LBOUND(so2_dryox_oh, dim=1) - kcs))
      enddo
      !NO!$acc loop gang vector independent
      do idx_29 = kcs, kce, 1
        delso2_drydep(idx_29) = zdryrt(idx_29 + (LBOUND(zdryrt, dim=1) - kcs),nn_so2) * y(idx_29 + (LBOUND(y, dim=1) - kcs),nn_so2) * dts
      enddo
      !NO!$acc loop gang vector independent
      do idx_30 = kcs, kce, 1
        delso2_wetdep(idx_30) = zwetrt(idx_30 + (LBOUND(zwetrt, dim=1) - kcs),nn_so2) * y(idx_30 + (LBOUND(y, dim=1) - kcs),nn_so2) * dts
      enddo
      !NO!$acc end kernels
      if (lflux) then
        !NO!$acc kernels
        !NO!$acc loop seq
        do i = 1, jpbk, 1
          !NO!$acc loop gang vector independent
          do idx_31 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
            prk(idx_31,ibimol(i)) = dflux(idx_31 + (LBOUND(dflux, dim=1) - LBOUND(prk, dim=1)),i)
          enddo
        enddo
        j = jpbk
        !NO!$acc loop seq
        do i = 1, jptk, 1
          !NO!$acc loop gang vector independent
          do idx_32 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
            prk(idx_32,itrimol(i)) = dflux(idx_32 + (LBOUND(dflux, dim=1) - LBOUND(prk, dim=1)),i + j)
          enddo
        enddo
        j = jpbk + jptk
        !NO!$acc loop seq
        do i = 1, jphk, 1
          !NO!$acc loop gang vector independent
          do idx_33 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
            prk(idx_33,ihetero(i)) = dflux(idx_33 + (LBOUND(dflux, dim=1) - LBOUND(prk, dim=1)),i + j)
          enddo
        enddo
        !NO!$acc end kernels
        call asad_chemical_diagnostics(row_length, rows, model_levels, theta_field_size, dpd_dummy, dpw_dummy, prk_dummy, y_dummy, ix, jy, k, volume, ierr)
      end if
    enddo
    if (ALLOCATED(rc)) then
      DEALLOCATE(rc)
    end if
    if (ALLOCATED(dflux)) then
      DEALLOCATE(dflux)
    end if
    if (ALLOCATED(ystore)) then
      DEALLOCATE(ystore)
    end if
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (tracer(:, n_h2o2) > h2o2_offline(:)) tracer(:, n_h2o2) = h2o2_offline(:)
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_chemistry_ctl_be
  subroutine ukca_chemco_be_offline(theta_field_size, zt, zp, zq, zfcloud, zclw, h_plus_1d_arr)
    use asad_bimol_mod, only : asad_bimol
    use asad_findreaction_mod, only : asad_findreaction
    use asad_hetero_mod, only : asad_hetero
    use asad_totnud_mod, only : asad_totnud
    use asad_trimol_mod, only : asad_trimol
    use ukca_chem_defs_mod, only : ratb_defs, rath_defs, ratt_defs
    use asad_mod, only : jpbk, jphk, jpspb, jpsph, jpspt, jptk, nbrkx, nhrkx, ntrkx, p, rk, spb, sph, spt, t, wp
    use ukca_missing_data_mod, only : imdi, rmdi
    use umPrintMgr, only : ummessage, umprint
    use ereport_mod, only : ereport
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: theta_field_size
    real, dimension(theta_field_size), intent(in) :: zt
    real, dimension(theta_field_size), intent(in) :: zp
    real, dimension(theta_field_size), intent(in) :: zq
    real, dimension(theta_field_size), intent(in) :: zfcloud
    real, dimension(theta_field_size), intent(in) :: zclw
    real, dimension(theta_field_size), intent(in) :: h_plus_1d_arr
    integer :: i
    integer :: j
    integer :: npnts
    integer :: errcode
    real, dimension(theta_field_size,2) :: dummy
    CHARACTER(LEN = 10) :: r1
    CHARACTER(LEN = 10) :: r2
    CHARACTER(LEN = 10) :: p1
    CHARACTER(LEN = 10) :: p2
    CHARACTER(LEN = 10) :: prods(2)
    CHARACTER(LEN = 72) :: cmessage
    logical, save :: first = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMCO_BE_OFFLINE'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    npnts = theta_field_size
    !NO!$acc loop gang vector independent
    do idx = 1, npnts, 1
      p(idx) = zp(idx)
    enddo
    !NO!$acc loop gang vector independent
    do idx_1 = 1, npnts, 1
      t(idx_1) = zt(idx_1)
    enddo
    !NO!$acc loop gang vector independent
    do idx_2 = 1, npnts, 1
      wp(idx_2) = zq(idx_2)
    enddo
    !NO!$acc end kernels
    call asad_totnud(npnts)
    call asad_bimol(npnts)
    call asad_trimol(npnts)
    call asad_hetero(npnts, zfcloud, zclw, dummy, h_plus_1d_arr)
    if (first) then
      if (.NOT.ALLOCATED(ibimol)) then
        ALLOCATE(ibimol(1:jpbk))
      end if
      if (.NOT.ALLOCATED(itrimol)) then
        ALLOCATE(itrimol(1:jptk))
      end if
      if (.NOT.ALLOCATED(ihetero)) then
        ALLOCATE(ihetero(1:jphk))
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_3 = LBOUND(ibimol, dim=1), UBOUND(ibimol, dim=1), 1
        ibimol(idx_3) = imdi
      enddo
      !NO!$acc loop gang vector independent
      do idx_4 = LBOUND(itrimol, dim=1), UBOUND(itrimol, dim=1), 1
        itrimol(idx_4) = imdi
      enddo
      !NO!$acc loop gang vector independent
      do idx_5 = LBOUND(ihetero, dim=1), UBOUND(ihetero, dim=1), 1
        ihetero(idx_5) = imdi
      enddo
      frac_dms_so2 = rmdi
      frac_dms_dmso = rmdi
      frac_dmso_so2 = rmdi
      frac_monoterp_sec_org(3) = rmdi
      !NO!$acc loop gang vector independent collapse(2)
      do idx_6 = LBOUND(rc, dim=2), UBOUND(rc, dim=2), 1
        do idx_7 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
          rc(idx_7,idx_6) = rmdi
        enddo
      enddo
      errcode = 0
      !NO!$acc end kernels
      do i = 1, jpbk, 1
        r1(:) = ratb_defs(i)%react1
        r2(:) = ratb_defs(i)%react2
        p1(:) = ratb_defs(i)%prod1
        p2(:) = ratb_defs(i)%prod2
        prods(:) = [p1, p2]
        ibimol(i) = asad_findreaction(r1, r2, prods, 2, spb, nbrkx, jpbk + 1, jpspb)
        if (ibimol(i) == imdi) then
          errcode = i
        end if
        if (p1(1:4) == 'SO2 ' .AND. p2(1:4) == 'DMSO') then
          frac_dms_so2 = ratb_defs(i)%pyield1
          frac_dms_dmso = ratb_defs(i)%pyield2
        end if
        if (r1(1:4) == 'DMSO') then
          frac_dmso_so2 = ratb_defs(i)%pyield1
        end if
        if (r1(1:8) == 'Monoterp' .AND. r2(1:4) == 'OH  ') then
          frac_monoterp_sec_org(1) = ratb_defs(i)%pyield1
        end if
        if (r1(1:8) == 'Monoterp' .AND. r2(1:4) == 'O3  ') then
          frac_monoterp_sec_org(2) = ratb_defs(i)%pyield1
        end if
        if (r1(1:8) == 'Monoterp' .AND. r2(1:4) == 'NO3 ') then
          frac_monoterp_sec_org(3) = ratb_defs(i)%pyield1
        end if
      enddo
      if (errcode /= 0) then
        cmessage(:) = 'Unidentified reaction: ' // ratb_defs(errcode) % react1 // ratb_defs(errcode) % react2
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      if (ABS(frac_dms_so2 - rmdi) < EPSILON(0.0) .OR. ABS(frac_dms_dmso - rmdi) < EPSILON(0.0) .OR. ABS(frac_dmso_so2 - rmdi) < EPSILON(0.0) .OR. ANY(ABS(frac_monoterp_sec_org - rmdi) < EPSILON(0.0))) then
        errcode = 1
        cmessage(:) = ' One or more product yields are undefined'
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      errcode = 0
      do i = 1, jptk, 1
        r1(:) = ratt_defs(i)%react1
        r2(:) = ratt_defs(i)%react2
        p1(:) = ratt_defs(i)%prod1
        p2(:) = ratt_defs(i)%prod2
        prods(:) = [p1, p2]
        itrimol(i) = asad_findreaction(r1, r2, prods, 2, spt, ntrkx, jptk + 1, jpspt)
        if (itrimol(i) == imdi) then
          errcode = i
        end if
      enddo
      if (errcode /= 0) then
        cmessage(:) = 'Unidentified reaction: ' // ratt_defs(errcode) % react1 // ratt_defs(errcode) % react2
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      errcode = 0
      do i = 1, jphk, 1
        r1(:) = rath_defs(i)%react1
        r2(:) = rath_defs(i)%react2
        p1(:) = rath_defs(i)%prod1
        p2(:) = rath_defs(i)%prod2
        prods(:) = [p1, p2]
        ihetero(i) = asad_findreaction(r1, r2, prods, 2, sph, nhrkx, jphk + 1, jpsph)
        if (ihetero(i) == imdi) then
          errcode = i
        end if
      enddo
      if (errcode /= 0) then
        cmessage(:) = 'Unidentified reaction: ' // rath_defs(errcode) % react1 // rath_defs(errcode) % react2
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(2)
    do i = 1, jpbk, 1
      do idx_8 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_8,i) = rk(idx_8 + (LBOUND(rk, dim=1) - LBOUND(rc, dim=1)),ibimol(i))
      enddo
    enddo
    j = jpbk
    !NO!$acc loop gang vector independent collapse(2)
    do i = 1, jptk, 1
      do idx_9 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_9,i + j) = rk(idx_9 + (LBOUND(rk, dim=1) - LBOUND(rc, dim=1)),itrimol(i))
      enddo
    enddo
    j = jpbk + jptk
    !NO!$acc loop gang vector independent collapse(2)
    do i = 1, jphk, 1
      do idx_10 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_10,i + j) = rk(idx_10 + (LBOUND(rk, dim=1) - LBOUND(rc, dim=1)),ihetero(i))
      enddo
    enddo
    !NO!$acc end kernels
    if (first .AND. ANY(ABS(rc - rmdi) < EPSILON(0.0))) then
      errcode = 1
      cmessage(:) = 'Missing rate coefficient data in array rc'
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A55,I2)') cmessage // ' at location: ', MAXLOC(rc)
      call umprint(ummessage, src='ukca_chemco_be_offline')
      call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
    end if
    if (first) then
      if (frac_dms_so2 > 1.0 .OR. frac_dms_dmso > 1.0 .OR. frac_dmso_so2 > 1.0 .OR. ANY(frac_monoterp_sec_org(:) > 1.0)) then
        errcode = 1
        cmessage(:) = 'A yield fraction > 1.0 was encountered'
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    end if
    first = .false.
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_chemco_be_offline
  subroutine ukca_deriv_offline(nr, n_be_calls, theta_field_size, dts, y, dryrt, wetrt, so2_wetox_h2o2, so2_wetox_o3, so2_dryox_oh, lflux, dflux)
    use asad_mod, only : jpspec
    use ukca_um_legacy_mod, only : oneover_v
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr
    integer, intent(in) :: n_be_calls
    integer, intent(in) :: theta_field_size
    real, intent(in) :: dts
    real, dimension(theta_field_size,jpspec), intent(in) :: dryrt
    real, dimension(theta_field_size,jpspec), intent(in) :: wetrt
    logical, intent(in) :: lflux
    real, dimension(theta_field_size), intent(out) :: so2_wetox_h2o2
    real, dimension(theta_field_size), intent(out) :: so2_wetox_o3
    real, dimension(theta_field_size), intent(out) :: so2_dryox_oh
    real, dimension(theta_field_size,nr), intent(out) :: dflux
    real, dimension(theta_field_size,jpspec), intent(inout) :: y
    integer :: j
    integer :: n
    integer :: npnts
    real, dimension(theta_field_size) :: p
    real, dimension(theta_field_size) :: l
    real, dimension(theta_field_size,jpspec) :: yp
    real, dimension(theta_field_size) :: tmp1
    real, dimension(theta_field_size) :: tmp2
    real :: f_inv
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DERIV_OFFLINE'
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
    integer :: idx_41
    integer :: idx_42
    integer :: idx_43
    integer :: idx_44
    integer :: idx_45
    integer :: idx_46
    integer :: idx_47
    integer :: idx_48
    integer :: idx_49
    integer :: idx_50

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    npnts = theta_field_size
    !NO!$acc loop gang vector independent
    do idx = LBOUND(so2_wetox_h2o2, dim=1), UBOUND(so2_wetox_h2o2, dim=1), 1
      so2_wetox_h2o2(idx) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_1 = LBOUND(so2_wetox_o3, dim=1), UBOUND(so2_wetox_o3, dim=1), 1
      so2_wetox_o3(idx_1) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_2 = LBOUND(so2_dryox_oh, dim=1), UBOUND(so2_dryox_oh, dim=1), 1
      so2_dryox_oh(idx_2) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_3 = LBOUND(dflux, dim=2), UBOUND(dflux, dim=2), 1
      do idx_4 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
        dflux(idx_4,idx_3) = 0.0
      enddo
    enddo
    !NO!$acc end kernels
    do n = 1, n_be_calls, 1
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do j = 1, jpspec, 1
        do idx_5 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
          yp(idx_5,j) = y(idx_5,j)
        enddo
      enddo
      !NO!$acc loop gang vector independent
      do idx_6 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_6) = 0.0
      enddo
      !NO!$acc loop gang vector independent
      do idx_7 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_7) = dryrt(idx_7,2) + wetrt(idx_7,2) + (rc(idx_7 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),1) + rc(idx_7 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),2)) * y(idx_7,9) + rc(idx_7 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),3) * y(idx_7,10)
      enddo
      !NO!$acc loop gang vector independent
      do idx_8 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_8,3) = (yp(idx_8,3) + dts * p(idx_8)) / (1.0 + dts * l(idx_8))
      enddo
      !NO!$acc loop gang vector independent
      do idx_9 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_9) = 1.0 + dts * l(idx_9)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_10 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_10,2) = (yp(idx_10,2) + dts * p(idx_10)) * tmp1(idx_10)
      enddo
      !NO!$acc loop gang vector independent
      do idx_11 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_11) = rc(idx_11 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),2) * y(idx_11,2) * y(idx_11,9) * frac_dms_dmso
      enddo
      !NO!$acc loop gang vector independent
      do idx_12 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_12) = dryrt(idx_12,5) + wetrt(idx_12,5) + rc(idx_12 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),4) * y(idx_12,9)
      enddo
      !NO!$acc loop gang vector independent
      do idx_13 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_13) = 1.0 + dts * l(idx_13)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_14 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_14,5) = (yp(idx_14,5) + dts * p(idx_14)) * tmp1(idx_14)
      enddo
      !NO!$acc loop gang vector independent
      do idx_15 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_15) = (rc(idx_15 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),1) + rc(idx_15 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),2) * frac_dms_so2) * y(idx_15,2) * y(idx_15,9) + rc(idx_15 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),3) * y(idx_15,2) * y(idx_15,10) + rc(idx_15 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),4) * y(idx_15,5) * y(idx_15,9) * frac_dmso_so2
      enddo
      !NO!$acc loop gang vector independent
      do idx_16 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_16) = dryrt(idx_16,3) + wetrt(idx_16,3) + rc(idx_16 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),10) * y(idx_16,9) + rc(idx_16 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),11) * y(idx_16,1) + (rc(idx_16 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),12) + rc(idx_16 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),13)) * y(idx_16,8)
      enddo
      !NO!$acc loop gang vector independent
      do idx_17 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_17) = 1.0 + dts * l(idx_17)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_18 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_18,3) = (yp(idx_18,3) + dts * p(idx_18)) * tmp1(idx_18)
      enddo
      !NO!$acc loop gang vector independent
      do idx_19 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_19) = rc(idx_19 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),10) * y(idx_19,3) * y(idx_19,9)
      enddo
      !NO!$acc loop gang vector independent
      do idx_20 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_20) = dryrt(idx_20,4) + wetrt(idx_20,4)
      enddo
      !NO!$acc loop gang vector independent
      do idx_21 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_21) = 1.0 + dts * l(idx_21)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_22 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_22,4) = (yp(idx_22,4) + dts * p(idx_22)) * tmp1(idx_22)
      enddo
      !NO!$acc loop gang vector independent
      do idx_23 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_23) = 0.0
      enddo
      !NO!$acc loop gang vector independent
      do idx_24 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_24) = dryrt(idx_24,6) + wetrt(idx_24,6) + rc(idx_24 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),5) * y(idx_24,9) + rc(idx_24 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),6) * y(idx_24,8) + rc(idx_24 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),7) * y(idx_24,10)
      enddo
      !NO!$acc loop gang vector independent
      do idx_25 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_25) = 1.0 + dts * l(idx_25)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_26 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_26,6) = (yp(idx_26,6) + dts * p(idx_26)) * tmp1(idx_26)
      enddo
      !NO!$acc loop gang vector independent
      do idx_27 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_27) = rc(idx_27 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),5) * y(idx_27,6) * y(idx_27,9) * frac_monoterp_sec_org(1) + rc(idx_27 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),6) * y(idx_27,6) * y(idx_27,8) * frac_monoterp_sec_org(2) + rc(idx_27 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),7) * y(idx_27,6) * y(idx_27,10) * frac_monoterp_sec_org(3)
      enddo
      !NO!$acc loop gang vector independent
      do idx_28 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_28) = dryrt(idx_28,7) + wetrt(idx_28,7)
      enddo
      !NO!$acc loop gang vector independent
      do idx_29 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_29) = 1.0 + dts * l(idx_29)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_30 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_30,7) = (yp(idx_30,7) + dts * p(idx_30)) * tmp1(idx_30)
      enddo
      !NO!$acc loop gang vector independent
      do idx_31 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_31) = rc(idx_31 + (LBOUND(rc, dim=1) - LBOUND(p, dim=1)),8) * y(idx_31,11) * y(idx_31,11)
      enddo
      !NO!$acc loop gang vector independent
      do idx_32 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_32) = dryrt(idx_32,1) + wetrt(idx_32,1) + rc(idx_32 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),9) * y(idx_32,9) + rc(idx_32 + (LBOUND(rc, dim=1) - LBOUND(l, dim=1)),11) * y(idx_32,3)
      enddo
      !NO!$acc loop gang vector independent
      do idx_33 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
        tmp2(idx_33) = 1.0 + dts * l(idx_33)
      enddo
      !NO!$acc end kernels
      call oneover_v(npnts, tmp2, tmp1)
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_34 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_34,1) = (yp(idx_34,1) + dts * p(idx_34)) * tmp1(idx_34)
      enddo
      !NO!$acc loop gang vector independent
      do idx_35 = LBOUND(so2_wetox_h2o2, dim=1), UBOUND(so2_wetox_h2o2, dim=1), 1
        so2_wetox_h2o2(idx_35) = so2_wetox_h2o2(idx_35) + rc(idx_35 + (LBOUND(rc, dim=1) - LBOUND(so2_wetox_h2o2, dim=1)),11) * y(idx_35,3) * y(idx_35,1) * dts
      enddo
      !NO!$acc loop gang vector independent
      do idx_36 = LBOUND(so2_wetox_o3, dim=1), UBOUND(so2_wetox_o3, dim=1), 1
        so2_wetox_o3(idx_36) = so2_wetox_o3(idx_36) + ((rc(idx_36 + (LBOUND(rc, dim=1) - LBOUND(so2_wetox_o3, dim=1)),12) + rc(idx_36 + (LBOUND(rc, dim=1) - LBOUND(so2_wetox_o3, dim=1)),13)) * y(idx_36,3) * y(idx_36,8)) * dts
      enddo
      !NO!$acc loop gang vector independent
      do idx_37 = LBOUND(so2_dryox_oh, dim=1), UBOUND(so2_dryox_oh, dim=1), 1
        so2_dryox_oh(idx_37) = so2_dryox_oh(idx_37) + rc(idx_37 + (LBOUND(rc, dim=1) - LBOUND(so2_dryox_oh, dim=1)),10) * y(idx_37,3) * y(idx_37,9) * dts
      enddo
      !NO!$acc end kernels
      if (lflux) then
        !NO!$acc kernels
        f_inv = 1.0 / n_be_calls
        !NO!$acc loop gang vector independent
        do idx_38 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_38,1) = dflux(idx_38,1) + rc(idx_38 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),1) * y(idx_38,2) * y(idx_38,9) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_39 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_39,2) = dflux(idx_39,1) + rc(idx_39 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),2) * y(idx_39,2) * y(idx_39,9) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_40 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_40,3) = dflux(idx_40,3) + rc(idx_40 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),3) * y(idx_40,2) * y(idx_40,10) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_41 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_41,4) = dflux(idx_41,4) + rc(idx_41 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),4) * y(idx_41,5) * y(idx_41,9) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_42 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_42,5) = dflux(idx_42,5) + rc(idx_42 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),5) * y(idx_42,6) * y(idx_42,9) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_43 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_43,6) = dflux(idx_43,6) + rc(idx_43 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),6) * y(idx_43,6) * y(idx_43,8) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_44 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_44,7) = dflux(idx_44,7) + rc(idx_44 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),7) * y(idx_44,6) * y(idx_44,10) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_45 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_45,8) = dflux(idx_45,8) + rc(idx_45 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),8) * y(idx_45,11) * y(idx_45,11) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_46 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_46,9) = dflux(idx_46,9) + rc(idx_46 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),9) * y(idx_46,1) * y(idx_46,9) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_47 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_47,10) = dflux(idx_47,10) + rc(idx_47 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),10) * y(idx_47,3) * y(idx_47,9) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_48 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_48,11) = dflux(idx_48,11) + rc(idx_48 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),11) * y(idx_48,3) * y(idx_48,1) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_49 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_49,12) = dflux(idx_49,12) + rc(idx_49 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),12) * y(idx_49,3) * y(idx_49,8) * f_inv
        enddo
        !NO!$acc loop gang vector independent
        do idx_50 = LBOUND(dflux, dim=1), UBOUND(dflux, dim=1), 1
          dflux(idx_50,13) = dflux(idx_50,13) + rc(idx_50 + (LBOUND(rc, dim=1) - LBOUND(dflux, dim=1)),13) * y(idx_50,3) * y(idx_50,8) * f_inv
        enddo
        !NO!$acc end kernels
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_deriv_offline

end module ukca_chemistry_ctl_be_mod
