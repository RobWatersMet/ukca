module ukca_chemistry_ctl_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_CTL_MOD'
  public

  contains
  subroutine ukca_chemistry_ctl(row_length, rows, model_levels, theta_field_size, tot_n_pnts, ntracers, istore_h2so4, pres, temp, q, qcf, qcl, tracer, all_ntp, cloud_frac, photol_rates, shno3_3d, volume, have_nat, uph2so4inaer, delso2_wet_h2o2, delso2_wet_o3, delh2so4_chem, so4_sa, atm_ch4_mol, atm_co_mol, atm_n2o_mol, atm_cf2cl2_mol, atm_cfcl3_mol, atm_mebr_mol, atm_h2_mol, h_plus, zdryrt, zwetrt, nlev_with_ddep, co2_interactive, l_stratosphere, firstcall)
    use asad_mod, only : advt, cdt, ctype, ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3, jpcspf, jpctr, jpdd, jpdw, jpnr, jppj, jpro2, jpspec, nadvt, nlnaro2, nprkx, o1d_in_ss, o3p_in_ss, rk, specf, speci, sph2o, sphno3, spro2, tnd, y, za
    use asad_chem_flux_diags, only : asad_chemical_diagnostics, asad_psc_diagnostic, l_asad_use_chem_diags, l_asad_use_drydep, l_asad_use_flux_rxns, l_asad_use_psc_diagnostic, l_asad_use_rxn_rates, l_asad_use_wetdep
    use asad_cdrive_mod, only : asad_cdrive
    use ukca_cspecies, only : c_na_species, c_species, n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2, n_h2so4, n_mebr, n_n2o, nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh, nn_so2
    use ukca_constants, only : c_co2, c_h2o, c_hono2, c_o1d, c_o3p
    use ukca_config_constants_mod, only : avogadro
    use ukca_config_specification_mod, only : ukca_config
    use ukca_ntp_mod, only : dim_ntp, name2ntpindex, ntp_type
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use ukca_missing_data_mod, only : rmdi
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: theta_field_size
    integer, intent(in) :: tot_n_pnts
    integer, intent(in) :: ntracers
    integer, intent(in) :: uph2so4inaer
    integer, intent(in) :: istore_h2so4
    integer, dimension(theta_field_size), intent(in) :: nlev_with_ddep
    real, dimension(tot_n_pnts), intent(in) :: pres
    real, dimension(tot_n_pnts), intent(in) :: temp
    real, dimension(tot_n_pnts), intent(in) :: volume
    real, dimension(tot_n_pnts), intent(in) :: h_plus
    real, dimension(tot_n_pnts), intent(in) :: qcf
    real, dimension(tot_n_pnts), intent(in) :: qcl
    real, dimension(tot_n_pnts), intent(in) :: cloud_frac
    real, dimension(tot_n_pnts), intent(in) :: so4_sa
    real, dimension(theta_field_size,jpdd), intent(in) :: zdryrt
    real, dimension(tot_n_pnts,jpdw), intent(in) :: zwetrt
    real, dimension(tot_n_pnts,jppj), intent(in) :: photol_rates
    real, dimension(tot_n_pnts), intent(in) :: co2_interactive
    real, dimension(tot_n_pnts), intent(out) :: shno3_3d
    real, dimension(tot_n_pnts), intent(inout) :: q
    real, dimension(tot_n_pnts,ntracers), intent(inout) :: tracer
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wet_h2o2
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wet_o3
    real, dimension(tot_n_pnts), intent(inout) :: delh2so4_chem
    real, dimension(tot_n_pnts), intent(inout) :: atm_ch4_mol
    real, dimension(tot_n_pnts), intent(inout) :: atm_co_mol
    real, dimension(tot_n_pnts), intent(inout) :: atm_n2o_mol
    real, dimension(tot_n_pnts), intent(inout) :: atm_cf2cl2_mol
    real, dimension(tot_n_pnts), intent(inout) :: atm_cfcl3_mol
    real, dimension(tot_n_pnts), intent(inout) :: atm_mebr_mol
    real, dimension(tot_n_pnts), intent(inout) :: atm_h2_mol
    type(ntp_type), dimension(dim_ntp), intent(inout) :: all_ntp
    logical, dimension(tot_n_pnts), intent(in) :: have_nat
    logical, dimension(tot_n_pnts), intent(in) :: l_stratosphere
    logical, intent(in) :: firstcall
    integer :: ix
    integer :: jy
    integer :: i
    integer :: j
    integer :: js
    integer :: jtr
    integer :: jro2
    integer :: jna
    integer :: jspf
    integer :: k
    integer :: l
    integer :: kcs
    integer :: kce
    integer :: ierr
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    CHARACTER(LEN = 10) :: prods(2)
    CHARACTER(LEN = 10) :: prods3(3)
    logical, dimension(theta_field_size) :: ddmask
    real, dimension(model_levels,jpspec) :: dpd_dummy
    real, dimension(model_levels,jpspec) :: dpw_dummy
    real, dimension(model_levels,jpnr) :: prk_dummy
    real, dimension(model_levels,jpspec) :: y_dummy
    real, dimension(model_levels) :: fpsc1_dummy
    real, dimension(model_levels) :: fpsc2_dummy
    real, allocatable, dimension(:) :: ystore
    real, dimension(theta_field_size,jpcspf) :: zftr
    real, dimension(theta_field_size,jpcspf) :: cdot
    real, dimension(theta_field_size) :: zq
    real, dimension(theta_field_size) :: co2_1d
    real, dimension(theta_field_size,jppj) :: zprt1d
    real, dimension(theta_field_size,jpdd) :: zdryrt2
    real, dimension(theta_field_size,2) :: rc_het
    real, dimension(tot_n_pnts,dim_ntp) :: ntp_data
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_CTL'
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
    integer :: widx1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    ix = 0
    jy = 0
    !$acc loop gang vector independent collapse(2)
    do idx = 1, jpspec, 1
      do idx_1 = 1, model_levels, 1
        dpd_dummy(idx_1,idx) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_2 = 1, jpspec, 1
      do idx_3 = 1, model_levels, 1
        dpw_dummy(idx_3,idx_2) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_4 = 1, model_levels, 1
      fpsc1_dummy(idx_4) = 0.0
    enddo
    !$acc loop gang vector independent
    do idx_5 = 1, model_levels, 1
      fpsc2_dummy(idx_5) = 0.0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_6 = 1, jpnr, 1
      do idx_7 = 1, model_levels, 1
        prk_dummy(idx_7,idx_6) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_8 = 1, jpspec, 1
      do idx_9 = 1, model_levels, 1
        y_dummy(idx_9,idx_8) = 0.0
      enddo
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_het_psc) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_10 = 1, tot_n_pnts, 1
        shno3_3d(idx_10) = 0.0
      enddo
      !$acc end kernels
    end if
    do l = 1, dim_ntp, 1
      if (all_ntp(l)%l_required) then
        ntp_data(:,l) = RESHAPE(all_ntp(l)%data_3d(:,:,:), [tot_n_pnts])
      end if
    enddo
    if (.NOT.ALLOCATED(ystore) .AND. uph2so4inaer == 1) then
      ALLOCATE(ystore(1:theta_field_size))
    end if
    do k = 1, model_levels, 1
      kcs = (k - 1) * theta_field_size + 1
      kce = k * theta_field_size
      if (ukca_config%l_ukca_het_psc) then
        if (k <= model_levels) then
          sph2o(:) = qcf(kcs:kce) / c_h2o
        else
          !$acc kernels
          !$acc loop gang vector independent
          do idx_11 = LBOUND(sph2o, dim=1), UBOUND(sph2o, dim=1), 1
            sph2o(idx_11) = 0.0
          enddo
          !$acc end kernels
        end if
      end if
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx_12 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
        do idx_13 = LBOUND(zdryrt2, dim=1), UBOUND(zdryrt2, dim=1), 1
          zdryrt2(idx_13,idx_12) = 0.0e0
        enddo
      enddo
      !$acc end kernels
      if (ukca_config%l_ukca_intdd) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_14 = LBOUND(ddmask, dim=1), UBOUND(ddmask, dim=1), 1
          ddmask(idx_14) = k <= nlev_with_ddep(idx_14)
        enddo
        !$acc loop gang vector independent collapse(2)
        do l = 1, jpdd, 1
          do widx1 = 1, theta_field_size, 1
            if (ddmask(widx1)) then
              zdryrt2(widx1,l) = zdryrt(widx1,l)
            end if
          enddo
        enddo
        !$acc end kernels
      else
        if (k == 1) then
          !$acc kernels
          !$acc loop gang vector independent collapse(2)
          do idx_15 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
            do idx_16 = LBOUND(zdryrt2, dim=1), UBOUND(zdryrt2, dim=1), 1
              zdryrt2(idx_16,idx_15) = zdryrt(idx_16,idx_15)
            enddo
          enddo
          !$acc end kernels
        end if
      end if
      !$acc kernels
      zq(:) = q(kcs:kce) / c_h2o
      !$acc end kernels
      if (ANY(speci(:) == 'CO2       ')) then
        if (ukca_config%l_chem_environ_co2_fld) then
          co2_1d(:) = co2_interactive(kcs:kce) / c_co2
        else
          co2_1d(:) = rmdi
        end if
      end if
      jspf = 0
      do js = 1, jpspec, 1
        do jtr = 1, jpctr, 1
          if (advt(jtr) == speci(js)) then
            jspf = jspf + 1
            zftr(:,jspf) = tracer(kcs:kce,jtr) / c_species(jtr)
          end if
        enddo
        if (ukca_config%l_ukca_ro2_ntp) then
          do jro2 = 1, jpro2, 1
            if (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') then
              l = name2ntpindex(spro2(jro2))
              jna = nlnaro2(jro2)
              if (nadvt(jna) == spro2(jro2)) then
                jspf = jspf + 1
                zftr(:,jspf) = ntp_data(kcs:kce,l) / c_na_species(jna)
              else
                errcode = jro2
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl'
                call umprint(ummessage, src='ukca_chemistry_ctl')
                cmessage(:) = 'ERROR: Indices for RO2 species do not match with nadvt'
                call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
              end if
            end if
          enddo
        end if
      enddo
      if (ukca_config%l_ukca_ro2_ntp) then
        if (jspf /= jpro2 + jpctr) then
          errcode = jspf
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl'
          call umprint(ummessage, src='ukca_chemistry_ctl')
          cmessage(:) = 'ERROR: Number of chemical active species /= jpro2+jpctr'
          call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
        end if
      else
        if (jspf /= jpctr) then
          errcode = jspf
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl'
          call umprint(ummessage, src='ukca_chemistry_ctl')
          cmessage(:) = 'ERROR: Number of chemical active species /= jpctr'
          call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
        end if
      end if
      if (ukca_config%l_ukca_offline) then
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx_17 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
          do idx_18 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
            zprt1d(idx_18,idx_17) = 0.0
          enddo
        enddo
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx_19 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
          do idx_20 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
            zprt1d(idx_20,idx_19) = photol_rates(idx_20 + (kcs - LBOUND(zprt1d, dim=1)),idx_19)
          enddo
        enddo
        !$acc end kernels
      end if
      if (ukca_config%l_ukca_trophet) then
        l = name2ntpindex('het_n2o5  ')
        !$acc kernels
        !$acc loop gang vector independent
        do idx_21 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
          rc_het(idx_21,1) = ntp_data(idx_21 + (kcs - LBOUND(rc_het, dim=1)),l)
        enddo
        !$acc end kernels
        l = name2ntpindex('het_ho2   ')
        !$acc kernels
        !$acc loop gang vector independent
        do idx_22 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
          rc_het(idx_22,2) = ntp_data(idx_22 + (kcs - LBOUND(rc_het, dim=1)),l)
        enddo
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx_23 = LBOUND(rc_het, dim=2), UBOUND(rc_het, dim=2), 1
          do idx_24 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
            rc_het(idx_24,idx_23) = 0.0
          enddo
        enddo
        !$acc end kernels
      end if
      !$acc kernels
      !$acc loop gang vector independent
      do idx_25 = LBOUND(za, dim=1), UBOUND(za, dim=1), 1
        za(idx_25) = so4_sa(idx_25 + (kcs - LBOUND(za, dim=1)))
      enddo
      !$acc end kernels
      if (uph2so4inaer == 1) then
        if (ukca_config%l_fix_ukca_h2so4_ystore) then
          !$acc kernels
          !$acc loop gang vector independent
          do idx_26 = LBOUND(ystore, dim=1), UBOUND(ystore, dim=1), 1
            ystore(idx_26) = zftr(idx_26 + (LBOUND(zftr, dim=1) - LBOUND(ystore, dim=1)),istore_h2so4)
          enddo
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop gang vector independent
          do idx_27 = LBOUND(ystore, dim=1), UBOUND(ystore, dim=1), 1
            ystore(idx_27) = y(idx_27 + (LBOUND(y, dim=1) - LBOUND(ystore, dim=1)),nn_h2so4)
          enddo
          !$acc end kernels
        end if
      end if
      call asad_cdrive(cdot, zftr, pres(kcs:kce), temp(kcs:kce), zq, co2_1d, cloud_frac(kcs:kce), qcl(kcs:kce), ix, jy, k, zdryrt2, zwetrt(kcs:kce,:), rc_het, zprt1d, theta_field_size, have_nat(kcs:kce), l_stratosphere(kcs:kce), h_plus(kcs:kce))
      if (ukca_config%l_ukca_het_psc) then
        if (ANY(sphno3(1:theta_field_size) > 0.0)) then
          shno3_3d(kcs:kce) = sphno3(:) / tnd(:) * c_hono2
        else
          !$acc kernels
          !$acc loop gang vector independent
          do idx_28 = kcs, kce, 1
            shno3_3d(idx_28) = 0.0
          enddo
          !$acc end kernels
        end if
      end if
      if (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) then
        if (ihso3_h2o2 > 0) then
          delso2_wet_h2o2(kcs:kce) = delso2_wet_h2o2(kcs:kce) + rk(:,ihso3_h2o2) * y(:,nn_so2) * y(:,nn_h2o2) * cdt
        end if
        if (ihso3_o3 > 0) then
          delso2_wet_o3(kcs:kce) = delso2_wet_o3(kcs:kce) + rk(:,ihso3_o3) * y(:,nn_so2) * y(:,nn_o3) * cdt
        end if
        if (iso3_o3 > 0) then
          delso2_wet_o3(kcs:kce) = delso2_wet_o3(kcs:kce) + rk(:,iso3_o3) * y(:,nn_so2) * y(:,nn_o3) * cdt
        end if
        if (iso2_oh > 0 .AND. ih2so4_hv > 0) then
          delh2so4_chem(kcs:kce) = delh2so4_chem(kcs:kce) + (rk(:,iso2_oh) * y(:,nn_so2) * y(:,nn_oh) - rk(:,ih2so4_hv) * y(:,nn_h2so4)) * cdt
        else
          if (iso2_oh > 0) then
            delh2so4_chem(kcs:kce) = delh2so4_chem(kcs:kce) + rk(:,iso2_oh) * y(:,nn_so2) * y(:,nn_oh) * cdt
          end if
        end if
        if (uph2so4inaer == 1) then
          if (ukca_config%l_fix_ukca_h2so4_ystore) then
            !$acc kernels
            delh2so4_chem(kcs:kce) = (zftr(:,istore_h2so4) - ystore(:)) / cdt
            !$acc loop gang vector independent
            do idx_29 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
              zftr(idx_29,istore_h2so4) = ystore(idx_29 + (LBOUND(ystore, dim=1) - LBOUND(zftr, dim=1)))
            enddo
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop gang vector independent
            do idx_30 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
              y(idx_30,nn_h2so4) = ystore(idx_30 + (LBOUND(ystore, dim=1) - LBOUND(y, dim=1)))
            enddo
            !$acc end kernels
          end if
        end if
      end if
      if (l_asad_use_chem_diags .AND. (l_asad_use_flux_rxns .OR. l_asad_use_rxn_rates .OR. (l_asad_use_wetdep .OR. l_asad_use_drydep))) then
        call asad_chemical_diagnostics(row_length, rows, model_levels, theta_field_size, dpd_dummy, dpw_dummy, prk_dummy, y_dummy, ix, jy, k, volume, ierr)
      end if
      if (l_asad_use_chem_diags .AND. l_asad_use_psc_diagnostic) then
        call asad_psc_diagnostic(row_length, rows, model_levels, theta_field_size, fpsc1_dummy, fpsc2_dummy, ix, jy, k, ierr)
      end if
      jspf = 0
      do js = 1, jpspec, 1
        do jtr = 1, jpctr, 1
          if (advt(jtr) == speci(js)) then
            jspf = jspf + 1
            tracer(kcs:kce,jtr) = zftr(:,jspf) * c_species(jtr)
          end if
        enddo
        if (ukca_config%l_ukca_ro2_ntp) then
          !$acc kernels
          !$acc loop seq
          do jro2 = 1, jpro2, 1
            if (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') then
              jspf = jspf + 1
              l = name2ntpindex(spro2(jro2))
              jna = nlnaro2(jro2)
              ntp_data(kcs:kce,l) = zftr(:,jspf) * c_na_species(jna)
            end if
          enddo
          !$acc end kernels
        end if
      enddo
      if (o1d_in_ss) then
        l = name2ntpindex('O(1D)     ')
        ntp_data(kcs:kce,l) = y(:,nn_o1d) / tnd(:) * c_o1d
      end if
      if (o3p_in_ss) then
        l = name2ntpindex('O(3P)     ')
        ntp_data(kcs:kce,l) = y(:,nn_o3p) / tnd(:) * c_o3p
      end if
      do jspf = 1, jpcspf, 1
        if (n_ch4 > 0) then
          if (specf(jspf) == advt(n_ch4)) then
            atm_ch4_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
        if (n_co > 0) then
          if (specf(jspf) == advt(n_co)) then
            atm_co_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
        if (n_n2o > 0) then
          if (specf(jspf) == advt(n_n2o)) then
            atm_n2o_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
        if (n_cf2cl2 > 0) then
          if (specf(jspf) == advt(n_cf2cl2)) then
            atm_cf2cl2_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
        if (n_cfcl3 > 0) then
          if (specf(jspf) == advt(n_cfcl3)) then
            atm_cfcl3_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
        if (n_mebr > 0) then
          if (specf(jspf) == advt(n_mebr)) then
            atm_mebr_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
        if (n_h2 > 0) then
          if (specf(jspf) == advt(n_h2)) then
            atm_h2_mol(kcs:kce) = zftr(:,jspf) * tnd(:) * volume(kcs:kce) * 1.0e6 / avogadro
          end if
        end if
      enddo
    enddo
    if (ALLOCATED(ystore)) then
      DEALLOCATE(ystore)
    end if
    do l = 1, dim_ntp, 1
      if (all_ntp(l)%l_required) then
        all_ntp(l)%data_3d(:,:,:) = RESHAPE(ntp_data(:,l), [row_length, rows, model_levels])
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_chemistry_ctl

end module ukca_chemistry_ctl_mod
