module ukca_chemistry_ctl_full_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_CTL_FULL_MOD'
  public

  contains
  subroutine ukca_chemistry_ctl_full(row_length, rows, model_levels, theta_field_size, tot_n_pnts, ntracers, istore_h2so4, pres, temp, q, qcf, qcl, tracer, all_ntp, cloud_frac, photol_rates, shno3, volume, have_nat, uph2so4inaer, delso2_wet_h2o2, delso2_wet_o3, delh2so4_chem, so4_sa, atm_ch4_mol, atm_co_mol, atm_n2o_mol, atm_cf2cl2_mol, atm_cfcl3_mol, atm_mebr_mol, atm_h2_mol, h_plus, zdryrt, zwetrt, nlev_with_ddep, l_stratosphere, co2_interactive, firstcall)
    use asad_mod, only : advt, cdt, ctype, dpd, dpw, fpsc1, fpsc2, ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3, jpcspf, jpctr, jpdd, jpdw, jpnr, jppj, jpro2, jpspec, nadvt, nlnaro2, nprkx, o1d_in_ss, o3p_in_ss, prk, rk, specf, speci, sph2o, sphno3, spro2, tnd, y, za
    use asad_cdrive_mod, only : asad_cdrive
    use asad_chem_flux_diags, only : asad_chemical_diagnostics, asad_psc_diagnostic, l_asad_use_chem_diags, l_asad_use_drydep, l_asad_use_flux_rxns, l_asad_use_psc_diagnostic, l_asad_use_rxn_rates, l_asad_use_wetdep
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
    real, dimension(tot_n_pnts), intent(out) :: shno3
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
    integer :: klevel
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
    real, allocatable, dimension(:) :: ystore
    real, dimension(tot_n_pnts,jpcspf) :: zftr
    real, dimension(tot_n_pnts,jpcspf) :: cdot
    real, dimension(tot_n_pnts) :: zq
    real, dimension(tot_n_pnts) :: co2_1d
    real, dimension(tot_n_pnts,jppj) :: zprt1d
    real, dimension(tot_n_pnts,jpdd) :: zdryrt2
    real, dimension(tot_n_pnts,2) :: rc_het
    logical, dimension(tot_n_pnts) :: stratflag
    real, dimension(tot_n_pnts,dim_ntp) :: ntp_data
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_CTL_FULL'
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
    integer :: widx1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    ix = 0
    jy = 0
    klevel = 0
    if (ukca_config%l_ukca_het_psc) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx = 1, tot_n_pnts, 1
        shno3(idx) = 0.0
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx_1 = LBOUND(stratflag, dim=1), UBOUND(stratflag, dim=1), 1
      stratflag(idx_1) = l_stratosphere(idx_1)
    enddo
    !$acc end kernels
    do l = 1, dim_ntp, 1
      if (all_ntp(l)%l_required) then
        ntp_data(:,l) = RESHAPE(all_ntp(l)%data_3d(:,:,:), [tot_n_pnts])
      end if
    enddo
    zq(:) = q(:) / c_h2o
    if (ukca_config%l_ukca_offline) then
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx_2 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
        do idx_3 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
          zprt1d(idx_3,idx_2) = 0.0
        enddo
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx_4 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
        do idx_5 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
          zprt1d(idx_5,idx_4) = photol_rates(idx_5,idx_4)
        enddo
      enddo
      !$acc end kernels
    end if
    if (ANY(speci(:) == 'CO2       ')) then
      if (ukca_config%l_chem_environ_co2_fld) then
        co2_1d(:) = co2_interactive(:) / c_co2
      else
        co2_1d(:) = rmdi
      end if
    end if
    if (ukca_config%l_ukca_trophet) then
      l = name2ntpindex('het_n2o5  ')
      !$acc kernels
      !$acc loop gang vector independent
      do idx_6 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
        rc_het(idx_6,1) = ntp_data(idx_6,l)
      enddo
      !$acc end kernels
      l = name2ntpindex('het_ho2   ')
      !$acc kernels
      !$acc loop gang vector independent
      do idx_7 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
        rc_het(idx_7,2) = ntp_data(idx_7,l)
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx_8 = LBOUND(rc_het, dim=2), UBOUND(rc_het, dim=2), 1
        do idx_9 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
          rc_het(idx_9,idx_8) = 0.0
        enddo
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_10 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
      do idx_11 = LBOUND(zdryrt2, dim=1), UBOUND(zdryrt2, dim=1), 1
        zdryrt2(idx_11,idx_10) = 0.0e0
      enddo
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_intdd) then
      !$acc kernels
      !$acc loop seq
      do k = 1, model_levels, 1
        kcs = (k - 1) * theta_field_size + 1
        kce = k * theta_field_size
        !$acc loop gang vector independent
        do idx_12 = LBOUND(ddmask, dim=1), UBOUND(ddmask, dim=1), 1
          ddmask(idx_12) = k <= nlev_with_ddep(idx_12)
        enddo
        !$acc loop gang vector independent collapse(2)
        do l = 1, jpdd, 1
          do widx1 = 1, theta_field_size, 1
            if (ddmask(widx1)) then
              zdryrt2(kcs + widx1 - 1,l) = zdryrt(widx1,l)
            end if
          enddo
        enddo
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx_13 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
        do idx_14 = 1, theta_field_size, 1
          zdryrt2(idx_14,idx_13) = zdryrt(idx_14,idx_13)
        enddo
      enddo
      !$acc end kernels
    end if
    jspf = 0
    do js = 1, jpspec, 1
      do jtr = 1, jpctr, 1
        if (advt(jtr) == speci(js)) then
          jspf = jspf + 1
          zftr(:,jspf) = tracer(:,jtr) / c_species(jtr)
        end if
      enddo
      if (ukca_config%l_ukca_ro2_ntp) then
        do jro2 = 1, jpro2, 1
          if (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') then
            l = name2ntpindex(spro2(jro2))
            jna = nlnaro2(jro2)
            if (nadvt(jna) == spro2(jro2)) then
              jspf = jspf + 1
              zftr(:,jspf) = ntp_data(:,l) / c_na_species(jna)
            else
              errcode = jro2
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl_full'
              call umprint(ummessage, src='ukca_chemistry_ctl_full')
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
        WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl_full'
        call umprint(ummessage, src='ukca_chemistry_ctl_full')
        cmessage(:) = 'ERROR: Number of chemical active species /= jpro2+jpctr'
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    else
      if (jspf /= jpctr) then
        errcode = jspf
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl_full'
        call umprint(ummessage, src='ukca_chemistry_ctl_full')
        cmessage(:) = 'ERROR: Number of chemical active species /= jpctr'
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    end if
    if (.NOT.ALLOCATED(ystore) .AND. uph2so4inaer == 1) then
      ALLOCATE(ystore(1:tot_n_pnts))
    end if
    if (ukca_config%l_ukca_het_psc) then
      sph2o(:) = qcf(:) / c_h2o
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx_15 = LBOUND(za, dim=1), UBOUND(za, dim=1), 1
      za(idx_15) = so4_sa(idx_15 + (LBOUND(so4_sa, dim=1) - LBOUND(za, dim=1)))
    enddo
    !$acc end kernels
    if (uph2so4inaer == 1) then
      if (ukca_config%l_fix_ukca_h2so4_ystore) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_16 = LBOUND(ystore, dim=1), UBOUND(ystore, dim=1), 1
          ystore(idx_16) = zftr(idx_16 + (LBOUND(zftr, dim=1) - LBOUND(ystore, dim=1)),istore_h2so4)
        enddo
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop gang vector independent
        do idx_17 = LBOUND(ystore, dim=1), UBOUND(ystore, dim=1), 1
          ystore(idx_17) = y(idx_17 + (LBOUND(y, dim=1) - LBOUND(ystore, dim=1)),nn_h2so4)
        enddo
        !$acc end kernels
      end if
    end if
    call asad_cdrive(cdot, zftr, pres, temp, zq, co2_1d, cloud_frac, qcl, ix, jy, k, zdryrt2, zwetrt, rc_het, zprt1d, tot_n_pnts, have_nat, stratflag, h_plus)
    if (ukca_config%l_ukca_het_psc) then
      if (ANY(sphno3(:) > 0.0)) then
        shno3(:) = sphno3(:) / tnd(:) * c_hono2
      else
        !$acc kernels
        !$acc loop gang vector independent
        do idx_18 = LBOUND(shno3, dim=1), UBOUND(shno3, dim=1), 1
          shno3(idx_18) = 0.0
        enddo
        !$acc end kernels
      end if
    end if
    if (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) then
      if (ihso3_h2o2 > 0) then
        delso2_wet_h2o2(:) = delso2_wet_h2o2(:) + rk(:,ihso3_h2o2) * y(:,nn_so2) * y(:,nn_h2o2) * cdt
      end if
      if (ihso3_o3 > 0) then
        delso2_wet_o3(:) = delso2_wet_o3(:) + rk(:,ihso3_o3) * y(:,nn_so2) * y(:,nn_o3) * cdt
      end if
      if (iso3_o3 > 0) then
        delso2_wet_o3(:) = delso2_wet_o3(:) + rk(:,iso3_o3) * y(:,nn_so2) * y(:,nn_o3) * cdt
      end if
      if (iso2_oh > 0 .AND. ih2so4_hv > 0) then
        delh2so4_chem(:) = delh2so4_chem(:) + (rk(:,iso2_oh) * y(:,nn_so2) * y(:,nn_oh) - rk(:,ih2so4_hv) * y(:,nn_h2so4)) * cdt
      else
        if (iso2_oh > 0) then
          delh2so4_chem(:) = delh2so4_chem(:) + rk(:,iso2_oh) * y(:,nn_so2) * y(:,nn_oh) * cdt
        end if
      end if
      if (uph2so4inaer == 1) then
        if (ukca_config%l_fix_ukca_h2so4_ystore) then
          !$acc kernels
          delh2so4_chem(:) = (zftr(:,istore_h2so4) - ystore(:)) / cdt
          !$acc loop gang vector independent
          do idx_19 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
            zftr(idx_19,istore_h2so4) = ystore(idx_19 + (LBOUND(ystore, dim=1) - LBOUND(zftr, dim=1)))
          enddo
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop gang vector independent
          do idx_20 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
            y(idx_20,nn_h2so4) = ystore(idx_20 + (LBOUND(ystore, dim=1) - LBOUND(y, dim=1)))
          enddo
          !$acc end kernels
        end if
      end if
    end if
    if (l_asad_use_chem_diags .AND. (l_asad_use_flux_rxns .OR. l_asad_use_rxn_rates .OR. (l_asad_use_wetdep .OR. l_asad_use_drydep))) then
      call asad_chemical_diagnostics(row_length, rows, model_levels, tot_n_pnts, dpd, dpw, prk, y, jy, ix, klevel, volume, ierr)
    end if
    if (l_asad_use_chem_diags .AND. l_asad_use_psc_diagnostic) then
      call asad_psc_diagnostic(row_length, rows, model_levels, tot_n_pnts, fpsc1, fpsc2, jy, ix, klevel, ierr)
    end if
    if (o1d_in_ss) then
      l = name2ntpindex('O(1D)     ')
      ntp_data(:,l) = y(:,nn_o1d) / tnd(:) * c_o1d
    end if
    if (o3p_in_ss) then
      l = name2ntpindex('O(3P)     ')
      ntp_data(:,l) = y(:,nn_o3p) / tnd(:) * c_o3p
    end if
    do jspf = 1, jpcspf, 1
      if (n_ch4 > 0) then
        if (specf(jspf) == advt(n_ch4)) then
          atm_ch4_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
      if (n_co > 0) then
        if (specf(jspf) == advt(n_co)) then
          atm_co_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
      if (n_n2o > 0) then
        if (specf(jspf) == advt(n_n2o)) then
          atm_n2o_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
      if (n_cf2cl2 > 0) then
        if (specf(jspf) == advt(n_cf2cl2)) then
          atm_cf2cl2_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
      if (n_cfcl3 > 0) then
        if (specf(jspf) == advt(n_cfcl3)) then
          atm_cfcl3_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
      if (n_mebr > 0) then
        if (specf(jspf) == advt(n_mebr)) then
          atm_mebr_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
      if (n_h2 > 0) then
        if (specf(jspf) == advt(n_h2)) then
          atm_h2_mol(:) = zftr(:,jspf) * tnd(:) * volume(:) * 1.0e6 / avogadro
        end if
      end if
    enddo
    if (ALLOCATED(ystore)) then
      DEALLOCATE(ystore)
    end if
    jspf = 0
    do js = 1, jpspec, 1
      do jtr = 1, jpctr, 1
        if (advt(jtr) == speci(js)) then
          jspf = jspf + 1
          tracer(:,jtr) = zftr(:,jspf) * c_species(jtr)
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
            ntp_data(:,l) = zftr(:,jspf) * c_na_species(jna)
          end if
        enddo
        !$acc end kernels
      end if
    enddo
    do l = 1, dim_ntp, 1
      if (all_ntp(l)%l_required) then
        all_ntp(l)%data_3d(:,:,:) = RESHAPE(ntp_data(:,l), [row_length, rows, model_levels])
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_chemistry_ctl_full

end module ukca_chemistry_ctl_full_mod
