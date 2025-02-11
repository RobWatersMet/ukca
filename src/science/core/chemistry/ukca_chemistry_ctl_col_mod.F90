module ukca_chemistry_ctl_col_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_CTL_COL_MOD'
  public

  contains
  subroutine ukca_chemistry_ctl_col(row_length, rows, model_levels, theta_field_size, ntracers, istore_h2so4, pres, temp, q, qcf, qcl, tracer, all_ntp, cloud_frac, photol_rates, shno3_3d, volume, have_nat3d, uph2so4inaer, delso2_wet_h2o2, delso2_wet_o3, delh2so4_chem, so4_sa, atm_ch4_mol, atm_co_mol, atm_n2o_mol, atm_cf2cl2_mol, atm_cfcl3_mol, atm_mebr_mol, atm_h2_mol, h_plus_3d_arr, zdryrt, zwetrt, nlev_with_ddep)
    use asad_mod, only : advt, cdt, ctype, dpd, dpw, fpsc1, fpsc2, ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3, jpcspf, jpctr, jpdd, jpdw, jpnr, jppj, jpro2, jpspec, nadvt, nlnaro2, nprkx, o1d_in_ss, o3p_in_ss, prk, rk, specf, speci, sph2o, sphno3, spro2, tnd, y, za
    use asad_chem_flux_diags, only : asad_chemical_diagnostics, asad_psc_diagnostic, l_asad_use_chem_diags, l_asad_use_drydep, l_asad_use_flux_rxns, l_asad_use_psc_diagnostic, l_asad_use_rxn_rates, l_asad_use_wetdep
    use ukca_cspecies, only : c_na_species, c_species, n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2, n_h2so4, n_mebr, n_n2o, nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh, nn_so2
    use ukca_tropopause, only : l_stratosphere
    use ukca_constants, only : c_co2, c_h2o, c_hono2, c_o1d, c_o3p
    use ukca_config_constants_mod, only : avogadro
    use ukca_config_specification_mod, only : ukca_config
    use ukca_ntp_mod, only : dim_ntp, name2ntpindex, ntp_type
    use ukca_environment_fields_mod, only : co2_interactive
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use ukca_missing_data_mod, only : rmdi
    use errormessagelength_mod, only : errormessagelength
    use asad_cdrive_mod, only : asad_cdrive
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: theta_field_size
    integer, intent(in) :: ntracers
    integer, intent(in) :: uph2so4inaer
    integer, intent(in) :: istore_h2so4
    integer, dimension(row_length,rows), intent(in) :: nlev_with_ddep
    real, dimension(row_length,rows,model_levels), intent(in) :: pres
    real, dimension(row_length,rows,model_levels), intent(in) :: temp
    real, dimension(row_length,rows,model_levels), intent(in) :: volume
    real, dimension(row_length,rows,model_levels), intent(in) :: h_plus_3d_arr
    real, dimension(row_length,rows,model_levels), intent(in) :: qcf
    real, dimension(row_length,rows,model_levels), intent(in) :: qcl
    real, dimension(row_length,rows,model_levels), intent(in) :: cloud_frac
    real, dimension(row_length,rows,model_levels), intent(in) :: so4_sa
    real, dimension(row_length,rows,jpdd), intent(in) :: zdryrt
    real, dimension(row_length,rows,model_levels,jpdw), intent(in) :: zwetrt
    real, dimension(row_length,rows,model_levels,jppj), intent(in) :: photol_rates
    real, dimension(row_length,rows,model_levels), intent(out) :: shno3_3d
    real, dimension(row_length,rows,model_levels), intent(inout) :: q
    real, dimension(row_length,rows,model_levels,ntracers), intent(inout) :: tracer
    real, dimension(row_length,rows,model_levels), intent(inout) :: delso2_wet_h2o2
    real, dimension(row_length,rows,model_levels), intent(inout) :: delso2_wet_o3
    real, dimension(row_length,rows,model_levels), intent(inout) :: delh2so4_chem
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_ch4_mol
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_co_mol
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_n2o_mol
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_cf2cl2_mol
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_cfcl3_mol
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_mebr_mol
    real, dimension(row_length,rows,model_levels), intent(inout) :: atm_h2_mol
    type(ntp_type), dimension(dim_ntp), intent(inout) :: all_ntp
    logical, dimension(row_length,rows,model_levels), intent(in) :: have_nat3d
    integer :: i
    integer :: j
    integer :: js
    integer :: jtr
    integer :: jro2
    integer :: jna
    integer :: jspf
    integer :: k
    integer :: klevel
    integer :: l
    integer :: kcs
    integer :: kce
    integer :: chunk_size
    integer :: ierr
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    CHARACTER(LEN = 10) :: prods(2)
    CHARACTER(LEN = 10) :: prods3(3)
    real, allocatable, dimension(:) :: ystore
    real, dimension(model_levels,jpcspf) :: zftr
    real, dimension(model_levels) :: zp
    real, dimension(model_levels) :: zt
    real, dimension(model_levels) :: zclw
    real, dimension(model_levels) :: zfcloud
    real, dimension(model_levels,jpcspf) :: cdot
    real, dimension(model_levels) :: zq
    real, dimension(model_levels) :: co2_1d
    real, dimension(model_levels,jppj) :: zprt1d
    real, dimension(model_levels,jpdd) :: zdryrt2
    real, dimension(model_levels,jpdw) :: zwetrt2
    real, dimension(model_levels,2) :: rc_het
    real, dimension(theta_field_size,model_levels) :: h_plus_2d_arr
    real, dimension(model_levels) :: h_plus_1d_arr
    real, dimension(model_levels,jpspec) :: dpd_full
    real, dimension(model_levels,jpspec) :: dpw_full
    real, dimension(model_levels) :: fpsc1_full
    real, dimension(model_levels) :: fpsc2_full
    real, dimension(model_levels,jpnr) :: prk_full
    real, dimension(model_levels,jpspec) :: y_full
    logical :: l_autotune_local
    logical, dimension(model_levels) :: stratflag
    logical, dimension(model_levels) :: have_nat1d
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_CTL_COL'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    klevel = 0
    if (ukca_config%l_ukca_het_psc) then
      !$acc kernels
      !$acc loop gang vector independent collapse(3)
      do idx = 1, model_levels, 1
        do idx_1 = 1, rows, 1
          do idx_2 = 1, row_length, 1
            shno3_3d(idx_2,idx_1,idx) = 0.0
          enddo
        enddo
      enddo
      !$acc end kernels
    end if
    l_autotune_local = .false.
    if (.NOT.ALLOCATED(ystore) .AND. uph2so4inaer == 1) then
      ALLOCATE(ystore(1:model_levels))
    end if
    if (l_autotune_local) then
      call ukca_reallocate_asad_arrays(ukca_config%ukca_chem_seg_size)
    end if
    do i = 1, rows, 1
      do j = 1, row_length, 1
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx_3 = 1, jpdd, 1
          do idx_4 = 1, model_levels, 1
            zdryrt2(idx_4,idx_3) = 0.0
          enddo
        enddo
        !$acc end kernels
        if (ukca_config%l_ukca_intdd) then
          !$acc kernels
          !$acc loop gang vector independent collapse(2)
          do l = 1, jpdd, 1
            do idx_5 = 1, nlev_with_ddep(j,i), 1
              zdryrt2(idx_5,l) = zdryrt(j,i,l)
            enddo
          enddo
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop gang vector independent
          do idx_6 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
            zdryrt2(1,idx_6) = zdryrt(j,i,idx_6)
          enddo
          !$acc end kernels
        end if
        !$acc kernels
        !$acc loop gang vector independent
        do idx_7 = LBOUND(zp, dim=1), UBOUND(zp, dim=1), 1
          zp(idx_7) = pres(j,i,idx_7)
        enddo
        !$acc loop gang vector independent
        do idx_8 = LBOUND(zt, dim=1), UBOUND(zt, dim=1), 1
          zt(idx_8) = temp(j,i,idx_8)
        enddo
        zq(:) = q(j,i,:) / c_h2o
        !$acc loop gang vector independent
        do idx_9 = LBOUND(h_plus_1d_arr, dim=1), UBOUND(h_plus_1d_arr, dim=1), 1
          h_plus_1d_arr(idx_9) = h_plus_3d_arr(j,i,idx_9)
        enddo
        !$acc end kernels
        if (ANY(speci(:) == 'CO2       ')) then
          if (ukca_config%l_chem_environ_co2_fld) then
            co2_1d(:) = co2_interactive(j,i,:) / c_co2
          else
            co2_1d(:) = rmdi
          end if
        end if
        !$acc kernels
        !$acc loop gang vector independent
        do idx_10 = LBOUND(zclw, dim=1), UBOUND(zclw, dim=1), 1
          zclw(idx_10) = qcl(j,i,idx_10)
        enddo
        !$acc loop gang vector independent
        do idx_11 = LBOUND(zfcloud, dim=1), UBOUND(zfcloud, dim=1), 1
          zfcloud(idx_11) = cloud_frac(j,i,idx_11)
        enddo
        jspf = 0
        !$acc end kernels
        do js = 1, jpspec, 1
          do jtr = 1, jpctr, 1
            if (advt(jtr) == speci(js)) then
              jspf = jspf + 1
              zftr(:,jspf) = tracer(j,i,:,jtr) / c_species(jtr)
            end if
          enddo
          if (ukca_config%l_ukca_ro2_ntp) then
            do jro2 = 1, jpro2, 1
              if (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') then
                l = name2ntpindex(spro2(jro2))
                jna = nlnaro2(jro2)
                if (nadvt(jna) == spro2(jro2)) then
                  jspf = jspf + 1
                  zftr(:,jspf) = all_ntp(l)%data_3d(j,i,:) / c_na_species(jna)
                else
                  errcode = jro2
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl'
                  call umprint(ummessage, src=routinename)
                  cmessage(:) = 'ERROR: Indices for RO2 species do not match w/ nadvt'
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
            call umprint(ummessage, src=routinename)
            cmessage(:) = 'ERROR: Number of chemical active species /= jpro2+jpctr'
            call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
          end if
        else
          if (jspf /= jpctr) then
            errcode = jspf
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl'
            call umprint(ummessage, src=routinename)
            cmessage(:) = 'ERROR: Number of chemical active species /= jpctr'
            call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
          end if
        end if
        if (ukca_config%l_ukca_offline) then
          !$acc kernels
          !$acc loop gang vector independent collapse(2)
          do idx_12 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
            do idx_13 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
              zprt1d(idx_13,idx_12) = 0.0
            enddo
          enddo
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop gang vector independent collapse(2)
          do idx_14 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
            do idx_15 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
              zprt1d(idx_15,idx_14) = photol_rates(j,i,idx_15,idx_14)
            enddo
          enddo
          !$acc end kernels
        end if
        if (.NOT.(ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero)) then
          if (ukca_config%l_ukca_trophet) then
            l = name2ntpindex('het_n2o5  ')
            !$acc kernels
            !$acc loop gang vector independent
            do idx_16 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
              rc_het(idx_16,1) = all_ntp(l)%data_3d(j,i,idx_16 + (LBOUND(all_ntp(l)%data_3d, dim=3) - LBOUND(rc_het, dim=1)))
            enddo
            !$acc end kernels
            l = name2ntpindex('het_ho2   ')
            !$acc kernels
            !$acc loop gang vector independent
            do idx_17 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
              rc_het(idx_17,2) = all_ntp(l)%data_3d(j,i,idx_17 + (LBOUND(all_ntp(l)%data_3d, dim=3) - LBOUND(rc_het, dim=1)))
            enddo
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop gang vector independent collapse(2)
            do idx_18 = LBOUND(rc_het, dim=2), UBOUND(rc_het, dim=2), 1
              do idx_19 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
                rc_het(idx_19,idx_18) = 0.0
              enddo
            enddo
            !$acc end kernels
          end if
          !$acc kernels
          !$acc loop gang vector independent
          do idx_20 = LBOUND(have_nat1d, dim=1), UBOUND(have_nat1d, dim=1), 1
            have_nat1d(idx_20) = have_nat3d(j,i,idx_20)
          enddo
          !$acc loop gang vector independent
          do idx_21 = LBOUND(stratflag, dim=1), UBOUND(stratflag, dim=1), 1
            stratflag(idx_21) = l_stratosphere(j,i,idx_21 + (LBOUND(l_stratosphere, dim=3) - LBOUND(stratflag, dim=1)))
          enddo
          !$acc loop gang vector independent collapse(2)
          do idx_22 = LBOUND(zwetrt2, dim=2), UBOUND(zwetrt2, dim=2), 1
            do idx_23 = LBOUND(zwetrt2, dim=1), UBOUND(zwetrt2, dim=1), 1
              zwetrt2(idx_23,idx_22) = zwetrt(j,i,idx_23,idx_22)
            enddo
          enddo
          !$acc end kernels
          do kcs = 1, model_levels, ukca_config%ukca_chem_seg_size
            kce = MIN(kcs + (ukca_config%ukca_chem_seg_size - 1), model_levels)
            chunk_size = kce + 1 - kcs
            if (chunk_size /= ukca_config%ukca_chem_seg_size) then
              call ukca_reallocate_asad_arrays(chunk_size)
            end if
            if (uph2so4inaer == 1) then
              if (ukca_config%l_fix_ukca_h2so4_ystore) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_24 = kcs, kce, 1
                  ystore(idx_24) = zftr(idx_24,istore_h2so4)
                enddo
                !$acc end kernels
              else
                !$acc kernels
                !$acc loop gang vector independent
                do idx_25 = kcs, kce, 1
                  ystore(idx_25) = y(idx_25 + (1 - kcs),nn_h2so4)
                enddo
                !$acc end kernels
              end if
            end if
            !$acc kernels
            !$acc loop gang vector independent
            do idx_26 = LBOUND(za, dim=1), UBOUND(za, dim=1), 1
              za(idx_26) = 0.0
            enddo
            !$acc loop gang vector independent
            do idx_27 = 1, chunk_size, 1
              za(idx_27) = so4_sa(j,i,idx_27 + (kcs - 1))
            enddo
            !$acc loop gang vector independent
            do idx_28 = LBOUND(sph2o, dim=1), UBOUND(sph2o, dim=1), 1
              sph2o(idx_28) = 0.0
            enddo
            !$acc end kernels
            if (ukca_config%l_ukca_het_psc) then
              sph2o(1:chunk_size) = qcf(j,i,kcs:kce) / c_h2o
            end if
            call asad_cdrive(cdot(kcs:kce,:), zftr(kcs:kce,:), zp(kcs:kce), zt(kcs:kce), zq(kcs:kce), co2_1d(kcs:kce), zfcloud(kcs:kce), zclw(kcs:kce), j, i, klevel, zdryrt2(kcs:kce,:), zwetrt2(kcs:kce,:), rc_het(kcs:kce,:), zprt1d(kcs:kce,:), chunk_size, have_nat1d(kcs:kce), stratflag(kcs:kce), h_plus_1d_arr(kcs:kce))
            !$acc kernels
            !$acc loop gang vector independent collapse(2)
            do idx_29 = LBOUND(dpd_full, dim=2), UBOUND(dpd_full, dim=2), 1
              do idx_30 = kcs, kce, 1
                dpd_full(idx_30,idx_29) = dpd(idx_30 + (1 - kcs),idx_29 + (LBOUND(dpd, dim=2) - LBOUND(dpd_full, dim=2)))
              enddo
            enddo
            !$acc loop gang vector independent collapse(2)
            do idx_31 = LBOUND(dpw_full, dim=2), UBOUND(dpw_full, dim=2), 1
              do idx_32 = kcs, kce, 1
                dpw_full(idx_32,idx_31) = dpw(idx_32 + (1 - kcs),idx_31 + (LBOUND(dpw, dim=2) - LBOUND(dpw_full, dim=2)))
              enddo
            enddo
            !$acc loop gang vector independent
            do idx_33 = kcs, kce, 1
              fpsc1_full(idx_33) = fpsc1(idx_33 + (1 - kcs))
            enddo
            !$acc loop gang vector independent
            do idx_34 = kcs, kce, 1
              fpsc2_full(idx_34) = fpsc2(idx_34 + (1 - kcs))
            enddo
            !$acc loop gang vector independent collapse(2)
            do idx_35 = LBOUND(prk_full, dim=2), UBOUND(prk_full, dim=2), 1
              do idx_36 = kcs, kce, 1
                prk_full(idx_36,idx_35) = prk(idx_36 + (1 - kcs),idx_35 + (LBOUND(prk, dim=2) - LBOUND(prk_full, dim=2)))
              enddo
            enddo
            !$acc loop gang vector independent collapse(2)
            do idx_37 = LBOUND(y_full, dim=2), UBOUND(y_full, dim=2), 1
              do idx_38 = kcs, kce, 1
                y_full(idx_38,idx_37) = y(idx_38 + (1 - kcs),idx_37 + (LBOUND(y, dim=2) - LBOUND(y_full, dim=2)))
              enddo
            enddo
            !$acc end kernels
            if (ukca_config%l_ukca_het_psc) then
              if (ANY(sphno3(:) > 0.0)) then
                shno3_3d(j,i,kcs:kce) = sphno3(:) / tnd(:) * c_hono2
              else
                !$acc kernels
                !$acc loop gang vector independent
                do idx_39 = kcs, kce, 1
                  shno3_3d(j,i,idx_39) = 0.0
                enddo
                !$acc end kernels
              end if
            end if
            if (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) then
              if (ihso3_h2o2 > 0) then
                delso2_wet_h2o2(j,i,kcs:kce) = delso2_wet_h2o2(j,i,kcs:kce) + rk(:,ihso3_h2o2) * y(:,nn_so2) * y(:,nn_h2o2) * cdt
              end if
              if (ihso3_o3 > 0) then
                delso2_wet_o3(j,i,kcs:kce) = delso2_wet_o3(j,i,kcs:kce) + rk(:,ihso3_o3) * y(:,nn_so2) * y(:,nn_o3) * cdt
              end if
              if (iso3_o3 > 0) then
                delso2_wet_o3(j,i,kcs:kce) = delso2_wet_o3(j,i,kcs:kce) + rk(:,iso3_o3) * y(:,nn_so2) * y(:,nn_o3) * cdt
              end if
              if (iso2_oh > 0 .AND. ih2so4_hv > 0) then
                delh2so4_chem(j,i,kcs:kce) = delh2so4_chem(j,i,kcs:kce) + (rk(:,iso2_oh) * y(:,nn_so2) * y(:,nn_oh) - rk(:,ih2so4_hv) * y(:,nn_h2so4)) * cdt
              else
                if (iso2_oh > 0) then
                  delh2so4_chem(j,i,kcs:kce) = delh2so4_chem(j,i,kcs:kce) + rk(:,iso2_oh) * y(:,nn_so2) * y(:,nn_oh) * cdt
                end if
              end if
              if (uph2so4inaer == 1) then
                if (ukca_config%l_fix_ukca_h2so4_ystore) then
                  !$acc kernels
                  delh2so4_chem(j,i,kcs:kce) = (zftr(kcs:kce,istore_h2so4) - ystore(kcs:kce)) / cdt
                  !$acc loop gang vector independent
                  do idx_40 = kcs, kce, 1
                    zftr(idx_40,istore_h2so4) = ystore(idx_40)
                  enddo
                  !$acc end kernels
                else
                  !$acc kernels
                  !$acc loop gang vector independent
                  do idx_41 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
                    y(idx_41,nn_h2so4) = ystore(idx_41 + (kcs - LBOUND(y, dim=1)))
                  enddo
                  !$acc end kernels
                end if
              end if
            end if
            jspf = 0
            do js = 1, jpspec, 1
              do jtr = 1, jpctr, 1
                if (advt(jtr) == speci(js)) then
                  jspf = jspf + 1
                  tracer(j,i,kcs:kce,jtr) = zftr(kcs:kce,jspf) * c_species(jtr)
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
                    all_ntp(l)%data_3d(j,i,kcs:kce) = zftr(kcs:kce,jspf) * c_na_species(jna)
                  end if
                enddo
                !$acc end kernels
              end if
            enddo
            if (o1d_in_ss) then
              l = name2ntpindex('O(1D)     ')
              all_ntp(l)%data_3d(j,i,kcs:kce) = y(:,nn_o1d) / tnd(:) * c_o1d
            end if
            if (o3p_in_ss) then
              l = name2ntpindex('O(3P)     ')
              all_ntp(l)%data_3d(j,i,kcs:kce) = y(:,nn_o3p) / tnd(:) * c_o3p
            end if
            do jspf = 1, jpcspf, 1
              if (n_ch4 > 0) then
                if (specf(jspf) == advt(n_ch4)) then
                  atm_ch4_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
              if (n_co > 0) then
                if (specf(jspf) == advt(n_co)) then
                  atm_co_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
              if (n_n2o > 0) then
                if (specf(jspf) == advt(n_n2o)) then
                  atm_n2o_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
              if (n_cf2cl2 > 0) then
                if (specf(jspf) == advt(n_cf2cl2)) then
                  atm_cf2cl2_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
              if (n_cfcl3 > 0) then
                if (specf(jspf) == advt(n_cfcl3)) then
                  atm_cfcl3_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
              if (n_mebr > 0) then
                if (specf(jspf) == advt(n_mebr)) then
                  atm_mebr_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
              if (n_h2 > 0) then
                if (specf(jspf) == advt(n_h2)) then
                  atm_h2_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf) * tnd(:) * volume(j,i,kcs:kce) * 1.0e6 / avogadro
                end if
              end if
            enddo
          enddo
          if (chunk_size /= ukca_config%ukca_chem_seg_size) then
            call ukca_reallocate_asad_arrays(ukca_config%ukca_chem_seg_size)
          end if
          if (l_asad_use_chem_diags .AND. (l_asad_use_flux_rxns .OR. l_asad_use_rxn_rates .OR. (l_asad_use_wetdep .OR. l_asad_use_drydep))) then
            call asad_chemical_diagnostics(row_length, rows, model_levels, model_levels, dpd_full, dpw_full, prk_full, y_full, j, i, klevel, volume, ierr)
          end if
          if (l_asad_use_chem_diags .AND. l_asad_use_psc_diagnostic) then
            call asad_psc_diagnostic(row_length, rows, model_levels, model_levels, fpsc1_full, fpsc2_full, j, i, klevel, ierr)
          end if
        else
          cmessage(:) = 'Column call is not available for Backward Euler schemes'
          errcode = 5
          call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
        end if
      enddo
    enddo
    if (ALLOCATED(ystore)) then
      DEALLOCATE(ystore)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_chemistry_ctl_col
  subroutine ukca_reallocate_asad_arrays(n_pnts)
    use asad_mod, only : co2, co3, deriv, dpd, dpw, ej, emr, f, fdot, fj, fpsc1, fpsc2, ftilde, ipa, jpcspf, jpnr, jpspec, lati, linfam, method, p, pd, pmintnd, prk, prod, qa, ratio, rk, sh2o, shno3, slos, spfj, spfjsize_max, sph2o, sphno3, t, t300, tnd, wp, y, ydot, za
    use ukca_config_specification_mod, only : int_method_nr, ukca_config
    integer, intent(in) :: n_pnts
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

    method = ukca_config%ukca_int_method
    if (method == int_method_nr) then
      if (ALLOCATED(spfj)) then
        DEALLOCATE(spfj)
      end if
    end if
    if (ALLOCATED(za)) then
      DEALLOCATE(za)
    end if
    if (ALLOCATED(ydot)) then
      DEALLOCATE(ydot)
    end if
    if (ALLOCATED(y)) then
      DEALLOCATE(y)
    end if
    if (ALLOCATED(co2)) then
      DEALLOCATE(co2)
    end if
    if (ALLOCATED(wp)) then
      DEALLOCATE(wp)
    end if
    if (ALLOCATED(tnd)) then
      DEALLOCATE(tnd)
    end if
    if (ALLOCATED(t300)) then
      DEALLOCATE(t300)
    end if
    if (ALLOCATED(t)) then
      DEALLOCATE(t)
    end if
    if (ALLOCATED(sphno3)) then
      DEALLOCATE(sphno3)
    end if
    if (ALLOCATED(sph2o)) then
      DEALLOCATE(sph2o)
    end if
    if (ALLOCATED(shno3)) then
      DEALLOCATE(shno3)
    end if
    if (ALLOCATED(sh2o)) then
      DEALLOCATE(sh2o)
    end if
    if (ALLOCATED(rk)) then
      DEALLOCATE(rk)
    end if
    if (ALLOCATED(ratio)) then
      DEALLOCATE(ratio)
    end if
    if (ALLOCATED(qa)) then
      DEALLOCATE(qa)
    end if
    if (ALLOCATED(prk)) then
      DEALLOCATE(prk)
    end if
    if (ALLOCATED(pmintnd)) then
      DEALLOCATE(pmintnd)
    end if
    if (ALLOCATED(p)) then
      DEALLOCATE(p)
    end if
    if (ALLOCATED(linfam)) then
      DEALLOCATE(linfam)
    end if
    if (ALLOCATED(lati)) then
      DEALLOCATE(lati)
    end if
    if (ALLOCATED(ipa)) then
      DEALLOCATE(ipa)
    end if
    if (ALLOCATED(ftilde)) then
      DEALLOCATE(ftilde)
    end if
    if (ALLOCATED(fpsc2)) then
      DEALLOCATE(fpsc2)
    end if
    if (ALLOCATED(fpsc1)) then
      DEALLOCATE(fpsc1)
    end if
    if (ALLOCATED(fj)) then
      DEALLOCATE(fj)
    end if
    if (ALLOCATED(fdot)) then
      DEALLOCATE(fdot)
    end if
    if (ALLOCATED(f)) then
      DEALLOCATE(f)
    end if
    if (ALLOCATED(emr)) then
      DEALLOCATE(emr)
    end if
    if (ALLOCATED(ej)) then
      DEALLOCATE(ej)
    end if
    if (ALLOCATED(dpw)) then
      DEALLOCATE(dpw)
    end if
    if (ALLOCATED(dpd)) then
      DEALLOCATE(dpd)
    end if
    if (ALLOCATED(deriv)) then
      DEALLOCATE(deriv)
    end if
    if (ALLOCATED(co3)) then
      DEALLOCATE(co3)
    end if
    if (ALLOCATED(pd)) then
      DEALLOCATE(pd)
    end if
    if (.NOT.ALLOCATED(pd)) then
      ALLOCATE(pd(1:n_pnts,1:2 * jpspec))
    end if
    if (.NOT.ALLOCATED(co3)) then
      ALLOCATE(co3(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(deriv)) then
      ALLOCATE(deriv(1:n_pnts,1:4,1:4))
    end if
    if (.NOT.ALLOCATED(dpd)) then
      ALLOCATE(dpd(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(dpw)) then
      ALLOCATE(dpw(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ej)) then
      ALLOCATE(ej(1:n_pnts,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(emr)) then
      ALLOCATE(emr(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(f)) then
      ALLOCATE(f(1:n_pnts,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fdot)) then
      ALLOCATE(fdot(1:n_pnts,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fj)) then
      ALLOCATE(fj(1:n_pnts,1:jpcspf,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fpsc1)) then
      ALLOCATE(fpsc1(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(fpsc2)) then
      ALLOCATE(fpsc2(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(ftilde)) then
      ALLOCATE(ftilde(1:n_pnts,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(ipa)) then
      ALLOCATE(ipa(1:n_pnts,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(lati)) then
      ALLOCATE(lati(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(linfam)) then
      ALLOCATE(linfam(1:n_pnts,0:jpcspf))
    end if
    if (.NOT.ALLOCATED(p)) then
      ALLOCATE(p(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(pmintnd)) then
      ALLOCATE(pmintnd(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(prk)) then
      ALLOCATE(prk(1:n_pnts,1:jpnr))
    end if
    if (.NOT.ALLOCATED(qa)) then
      ALLOCATE(qa(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ratio)) then
      ALLOCATE(ratio(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(rk)) then
      ALLOCATE(rk(1:n_pnts,1:jpnr))
    end if
    if (.NOT.ALLOCATED(sh2o)) then
      ALLOCATE(sh2o(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(shno3)) then
      ALLOCATE(shno3(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(sph2o)) then
      ALLOCATE(sph2o(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(sphno3)) then
      ALLOCATE(sphno3(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(t)) then
      ALLOCATE(t(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(t300)) then
      ALLOCATE(t300(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(tnd)) then
      ALLOCATE(tnd(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(wp)) then
      ALLOCATE(wp(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(co2)) then
      ALLOCATE(co2(1:n_pnts))
    end if
    if (.NOT.ALLOCATED(y)) then
      ALLOCATE(y(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ydot)) then
      ALLOCATE(ydot(1:n_pnts,1:jpspec))
    end if
    if (.NOT.ALLOCATED(za)) then
      ALLOCATE(za(1:n_pnts))
    end if
    if (method == int_method_nr) then
      if (.NOT.ALLOCATED(spfj)) then
        ALLOCATE(spfj(1:n_pnts,1:spfjsize_max))
      end if
    end if
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - Unsupported statement: Nullify_Stmt
    !  - Unsupported statement: Nullify_Stmt
    NULLIFY(prod)
    NULLIFY(slos)
    !$acc kernels
    prod => pd(:,1:jpspec)
    slos => pd(:,jpspec + 1:2 * jpspec)
    !$acc loop gang vector independent collapse(3)
    do idx = LBOUND(deriv, dim=3), UBOUND(deriv, dim=3), 1
      do idx_1 = LBOUND(deriv, dim=2), UBOUND(deriv, dim=2), 1
        do idx_2 = LBOUND(deriv, dim=1), UBOUND(deriv, dim=1), 1
          deriv(idx_2,idx_1,idx) = 1.0
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_3 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
      do idx_4 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_4,idx_3) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_5 = LBOUND(fdot, dim=2), UBOUND(fdot, dim=2), 1
      do idx_6 = LBOUND(fdot, dim=1), UBOUND(fdot, dim=1), 1
        fdot(idx_6,idx_5) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_7 = LBOUND(ej, dim=2), UBOUND(ej, dim=2), 1
      do idx_8 = LBOUND(ej, dim=1), UBOUND(ej, dim=1), 1
        ej(idx_8,idx_7) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_9 = LBOUND(linfam, dim=2), UBOUND(linfam, dim=2), 1
      do idx_10 = LBOUND(linfam, dim=1), UBOUND(linfam, dim=1), 1
        linfam(idx_10,idx_9) = .false.
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_11 = LBOUND(y, dim=2), UBOUND(y, dim=2), 1
      do idx_12 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_12,idx_11) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_13 = LBOUND(ydot, dim=2), UBOUND(ydot, dim=2), 1
      do idx_14 = LBOUND(ydot, dim=1), UBOUND(ydot, dim=1), 1
        ydot(idx_14,idx_13) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_15 = LBOUND(prod, dim=2), UBOUND(prod, dim=2), 1
      do idx_16 = LBOUND(prod, dim=1), UBOUND(prod, dim=1), 1
        prod(idx_16,idx_15) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_17 = LBOUND(slos, dim=2), UBOUND(slos, dim=2), 1
      do idx_18 = LBOUND(slos, dim=1), UBOUND(slos, dim=1), 1
        slos(idx_18,idx_17) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_19 = LBOUND(dpd, dim=2), UBOUND(dpd, dim=2), 1
      do idx_20 = LBOUND(dpd, dim=1), UBOUND(dpd, dim=1), 1
        dpd(idx_20,idx_19) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_21 = LBOUND(dpw, dim=2), UBOUND(dpw, dim=2), 1
      do idx_22 = LBOUND(dpw, dim=1), UBOUND(dpw, dim=1), 1
        dpw(idx_22,idx_21) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_23 = LBOUND(emr, dim=2), UBOUND(emr, dim=2), 1
      do idx_24 = LBOUND(emr, dim=1), UBOUND(emr, dim=1), 1
        emr(idx_24,idx_23) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_25 = LBOUND(rk, dim=2), UBOUND(rk, dim=2), 1
      do idx_26 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
        rk(idx_26,idx_25) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_27 = LBOUND(prk, dim=2), UBOUND(prk, dim=2), 1
      do idx_28 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
        prk(idx_28,idx_27) = 0.0
      enddo
    enddo
    !$acc end kernels
    return

  end subroutine ukca_reallocate_asad_arrays

end module ukca_chemistry_ctl_col_mod
