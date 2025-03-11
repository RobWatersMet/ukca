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
    use nvtx
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
    integer :: idx_42
    integer :: idx_43
    integer :: idx_44
    integer :: idx_45
    integer :: idx_46
    integer :: idx_47
    integer :: idx_48
    integer :: idx_49
    integer :: idx_50
    integer :: idx_51
    integer :: idx_52
    integer :: idx_53
    integer :: idx_54
    integer :: idx_55
    integer :: idx_56
    integer :: idx_57
    integer :: idx_58
    integer :: idx_59
    integer :: idx_60
    integer :: idx_61
    integer :: idx_62
    integer :: idx_63
    integer :: idx_64
    integer :: idx_65

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    klevel = 0
    if (ukca_config%l_ukca_het_psc) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(3)
      do idx = 1, model_levels, 1
        do idx_1 = 1, rows, 1
          do idx_2 = 1, row_length, 1
            shno3_3d(idx_2,idx_1,idx) = 0.0
          enddo
        enddo
      enddo
      !NO!$acc end kernels
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
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(2)
        do idx_3 = 1, jpdd, 1
          do idx_4 = 1, model_levels, 1
            zdryrt2(idx_4,idx_3) = 0.0
          enddo
        enddo
        !NO!$acc end kernels
        if (ukca_config%l_ukca_intdd) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do l = 1, jpdd, 1
            do idx_5 = 1, nlev_with_ddep(j,i), 1
              zdryrt2(idx_5,l) = zdryrt(j,i,l)
            enddo
          enddo
          !NO!$acc end kernels
        else
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_6 = LBOUND(zdryrt2, dim=2), UBOUND(zdryrt2, dim=2), 1
            zdryrt2(1,idx_6) = zdryrt(j,i,idx_6)
          enddo
          !NO!$acc end kernels
        end if
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_7 = LBOUND(zp, dim=1), UBOUND(zp, dim=1), 1
          zp(idx_7) = pres(j,i,idx_7)
        enddo
        !NO!$acc loop gang vector independent
        do idx_8 = LBOUND(zt, dim=1), UBOUND(zt, dim=1), 1
          zt(idx_8) = temp(j,i,idx_8)
        enddo
        !NO!$acc loop gang vector independent
        do idx_9 = LBOUND(zq, dim=1), UBOUND(zq, dim=1), 1
          zq(idx_9) = q(j,i,idx_9) / c_h2o
        enddo
        !NO!$acc loop gang vector independent
        do idx_10 = LBOUND(h_plus_1d_arr, dim=1), UBOUND(h_plus_1d_arr, dim=1), 1
          h_plus_1d_arr(idx_10) = h_plus_3d_arr(j,i,idx_10)
        enddo
        !NO!$acc end kernels
        if (ANY(speci(:) == 'CO2       ')) then
          if (ukca_config%l_chem_environ_co2_fld) then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_11 = LBOUND(co2_1d, dim=1), UBOUND(co2_1d, dim=1), 1
              co2_1d(idx_11) = co2_interactive(j,i,idx_11 + (LBOUND(co2_interactive, dim=3) - LBOUND(co2_1d, dim=1))) / c_co2
            enddo
            !NO!$acc end kernels
          else
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_12 = LBOUND(co2_1d, dim=1), UBOUND(co2_1d, dim=1), 1
              co2_1d(idx_12) = rmdi
            enddo
            !NO!$acc end kernels
          end if
        end if
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_13 = LBOUND(zclw, dim=1), UBOUND(zclw, dim=1), 1
          zclw(idx_13) = qcl(j,i,idx_13)
        enddo
        !NO!$acc loop gang vector independent
        do idx_14 = LBOUND(zfcloud, dim=1), UBOUND(zfcloud, dim=1), 1
          zfcloud(idx_14) = cloud_frac(j,i,idx_14)
        enddo
        jspf = 0
        !NO!$acc end kernels
        do js = 1, jpspec, 1
          do jtr = 1, jpctr, 1
            if (advt(jtr) == speci(js)) then
              !NO!$acc kernels
              jspf = jspf + 1
              !NO!$acc loop gang vector independent
              do idx_15 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
                zftr(idx_15,jspf) = tracer(j,i,idx_15,jtr) / c_species(jtr)
              enddo
              !NO!$acc end kernels
            end if
          enddo
          if (ukca_config%l_ukca_ro2_ntp) then
            do jro2 = 1, jpro2, 1
              if (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') then
                l = name2ntpindex(spro2(jro2))
                jna = nlnaro2(jro2)
                if (nadvt(jna) == spro2(jro2)) then
                  !NO!$acc kernels
                  jspf = jspf + 1
                  !NO!$acc loop gang vector independent
                  do idx_16 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
                    zftr(idx_16,jspf) = all_ntp(l)%data_3d(j,i,idx_16 + (LBOUND(all_ntp(l)%data_3d, dim=3) - LBOUND(zftr, dim=1))) / c_na_species(jna)
                  enddo
                  !NO!$acc end kernels
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
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do idx_17 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
            do idx_18 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
              zprt1d(idx_18,idx_17) = 0.0
            enddo
          enddo
          !NO!$acc end kernels
        else
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do idx_19 = LBOUND(zprt1d, dim=2), UBOUND(zprt1d, dim=2), 1
            do idx_20 = LBOUND(zprt1d, dim=1), UBOUND(zprt1d, dim=1), 1
              zprt1d(idx_20,idx_19) = photol_rates(j,i,idx_20,idx_19)
            enddo
          enddo
          !NO!$acc end kernels
        end if
        if (.NOT.(ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero)) then
          if (ukca_config%l_ukca_trophet) then
            l = name2ntpindex('het_n2o5  ')
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_21 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
              rc_het(idx_21,1) = all_ntp(l)%data_3d(j,i,idx_21 + (LBOUND(all_ntp(l)%data_3d, dim=3) - LBOUND(rc_het, dim=1)))
            enddo
            !NO!$acc end kernels
            l = name2ntpindex('het_ho2   ')
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_22 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
              rc_het(idx_22,2) = all_ntp(l)%data_3d(j,i,idx_22 + (LBOUND(all_ntp(l)%data_3d, dim=3) - LBOUND(rc_het, dim=1)))
            enddo
            !NO!$acc end kernels
          else
            !NO!$acc kernels
            !NO!$acc loop gang vector independent collapse(2)
            do idx_23 = LBOUND(rc_het, dim=2), UBOUND(rc_het, dim=2), 1
              do idx_24 = LBOUND(rc_het, dim=1), UBOUND(rc_het, dim=1), 1
                rc_het(idx_24,idx_23) = 0.0
              enddo
            enddo
            !NO!$acc end kernels
          end if
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_25 = LBOUND(have_nat1d, dim=1), UBOUND(have_nat1d, dim=1), 1
            have_nat1d(idx_25) = have_nat3d(j,i,idx_25)
          enddo
          !NO!$acc loop gang vector independent
          do idx_26 = LBOUND(stratflag, dim=1), UBOUND(stratflag, dim=1), 1
            stratflag(idx_26) = l_stratosphere(j,i,idx_26 + (LBOUND(l_stratosphere, dim=3) - LBOUND(stratflag, dim=1)))
          enddo
          !NO!$acc loop gang vector independent collapse(2)
          do idx_27 = LBOUND(zwetrt2, dim=2), UBOUND(zwetrt2, dim=2), 1
            do idx_28 = LBOUND(zwetrt2, dim=1), UBOUND(zwetrt2, dim=1), 1
              zwetrt2(idx_28,idx_27) = zwetrt(j,i,idx_28,idx_27)
            enddo
          enddo
          !NO!$acc end kernels
          do kcs = 1, model_levels, ukca_config%ukca_chem_seg_size
            kce = MIN(kcs + (ukca_config%ukca_chem_seg_size - 1), model_levels)
            chunk_size = (kce + 1) - kcs
            if (chunk_size /= ukca_config%ukca_chem_seg_size) then
              call ukca_reallocate_asad_arrays(chunk_size)
            end if
            if (uph2so4inaer == 1) then
              if (ukca_config%l_fix_ukca_h2so4_ystore) then
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_29 = kcs, kce, 1
                  ystore(idx_29) = zftr(idx_29,istore_h2so4)
                enddo
                !NO!$acc end kernels
              else
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_30 = kcs, kce, 1
                  ystore(idx_30) = y(idx_30 + (1 - kcs),nn_h2so4)
                enddo
                !NO!$acc end kernels
              end if
            end if
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_31 = LBOUND(za, dim=1), UBOUND(za, dim=1), 1
              za(idx_31) = 0.0
            enddo
            !NO!$acc loop gang vector independent
            do idx_32 = 1, chunk_size, 1
              za(idx_32) = so4_sa(j,i,idx_32 + (kcs - 1))
            enddo
            !NO!$acc loop gang vector independent
            do idx_33 = LBOUND(sph2o, dim=1), UBOUND(sph2o, dim=1), 1
              sph2o(idx_33) = 0.0
            enddo
            !NO!$acc end kernels
            if (ukca_config%l_ukca_het_psc) then
              !NO!$acc kernels
              !NO!$acc loop gang vector independent
              do idx_34 = 1, chunk_size, 1
                sph2o(idx_34) = qcf(j,i,idx_34 + (kcs - 1)) / c_h2o
              enddo
              !NO!$acc end kernels
            end if
            call asad_cdrive(cdot(kcs:kce,:), zftr(kcs:kce,:), zp(kcs:kce), zt(kcs:kce), zq(kcs:kce), co2_1d(kcs:kce), zfcloud(kcs:kce), zclw(kcs:kce), j, i, klevel, zdryrt2(kcs:kce,:), zwetrt2(kcs:kce,:), rc_het(kcs:kce,:), zprt1d(kcs:kce,:), chunk_size, have_nat1d(kcs:kce), stratflag(kcs:kce), h_plus_1d_arr(kcs:kce))
            !NO!$acc kernels
            !NO!$acc loop gang vector independent collapse(2)
            do idx_35 = LBOUND(dpd_full, dim=2), UBOUND(dpd_full, dim=2), 1
              do idx_36 = kcs, kce, 1
                dpd_full(idx_36,idx_35) = dpd(idx_36 + (1 - kcs),idx_35 + (LBOUND(dpd, dim=2) - LBOUND(dpd_full, dim=2)))
              enddo
            enddo
            !NO!$acc loop gang vector independent collapse(2)
            do idx_37 = LBOUND(dpw_full, dim=2), UBOUND(dpw_full, dim=2), 1
              do idx_38 = kcs, kce, 1
                dpw_full(idx_38,idx_37) = dpw(idx_38 + (1 - kcs),idx_37 + (LBOUND(dpw, dim=2) - LBOUND(dpw_full, dim=2)))
              enddo
            enddo
            !NO!$acc loop gang vector independent
            do idx_39 = kcs, kce, 1
              fpsc1_full(idx_39) = fpsc1(idx_39 + (1 - kcs))
            enddo
            !NO!$acc loop gang vector independent
            do idx_40 = kcs, kce, 1
              fpsc2_full(idx_40) = fpsc2(idx_40 + (1 - kcs))
            enddo
            !NO!$acc loop gang vector independent collapse(2)
            do idx_41 = LBOUND(prk_full, dim=2), UBOUND(prk_full, dim=2), 1
              do idx_42 = kcs, kce, 1
                prk_full(idx_42,idx_41) = prk(idx_42 + (1 - kcs),idx_41 + (LBOUND(prk, dim=2) - LBOUND(prk_full, dim=2)))
              enddo
            enddo
            !NO!$acc loop gang vector independent collapse(2)
            do idx_43 = LBOUND(y_full, dim=2), UBOUND(y_full, dim=2), 1
              do idx_44 = kcs, kce, 1
                y_full(idx_44,idx_43) = y(idx_44 + (1 - kcs),idx_43 + (LBOUND(y, dim=2) - LBOUND(y_full, dim=2)))
              enddo
            enddo
            !NO!$acc end kernels
            if (ukca_config%l_ukca_het_psc) then
              if (ANY(sphno3(:) > 0.0)) then
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_45 = kcs, kce, 1
                  shno3_3d(j,i,idx_45) = (sphno3(idx_45 + (LBOUND(sphno3, dim=1) - kcs)) / tnd(idx_45 + (LBOUND(tnd, dim=1) - kcs))) * c_hono2
                enddo
                !NO!$acc end kernels
              else
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_46 = kcs, kce, 1
                  shno3_3d(j,i,idx_46) = 0.0
                enddo
                !NO!$acc end kernels
              end if
            end if
            if (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) then
              if (ihso3_h2o2 > 0) then
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_47 = kcs, kce, 1
                  delso2_wet_h2o2(j,i,idx_47) = delso2_wet_h2o2(j,i,idx_47) + (rk(idx_47 + (LBOUND(rk, dim=1) - kcs),ihso3_h2o2) * y(idx_47 + (LBOUND(y, dim=1) - kcs),nn_so2) * y(idx_47 + (LBOUND(y, dim=1) - kcs),nn_h2o2)) * cdt
                enddo
                !NO!$acc end kernels
              end if
              if (ihso3_o3 > 0) then
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_48 = kcs, kce, 1
                  delso2_wet_o3(j,i,idx_48) = delso2_wet_o3(j,i,idx_48) + (rk(idx_48 + (LBOUND(rk, dim=1) - kcs),ihso3_o3) * y(idx_48 + (LBOUND(y, dim=1) - kcs),nn_so2) * y(idx_48 + (LBOUND(y, dim=1) - kcs),nn_o3)) * cdt
                enddo
                !NO!$acc end kernels
              end if
              if (iso3_o3 > 0) then
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_49 = kcs, kce, 1
                  delso2_wet_o3(j,i,idx_49) = delso2_wet_o3(j,i,idx_49) + (rk(idx_49 + (LBOUND(rk, dim=1) - kcs),iso3_o3) * y(idx_49 + (LBOUND(y, dim=1) - kcs),nn_so2) * y(idx_49 + (LBOUND(y, dim=1) - kcs),nn_o3)) * cdt
                enddo
                !NO!$acc end kernels
              end if
              if (iso2_oh > 0 .AND. ih2so4_hv > 0) then
                !NO!$acc kernels
                !NO!$acc loop gang vector independent
                do idx_50 = kcs, kce, 1
                  delh2so4_chem(j,i,idx_50) = delh2so4_chem(j,i,idx_50) + (rk(idx_50 + (LBOUND(rk, dim=1) - kcs),iso2_oh) * y(idx_50 + (LBOUND(y, dim=1) - kcs),nn_so2) * y(idx_50 + (LBOUND(y, dim=1) - kcs),nn_oh) - rk(idx_50 + (LBOUND(rk, dim=1) - kcs),ih2so4_hv) * y(idx_50 + (LBOUND(y, dim=1) - kcs),nn_h2so4)) * cdt
                enddo
                !NO!$acc end kernels
              else
                if (iso2_oh > 0) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_51 = kcs, kce, 1
                    delh2so4_chem(j,i,idx_51) = delh2so4_chem(j,i,idx_51) + (rk(idx_51 + (LBOUND(rk, dim=1) - kcs),iso2_oh) * y(idx_51 + (LBOUND(y, dim=1) - kcs),nn_so2) * y(idx_51 + (LBOUND(y, dim=1) - kcs),nn_oh)) * cdt
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (uph2so4inaer == 1) then
                if (ukca_config%l_fix_ukca_h2so4_ystore) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_52 = kcs, kce, 1
                    delh2so4_chem(j,i,idx_52) = (zftr(idx_52,istore_h2so4) - ystore(idx_52)) / cdt
                  enddo
                  !NO!$acc loop gang vector independent
                  do idx_53 = kcs, kce, 1
                    zftr(idx_53,istore_h2so4) = ystore(idx_53)
                  enddo
                  !NO!$acc end kernels
                else
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_54 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
                    y(idx_54,nn_h2so4) = ystore(idx_54 + (kcs - LBOUND(y, dim=1)))
                  enddo
                  !NO!$acc end kernels
                end if
              end if
            end if
            jspf = 0
            do js = 1, jpspec, 1
              do jtr = 1, jpctr, 1
                if (advt(jtr) == speci(js)) then
                  !NO!$acc kernels
                  jspf = jspf + 1
                  !NO!$acc loop gang vector independent
                  do idx_55 = kcs, kce, 1
                    tracer(j,i,idx_55,jtr) = zftr(idx_55,jspf) * c_species(jtr)
                  enddo
                  !NO!$acc end kernels
                end if
              enddo
              if (ukca_config%l_ukca_ro2_ntp) then
                do jro2 = 1, jpro2, 1
                  if (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') then
                    jspf = jspf + 1
                    l = name2ntpindex(spro2(jro2))
                    !NO!$acc kernels
                    jna = nlnaro2(jro2)
                    !NO!$acc loop gang vector independent
                    do idx_56 = kcs, kce, 1
                      all_ntp(l)%data_3d(j,i,idx_56) = zftr(idx_56,jspf) * c_na_species(jna)
                    enddo
                    !NO!$acc end kernels
                  end if
                enddo
              end if
            enddo
            if (o1d_in_ss) then
              l = name2ntpindex('O(1D)     ')
              !NO!$acc kernels
              !NO!$acc loop gang vector independent
              do idx_57 = kcs, kce, 1
                all_ntp(l)%data_3d(j,i,idx_57) = (y(idx_57 + (LBOUND(y, dim=1) - kcs),nn_o1d) / tnd(idx_57 + (LBOUND(tnd, dim=1) - kcs))) * c_o1d
              enddo
              !NO!$acc end kernels
            end if
            if (o3p_in_ss) then
              l = name2ntpindex('O(3P)     ')
              !NO!$acc kernels
              !NO!$acc loop gang vector independent
              do idx_58 = kcs, kce, 1
                all_ntp(l)%data_3d(j,i,idx_58) = (y(idx_58 + (LBOUND(y, dim=1) - kcs),nn_o3p) / tnd(idx_58 + (LBOUND(tnd, dim=1) - kcs))) * c_o3p
              enddo
              !NO!$acc end kernels
            end if
            do jspf = 1, jpcspf, 1
              if (n_ch4 > 0) then
                if (specf(jspf) == advt(n_ch4)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_59 = kcs, kce, 1
                    atm_ch4_mol(j,i,idx_59) = zftr(idx_59,jspf) * tnd(idx_59 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_59) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (n_co > 0) then
                if (specf(jspf) == advt(n_co)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_60 = kcs, kce, 1
                    atm_co_mol(j,i,idx_60) = zftr(idx_60,jspf) * tnd(idx_60 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_60) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (n_n2o > 0) then
                if (specf(jspf) == advt(n_n2o)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_61 = kcs, kce, 1
                    atm_n2o_mol(j,i,idx_61) = zftr(idx_61,jspf) * tnd(idx_61 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_61) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (n_cf2cl2 > 0) then
                if (specf(jspf) == advt(n_cf2cl2)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_62 = kcs, kce, 1
                    atm_cf2cl2_mol(j,i,idx_62) = zftr(idx_62,jspf) * tnd(idx_62 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_62) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (n_cfcl3 > 0) then
                if (specf(jspf) == advt(n_cfcl3)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_63 = kcs, kce, 1
                    atm_cfcl3_mol(j,i,idx_63) = zftr(idx_63,jspf) * tnd(idx_63 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_63) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (n_mebr > 0) then
                if (specf(jspf) == advt(n_mebr)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_64 = kcs, kce, 1
                    atm_mebr_mol(j,i,idx_64) = zftr(idx_64,jspf) * tnd(idx_64 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_64) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
              if (n_h2 > 0) then
                if (specf(jspf) == advt(n_h2)) then
                  !NO!$acc kernels
                  !NO!$acc loop gang vector independent
                  do idx_65 = kcs, kce, 1
                    atm_h2_mol(j,i,idx_65) = zftr(idx_65,jspf) * tnd(idx_65 + (LBOUND(tnd, dim=1) - kcs)) * volume(j,i,idx_65) * 1.0e6 / avogadro
                  enddo
                  !NO!$acc end kernels
                end if
              end if
            enddo
          enddo
          if (chunk_size /= ukca_config%ukca_chem_seg_size) then
            call ukca_reallocate_asad_arrays(ukca_config%ukca_chem_seg_size)
          end if
          if (l_asad_use_chem_diags .AND. ((l_asad_use_flux_rxns .OR. l_asad_use_rxn_rates) .OR. (l_asad_use_wetdep .OR. l_asad_use_drydep))) then
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
    call nvtxEndRange()
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
    !NO!$acc kernels
    prod => pd(:,1:jpspec)
    slos => pd(:,jpspec + 1:2 * jpspec)
    !NO!$acc loop gang vector independent collapse(3)
    do idx = LBOUND(deriv, dim=3), UBOUND(deriv, dim=3), 1
      do idx_1 = LBOUND(deriv, dim=2), UBOUND(deriv, dim=2), 1
        do idx_2 = LBOUND(deriv, dim=1), UBOUND(deriv, dim=1), 1
          deriv(idx_2,idx_1,idx) = 1.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_3 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
      do idx_4 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_4,idx_3) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_5 = LBOUND(fdot, dim=2), UBOUND(fdot, dim=2), 1
      do idx_6 = LBOUND(fdot, dim=1), UBOUND(fdot, dim=1), 1
        fdot(idx_6,idx_5) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_7 = LBOUND(ej, dim=2), UBOUND(ej, dim=2), 1
      do idx_8 = LBOUND(ej, dim=1), UBOUND(ej, dim=1), 1
        ej(idx_8,idx_7) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_9 = LBOUND(linfam, dim=2), UBOUND(linfam, dim=2), 1
      do idx_10 = LBOUND(linfam, dim=1), UBOUND(linfam, dim=1), 1
        linfam(idx_10,idx_9) = .false.
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_11 = LBOUND(y, dim=2), UBOUND(y, dim=2), 1
      do idx_12 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_12,idx_11) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_13 = LBOUND(ydot, dim=2), UBOUND(ydot, dim=2), 1
      do idx_14 = LBOUND(ydot, dim=1), UBOUND(ydot, dim=1), 1
        ydot(idx_14,idx_13) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_15 = LBOUND(prod, dim=2), UBOUND(prod, dim=2), 1
      do idx_16 = LBOUND(prod, dim=1), UBOUND(prod, dim=1), 1
        prod(idx_16,idx_15) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_17 = LBOUND(slos, dim=2), UBOUND(slos, dim=2), 1
      do idx_18 = LBOUND(slos, dim=1), UBOUND(slos, dim=1), 1
        slos(idx_18,idx_17) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_19 = LBOUND(dpd, dim=2), UBOUND(dpd, dim=2), 1
      do idx_20 = LBOUND(dpd, dim=1), UBOUND(dpd, dim=1), 1
        dpd(idx_20,idx_19) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_21 = LBOUND(dpw, dim=2), UBOUND(dpw, dim=2), 1
      do idx_22 = LBOUND(dpw, dim=1), UBOUND(dpw, dim=1), 1
        dpw(idx_22,idx_21) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_23 = LBOUND(emr, dim=2), UBOUND(emr, dim=2), 1
      do idx_24 = LBOUND(emr, dim=1), UBOUND(emr, dim=1), 1
        emr(idx_24,idx_23) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_25 = LBOUND(rk, dim=2), UBOUND(rk, dim=2), 1
      do idx_26 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
        rk(idx_26,idx_25) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_27 = LBOUND(prk, dim=2), UBOUND(prk, dim=2), 1
      do idx_28 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
        prk(idx_28,idx_27) = 0.0
      enddo
    enddo
    !NO!$acc end kernels
    return

  end subroutine ukca_reallocate_asad_arrays

end module ukca_chemistry_ctl_col_mod
