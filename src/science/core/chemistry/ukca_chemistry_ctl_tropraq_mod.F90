module ukca_chemistry_ctl_tropraq_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_CTL_TROPRAQ_MOD'
  public

  contains
  subroutine ukca_chemistry_ctl_tropraq(row_length, rows, model_levels, theta_field_size, tot_n_pnts, ntracers, secs_per_step, pres, temp, q, qcf, qcl, rh, tracer, all_ntp, cloud_frac, photol_rates, volume, so4_aitken, so4_accum, soot_fresh, soot_aged, ocff_fresh, ocff_aged, biogenic, sea_salt_film, sea_salt_jet, uph2so4inaer, delso2_wet_h2o2, delso2_wet_o3, delh2so4_chem, delso2_drydep, delso2_wetdep, trop_ch4_mol, trop_o3_mol, trop_oh_mol, strat_ch4_mol, strat_ch4loss, len_stashwork, stashwork, h_plus, zdryrt, zwetrt, nlev_with_ddep, l_stratosphere, firstcall)
    use ukca_um_legacy_mod, only : rgas=>r
    use asad_mod, only : advt, cdt, ctype, jpcspf, jpctr, jpdd, jpdw, jpeq, jphk, jpnr, jppj, jpro2, jpspec, jpspj, jpspt, jptk, ldepd, ldepw, nadvt, nlnaro2, nnaf, nprkx, ntrkx, speci, sph2o, spj, spro2, spt, y
    use ukca_config_defs_mod, only : nr_phot, nr_therm
    use ukca_cspecies, only : c_na_species, c_species, n_ch4, n_h2o, n_h2so4, n_hono2, n_o3, nn_ch4, nn_cl, nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh, nn_so2
    use ukca_constants, only : c_h2o, c_hono2
    use ukca_config_constants_mod, only : avogadro, boltzmann
    use ukca_config_specification_mod, only : ukca_config
    use ukca_raq_diags_mod, only : ukca_raq_diags
    use ukca_ntp_mod, only : dim_ntp, name2ntpindex, ntp_type
    use ukca_chemco_raq_mod, only : ukca_chemco_raq
    use ukca_deriv_raqaero_mod, only : ukca_deriv_raqaero
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use ukca_missing_data_mod, only : rmdi
    use errormessagelength_mod, only : errormessagelength
    use ukca_be_drydep_mod, only : ukca_be_drydep
    use ukca_be_wetdep_mod, only : ukca_be_wetdep
    use ukca_ch4_stratloss_mod, only : ukca_ch4_stratloss
    use ukca_chemco_mod, only : ukca_chemco
    use ukca_deriv_mod, only : ukca_deriv
    use ukca_deriv_aero_mod, only : ukca_deriv_aero
    use ukca_deriv_raq_mod, only : ukca_deriv_raq
    use ukca_fracdiss_mod, only : ukca_fracdiss
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: theta_field_size
    integer, intent(in) :: tot_n_pnts
    integer, intent(in) :: ntracers
    integer, intent(in) :: uph2so4inaer
    integer, dimension(theta_field_size), intent(in) :: nlev_with_ddep
    real, intent(in) :: secs_per_step
    real, dimension(tot_n_pnts), intent(in) :: pres
    real, dimension(tot_n_pnts), intent(in) :: temp
    real, dimension(tot_n_pnts), intent(in) :: volume
    real, dimension(tot_n_pnts), intent(in) :: h_plus
    real, dimension(tot_n_pnts), intent(in) :: so4_aitken
    real, dimension(tot_n_pnts), intent(in) :: so4_accum
    real, dimension(tot_n_pnts), intent(in) :: soot_fresh
    real, dimension(tot_n_pnts), intent(in) :: soot_aged
    real, dimension(tot_n_pnts), intent(in) :: ocff_fresh
    real, dimension(tot_n_pnts), intent(in) :: ocff_aged
    real, dimension(tot_n_pnts), intent(in) :: biogenic
    real, dimension(tot_n_pnts), intent(in) :: sea_salt_film
    real, dimension(tot_n_pnts), intent(in) :: sea_salt_jet
    real, dimension(tot_n_pnts), intent(in) :: qcf
    real, dimension(tot_n_pnts), intent(in) :: qcl
    real, dimension(tot_n_pnts), intent(in) :: rh
    real, dimension(tot_n_pnts), intent(in) :: cloud_frac
    real, dimension(theta_field_size,jpdd), intent(in) :: zdryrt
    real, dimension(theta_field_size,model_levels,jpdw), intent(in) :: zwetrt
    real, dimension(tot_n_pnts,jppj), intent(in) :: photol_rates
    real, dimension(tot_n_pnts), intent(inout) :: q
    real, dimension(tot_n_pnts,ntracers), intent(inout) :: tracer
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wet_h2o2
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wet_o3
    real, dimension(tot_n_pnts), intent(inout) :: delh2so4_chem
    real, dimension(tot_n_pnts), intent(inout) :: delso2_drydep
    real, dimension(tot_n_pnts), intent(inout) :: delso2_wetdep
    real, dimension(tot_n_pnts), intent(inout) :: trop_ch4_mol
    real, dimension(tot_n_pnts), intent(inout) :: trop_o3_mol
    real, dimension(tot_n_pnts), intent(inout) :: trop_oh_mol
    real, dimension(tot_n_pnts), intent(inout) :: strat_ch4_mol
    real, dimension(tot_n_pnts), intent(inout) :: strat_ch4loss
    type(ntp_type), dimension(dim_ntp), intent(inout) :: all_ntp
    integer, intent(in) :: len_stashwork
    real, dimension(len_stashwork), intent(inout) :: stashwork
    logical, dimension(tot_n_pnts), intent(in) :: l_stratosphere
    logical, intent(in) :: firstcall
    integer :: nlev_ch4_stratloss
    integer, save :: nr
    integer, save :: n_be_calls
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
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    real, save :: dts
    real, dimension(theta_field_size) :: so2_wetox_h2o2
    real, dimension(theta_field_size) :: so2_wetox_o3
    real, dimension(theta_field_size) :: so2_dryox_oh
    real, allocatable, dimension(:,:) :: be_rc
    real, allocatable, dimension(:,:) :: be_hrc
    real, allocatable, dimension(:,:) :: zfnatr
    real, allocatable, dimension(:) :: ystore
    real, dimension(theta_field_size,jpcspf) :: zftr
    real, dimension(theta_field_size,model_levels,jpdw,jpeq + 1) :: zfrdiss
    real, dimension(row_length,rows,model_levels) :: kp_nh
    real, dimension(theta_field_size) :: be_tnd
    real, dimension(theta_field_size) :: be_h2o
    real, dimension(theta_field_size) :: be_o2
    real, dimension(theta_field_size) :: be_vol
    real, dimension(theta_field_size) :: be_rho
    real, dimension(theta_field_size) :: be_rh_frac
    real, dimension(theta_field_size) :: be_so4_aitken
    real, dimension(theta_field_size) :: be_so4_accum
    real, dimension(theta_field_size) :: be_soot_fresh
    real, dimension(theta_field_size) :: be_soot_aged
    real, dimension(theta_field_size) :: be_ocff_fresh
    real, dimension(theta_field_size) :: be_ocff_aged
    real, dimension(theta_field_size) :: be_biogenic
    real, dimension(theta_field_size) :: be_sea_salt_film
    real, dimension(theta_field_size) :: be_sea_salt_jet
    real, dimension(theta_field_size,jpspec) :: be_wetrt
    real, dimension(theta_field_size,jpspec) :: be_dryrt
    real, dimension(theta_field_size,jpspec,jpeq + 1) :: be_frdiss
    real, dimension(theta_field_size,jpspec) :: be_y
    real, dimension(theta_field_size,5) :: k_dms
    real, dimension(theta_field_size,model_levels,jpdd) :: dry_dep_3d
    real, dimension(theta_field_size,model_levels,jpdw) :: wet_dep_3d
    real, dimension(tot_n_pnts,dim_ntp) :: ntp_data
    logical, dimension(theta_field_size) :: stratflag
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_CTL_TROPRAQ'
    integer :: widx1
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    nr = nr_therm + nr_phot
    if (firstcall) then
      n_be_calls = INT(secs_per_step / ukca_config%dts0)
      if (n_be_calls == 0) then
        n_be_calls = 1
      end if
      dts = secs_per_step / n_be_calls
      if (printstatus >= prstatus_oper) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I0,E12.4)') 'n_be_calls, dts= ', n_be_calls, dts
        call umprint(ummessage, src='ukca_chemistry_ctl_tropraq')
      end if
    end if
    if (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) then
      call ukca_fracdiss(row_length, rows, model_levels, temp, pres, rh, qcl, zfrdiss, kp_nh, h_plus)
    end if
    if (ukca_config%nlev_ch4_stratloss >= 0) then
      nlev_ch4_stratloss = ukca_config%nlev_ch4_stratloss
    else
      nlev_ch4_stratloss = 3
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
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx = LBOUND(sph2o, dim=1), UBOUND(sph2o, dim=1), 1
            sph2o(idx) = qcf(idx + (kcs - LBOUND(sph2o, dim=1))) / c_h2o
          enddo
          !NO!$acc end kernels
        else
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_1 = LBOUND(sph2o, dim=1), UBOUND(sph2o, dim=1), 1
            sph2o(idx_1) = 0.0
          enddo
          !NO!$acc end kernels
        end if
      end if
      jspf = 0
      do js = 1, jpspec, 1
        do jtr = 1, jpctr, 1
          if (advt(jtr) == speci(js)) then
            !NO!$acc kernels
            jspf = jspf + 1
            !NO!$acc loop gang vector independent
            do idx_2 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
              zftr(idx_2,jspf) = tracer(idx_2 + (kcs - LBOUND(zftr, dim=1)),jtr) / c_species(jtr)
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
                do idx_3 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
                  zftr(idx_3,jspf) = ntp_data(idx_3 + (kcs - LBOUND(zftr, dim=1)),l) / c_na_species(jna)
                enddo
                !NO!$acc end kernels
              else
                errcode = jro2
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl_tropraq'
                call umprint(ummessage, src='ukca_chemistry_ctl_tropraq')
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
          WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl_tropraq'
          call umprint(ummessage, src='ukca_chemistry_ctl_tropraq')
          cmessage(:) = 'ERROR: Number of chemical active species /= jpro2+jpctr'
          call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
        end if
      else
        if (jspf /= jpctr) then
          errcode = jspf
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') '** ERROR in ukca_chemistry_ctl_tropraq'
          call umprint(ummessage, src='ukca_chemistry_ctl_tropraq')
          cmessage(:) = 'ERROR: Number of chemical active species /= jpctr'
          call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
        end if
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_4 = LBOUND(be_tnd, dim=1), UBOUND(be_tnd, dim=1), 1
        be_tnd(idx_4) = pres(idx_4 + (kcs - LBOUND(be_tnd, dim=1))) / (boltzmann * 1.0e6 * temp(idx_4 + (kcs - LBOUND(be_tnd, dim=1))))
      enddo
      !NO!$acc loop gang vector independent
      do idx_5 = LBOUND(be_o2, dim=1), UBOUND(be_o2, dim=1), 1
        be_o2(idx_5) = 0.2095 * be_tnd(idx_5)
      enddo
      !NO!$acc loop gang vector independent
      do idx_6 = LBOUND(be_h2o, dim=1), UBOUND(be_h2o, dim=1), 1
        be_h2o(idx_6) = q(idx_6 + (kcs - LBOUND(be_h2o, dim=1))) * be_tnd(idx_6)
      enddo
      !NO!$acc loop gang vector independent
      do idx_7 = LBOUND(be_vol, dim=1), UBOUND(be_vol, dim=1), 1
        be_vol(idx_7) = volume(idx_7 + (kcs - LBOUND(be_vol, dim=1))) * 1.0e6
      enddo
      be_rho(:) = pres(kcs:kce) / (rgas * temp(kcs:kce))
      !NO!$acc loop gang vector independent collapse(2)
      do l = 1, jpctr, 1
        do idx_8 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
          zftr(idx_8,l) = zftr(idx_8,l) * be_tnd(idx_8)
        enddo
      enddo
      !NO!$acc end kernels
      if (.NOT.ALLOCATED(zfnatr)) then
        ALLOCATE(zfnatr(1:theta_field_size,1:nnaf))
      end if
      do js = 1, nnaf, 1
        l = name2ntpindex(nadvt(js))
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_9 = LBOUND(zfnatr, dim=1), UBOUND(zfnatr, dim=1), 1
          zfnatr(idx_9,js) = ntp_data(idx_9 + (kcs - LBOUND(zfnatr, dim=1)),l) / c_na_species(js)
        enddo
        !NO!$acc end kernels
      enddo
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do l = 1, nnaf, 1
        do idx_10 = LBOUND(zfnatr, dim=1), UBOUND(zfnatr, dim=1), 1
          zfnatr(idx_10,l) = zfnatr(idx_10,l) * be_tnd(idx_10 + (LBOUND(be_tnd, dim=1) - LBOUND(zfnatr, dim=1)))
        enddo
      enddo
      !NO!$acc end kernels
      call ukca_be_wetdep(theta_field_size, zwetrt(:,k,:), be_wetrt)
      call ukca_be_drydep(k, theta_field_size, nlev_with_ddep, zdryrt, be_dryrt)
      if (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) then
        do i = 1, jpeq + 1, 1
          call ukca_be_wetdep(theta_field_size, zfrdiss(:,k,:,i), be_frdiss(:,:,i))
        enddo
      end if
      if (.NOT.ALLOCATED(be_rc)) then
        ALLOCATE(be_rc(1:theta_field_size,1:nr_therm))
      end if
      if (ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_11 = LBOUND(be_rh_frac, dim=1), UBOUND(be_rh_frac, dim=1), 1
          be_rh_frac(idx_11) = rh(idx_11 + (kcs - LBOUND(be_rh_frac, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_12 = LBOUND(be_so4_aitken, dim=1), UBOUND(be_so4_aitken, dim=1), 1
          be_so4_aitken(idx_12) = so4_aitken(idx_12 + (kcs - LBOUND(be_so4_aitken, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_13 = LBOUND(be_so4_accum, dim=1), UBOUND(be_so4_accum, dim=1), 1
          be_so4_accum(idx_13) = so4_accum(idx_13 + (kcs - LBOUND(be_so4_accum, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_14 = LBOUND(be_soot_fresh, dim=1), UBOUND(be_soot_fresh, dim=1), 1
          be_soot_fresh(idx_14) = soot_fresh(idx_14 + (kcs - LBOUND(be_soot_fresh, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_15 = LBOUND(be_soot_aged, dim=1), UBOUND(be_soot_aged, dim=1), 1
          be_soot_aged(idx_15) = soot_aged(idx_15 + (kcs - LBOUND(be_soot_aged, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_16 = LBOUND(be_ocff_fresh, dim=1), UBOUND(be_ocff_fresh, dim=1), 1
          be_ocff_fresh(idx_16) = ocff_fresh(idx_16 + (kcs - LBOUND(be_ocff_fresh, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_17 = LBOUND(be_ocff_aged, dim=1), UBOUND(be_ocff_aged, dim=1), 1
          be_ocff_aged(idx_17) = ocff_aged(idx_17 + (kcs - LBOUND(be_ocff_aged, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_18 = LBOUND(be_biogenic, dim=1), UBOUND(be_biogenic, dim=1), 1
          be_biogenic(idx_18) = biogenic(idx_18 + (kcs - LBOUND(be_biogenic, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_19 = LBOUND(be_sea_salt_film, dim=1), UBOUND(be_sea_salt_film, dim=1), 1
          be_sea_salt_film(idx_19) = sea_salt_film(idx_19 + (kcs - LBOUND(be_sea_salt_film, dim=1)))
        enddo
        !NO!$acc loop gang vector independent
        do idx_20 = LBOUND(be_sea_salt_jet, dim=1), UBOUND(be_sea_salt_jet, dim=1), 1
          be_sea_salt_jet(idx_20) = sea_salt_jet(idx_20 + (kcs - LBOUND(be_sea_salt_jet, dim=1)))
        enddo
        !NO!$acc end kernels
        if (.NOT.ALLOCATED(be_hrc)) then
          ALLOCATE(be_hrc(1:theta_field_size,1:jphk))
        end if
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_21 = LBOUND(stratflag, dim=1), UBOUND(stratflag, dim=1), 1
          stratflag(idx_21) = l_stratosphere(idx_21 + (kcs - LBOUND(stratflag, dim=1)))
        enddo
        !NO!$acc end kernels
        call ukca_chemco_raq(nr_therm, theta_field_size, temp(kcs:kce), be_tnd, be_h2o, be_o2, qcl(kcs:kce), cloud_frac(kcs:kce), be_frdiss, be_rho, be_rh_frac, be_so4_aitken, be_so4_accum, be_soot_fresh, be_soot_aged, be_ocff_fresh, be_ocff_aged, be_biogenic, be_sea_salt_film, be_sea_salt_jet, stratflag, be_rc, be_hrc, h_plus(kcs:kce))
      else
        call ukca_chemco(nr_therm, theta_field_size, temp(kcs:kce), be_tnd, be_h2o, be_o2, qcl(kcs:kce), cloud_frac(kcs:kce), be_frdiss, k_dms, be_rc, h_plus(kcs:kce))
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx_22 = LBOUND(be_y, dim=2), UBOUND(be_y, dim=2), 1
        do idx_23 = LBOUND(be_y, dim=1), UBOUND(be_y, dim=1), 1
          be_y(idx_23,idx_22) = 0.0
        enddo
      enddo
      !NO!$acc end kernels
      do i = 1, jpspec, 1
        do j = 1, jpctr, 1
          if (speci(i) == advt(j)) then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_24 = LBOUND(be_y, dim=1), UBOUND(be_y, dim=1), 1
              be_y(idx_24,i) = zftr(idx_24,j)
            enddo
            !NO!$acc end kernels
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
      enddo
      do i = 1, jpspec, 1
        do j = 1, nnaf, 1
          if (speci(i) == nadvt(j)) then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_25 = LBOUND(be_y, dim=1), UBOUND(be_y, dim=1), 1
              be_y(idx_25,i) = zfnatr(idx_25 + (LBOUND(zfnatr, dim=1) - LBOUND(be_y, dim=1)),j)
            enddo
            !NO!$acc end kernels
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
      enddo
      if (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_26 = LBOUND(so2_wetox_h2o2, dim=1), UBOUND(so2_wetox_h2o2, dim=1), 1
          so2_wetox_h2o2(idx_26) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_27 = LBOUND(so2_dryox_oh, dim=1), UBOUND(so2_dryox_oh, dim=1), 1
          so2_dryox_oh(idx_27) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_28 = LBOUND(so2_wetox_o3, dim=1), UBOUND(so2_wetox_o3, dim=1), 1
          so2_wetox_o3(idx_28) = 0.0
        enddo
        !NO!$acc end kernels
      else
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_29 = kcs, kce, 1
          delso2_wet_h2o2(idx_29) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_30 = kcs, kce, 1
          delso2_wet_o3(idx_30) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_31 = kcs, kce, 1
          delh2so4_chem(idx_31) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_32 = kcs, kce, 1
          delso2_drydep(idx_32) = 0.0
        enddo
        !NO!$acc loop gang vector independent
        do idx_33 = kcs, kce, 1
          delso2_wetdep(idx_33) = 0.0
        enddo
        !NO!$acc end kernels
      end if
      if (ukca_config%l_ukca_aerchem) then
        call ukca_deriv_aero(nr_therm, n_be_calls, theta_field_size, be_rc, be_wetrt, be_dryrt, photol_rates(kcs:kce,:), k_dms, be_h2o, be_tnd, be_o2, dts, be_y, so2_wetox_h2o2, so2_wetox_o3, so2_dryox_oh)
      else
        if (ukca_config%l_ukca_raq) then
          call ukca_deriv_raq(nr_therm, n_be_calls, theta_field_size, dts, be_h2o, be_vol, be_rc, be_hrc, be_dryrt, be_wetrt, photol_rates(kcs:kce,:), ldepd, ldepw, be_y, dry_dep_3d(:,k,:), wet_dep_3d(:,k,:))
        else
          if (ukca_config%l_ukca_raqaero) then
            if (uph2so4inaer == 1) then
              !NO!$acc kernels
              !NO!$acc loop gang vector independent
              do idx_34 = LBOUND(ystore, dim=1), UBOUND(ystore, dim=1), 1
                ystore(idx_34) = be_y(idx_34 + (LBOUND(be_y, dim=1) - LBOUND(ystore, dim=1)),nn_h2so4)
              enddo
              !NO!$acc end kernels
            end if
            call ukca_deriv_raqaero(nr_therm, n_be_calls, theta_field_size, dts, be_rc, be_wetrt, be_dryrt, photol_rates(kcs:kce,:), be_h2o, be_y, so2_wetox_h2o2, so2_wetox_o3, so2_dryox_oh)
            if (uph2so4inaer == 1) then
              !NO!$acc kernels
              !NO!$acc loop gang vector independent
              do idx_35 = LBOUND(be_y, dim=1), UBOUND(be_y, dim=1), 1
                be_y(idx_35,nn_h2so4) = ystore(idx_35 + (LBOUND(ystore, dim=1) - LBOUND(be_y, dim=1)))
              enddo
              !NO!$acc end kernels
            end if
          else
            call ukca_deriv(nr, n_be_calls, theta_field_size, be_rc, be_wetrt, be_dryrt, photol_rates(kcs:kce,:), be_h2o, be_tnd, be_o2, dts, be_y)
          end if
        end if
      end if
      if (ALLOCATED(be_hrc)) then
        DEALLOCATE(be_hrc)
      end if
      if (ALLOCATED(be_rc)) then
        DEALLOCATE(be_rc)
      end if
      if (k >= model_levels + 1 - nlev_ch4_stratloss) then
        call ukca_ch4_stratloss(n_be_calls, theta_field_size, be_vol, dts, be_y(:,nn_ch4), strat_ch4loss(kcs:kce))
      end if
      do j = 1, jpctr, 1
        do i = 1, jpspec, 1
          if (advt(j) == speci(i)) then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_36 = LBOUND(zftr, dim=1), UBOUND(zftr, dim=1), 1
              zftr(idx_36,j) = be_y(idx_36,i)
            enddo
            !NO!$acc end kernels
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
      enddo
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do js = 1, jpctr, 1
        do idx_37 = kcs, kce, 1
          tracer(idx_37,js) = zftr(idx_37 + (LBOUND(zftr, dim=1) - kcs),js) / be_tnd(idx_37 + (LBOUND(be_tnd, dim=1) - kcs)) * c_species(js)
        enddo
      enddo
      !NO!$acc end kernels
      do j = 1, nnaf, 1
        do i = 1, jpspec, 1
          if (nadvt(j) == speci(i)) then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_38 = LBOUND(zfnatr, dim=1), UBOUND(zfnatr, dim=1), 1
              zfnatr(idx_38,j) = be_y(idx_38 + (LBOUND(be_y, dim=1) - LBOUND(zfnatr, dim=1)),i) / be_tnd(idx_38 + (LBOUND(be_tnd, dim=1) - LBOUND(zfnatr, dim=1)))
            enddo
            !NO!$acc end kernels
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
      enddo
      do js = 1, nnaf, 1
        l = name2ntpindex(nadvt(js))
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_39 = kcs, kce, 1
          ntp_data(idx_39,l) = zfnatr(idx_39 + (LBOUND(zfnatr, dim=1) - kcs),js) * c_na_species(js)
        enddo
        !NO!$acc end kernels
      enddo
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_40 = kcs, kce, 1
        trop_ch4_mol(idx_40) = zftr(idx_40 + (LBOUND(zftr, dim=1) - kcs),n_ch4) * volume(idx_40) * 1.0e6 / avogadro
      enddo
      !NO!$acc loop gang vector independent
      do idx_41 = kcs, kce, 1
        strat_ch4_mol(idx_41) = trop_ch4_mol(idx_41)
      enddo
      !NO!$acc loop gang vector independent
      do idx_42 = kcs, kce, 1
        trop_o3_mol(idx_42) = zftr(idx_42 + (LBOUND(zftr, dim=1) - kcs),n_o3) * volume(idx_42) * 1.0e6 / avogadro
      enddo
      !NO!$acc loop gang vector independent
      do idx_43 = kcs, kce, 1
        trop_oh_mol(idx_43) = be_y(idx_43 + (LBOUND(be_y, dim=1) - kcs),nn_oh) * volume(idx_43) * 1.0e6 / avogadro
      enddo
      !NO!$acc end kernels
      if (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_44 = kcs, kce, 1
          delso2_wet_h2o2(idx_44) = so2_wetox_h2o2(idx_44 + (LBOUND(so2_wetox_h2o2, dim=1) - kcs))
        enddo
        !NO!$acc loop gang vector independent
        do idx_45 = kcs, kce, 1
          delso2_wet_o3(idx_45) = so2_wetox_o3(idx_45 + (LBOUND(so2_wetox_o3, dim=1) - kcs))
        enddo
        !NO!$acc loop gang vector independent
        do idx_46 = kcs, kce, 1
          delh2so4_chem(idx_46) = so2_dryox_oh(idx_46 + (LBOUND(so2_dryox_oh, dim=1) - kcs))
        enddo
        !NO!$acc loop gang vector independent
        do idx_47 = kcs, kce, 1
          delso2_drydep(idx_47) = be_dryrt(idx_47 + (LBOUND(be_dryrt, dim=1) - kcs),nn_so2) * be_y(idx_47 + (LBOUND(be_y, dim=1) - kcs),nn_so2) * dts
        enddo
        !NO!$acc loop gang vector independent
        do idx_48 = kcs, kce, 1
          delso2_wetdep(idx_48) = be_wetrt(idx_48 + (LBOUND(be_wetrt, dim=1) - kcs),nn_so2) * be_y(idx_48 + (LBOUND(be_y, dim=1) - kcs),nn_so2) * dts
        enddo
        !NO!$acc end kernels
      end if
    enddo
    if (ALLOCATED(zfnatr)) then
      DEALLOCATE(zfnatr)
    end if
    if (ALLOCATED(ystore)) then
      DEALLOCATE(ystore)
    end if
    do l = 1, dim_ntp, 1
      if (all_ntp(l)%l_required) then
        all_ntp(l)%data_3d(:,:,:) = RESHAPE(ntp_data(:,l), [row_length, rows, model_levels])
      end if
    enddo
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do widx1 = 1, tot_n_pnts, 1
      if (l_stratosphere(widx1)) then
        trop_ch4_mol(widx1) = 0.0
        trop_o3_mol(widx1) = 0.0
        trop_oh_mol(widx1) = 0.0
      else
        strat_ch4_mol(widx1) = 0.0
      end if
    enddo
    !NO!$acc end kernels
    if (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_raq) then
      call ukca_raq_diags(row_length, rows, model_levels, ntracers, dry_dep_3d, wet_dep_3d, tracer, pres, temp, len_stashwork, stashwork)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_chemistry_ctl_tropraq

end module ukca_chemistry_ctl_tropraq_mod
