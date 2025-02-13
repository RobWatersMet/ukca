module ukca_chemistry_setup_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_SETUP_MOD'
  public

  contains
  subroutine ukca_chemistry_setup(row_length, rows, model_levels, theta_field_size, bl_levels, ntracers, ntype, npft, len_stashwork, i_month, i_day_number, i_hour, land_points, land_index, tile_pts, tile_index, istore_h2so4, nlev_with_ddep, r_minute, secs_per_step, latitude, longitude, sinlat, tanlat, temp, pres, rh, p_layer_boundaries, t_surf, dzl, u_s, z0m, drain, crain, frac_types, seaice_frac, stcon, surf_hf, soilmc_lp, fland, laift_lp, canhtft_lp, z0tile_lp, t0tile_lp, zbl, h_plus, zdryrt, zwetrt, tracer, stashwork, firstcall)
    use asad_findreaction_mod, only : asad_findreaction
    use asad_mod, only : advt, ih2so4_hv, iso2_oh, jpcspf, jpdd, jpdw, jppj, jpspec, jpspj, jpspt, jptk, ndepd, ndepw, nldepd, nprkx, ntrkx, specf, speci, spj, spt
    use ukca_config_specification_mod, only : ukca_config
    use ukca_conserve_mod, only : ukca_conserve
    use ukca_cspecies, only : n_h2o, n_h2o2, n_h2so4, nn_cl
    use ukca_ddepctl_mod, only : ukca_ddepctl
    use ukca_ddeprt_mod, only : ukca_ddeprt
    use ukca_environment_fields_mod, only : h2o2_offline, surf_wetness
    use ukca_wdeprt_mod, only : ukca_wdeprt
    use ukca_um_legacy_mod, only : deposition_from_ukca_chemistry
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    use parkind1, only : jpim, jprb
    use umPrintMgr, only : ummessage, umprint
    use yomhook, only : dr_hook, lhook
    use nvtx
    logical, parameter :: before_chem = .true.
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: theta_field_size
    integer, intent(in) :: bl_levels
    integer, intent(in) :: ntracers
    integer, intent(in) :: ntype
    integer, intent(in) :: npft
    integer, intent(in) :: len_stashwork
    integer, intent(in) :: i_month
    integer, intent(in) :: i_day_number
    integer, intent(in) :: i_hour
    integer, intent(in) :: land_points
    integer, dimension(land_points), intent(in) :: land_index
    integer, dimension(ntype), intent(in) :: tile_pts
    integer, dimension(land_points,ntype), intent(in) :: tile_index
    integer, intent(out) :: istore_h2so4
    integer, dimension(row_length,rows), intent(out) :: nlev_with_ddep
    real, intent(in) :: r_minute
    real, intent(in) :: secs_per_step
    real, dimension(row_length,rows), intent(in) :: latitude
    real, dimension(row_length,rows), intent(in) :: longitude
    real, dimension(row_length,rows), intent(in) :: sinlat
    real, dimension(row_length,rows), intent(in) :: tanlat
    real, dimension(row_length,rows,model_levels), intent(in) :: temp
    real, dimension(row_length,rows,model_levels), intent(in) :: pres
    real, dimension(row_length,rows,model_levels), intent(in) :: rh
    real, dimension(row_length,rows,0:model_levels), intent(in) :: p_layer_boundaries
    real, dimension(row_length,rows), intent(in) :: t_surf
    real, dimension(row_length,rows,bl_levels), intent(in) :: dzl
    real, dimension(row_length,rows), intent(in) :: u_s
    real, dimension(row_length,rows), intent(in) :: z0m
    real, dimension(row_length,rows,model_levels), intent(in) :: drain
    real, dimension(row_length,rows,model_levels), intent(in) :: crain
    real, dimension(land_points,ntype), intent(in) :: frac_types
    real, dimension(row_length,rows), intent(in) :: seaice_frac
    real, dimension(row_length,rows,npft), intent(in) :: stcon
    real, dimension(row_length,rows), intent(in) :: surf_hf
    real, dimension(land_points), intent(in) :: soilmc_lp
    real, dimension(land_points), intent(in) :: fland
    real, dimension(land_points,npft), intent(in) :: laift_lp
    real, dimension(land_points,npft), intent(in) :: canhtft_lp
    real, dimension(land_points,ntype), intent(in) :: z0tile_lp
    real, dimension(land_points,ntype), intent(in) :: t0tile_lp
    real, dimension(row_length,rows), intent(in) :: zbl
    real, dimension(row_length,rows,model_levels), intent(in) :: h_plus
    real, dimension(row_length,rows,jpdd), intent(out) :: zdryrt
    real, dimension(row_length,rows,model_levels,jpdw), intent(out) :: zwetrt
    real, dimension(row_length,rows,model_levels,ntracers), intent(inout) :: tracer
    real, dimension(len_stashwork), intent(inout) :: stashwork
    logical, intent(in) :: firstcall
    integer :: jspf
    logical :: tropraq
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    CHARACTER(LEN = 10) :: prods(2)
    CHARACTER(LEN = 10) :: prods3(3)
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_SETUP'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    tropraq = ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero
    if (firstcall .AND. (.NOT.ukca_config%l_ukca_offline_be)) then
      if (ukca_config%l_ukca_advh2o .AND. n_h2o == 0) then
        errcode = 4
        cmessage(:) = 'No tracer for advected water vapour'
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      if (.NOT.tropraq) then
        iso2_oh = 0
        ih2so4_hv = 0
        istore_h2so4 = 0
        if (ukca_config%l_ukca_nr_aqchem) then
          do jspf = 1, jpcspf, 1
            if (specf(jspf) == advt(n_h2so4)) then
              istore_h2so4 = jspf
            end if
          enddo
          if (ukca_config%l_ukca_offline) then
            prods(:) = ['H2SO4     ', '          ']
          else
            prods(:) = ['H2SO4     ', 'HO2       ']
          end if
          iso2_oh = asad_findreaction('SO2       ', 'OH        ', prods, 2, spt, ntrkx, jptk + 1, jpspt)
          if (iso2_oh == 0) then
            prods(:) = ['SO3       ', 'HO2       ']
            iso2_oh = asad_findreaction('SO2       ', 'OH        ', prods, 2, spt, ntrkx, jptk + 1, jpspt)
            if (ukca_config%i_ukca_chem_version >= 121) then
              prods3(:) = ['SO3       ', 'OH        ', 'H         ']
              ih2so4_hv = asad_findreaction('H2SO4     ', 'PHOTON    ', prods3, 3, spj, nprkx, jppj + 1, jpspj)
            else
              prods(:) = ['SO3       ', 'OH        ']
              ih2so4_hv = asad_findreaction('H2SO4     ', 'PHOTON    ', prods, 2, spj, nprkx, jppj + 1, jpspj)
            end if
          end if
          if (iso2_oh == 0 .AND. ih2so4_hv == 0) then
            errcode = 1
            cmessage(:) = ' Sulphur chemistry reactions not found'
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') cmessage
            call umprint(ummessage, src='ukca_chemistry_ctl')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A,I0,A,I0)') 'iso2_oh: ', iso2_oh, ' ih2so4_hv: ', ih2so4_hv
            call umprint(ummessage, src='ukca_chemistry_ctl')
            call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
          end if
        end if
      end if
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx = 1, jpdd, 1
      do idx_1 = 1, rows, 1
        do idx_2 = 1, row_length, 1
          zdryrt(idx_2,idx_1,idx) = 0.0
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ndepd /= 0 .AND. (.NOT.ukca_config%l_ukca_drydep_off)) then
      if (ukca_config%l_ukca_intdd) then
        if (ukca_config%l_deposition_jules) then
          call deposition_from_ukca_chemistry(secs_per_step, bl_levels, row_length, rows, ntype, npft, jpspec, ndepd, nldepd, speci, land_points, land_index, tile_pts, tile_index, seaice_frac, fland, sinlat, p_layer_boundaries(:,:,0), rh(:,:,1), t_surf, surf_hf, surf_wetness, z0tile_lp, stcon, laift_lp, canhtft_lp, t0tile_lp, soilmc_lp, zbl, dzl, frac_types, u_s, zdryrt, nlev_with_ddep, len_stashwork, stashwork)
        else
          if (tropraq .AND. (.NOT.ukca_config%l_ukca_offline_be)) then
            call ukca_ddepctl(row_length, rows, bl_levels, ntype, npft, land_points, land_index, tile_pts, tile_index, secs_per_step, sinlat, frac_types, t_surf, p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s, rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp, canhtft_lp, z0tile_lp, t0tile_lp, nlev_with_ddep, zdryrt, len_stashwork, stashwork)
          else
            call ukca_ddepctl(row_length, rows, bl_levels, ntype, npft, land_points, land_index, tile_pts, tile_index, secs_per_step, sinlat, frac_types, t_surf, p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s, rh(:,:,1), stcon, soilmc_lp, fland, seaice_frac, laift_lp, canhtft_lp, z0tile_lp, t0tile_lp, nlev_with_ddep, zdryrt, len_stashwork, stashwork)
          end if
        end if
      else
        call ukca_ddeprt(theta_field_size, bl_levels, i_day_number, i_month, i_hour, r_minute, secs_per_step, longitude, latitude, tanlat, dzl, z0m, u_s, t_surf, zdryrt)
      end if
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(4)
    do idx_3 = 1, jpdw, 1
      do idx_4 = 1, model_levels, 1
        do idx_5 = 1, rows, 1
          do idx_6 = 1, row_length, 1
            zwetrt(idx_6,idx_5,idx_4,idx_3) = 0.0
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ndepw /= 0 .AND. (.NOT.ukca_config%l_ukca_wetdep_off)) then
      call ukca_wdeprt(theta_field_size, model_levels, drain, crain, temp, latitude, secs_per_step, zwetrt, h_plus)
    end if
    if (ukca_config%l_ukca_offline_be) then
      return
    end if
    if (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR. ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) then
      if (nn_cl > 0) then
        call ukca_conserve(row_length, rows, model_levels, ntracers, tracer, pres, drain, crain, before_chem)
      end if
    end if
    if (ukca_config%l_ukca_offline) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (tracer(:, :, :, n_h2o2) > h2o2_offline(:, :, :)) tracer(:, :, :, n_h2o2) = h2o2_offline(:, :, :)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine ukca_chemistry_setup

end module ukca_chemistry_setup_mod
