module ukca_ddepctl_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DDEPCTL_MOD'
  public

  contains
  subroutine ukca_ddepctl(row_length, rows, bl_levels, ntype, npft, land_points, land_index, tile_pts, tile_index, timestep, sinlat, frac_types, t_surf, p_surf, dzl, zbl, surf_hf, u_s, rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp, canhtft_lp, z0tile_lp, t0tile_lp, nlev_with_ddep, zdryrt, len_stashwork, stashwork)
    use ukca_config_specification_mod, only : ukca_config
    use ukca_config_constants_mod, only : rho_water, tfs
    use asad_mod, only : jpdd
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    use ukca_aerod_mod, only : ukca_aerod
    use ukca_ddcalc_mod, only : ukca_ddcalc
    use ukca_surfddr_mod, only : ukca_surfddr
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: bl_levels
    integer, intent(in) :: land_points
    integer, intent(in) :: ntype
    integer, intent(in) :: npft
    integer, dimension(land_points), intent(in) :: land_index
    integer, dimension(ntype), intent(in) :: tile_pts
    integer, dimension(land_points,ntype), intent(in) :: tile_index
    real, intent(in) :: timestep
    real, dimension(row_length,rows), intent(in) :: sinlat
    real, dimension(:,:), intent(in) :: frac_types
    real, dimension(row_length,rows), intent(in) :: t_surf
    real, dimension(row_length,rows), intent(in) :: p_surf
    real, dimension(row_length,rows,bl_levels), intent(in) :: dzl
    real, dimension(row_length,rows), intent(in) :: zbl
    real, dimension(row_length,rows), intent(in) :: surf_hf
    real, dimension(row_length,rows), intent(in) :: u_s
    real, dimension(row_length,rows), intent(in) :: rh
    real, dimension(row_length,rows), intent(in) :: seaice_frac
    real, dimension(row_length,rows,npft), intent(in) :: stcon
    real, dimension(land_points), intent(in) :: soilmc_lp
    real, dimension(land_points), intent(in) :: fland
    real, dimension(land_points,npft), intent(in) :: laift_lp
    real, dimension(land_points,npft), intent(in) :: canhtft_lp
    real, dimension(land_points,ntype), intent(in) :: z0tile_lp
    real, dimension(land_points,ntype), intent(in) :: t0tile_lp
    integer, intent(in) :: len_stashwork
    real, dimension(len_stashwork), intent(inout) :: stashwork
    integer, dimension(row_length,rows), intent(out) :: nlev_with_ddep
    real, dimension(row_length,rows,jpdd), intent(out) :: zdryrt
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    integer :: lake
    real :: seafrac
    real :: lftotal
    real, dimension(row_length,rows) :: smr
    real, dimension(row_length,rows) :: land_fraction
    real, dimension(row_length,rows,npft) :: lai_ft
    real, dimension(row_length,rows,npft) :: canht_ft
    real, dimension(row_length,rows,ntype) :: z0tile
    real, dimension(row_length,rows,ntype) :: t0tile
    real, dimension(row_length,rows) :: o3_stom_frac
    real, dimension(row_length,rows,ntype) :: resa
    real, dimension(row_length,rows,jpdd) :: rb
    real, dimension(row_length,rows,ntype,jpdd) :: rc
    real, dimension(row_length,rows,ntype) :: gsf
    real, dimension(row_length,rows) :: so4_vd
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DDEPCTL'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(2)
    do idx = 1, rows, 1
      do idx_1 = 1, row_length, 1
        land_fraction(idx_1,idx) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_2 = 1, rows, 1
      do idx_3 = 1, row_length, 1
        smr(idx_3,idx_2) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_4 = 1, npft, 1
      do idx_5 = 1, rows, 1
        do idx_6 = 1, row_length, 1
          lai_ft(idx_6,idx_5,idx_4) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_7 = 1, npft, 1
      do idx_8 = 1, rows, 1
        do idx_9 = 1, row_length, 1
          canht_ft(idx_9,idx_8,idx_7) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_10 = 1, ntype, 1
      do idx_11 = 1, rows, 1
        do idx_12 = 1, row_length, 1
          z0tile(idx_12,idx_11,idx_10) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_13 = 1, ntype, 1
      do idx_14 = 1, rows, 1
        do idx_15 = 1, row_length, 1
          t0tile(idx_15,idx_14,idx_13) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_16 = 1, ntype, 1
      do idx_17 = 1, rows, 1
        do idx_18 = 1, row_length, 1
          gsf(idx_18,idx_17,idx_16) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop seq
    do l = 1, land_points, 1
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j - 1) * row_length
      smr(i,j) = soilmc_lp(l) / (ukca_config%dzsoil_layer1 * rho_water)
      land_fraction(i,j) = fland(l)
    enddo
    !NO!$acc loop gang vector independent collapse(1)
    do m = 1, npft, 1
      !NO!$acc loop seq
      do n = 1, tile_pts(m), 1
        l = tile_index(n,m)
        j = (land_index(l) - 1) / row_length + 1
        i = land_index(l) - (j - 1) * row_length
        lai_ft(i,j,m) = laift_lp(l,m)
        canht_ft(i,j,m) = canhtft_lp(l,m)
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(1)
    do m = 1, ntype, 1
      !NO!$acc loop seq
      do n = 1, tile_pts(m), 1
        l = tile_index(n,m)
        j = (land_index(l) - 1) / row_length + 1
        i = land_index(l) - (j - 1) * row_length
        z0tile(i,j,m) = z0tile_lp(l,m)
        t0tile(i,j,m) = t0tile_lp(l,m)
        gsf(i,j,m) = frac_types(l,m)
      enddo
    enddo
    lake = ukca_config%i_lake
    !NO!$acc end kernels
    do k = 1, rows, 1
      do i = 1, row_length, 1
        if (gsf(i,k,lake) == 1.0) then
          land_fraction(i,k) = 0.0
        end if
        seafrac = 1.0 - land_fraction(i,k)
        if (land_fraction(i,k) < 1.0 .AND. gsf(i,k,lake) < 1.0) then
          gsf(i,k,lake) = seafrac
          if (land_fraction(i,k) > 0.0) then
            lftotal = 0.0
            n = 0
            do j = 1, ntype, 1
              n = n + 1
              if (j == lake) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Cycle_Stmt
                CYCLE
              end if
              lftotal = lftotal + gsf(i,k,n)
            enddo
            n = 0
            do j = 1, ntype, 1
              n = n + 1
              if (j == lake) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Cycle_Stmt
                CYCLE
              end if
              gsf(i,k,n) = gsf(i,k,n) * land_fraction(i,k) / lftotal
            enddo
          end if
        end if
        if (seaice_frac(i,k) > 0.0) then
          gsf(i,k,lake) = (1.0 - seaice_frac(i,k)) * seafrac
          gsf(i,k,ntype) = gsf(i,k,ntype) + seaice_frac(i,k) * seafrac
        end if
      enddo
    enddo
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(3)
    do n = 1, ntype, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          if (t0tile(i,k,n) < 100.0) then
            t0tile(i,k,n) = t_surf(i,k)
          end if
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do k = 1, rows, 1
      do i = 1, row_length, 1
        if (seaice_frac(i,k) > 0.0) then
          t0tile(i,k,lake) = tfs
        end if
      enddo
    enddo
    !NO!$acc end kernels
    call ukca_aerod(row_length, rows, ntype, npft, lake, t_surf, p_surf, surf_hf, u_s, canht_ft, gsf, zbl, z0tile, resa, rb, so4_vd)
    call ukca_surfddr(row_length, rows, ntype, npft, sinlat, t_surf, p_surf, rh, smr, u_s, gsf, stcon, t0tile, lai_ft, so4_vd, rc, o3_stom_frac)
    call ukca_ddcalc(row_length, rows, bl_levels, ntype, npft, timestep, dzl, zbl, gsf, resa, rb, rc, nlev_with_ddep, zdryrt, len_stashwork, stashwork)
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_ddepctl

end module ukca_ddepctl_mod
