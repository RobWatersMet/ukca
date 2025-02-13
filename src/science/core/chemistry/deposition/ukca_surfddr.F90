module ukca_surfddr_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_SURFDDR_MOD'
  public

  contains
  subroutine ukca_surfddr(row_length, rows, ntype, npft, sinlat, t0, p0, rh, smr, u_s, gsf, stcon, t0tile, lai_ft, so4_vd, rc, o3_stom_frac)
    use asad_mod, only : jpdd, ndepd, nldepd, speci
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    use ukca_ddepo3_ocean_mod, only : ukca_ddepo3_ocean
    use ukca_environment_fields_mod, only : surf_wetness
    use parkind1, only : jpim, jprb
    use ukca_config_constants_mod, only : rmol
    use ukca_config_specification_mod, only : ukca_config
    use umPrintMgr, only : printstatus, prstatus_diag, ummessage, umprint
    use yomhook, only : dr_hook, lhook
    use nvtx
    real, parameter :: r_null = 1.0e50
    integer, parameter :: n_elev_ice = 10
    integer, parameter :: first_term = 1
    integer, parameter :: second_term = 2
    integer, parameter :: third_term = 3
    integer, parameter :: fourth_term = 4
    integer, parameter :: r_stom_spec = 5
    integer, parameter :: r_stom_no2 = 1
    integer, parameter :: r_stom_o3 = 2
    integer, parameter :: r_stom_pan = 3
    integer, parameter :: r_stom_so2 = 4
    integer, parameter :: r_stom_nh3 = 5
    real, parameter :: tol = 1.0e-10
    real, parameter :: soil_moist_thresh = 0.3
    real, parameter :: rh_thresh = 0.813
    real, parameter :: minus1degc = 272.15
    real, parameter :: minus5degc = 268.15
    real, parameter :: r_wet_so2 = 1.0
    real, parameter :: r_dry_so2 = 2000.0
    real, parameter :: r_cut_so2_5degc = 500.0
    real, parameter :: r_cut_so2_5to_1 = 200.0
    real, parameter :: so2con1 = 2.5e4
    real, parameter :: so2con2 = 0.58e12
    real, parameter :: so2con3 = -6.93
    real, parameter :: so2con4 = -27.8
    real, parameter :: smfrac1 = 0.16
    real, parameter :: smfrac2 = 0.30
    real, parameter :: smfrac3 = 0.50
    real, parameter :: smfrac4 = 0.20
    real, parameter :: tar_scaling = 15.0
    real, parameter :: glmin = 1.0e-6
    integer, parameter :: r_tundra_spec = 5
    integer, parameter :: r_tundra_no2 = 1
    integer, parameter :: r_tundra_co = 2
    integer, parameter :: r_tundra_o3 = 3
    integer, parameter :: r_tundra_pan = 5
    real, parameter :: min_tile_frac = 0.0
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: ntype
    integer, intent(in) :: npft
    real, dimension(row_length,rows), intent(in) :: sinlat
    real, dimension(row_length,rows), intent(in) :: t0
    real, dimension(row_length,rows), intent(in) :: p0
    real, dimension(row_length,rows), intent(in) :: rh
    real, dimension(row_length,rows), intent(in) :: smr
    real, dimension(row_length,rows), intent(in) :: u_s
    real, dimension(row_length,rows,ntype), intent(in) :: gsf
    real, dimension(row_length,rows,ntype), intent(in) :: t0tile
    real, dimension(row_length,rows,npft), intent(in) :: stcon
    real, dimension(row_length,rows,npft), intent(in) :: lai_ft
    real, dimension(row_length,rows), intent(in) :: so4_vd
    real, dimension(row_length,rows,ntype,jpdd), intent(out) :: rc
    real, dimension(row_length,rows), intent(out) :: o3_stom_frac
    integer :: errcode
    logical, save :: first = .true.
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer, save :: brd_leaf
    integer, save :: brd_leaf_dec
    integer, save :: brd_leaf_eg_trop
    integer, save :: brd_leaf_eg_temp
    integer, save :: ndl_leaf
    integer, save :: ndl_leaf_dec
    integer, save :: ndl_leaf_eg
    integer, save :: c3_grass
    integer, save :: c3_crop
    integer, save :: c3_pasture
    integer, save :: c4_grass
    integer, save :: c4_crop
    integer, save :: c4_pasture
    integer, save :: shrub
    integer, save :: shrub_dec
    integer, save :: shrub_eg
    integer, save :: urban
    integer, save :: lake
    integer, save :: soil
    integer, save :: ice
    integer, dimension(n_elev_ice), save :: elev_ice
    integer :: i
    integer :: k
    integer :: j
    integer :: n
    real :: sm
    real :: rr
    real :: ts
    real :: f_ch4_uptk
    real :: t_surf1
    real :: p_surf1
    real :: ustar
    real :: sst
    real :: rc_ocean
    real, save :: mml = 1.008e5
    real, save :: r_wet_o3 = 500.0
    real, save :: cuticle_o3 = 5000.0
    real, save :: tundra_s_limit = 0.866
    real, dimension(row_length,npft) :: r_cut_o3
    real :: r_cut_so2
    real, dimension(row_length,rows,npft,r_stom_spec) :: r_stom
    real, allocatable, dimension(:,:), save :: rsurf
    real, allocatable, dimension(:), save :: r_null_tile
    real, allocatable, dimension(:), save :: rooh
    real, allocatable, dimension(:), save :: aerosol
    real, allocatable, dimension(:), save :: tenpointzero
    real, allocatable, dimension(:), save :: ch4_up_flux
    real, allocatable, dimension(:), save :: h2dd_c
    real, allocatable, dimension(:), save :: h2dd_m
    real, allocatable, dimension(:), save :: h2dd_q
    real, dimension(r_tundra_spec), save :: r_tundra = [1200.0, 25000.0, 800.0, 3850.0, 1100.0]
    real, dimension(4), save :: ch4dd_tun = [- 4.757E-6, 4.0288E-3, - 1.13592, 106.636]
    real, dimension(3), save :: hno3dd_ice = [- 13.57, 6841.9, - 857410.6]
    real, dimension(5), save :: dif = [1.6, 1.6, 2.6, 1.9, 0.97]
    logical, dimension(row_length,rows,ntype) :: todo
    logical, dimension(row_length,rows) :: microb
    integer :: n_nosurf
    logical, save :: l_fix_improve_drydep
    logical, save :: l_fix_ukca_h2dd_x
    logical, save :: l_fix_drydep_so2_water
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SURFDDR'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first) then
      l_fix_improve_drydep = ukca_config%l_fix_improve_drydep
      l_fix_ukca_h2dd_x = ukca_config%l_fix_ukca_h2dd_x
      l_fix_drydep_so2_water = ukca_config%l_fix_drydep_so2_water
      brd_leaf = ukca_config%i_brd_leaf
      brd_leaf_dec = ukca_config%i_brd_leaf_dec
      brd_leaf_eg_trop = ukca_config%i_brd_leaf_eg_trop
      brd_leaf_eg_temp = ukca_config%i_brd_leaf_eg_temp
      ndl_leaf = ukca_config%i_ndl_leaf
      ndl_leaf_dec = ukca_config%i_ndl_leaf_dec
      ndl_leaf_eg = ukca_config%i_ndl_leaf_eg
      c3_grass = ukca_config%i_c3_grass
      c3_crop = ukca_config%i_c3_crop
      c3_pasture = ukca_config%i_c3_pasture
      c4_grass = ukca_config%i_c4_grass
      c4_crop = ukca_config%i_c4_crop
      c4_pasture = ukca_config%i_c4_pasture
      shrub = ukca_config%i_shrub
      shrub_dec = ukca_config%i_shrub_dec
      shrub_eg = ukca_config%i_shrub_eg
      urban = ukca_config%i_urban
      lake = ukca_config%i_lake
      soil = ukca_config%i_soil
      ice = ukca_config%i_ice
      if (ntype == 27) then
        if (SIZE(ukca_config%i_elev_ice) /= n_elev_ice) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(cmessage, '(A,I0,A,I0,A)') 'No of elevated ice surface types (', SIZE(ukca_config % i_elev_ice), ') does not match expected (', n_elev_ice, ')'
          errcode = 1
          call ereport(routinename, errcode, cmessage)
        end if
        elev_ice(:) = ukca_config%i_elev_ice
      end if
      ALLOCATE(rsurf(1:ntype,1:jpdd))
      ALLOCATE(r_null_tile(1:ntype))
      ALLOCATE(rooh(1:ntype))
      ALLOCATE(aerosol(1:ntype))
      ALLOCATE(tenpointzero(1:ntype))
      ALLOCATE(ch4_up_flux(1:ntype))
      ALLOCATE(h2dd_c(1:ntype))
      ALLOCATE(h2dd_m(1:ntype))
      ALLOCATE(h2dd_q(1:ntype))
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx = LBOUND(rsurf, dim=2), UBOUND(rsurf, dim=2), 1
        do idx_1 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
          rsurf(idx_1,idx) = r_null
        enddo
      enddo
      !$acc end kernels
      !$acc kernels
      !$acc loop gang vector independent
      do idx_2 = LBOUND(r_null_tile, dim=1), UBOUND(r_null_tile, dim=1), 1
        r_null_tile(idx_2) = r_null
      enddo
      !$acc end kernels
      !$acc kernels
      !$acc loop gang vector independent
      do idx_3 = LBOUND(tenpointzero, dim=1), UBOUND(tenpointzero, dim=1), 1
        tenpointzero(idx_3) = 10.0
      enddo
      !$acc end kernels
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
      if (ntype == 9) then
        if (brd_leaf /= 1 .OR. ndl_leaf /= 2 .OR. c3_grass /= 3 .OR. c4_grass /= 4 .OR. shrub /= 5 .OR. urban /= 6 .OR. lake /= 7 .OR. soil /= 8 .OR. ice /= 9) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') 'Nine tile order changed from standard setup'
          call umprint(ummessage, src=routinename)
          errcode = 1009
          call ereport(routinename, errcode, cmessage)
        end if
      else
        if (ntype == 13 .OR. (ntype == 17 .OR. ntype == 27)) then
          if (brd_leaf_dec /= 1 .OR. brd_leaf_eg_trop /= 2 .OR. brd_leaf_eg_temp /= 3 .OR. ndl_leaf_dec /= 4 .OR. ndl_leaf_eg /= 5 .OR. c3_grass /= 6) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') 'Tile order changed from standard setup'
            call umprint(ummessage, src=routinename)
            errcode = 1001
            call ereport(routinename, errcode, cmessage)
          end if
        else
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,I0)') 'ntype must be set to 9, 13, 17, 27.' // 'ntype = ', ntype
          call umprint(ummessage, src=routinename)
          errcode = 1002
          call ereport(routinename, errcode, cmessage)
        end if
      end if
      if (ntype == 13) then
        if (c4_grass /= 7 .OR. shrub_dec /= 8 .OR. shrub_eg /= 9 .OR. urban /= 10 .OR. lake /= 11 .OR. soil /= 12 .OR. ice /= 13) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') 'Thirteen tile order changed from standard setup'
          call umprint(ummessage, src=routinename)
          errcode = 1013
          call ereport(routinename, errcode, cmessage)
        end if
      else
        if (ntype == 17 .OR. ntype == 27) then
          if (c3_crop /= 7 .OR. c3_pasture /= 8 .OR. c4_grass /= 9 .OR. c4_crop /= 10 .OR. c4_pasture /= 11 .OR. shrub_dec /= 12 .OR. shrub_eg /= 13 .OR. urban /= 14 .OR. lake /= 15 .OR. soil /= 16 .OR. ice /= 17) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') 'Tile order changed from standard setup'
            call umprint(ummessage, src=routinename)
            errcode = 1727
            call ereport(routinename, errcode, cmessage)
          end if
        end if
      end if
      if (ntype == 27) then
        if (elev_ice(1) /= 18 .OR. elev_ice(2) /= 19 .OR. elev_ice(3) /= 20 .OR. elev_ice(4) /= 21 .OR. elev_ice(5) /= 22 .OR. elev_ice(6) /= 23 .OR. elev_ice(7) /= 24 .OR. elev_ice(8) /= 25 .OR. elev_ice(9) /= 26 .OR. elev_ice(10) /= 27) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') 'Tile order changed from standard setup'
          call umprint(ummessage, src=routinename)
          errcode = 1027
          call ereport(routinename, errcode, cmessage)
        end if
      end if
      if (ntype == 9) then
        if (l_fix_improve_drydep) then
          rooh(:) = [279.2, 238.2, 366.3, 322.9, 362.5, 424.9, 933.1, 585.4, 1156.1]
        else
          rooh(:) = [30.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        end if
        aerosol(:) = [r_null, r_null, r_null, r_null, r_null, r_null, 1000.0, r_null, 20000.0]
        if (l_fix_improve_drydep) then
          ch4_up_flux(:) = [39.5, 50.0, 30.0, 37.0, 27.5, 0.0, 0.0, 30.0, 0.0]
        else
          ch4_up_flux(:) = [39.5, 50.0, 30.0, 37.0, 27.5, 0.0, 0.0, 27.5, 0.0]
        end if
        if (l_fix_ukca_h2dd_x) then
          h2dd_c(:) = [19.70, 19.70, 17.70, 1.235, 1.000, 0.000, 0.000, 0.000, 0.000]
          h2dd_m(:) = [- 41.90, - 41.90, - 41.39, - 0.472, 0.000, 0.000, 0.000, - 0.000, 0.000]
          h2dd_q(:) = [0.00, 0.00, 0.00, 0.27, 0.00, 0.00, 0.00, 0.00, 0.00]
        else
          h2dd_c(:) = [0.00197, 0.00197, 0.00177, 1.23460, 0.00010, 0.00000, 0.00000, 0.00000, 0.00000]
          h2dd_m(:) = [- 0.00419, - 0.00419, - 0.00414, - 0.47200, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
          h2dd_q(:) = [0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27]
        end if
      else
        if (ntype == 13 .OR. (ntype == 17 .OR. ntype == 27)) then
          rooh(1:6) = [300.3, 270.3, 266.9, 238.0, 238.5, 366.3]
          aerosol(1:6) = [r_null, r_null, r_null, r_null, r_null, r_null]
          ch4_up_flux(1:6) = [39.5, 39.5, 39.5, 50.0, 50.0, 30.0]
          h2dd_c(1:6) = [19.70, 19.70, 19.70, 19.70, 19.70, 17.70]
          h2dd_m(1:6) = [- 41.90, - 41.90, - 41.90, - 41.90, - 41.90, - 41.39]
          if (l_fix_ukca_h2dd_x) then
            h2dd_q(1:6) = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
          else
            h2dd_q(1:6) = [0.27, 0.27, 0.27, 0.27, 0.27, 0.27]
          end if
        else
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,I0)') 'ntype = ', ntype
          call umprint(ummessage, src=routinename)
          errcode = 1
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
          call ereport(routinename, errcode, cmessage)
        end if
      end if
      if (ntype == 13) then
        rooh(7:13) = [322.9, 332.8, 392.2, 424.9, 933.1, 585.4, 1156.1]
        aerosol(7:13) = [r_null, r_null, r_null, r_null, 1000.0, r_null, 20000.0]
        if (l_fix_ukca_h2dd_x) then
          ch4_up_flux(7:13) = [37.0, 27.5, 27.5, 0.0, 0.0, 30.0, 0.0]
          h2dd_c(7:13) = [1.235, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
          h2dd_m(7:13) = [- 0.472, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
          h2dd_q(7:13) = [0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        else
          ch4_up_flux(7:13) = [37.0, 27.5, 27.5, 0.0, 0.0, 27.5, 0.0]
          h2dd_c(7:13) = [1.2346, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0]
          h2dd_m(7:13) = [- 0.472, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
          h2dd_q(7:13) = [0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27]
        end if
      else
        if (ntype == 17 .OR. ntype == 27) then
          rooh(7:17) = [366.3, 366.3, 322.9, 322.9, 322.9, 332.8, 392.2, 424.9, 933.1, 585.4, 1156.1]
          aerosol(7:17) = [r_null, r_null, r_null, r_null, r_null, r_null, r_null, r_null, 1000.0, r_null, 20000.0]
          if (l_fix_ukca_h2dd_x) then
            ch4_up_flux(7:17) = [30.0, 30.0, 37.0, 37.0, 37.0, 27.5, 27.5, 0.0, 0.0, 30.0, 0.0]
            h2dd_c(7:17) = [17.70, 17.70, 1.235, 1.235, 1.235, 1.000, 1.000, 0.000, 0.000, 0.000, 0.0]
            h2dd_m(7:17) = [- 41.39, - 41.39, - 0.472, - 0.472, - 0.472, 0.000, 0.000, 0.000, 0.000, - 41.39, 0.0]
            h2dd_q(7:17) = [0.00, 0.00, 0.27, 0.27, 0.27, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0]
          else
            ch4_up_flux(7:17) = [30.0, 30.0, 37.0, 37.0, 37.0, 27.5, 27.5, 0.0, 0.0, 27.5, 0.0]
            h2dd_c(7:17) = [0.00177, 0.00177, 1.23460, 1.23460, 1.23460, 0.00010, 0.00010, 0.00000, 0.00000, 0.00000, 0.0]
            h2dd_m(7:17) = [- 0.00414, - 0.00414, - 0.47200, - 0.47200, - 0.47200, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.0]
            h2dd_q(7:17) = [0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27]
          end if
        end if
      end if
      if (ntype == 27) then
        rooh(18:27) = [1156.1, 1156.1, 1156.1, 1156.1, 1156.1, 1156.1, 1156.1, 1156.1, 1156.1, 1156.1]
        aerosol(18:27) = [20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0]
        ch4_up_flux(18:27) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        h2dd_c(18:27) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        h2dd_m(18:27) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        if (l_fix_ukca_h2dd_x) then
          h2dd_q(18:27) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        else
          h2dd_q(18:27) = [0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27]
        end if
      end if
      if (ntype == 9) then
        n_nosurf = 0
        do n = 1, ndepd, 1
          if (speci(nldepd(n)) == 'O3        ' .OR. speci(nldepd(n)) == 'O3S       ') then
            if (l_fix_improve_drydep) then
              rsurf(:,n) = [219.3, 233.0, 355.0, 309.3, 358.2, 444.4, 2000.0, 645.2, 2000.0]
            else
              rsurf(:,n) = [200.0, 200.0, 200.0, 200.0, 400.0, 800.0, 2200.0, 800.0, 2500.0]
            end if
          else
            if (speci(nldepd(n)) == 'NO2       ' .OR. speci(nldepd(n)) == 'NO3       ') then
              if (l_fix_improve_drydep) then
                rsurf(:,n) = [364.1, 291.3, 443.8, 386.6, 447.8, 555.6, 2500.0, 806.5, 2500.0]
              else
                rsurf(:,n) = [225.0, 225.0, 400.0, 400.0, 600.0, 1200.0, 2600.0, 1200.0, 3500.0]
              end if
            else
              if (speci(nldepd(n)) == 'NO        ') then
                if (l_fix_improve_drydep) then
                  rsurf(:,n) = [2184.5, 1747.6, 3662.7, 2319.6, 2686.8, 3333.3, 15000.0, 4838.7, 15000.0]
                else
                  rsurf(:,n) = [1350.0, 1350.0, 2400.0, 2400.0, 3600.0, 72000.0, r_null, 72000.0, 21000.0]
                end if
              else
                if (speci(nldepd(n)) == 'HNO3      ' .OR. (speci(nldepd(n)) == 'HONO2     ' .OR. (speci(nldepd(n)) == 'B2ndry    ' .OR. (speci(nldepd(n)) == 'A2ndry    ' .OR. speci(nldepd(n)) == 'N2O5      ')))) then
                  if (l_fix_improve_drydep) then
                    rsurf(:,n) = [8.5, 8.4, 13.2, 12.0, 62.3, 12.8, 13.9, 16.0, 19.4]
                  else
                    !$acc kernels
                    !$acc loop gang vector independent
                    do idx_4 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                      rsurf(idx_4,n) = tenpointzero(idx_4 + (LBOUND(tenpointzero, dim=1) - LBOUND(rsurf, dim=1)))
                    enddo
                    !$acc end kernels
                  end if
                else
                  if (speci(nldepd(n)) == 'HNO4      ' .OR. speci(nldepd(n)) == 'HO2NO2    ') then
                    if (l_fix_improve_drydep) then
                      rsurf(:,n) = [17.0, 16.8, 26.4, 24.0, 24.9, 25.7, 27.7, 32.1, 38.8]
                    else
                      !$acc kernels
                      !$acc loop gang vector independent
                      do idx_5 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                        rsurf(idx_5,n) = tenpointzero(idx_5 + (LBOUND(tenpointzero, dim=1) - LBOUND(rsurf, dim=1)))
                      enddo
                      !$acc end kernels
                    end if
                  else
                    if (speci(nldepd(n)) == 'ISON      ' .OR. (speci(nldepd(n)) == 'HOC2H4NO3 ' .OR. (speci(nldepd(n)) == 'RTX24NO3  ' .OR. (speci(nldepd(n)) == 'RN9NO3    ' .OR. (speci(nldepd(n)) == 'RN12NO3   ' .OR. (speci(nldepd(n)) == 'RN15NO3   ' .OR. (speci(nldepd(n)) == 'RN18NO3   ' .OR. (speci(nldepd(n)) == 'RU14NO3   ' .OR. (speci(nldepd(n)) == 'RTN28NO3  ' .OR. (speci(nldepd(n)) == 'RTN25NO3  ' .OR. (speci(nldepd(n)) == 'RTX28NO3  ' .OR. (speci(nldepd(n)) == 'RTX22NO3  ' .OR. (speci(nldepd(n)) == 'RA22NO3   ' .OR. (speci(nldepd(n)) == 'RA25NO3   ' .OR. (speci(nldepd(n)) == 'RTN23NO3  ' .OR. (speci(nldepd(n)) == 'RU12NO3   ' .OR. (speci(nldepd(n)) == 'RU10NO3   ' .OR. (speci(nldepd(n)) == 'RA13NO3   ' .OR. (speci(nldepd(n)) == 'RA16NO3   ' .OR. speci(nldepd(n)) == 'RA19NO3   '))))))))))))))))))) then
                      if (l_fix_improve_drydep) then
                        rsurf(:,n) = [582.5, 466.0, 710.1, 618.6, 716.5, 888.9, 4000.0, 1290.3, 4000.0]
                      else
                        !$acc kernels
                        !$acc loop gang vector independent
                        do idx_6 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                          rsurf(idx_6,n) = tenpointzero(idx_6 + (LBOUND(tenpointzero, dim=1) - LBOUND(rsurf, dim=1)))
                        enddo
                        !$acc end kernels
                      end if
                    else
                      if (speci(nldepd(n)) == 'HCl       ' .OR. (speci(nldepd(n)) == 'HOCl      ' .OR. (speci(nldepd(n)) == 'HBr       ' .OR. speci(nldepd(n)) == 'HOBr      '))) then
                        if (l_fix_improve_drydep) then
                          rsurf(:,n) = [8.5, 8.4, 13.2, 12.0, 62.3, 12.8, 13.9, 16.0, 19.4]
                        else
                          !$acc kernels
                          !$acc loop gang vector independent
                          do idx_7 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                            rsurf(idx_7,n) = r_null_tile(idx_7 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                          enddo
                          !$acc end kernels
                        end if
                      else
                        if (speci(nldepd(n)) == 'H2SO4     ') then
                          if (l_fix_improve_drydep) then
                            rsurf(:,n) = [84.9, 83.8, 131.9, 120.0, 124.4, 128.5, 138.6, 160.4, 194.2]
                          else
                            !$acc kernels
                            !$acc loop gang vector independent
                            do idx_8 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                              rsurf(idx_8,n) = r_null_tile(idx_8 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                            enddo
                            !$acc end kernels
                          end if
                        else
                          if (speci(nldepd(n)) == 'H2O2      ' .OR. speci(nldepd(n)) == 'HOOH      ') then
                            if (l_fix_improve_drydep) then
                              rsurf(:,n) = [84.9, 83.8, 131.9, 120.0, 124.4, 128.5, 138.6, 160.4, 194.2]
                            else
                              !$acc kernels
                              !$acc loop gang vector independent
                              do idx_9 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                rsurf(idx_9,n) = tenpointzero(idx_9 + (LBOUND(tenpointzero, dim=1) - LBOUND(rsurf, dim=1)))
                              enddo
                              !$acc end kernels
                            end if
                          else
                            if (speci(nldepd(n)) == 'CH3OOH    ' .OR. (speci(nldepd(n)) == 'MeOOH     ' .OR. (speci(nldepd(n)) == 'C2H5OOH   ' .OR. (speci(nldepd(n)) == 'EtOOH     ' .OR. (speci(nldepd(n)) == 'n_C3H7OOH ' .OR. (speci(nldepd(n)) == 'i_C3H7OOH ' .OR. (speci(nldepd(n)) == 'n-PrOOH   ' .OR. (speci(nldepd(n)) == 'i-PrOOH   ' .OR. (speci(nldepd(n)) == 'MeCOCH2OOH' .OR. (speci(nldepd(n)) == 'ISOOH     ' .OR. (speci(nldepd(n)) == 'MACROOH   ' .OR. (speci(nldepd(n)) == 'MeCO3H    ' .OR. (speci(nldepd(n)) == 'MeCO2H    ' .OR. (speci(nldepd(n)) == 'HCOOH     ' .OR. (speci(nldepd(n)) == 'PropeOOH  ' .OR. (speci(nldepd(n)) == 'MEKOOH    ' .OR. (speci(nldepd(n)) == 'ALKAOOH   ' .OR. (speci(nldepd(n)) == 'AROMOOH   ' .OR. (speci(nldepd(n)) == 'BSVOC1    ' .OR. (speci(nldepd(n)) == 'BSVOC2    ' .OR. (speci(nldepd(n)) == 'ASVOC1    ' .OR. (speci(nldepd(n)) == 'ASVOC2    ' .OR. (speci(nldepd(n)) == 'ISOSVOC1  ' .OR. (speci(nldepd(n)) == 'ISOSVOC2  ' .OR. (speci(nldepd(n)) == 's-BuOOH   ' .OR. (speci(nldepd(n)) == 'MVKOOH    ' .OR. (speci(nldepd(n)) == 'HACET     ' .OR. (speci(nldepd(n)) == 'EtCO3H    ' .OR. (speci(nldepd(n)) == 'HOCH2CO3H ' .OR. (speci(nldepd(n)) == 'RCOOH25   ' .OR. (speci(nldepd(n)) == 'HOC2H4OOH ' .OR. (speci(nldepd(n)) == 'RN10OOH   ' .OR. (speci(nldepd(n)) == 'RN13OOH   ' .OR. (speci(nldepd(n)) == 'RN16OOH   ' .OR. (speci(nldepd(n)) == 'RN19OOH   ' .OR. (speci(nldepd(n)) == 'RN8OOH    ' .OR. (speci(nldepd(n)) == 'RN11OOH   ' .OR. (speci(nldepd(n)) == 'RN14OOH   ' .OR. (speci(nldepd(n)) == 'RN17OOH   ' .OR. (speci(nldepd(n)) == 'RU14OOH   ' .OR. (speci(nldepd(n)) == 'RU12OOH   ' .OR. (speci(nldepd(n)) == 'RU10OOH   ' .OR. (speci(nldepd(n)) == 'NRU14OOH  ' .OR. (speci(nldepd(n)) == 'NRU12OOH  ' .OR. (speci(nldepd(n)) == 'RN9OOH    ' .OR. (speci(nldepd(n)) == 'RN12OOH   ' .OR. (speci(nldepd(n)) == 'RN15OOH   ' .OR. (speci(nldepd(n)) == 'RN18OOH   ' .OR. (speci(nldepd(n)) == 'NRN6OOH   ' .OR. (speci(nldepd(n)) == 'NRN9OOH   ' .OR. (speci(nldepd(n)) == 'NRN12OOH  ' .OR. (speci(nldepd(n)) == 'RA13OOH   ' .OR. (speci(nldepd(n)) == 'RA16OOH   ' .OR. (speci(nldepd(n)) == 'RA19OOH   ' .OR. (speci(nldepd(n)) == 'RTN28OOH  ' .OR. (speci(nldepd(n)) == 'NRTN28OOH ' .OR. (speci(nldepd(n)) == 'RTN26OOH  ' .OR. (speci(nldepd(n)) == 'RTN25OOH  ' .OR. (speci(nldepd(n)) == 'RTN24OOH  ' .OR. (speci(nldepd(n)) == 'RTN23OOH  ' .OR. (speci(nldepd(n)) == 'RTN14OOH  ' .OR. (speci(nldepd(n)) == 'RTN10OOH  ' .OR. (speci(nldepd(n)) == 'RTX28OOH  ' .OR. (speci(nldepd(n)) == 'RTX24OOH  ' .OR. (speci(nldepd(n)) == 'RTX22OOH  ' .OR. (speci(nldepd(n)) == 'NRTX28OOH ' .OR. (speci(nldepd(n)) == 'CARB7     ' .OR. (speci(nldepd(n)) == 'CARB10    ' .OR. (speci(nldepd(n)) == 'CARB13    ' .OR. (speci(nldepd(n)) == 'CARB16    ' .OR. (speci(nldepd(n)) == 'RA22OOH   ' .OR. (speci(nldepd(n)) == 'RA25OOH   ' .OR. (speci(nldepd(n)) == 'HPUCARB12 ' .OR. (speci(nldepd(n)) == 'DHPCARB9  ' .OR. speci(nldepd(n)) == 'DHPR12OOH')))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))) then
                              !$acc kernels
                              !$acc loop gang vector independent
                              do idx_10 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                rsurf(idx_10,n) = rooh(idx_10 + (LBOUND(rooh, dim=1) - LBOUND(rsurf, dim=1)))
                              enddo
                              !$acc end kernels
                            else
                              if (speci(nldepd(n)) == 'PAN       ' .OR. (speci(nldepd(n)) == 'PPAN      ' .OR. (speci(nldepd(n)) == 'PHAN      ' .OR. (speci(nldepd(n)) == 'RU12PAN   ' .OR. speci(nldepd(n)) == 'RTN26PAN  ')))) then
                                if (l_fix_improve_drydep) then
                                  rsurf(:,n) = [485.4, 388.4, 591.7, 515.5, 597.1, 740.7, 3333.3, 1075.3, 3333.3]
                                else
                                  rsurf(:,n) = [500.0, 500.0, 500.0, 500.0, 500.0, r_null, 12500.0, 500.0, 12500.0]
                                end if
                              else
                                if (speci(nldepd(n)) == 'MPAN      ') then
                                  if (l_fix_improve_drydep) then
                                    rsurf(:,n) = [970.9, 776.7, 1183.4, 1030.9, 1194.2, 1481.5, 6666.7, 2150.5, 6666.7]
                                  else
                                    rsurf(:,n) = [500.0, 500.0, 500.0, 500.0, 500.0, r_null, 12500.0, 500.0, 12500.0]
                                  end if
                                else
                                  if (speci(nldepd(n)) == 'OnitU     ') then
                                    if (l_fix_improve_drydep) then
                                      rsurf(:,n) = [582.5, 466.0, 710.1, 618.6, 716.5, 888.9, 4000.0, 1290.3, 4000.0]
                                    else
                                      rsurf(:,n) = [500.0, 500.0, 500.0, 500.0, 500.0, r_null, 12500.0, 500.0, 12500.0]
                                    end if
                                  else
                                    if (speci(nldepd(n)) == 'NH3       ') then
                                      if (l_fix_improve_drydep) then
                                        rsurf(:,n) = [120.0, 130.9, 209.8, 196.1, 191.0, 180.7, 148.9, 213.5, 215.1]
                                      else
                                        !$acc kernels
                                        !$acc loop gang vector independent
                                        do idx_11 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                          rsurf(idx_11,n) = tenpointzero(idx_11 + (LBOUND(tenpointzero, dim=1) - LBOUND(rsurf, dim=1)))
                                        enddo
                                        !$acc end kernels
                                      end if
                                    else
                                      if (speci(nldepd(n)) == 'CO        ') then
                                        if (l_fix_improve_drydep) then
                                          rsurf(:,n) = [3700.0, 7300.0, 4550.0, 1960.0, 4550.0, r_null, r_null, 4550.0, r_null]
                                        else
                                          rsurf(:,n) = [3700.0, 7300.0, 4550.0, 1960.0, 4550.0, r_null, r_null, 4550.0, r_null]
                                        end if
                                      else
                                        if (speci(nldepd(n)) == 'CH4       ') then
                                          if (ukca_config%l_ukca_emsdrvn_ch4) then
                                            !$acc kernels
                                            !$acc loop gang vector independent
                                            do idx_12 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                              rsurf(idx_12,n) = 1.0 / r_null_tile(idx_12 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                            enddo
                                            !$acc end kernels
                                          else
                                            !$acc kernels
                                            !$acc loop gang vector independent
                                            do idx_13 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                              rsurf(idx_13,n) = r_null_tile(idx_13 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                            enddo
                                            !$acc end kernels
                                          end if
                                        else
                                          if (speci(nldepd(n)) == 'HONO      ') then
                                            !$acc kernels
                                            !$acc loop gang vector independent
                                            do idx_14 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                              rsurf(idx_14,n) = r_null_tile(idx_14 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                            enddo
                                            !$acc end kernels
                                          else
                                            if (speci(nldepd(n)) == 'H2        ') then
                                              !$acc kernels
                                              !$acc loop gang vector independent
                                              do idx_15 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                rsurf(idx_15,n) = r_null_tile(idx_15 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                              enddo
                                              !$acc end kernels
                                            else
                                              if (speci(nldepd(n)) == 'SO2       ') then
                                                if (ukca_config%l_ukca_dry_dep_so2wet) then
                                                  rsurf(:,n) = [120.0, 130.9, 209.8, 196.1, 191.0, 180.7, 1.0, 213.5, 215.1]
                                                else
                                                  if (l_fix_improve_drydep) then
                                                    rsurf(:,n) = [120.0, 130.9, 209.8, 196.1, 191.0, 180.7, 10.0, 213.5, 215.1]
                                                  else
                                                    rsurf(:,n) = [100.0, 100.0, 150.0, 350.0, 400.0, 400.0, 10.0, 700.0, r_null]
                                                  end if
                                                end if
                                              else
                                                if (speci(nldepd(n)) == 'BSOA      ' .OR. (speci(nldepd(n)) == 'ASOA      ' .OR. speci(nldepd(n)) == 'ISOSOA    ')) then
                                                  !$acc kernels
                                                  !$acc loop gang vector independent
                                                  do idx_16 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                    rsurf(idx_16,n) = aerosol(idx_16 + (LBOUND(aerosol, dim=1) - LBOUND(rsurf, dim=1)))
                                                  enddo
                                                  !$acc end kernels
                                                else
                                                  if (speci(nldepd(n)) == 'ORGNIT    ') then
                                                    if (l_fix_improve_drydep) then
                                                      rsurf(:,n) = [582.5, 466.0, 710.1, 618.6, 716.5, 888.9, 4000.0, 1290.3, 4000.0]
                                                    else
                                                      !$acc kernels
                                                      !$acc loop gang vector independent
                                                      do idx_17 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                        rsurf(idx_17,n) = aerosol(idx_17 + (LBOUND(aerosol, dim=1) - LBOUND(rsurf, dim=1)))
                                                      enddo
                                                      !$acc end kernels
                                                    end if
                                                  else
                                                    if (speci(nldepd(n)) == 'Sec_Org   ' .OR. speci(nldepd(n)) == 'SEC_ORG_I ') then
                                                      if (l_fix_improve_drydep) then
                                                        !$acc kernels
                                                        !$acc loop gang vector independent
                                                        do idx_18 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                          rsurf(idx_18,n) = aerosol(idx_18 + (LBOUND(aerosol, dim=1) - LBOUND(rsurf, dim=1)))
                                                        enddo
                                                        !$acc end kernels
                                                      else
                                                        !$acc kernels
                                                        !$acc loop gang vector independent
                                                        do idx_19 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                          rsurf(idx_19,n) = r_null_tile(idx_19 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                        enddo
                                                        !$acc end kernels
                                                      end if
                                                    else
                                                      if (speci(nldepd(n)) == 'HCHO      ') then
                                                        if (l_fix_improve_drydep) then
                                                          rsurf(:,n) = [136.0, 143.5, 228.5, 211.6, 210.5, 205.1, 182.7, 246.5, 261.8]
                                                        else
                                                          rsurf(:,n) = [100.0, 100.0, 150.0, 350.0, 600.0, 400.0, 200.0, 700.0, 700.0]
                                                        end if
                                                      else
                                                        if (speci(nldepd(n)) == 'MeOH      ' .OR. (speci(nldepd(n)) == 'EtOH      ' .OR. (speci(nldepd(n)) == 'i-PrOH    ' .OR. (speci(nldepd(n)) == 'n-PrOH    ' .OR. (speci(nldepd(n)) == 'AROH14    ' .OR. (speci(nldepd(n)) == 'ARNOH14   ' .OR. (speci(nldepd(n)) == 'AROH17    ' .OR. (speci(nldepd(n)) == 'ARNOH17   ' .OR. (speci(nldepd(n)) == 'IEPOX     ' .OR. (speci(nldepd(n)) == 'HMML      ' .OR. (speci(nldepd(n)) == 'HUCARB9   ' .OR. speci(nldepd(n)) == 'DHCARB9   '))))))))))) then
                                                          if (l_fix_improve_drydep) then
                                                            rsurf(:,n) = [187.1, 199.4, 318.3, 295.6, 292.2, 282.1, 245.1, 337.2, 352.2]
                                                          else
                                                            !$acc kernels
                                                            !$acc loop gang vector independent
                                                            do idx_20 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                              rsurf(idx_20,n) = r_null_tile(idx_20 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                            enddo
                                                            !$acc end kernels
                                                          end if
                                                        else
                                                          if (speci(nldepd(n)) == 'MeCHO     ' .OR. (speci(nldepd(n)) == 'EtCHO     ' .OR. (speci(nldepd(n)) == 'MACR      ' .OR. (speci(nldepd(n)) == 'NALD      ' .OR. (speci(nldepd(n)) == 'HOCH2CHO  ' .OR. (speci(nldepd(n)) == 'CARB14    ' .OR. (speci(nldepd(n)) == 'CARB17    ' .OR. (speci(nldepd(n)) == 'CARB11A   ' .OR. (speci(nldepd(n)) == 'UCARB10   ' .OR. (speci(nldepd(n)) == 'CARB15    ' .OR. (speci(nldepd(n)) == 'UCARB12   ' .OR. (speci(nldepd(n)) == 'UDCARB8   ' .OR. (speci(nldepd(n)) == 'UDCARB11  ' .OR. (speci(nldepd(n)) == 'UDCARB14  ' .OR. (speci(nldepd(n)) == 'TNCARB26  ' .OR. (speci(nldepd(n)) == 'TNCARB10  ' .OR. (speci(nldepd(n)) == 'TNCARB12  ' .OR. (speci(nldepd(n)) == 'TNCARB11  ' .OR. (speci(nldepd(n)) == 'CCARB12   ' .OR. (speci(nldepd(n)) == 'TNCARB15  ' .OR. (speci(nldepd(n)) == 'TXCARB24  ' .OR. (speci(nldepd(n)) == 'TXCARB22  ' .OR. (speci(nldepd(n)) == 'UDCARB17  ' .OR. (speci(nldepd(n)) == 'NOA       ' .OR. (speci(nldepd(n)) == 'NUCARB12  ' .OR. speci(nldepd(n)) == 'ANHY      '))))))))))))))))))))))))) then
                                                            if (l_fix_improve_drydep) then
                                                              rsurf(:,n) = [5825.2, 4660.3, 7100.6, 6185.6, 7164.8, 8888.9, 40000.0, 12903.2, 40000.0]
                                                            else
                                                              rsurf(:,n) = [1200.0, 1200.0, 1200.0, 1200.0, 1000.0, 2400.0, r_null, r_null, r_null]
                                                            end if
                                                          else
                                                            if (speci(nldepd(n)) == 'MGLY      ' .OR. (speci(nldepd(n)) == 'CARB3     ' .OR. (speci(nldepd(n)) == 'CARB6     ' .OR. (speci(nldepd(n)) == 'CARB9     ' .OR. speci(nldepd(n)) == 'CARB12    ')))) then
                                                              if (l_fix_improve_drydep) then
                                                                rsurf(:,n) = [12001.2, 13086.3, 20979.0, 19607.8, 19091.9, 18072.3, 14888.3, 21352.3, 21505.4]
                                                              else
                                                                rsurf(:,n) = [1200.0, 1200.0, 1200.0, 1200.0, 1000.0, 2400.0, r_null, r_null, r_null]
                                                              end if
                                                            else
                                                              if (speci(nldepd(n)) == 'DMSO      ' .OR. (speci(nldepd(n)) == 'Monoterp  ' .OR. (speci(nldepd(n)) == 'APINENE   ' .OR. speci(nldepd(n)) == 'BPINENE   '))) then
                                                                !$acc kernels
                                                                !$acc loop gang vector independent
                                                                do idx_21 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                                  rsurf(idx_21,n) = r_null_tile(idx_21 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                                enddo
                                                                !$acc end kernels
                                                              else
                                                                n_nosurf = n_nosurf + 1
                                                                if (first .AND. printstatus == prstatus_diag) then
                                                                  ! PSyclone CodeBlock (unsupported code) reason:
                                                                  !  - Unsupported statement: Write_Stmt
                                                                  WRITE(umMessage, '(A)') 'Surface resistance values not set for ' // speci(nldepd(n))
                                                                  call umprint(ummessage, src=routinename)
                                                                end if
                                                              end if
                                                            end if
                                                          end if
                                                        end if
                                                      end if
                                                    end if
                                                  end if
                                                end if
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        enddo
        if (first .AND. n_nosurf > 0) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(cmessage, '(A,I0,A)') ' Surface resistance values not found for ', n_nosurf, ' species.'
          errcode = ABS(n_nosurf)
          call ereport(routinename, errcode, cmessage)
        end if
      else
        if (ntype == 13 .OR. (ntype == 17 .OR. ntype == 27)) then
          n_nosurf = 0
          do n = 1, ndepd, 1
            if (speci(nldepd(n)) == 'O3        ' .OR. speci(nldepd(n)) == 'O3S       ') then
              rsurf(1:6,n) = [307.7, 285.7, 280.4, 232.6, 233.5, 355.0]
            else
              if (speci(nldepd(n)) == 'SO2       ') then
                rsurf(1:6,n) = [137.0, 111.1, 111.9, 131.3, 130.4, 209.8]
              else
                if (speci(nldepd(n)) == 'NO2       ' .OR. speci(nldepd(n)) == 'NO3       ') then
                  rsurf(1:6,n) = [384.6, 357.1, 350.5, 290.7, 291.8, 443.8]
                else
                  if (speci(nldepd(n)) == 'NO        ') then
                    rsurf(1:6,n) = [2307.7, 2142.9, 2102.8, 1744.2, 1751.0, 3662.7]
                  else
                    if (speci(nldepd(n)) == 'HNO3      ' .OR. (speci(nldepd(n)) == 'HONO2     ' .OR. (speci(nldepd(n)) == 'B2ndry    ' .OR. (speci(nldepd(n)) == 'A2ndry    ' .OR. speci(nldepd(n)) == 'N2O5      ')))) then
                      rsurf(1:6,n) = [9.5, 8.0, 8.0, 8.4, 8.4, 13.2]
                    else
                      if (speci(nldepd(n)) == 'HCl       ' .OR. (speci(nldepd(n)) == 'HOCl      ' .OR. (speci(nldepd(n)) == 'HBr       ' .OR. speci(nldepd(n)) == 'HOBr      '))) then
                        if (l_fix_improve_drydep) then
                          rsurf(1:6,n) = [9.5, 8.0, 8.0, 8.4, 8.4, 13.2]
                        else
                          !$acc kernels
                          !$acc loop gang vector independent
                          do idx_22 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                            rsurf(idx_22,n) = r_null_tile(idx_22 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                          enddo
                          !$acc end kernels
                        end if
                      else
                        if (speci(nldepd(n)) == 'H2SO4     ') then
                          if (l_fix_improve_drydep) then
                            rsurf(1:6,n) = [94.8, 80.0, 80.0, 83.9, 83.7, 131.9]
                          else
                            !$acc kernels
                            !$acc loop gang vector independent
                            do idx_23 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                              rsurf(idx_23,n) = r_null_tile(idx_23 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                            enddo
                            !$acc end kernels
                          end if
                        else
                          if (speci(nldepd(n)) == 'HNO4      ' .OR. speci(nldepd(n)) == 'HO2NO2    ') then
                            rsurf(1:6,n) = [19.0, 16.0, 16.0, 16.8, 16.7, 26.4]
                          else
                            if (speci(nldepd(n)) == 'HONO      ') then
                              rsurf(1:6,n) = [47.4, 40.0, 40.0, 42.0, 41.8, 65.9]
                            else
                              if (speci(nldepd(n)) == 'H2O2      ' .OR. speci(nldepd(n)) == 'HOOH      ') then
                                rsurf(1:6,n) = [94.8, 80.0, 80.0, 83.9, 83.7, 131.9]
                              else
                                if (speci(nldepd(n)) == 'PAN       ' .OR. (speci(nldepd(n)) == 'PPAN      ' .OR. (speci(nldepd(n)) == 'PHAN      ' .OR. (speci(nldepd(n)) == 'RU12PAN   ' .OR. speci(nldepd(n)) == 'RTN26PAN  ')))) then
                                  rsurf(1:6,n) = [512.8, 476.2, 467.3, 387.6, 389.1, 591.7]
                                else
                                  if (speci(nldepd(n)) == 'MPAN      ') then
                                    rsurf(1:6,n) = [1025.6, 952.4, 934.6, 775.2, 778.2, 1183.4]
                                  else
                                    if (speci(nldepd(n)) == 'OnitU     ' .OR. (speci(nldepd(n)) == 'ISON      ' .OR. (speci(nldepd(n)) == 'ORGNIT    ' .OR. (speci(nldepd(n)) == 'HOC2H4NO3 ' .OR. (speci(nldepd(n)) == 'RN9NO3    ' .OR. (speci(nldepd(n)) == 'RN12NO3   ' .OR. (speci(nldepd(n)) == 'RTN25NO3  ' .OR. (speci(nldepd(n)) == 'RN15NO3   ' .OR. (speci(nldepd(n)) == 'RN18NO3   ' .OR. (speci(nldepd(n)) == 'RU14NO3   ' .OR. (speci(nldepd(n)) == 'RTN28NO3  ' .OR. (speci(nldepd(n)) == 'RTX28NO3  ' .OR. (speci(nldepd(n)) == 'RA22NO3   ' .OR. (speci(nldepd(n)) == 'RA25NO3   ' .OR. (speci(nldepd(n)) == 'RU12NO3   ' .OR. (speci(nldepd(n)) == 'RU10NO3   ' .OR. (speci(nldepd(n)) == 'RTX22NO3  ' .OR. (speci(nldepd(n)) == 'RTN23NO3  ' .OR. (speci(nldepd(n)) == 'RTX24NO3  ' .OR. (speci(nldepd(n)) == 'RA13NO3   ' .OR. (speci(nldepd(n)) == 'RA16NO3   ' .OR. speci(nldepd(n)) == 'RA19NO3   '))))))))))))))))))))) then
                                      rsurf(1:6,n) = [615.4, 571.4, 560.7, 465.1, 466.9, 710.1]
                                    else
                                      if (speci(nldepd(n)) == 'CH3OOH    ' .OR. (speci(nldepd(n)) == 'MeOOH     ' .OR. (speci(nldepd(n)) == 'C2H5OOH   ' .OR. (speci(nldepd(n)) == 'EtOOH     ' .OR. (speci(nldepd(n)) == 'n_C3H7OOH ' .OR. (speci(nldepd(n)) == 'i_C3H7OOH ' .OR. (speci(nldepd(n)) == 'n-PrOOH   ' .OR. (speci(nldepd(n)) == 'i-PrOOH   ' .OR. (speci(nldepd(n)) == 'MeCOCH2OOH' .OR. (speci(nldepd(n)) == 'ISOOH     ' .OR. (speci(nldepd(n)) == 'MACROOH   ' .OR. (speci(nldepd(n)) == 'MeCO3H    ' .OR. (speci(nldepd(n)) == 'MeCO2H    ' .OR. (speci(nldepd(n)) == 'HCOOH     ' .OR. (speci(nldepd(n)) == 'PropeOOH  ' .OR. (speci(nldepd(n)) == 'MEKOOH    ' .OR. (speci(nldepd(n)) == 'ALKAOOH   ' .OR. (speci(nldepd(n)) == 'AROMOOH   ' .OR. (speci(nldepd(n)) == 'BSVOC1    ' .OR. (speci(nldepd(n)) == 'BSVOC2    ' .OR. (speci(nldepd(n)) == 'ASVOC1    ' .OR. (speci(nldepd(n)) == 'ASVOC2    ' .OR. (speci(nldepd(n)) == 'ISOSVOC1  ' .OR. (speci(nldepd(n)) == 'ISOSVOC2  ' .OR. (speci(nldepd(n)) == 's-BuOOH   ' .OR. (speci(nldepd(n)) == 'MVKOOH    ' .OR. (speci(nldepd(n)) == 'HACET     ' .OR. (speci(nldepd(n)) == 'EtCO3H    ' .OR. (speci(nldepd(n)) == 'HOCH2CO3H ' .OR. (speci(nldepd(n)) == 'RCOOH25   ' .OR. (speci(nldepd(n)) == 'HOC2H4OOH ' .OR. (speci(nldepd(n)) == 'RN10OOH   ' .OR. (speci(nldepd(n)) == 'RN13OOH   ' .OR. (speci(nldepd(n)) == 'RN16OOH   ' .OR. (speci(nldepd(n)) == 'RN19OOH   ' .OR. (speci(nldepd(n)) == 'RN8OOH    ' .OR. (speci(nldepd(n)) == 'RN11OOH   ' .OR. (speci(nldepd(n)) == 'RN14OOH   ' .OR. (speci(nldepd(n)) == 'RN17OOH   ' .OR. (speci(nldepd(n)) == 'RU14OOH   ' .OR. (speci(nldepd(n)) == 'RU12OOH   ' .OR. (speci(nldepd(n)) == 'RU10OOH   ' .OR. (speci(nldepd(n)) == 'NRU14OOH  ' .OR. (speci(nldepd(n)) == 'NRU12OOH  ' .OR. (speci(nldepd(n)) == 'RN9OOH    ' .OR. (speci(nldepd(n)) == 'RN12OOH   ' .OR. (speci(nldepd(n)) == 'RN15OOH   ' .OR. (speci(nldepd(n)) == 'RN18OOH   ' .OR. (speci(nldepd(n)) == 'NRN6OOH   ' .OR. (speci(nldepd(n)) == 'NRN9OOH   ' .OR. (speci(nldepd(n)) == 'NRN12OOH  ' .OR. (speci(nldepd(n)) == 'RA13OOH   ' .OR. (speci(nldepd(n)) == 'RA16OOH   ' .OR. (speci(nldepd(n)) == 'RA19OOH   ' .OR. (speci(nldepd(n)) == 'RTN28OOH  ' .OR. (speci(nldepd(n)) == 'NRTN28OOH ' .OR. (speci(nldepd(n)) == 'RTN26OOH  ' .OR. (speci(nldepd(n)) == 'RTN25OOH  ' .OR. (speci(nldepd(n)) == 'RTN24OOH  ' .OR. (speci(nldepd(n)) == 'RTN23OOH  ' .OR. (speci(nldepd(n)) == 'RTN14OOH  ' .OR. (speci(nldepd(n)) == 'RTN10OOH  ' .OR. (speci(nldepd(n)) == 'RTX28OOH  ' .OR. (speci(nldepd(n)) == 'RTX24OOH  ' .OR. (speci(nldepd(n)) == 'RTX22OOH  ' .OR. (speci(nldepd(n)) == 'NRTX28OOH ' .OR. (speci(nldepd(n)) == 'CARB7     ' .OR. (speci(nldepd(n)) == 'CARB10    ' .OR. (speci(nldepd(n)) == 'CARB13    ' .OR. (speci(nldepd(n)) == 'CARB16    ' .OR. (speci(nldepd(n)) == 'RA22OOH   ' .OR. (speci(nldepd(n)) == 'RA25OOH   ' .OR. (speci(nldepd(n)) == 'HPUCARB12 ' .OR. (speci(nldepd(n)) == 'DHPCARB9  ' .OR. speci(nldepd(n)) == 'DHPR12OOH ')))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))) then
                                        !$acc kernels
                                        !$acc loop gang vector independent
                                        do idx_24 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                          rsurf(idx_24,n) = rooh(idx_24 + (LBOUND(rooh, dim=1) - LBOUND(rsurf, dim=1)))
                                        enddo
                                        !$acc end kernels
                                      else
                                        if (speci(nldepd(n)) == 'NH3       ') then
                                          rsurf(1:6,n) = [137.0, 111.1, 111.9, 131.3, 130.4, 209.8]
                                        else
                                          if (speci(nldepd(n)) == 'CO        ') then
                                            rsurf(1:6,n) = [3700.0, 3700.0, 3700.0, 7300.0, 7300.0, 4550.0]
                                          else
                                            if (speci(nldepd(n)) == 'CH4       ') then
                                              if (ukca_config%l_ukca_emsdrvn_ch4) then
                                                !$acc kernels
                                                !$acc loop gang vector independent
                                                do idx_25 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                  rsurf(idx_25,n) = 1.0 / r_null_tile(idx_25 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                enddo
                                                !$acc end kernels
                                              else
                                                !$acc kernels
                                                !$acc loop gang vector independent
                                                do idx_26 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                  rsurf(idx_26,n) = r_null_tile(idx_26 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                enddo
                                                !$acc end kernels
                                              end if
                                            else
                                              if (speci(nldepd(n)) == 'H2        ') then
                                                !$acc kernels
                                                !$acc loop gang vector independent
                                                do idx_27 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                  rsurf(idx_27,n) = r_null_tile(idx_27 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                enddo
                                                !$acc end kernels
                                              else
                                                if (speci(nldepd(n)) == 'HCHO      ') then
                                                  rsurf(1:6,n) = [154.1, 126.6, 127.2, 143.8, 143.1, 228.5]
                                                else
                                                  if (speci(nldepd(n)) == 'MeOH      ' .OR. (speci(nldepd(n)) == 'EtOH      ' .OR. (speci(nldepd(n)) == 'i-PrOH    ' .OR. (speci(nldepd(n)) == 'n-PrOH    ' .OR. (speci(nldepd(n)) == 'AROH14    ' .OR. (speci(nldepd(n)) == 'ARNOH14   ' .OR. (speci(nldepd(n)) == 'AROH17    ' .OR. (speci(nldepd(n)) == 'ARNOH17   ' .OR. (speci(nldepd(n)) == 'IEPOX     ' .OR. (speci(nldepd(n)) == 'HMML      ' .OR. (speci(nldepd(n)) == 'HUCARB9   ' .OR. speci(nldepd(n)) == 'DHCARB9   '))))))))))) then
                                                    if (l_fix_improve_drydep) then
                                                      rsurf(1:6,n) = [212.6, 173.9, 174.9, 200.0, 198.0, 318.3]
                                                    else
                                                      !$acc kernels
                                                      !$acc loop gang vector independent
                                                      do idx_28 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                        rsurf(idx_28,n) = r_null_tile(idx_28 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                      enddo
                                                      !$acc end kernels
                                                    end if
                                                  else
                                                    if (speci(nldepd(n)) == 'MeCHO     ' .OR. (speci(nldepd(n)) == 'EtCHO     ' .OR. (speci(nldepd(n)) == 'MACR      ' .OR. (speci(nldepd(n)) == 'NALD      ' .OR. (speci(nldepd(n)) == 'HOCH2CHO  ' .OR. (speci(nldepd(n)) == 'CARB14    ' .OR. (speci(nldepd(n)) == 'CARB17    ' .OR. (speci(nldepd(n)) == 'CARB11A   ' .OR. (speci(nldepd(n)) == 'UCARB10   ' .OR. (speci(nldepd(n)) == 'UCARB12   ' .OR. (speci(nldepd(n)) == 'UDCARB8   ' .OR. (speci(nldepd(n)) == 'UDCARB11  ' .OR. (speci(nldepd(n)) == 'UDCARB14  ' .OR. (speci(nldepd(n)) == 'TNCARB26  ' .OR. (speci(nldepd(n)) == 'TNCARB10  ' .OR. (speci(nldepd(n)) == 'TNCARB12  ' .OR. (speci(nldepd(n)) == 'TNCARB11  ' .OR. (speci(nldepd(n)) == 'CCARB12   ' .OR. (speci(nldepd(n)) == 'TNCARB15  ' .OR. (speci(nldepd(n)) == 'TXCARB24  ' .OR. (speci(nldepd(n)) == 'TXCARB22  ' .OR. (speci(nldepd(n)) == 'UDCARB17  ' .OR. (speci(nldepd(n)) == 'NOA       ' .OR. (speci(nldepd(n)) == 'NUCARB12  ' .OR. speci(nldepd(n)) == 'ANHY      ')))))))))))))))))))))))) then
                                                      rsurf(1:6,n) = [6153.8, 5714.3, 5607.5, 4651.2, 4669.3, 7100.6]
                                                    else
                                                      if (speci(nldepd(n)) == 'MGLY      ' .OR. (speci(nldepd(n)) == 'CARB3     ' .OR. (speci(nldepd(n)) == 'CARB6     ' .OR. (speci(nldepd(n)) == 'CARB9     ' .OR. (speci(nldepd(n)) == 'CARB12    ' .OR. speci(nldepd(n)) == 'CARB15    '))))) then
                                                        rsurf(1:6,n) = [13698.6, 11111.1, 11194.0, 13129.1, 13043.5, 20979.0]
                                                      else
                                                        if (speci(nldepd(n)) == 'BSOA      ' .OR. (speci(nldepd(n)) == 'ASOA      ' .OR. speci(nldepd(n)) == 'ISOSOA    ')) then
                                                          !$acc kernels
                                                          !$acc loop gang vector independent
                                                          do idx_29 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                            rsurf(idx_29,n) = aerosol(idx_29 + (LBOUND(aerosol, dim=1) - LBOUND(rsurf, dim=1)))
                                                          enddo
                                                          !$acc end kernels
                                                        else
                                                          if (speci(nldepd(n)) == 'Sec_Org   ' .OR. speci(nldepd(n)) == 'SEC_ORG_I ') then
                                                            if (l_fix_improve_drydep) then
                                                              !$acc kernels
                                                              !$acc loop gang vector independent
                                                              do idx_30 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                                rsurf(idx_30,n) = aerosol(idx_30 + (LBOUND(aerosol, dim=1) - LBOUND(rsurf, dim=1)))
                                                              enddo
                                                              !$acc end kernels
                                                            else
                                                              !$acc kernels
                                                              !$acc loop gang vector independent
                                                              do idx_31 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                                rsurf(idx_31,n) = r_null_tile(idx_31 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                              enddo
                                                              !$acc end kernels
                                                            end if
                                                          else
                                                            if (speci(nldepd(n)) == 'DMSO      ' .OR. (speci(nldepd(n)) == 'Monoterp  ' .OR. (speci(nldepd(n)) == 'APINENE   ' .OR. speci(nldepd(n)) == 'BPINENE   '))) then
                                                              !$acc kernels
                                                              !$acc loop gang vector independent
                                                              do idx_32 = LBOUND(rsurf, dim=1), UBOUND(rsurf, dim=1), 1
                                                                rsurf(idx_32,n) = r_null_tile(idx_32 + (LBOUND(r_null_tile, dim=1) - LBOUND(rsurf, dim=1)))
                                                              enddo
                                                              !$acc end kernels
                                                            else
                                                              n_nosurf = n_nosurf + 1
                                                              if (first .AND. printstatus == prstatus_diag) then
                                                                ummessage(:) = ' Surface resistance values not set for ' // speci(nldepd(n))
                                                                call umprint(ummessage, src=routinename)
                                                              end if
                                                            end if
                                                          end if
                                                        end if
                                                      end if
                                                    end if
                                                  end if
                                                end if
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          enddo
          if (first .AND. n_nosurf > 0) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(cmessage, '(A,I0,A)') ' Surface resistance values not found for ', n_nosurf, ' species.'
            errcode = ABS(n_nosurf)
            call ereport(routinename, errcode, cmessage)
          end if
        else
          errcode = 319
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
          call ereport(routinename, errcode, cmessage)
        end if
      end if
      if (ntype == 13) then
        do n = 1, ndepd, 1
          if (speci(nldepd(n)) == 'O3        ' .OR. speci(nldepd(n)) == 'O3S       ') then
            rsurf(7:13,n) = [309.3, 324.3, 392.2, 444.4, 2000.0, 645.2, 2000.0]
          else
            if (speci(nldepd(n)) == 'SO2       ') then
              if (ukca_config%l_ukca_dry_dep_so2wet) then
                rsurf(7:13,n) = [196.1, 185.8, 196.1, 180.7, 1.0, 213.5, 215.1]
              else
                if (l_fix_drydep_so2_water) then
                  rsurf(7:13,n) = [196.1, 185.8, 196.1, 180.7, 10.0, 213.5, 215.1]
                else
                  rsurf(7:13,n) = [196.1, 185.8, 196.1, 180.7, 148.9, 213.5, 215.1]
                end if
              end if
            else
              if (speci(nldepd(n)) == 'NO2       ' .OR. speci(nldepd(n)) == 'NO3       ') then
                rsurf(7:13,n) = [386.6, 405.4, 490.2, 555.6, 2500.0, 806.5, 2500.0]
              else
                if (speci(nldepd(n)) == 'NO        ') then
                  rsurf(7:13,n) = [2319.6, 2432.4, 2941.2, 3333.3, 15000.0, 4838.7, 15000.0]
                else
                  if (speci(nldepd(n)) == 'HNO3      ' .OR. (speci(nldepd(n)) == 'HONO2     ' .OR. (speci(nldepd(n)) == 'B2ndry    ' .OR. (speci(nldepd(n)) == 'A2ndry    ' .OR. speci(nldepd(n)) == 'N2O5      ')))) then
                    rsurf(7:13,n) = [12.0, 59.1, 65.4, 12.8, 13.9, 16.0, 19.4]
                  else
                    if (speci(nldepd(n)) == 'HCl       ' .OR. (speci(nldepd(n)) == 'HOCl      ' .OR. (speci(nldepd(n)) == 'HBr       ' .OR. speci(nldepd(n)) == 'HOBr      '))) then
                      if (l_fix_improve_drydep) then
                        rsurf(7:13,n) = [12.0, 59.1, 65.4, 12.8, 13.9, 16.0, 19.4]
                      end if
                    else
                      if (speci(nldepd(n)) == 'H2SO4     ') then
                        if (l_fix_improve_drydep) then
                          rsurf(7:13,n) = [120.0, 118.1, 130.7, 128.5, 138.6, 160.4, 194.2]
                        end if
                      else
                        if (speci(nldepd(n)) == 'HNO4      ' .OR. speci(nldepd(n)) == 'HO2NO2    ') then
                          rsurf(7:13,n) = [24.0, 23.6, 26.1, 25.7, 27.7, 32.1, 38.8]
                        else
                          if (speci(nldepd(n)) == 'HONO      ') then
                            rsurf(7:13,n) = [60.0, 59.1, 65.4, 64.2, 69.3, 80.2, 97.1]
                          else
                            if (speci(nldepd(n)) == 'H2O2      ' .OR. speci(nldepd(n)) == 'HOOH      ') then
                              rsurf(7:13,n) = [120.0, 118.1, 130.7, 128.5, 138.6, 160.4, 194.2]
                            else
                              if (speci(nldepd(n)) == 'PAN       ' .OR. (speci(nldepd(n)) == 'PPAN      ' .OR. (speci(nldepd(n)) == 'PHAN      ' .OR. (speci(nldepd(n)) == 'RU12PAN   ' .OR. speci(nldepd(n)) == 'RTN26PAN  ')))) then
                                rsurf(7:13,n) = [515.5, 540.5, 653.6, 740.7, 3333.3, 1075.3, 3333.3]
                              else
                                if (speci(nldepd(n)) == 'MPAN      ') then
                                  rsurf(7:13,n) = [1030.9, 1081.1, 1307.2, 1481.5, 6666.7, 2150.5, 6666.7]
                                else
                                  if (speci(nldepd(n)) == 'OnitU     ' .OR. (speci(nldepd(n)) == 'ISON      ' .OR. (speci(nldepd(n)) == 'ORGNIT    ' .OR. (speci(nldepd(n)) == 'HOC2H4NO3 ' .OR. (speci(nldepd(n)) == 'RTX24NO3  ' .OR. (speci(nldepd(n)) == 'RN9NO3    ' .OR. (speci(nldepd(n)) == 'RN12NO3   ' .OR. (speci(nldepd(n)) == 'RN15NO3   ' .OR. (speci(nldepd(n)) == 'RN18NO3   ' .OR. (speci(nldepd(n)) == 'RU14NO3   ' .OR. (speci(nldepd(n)) == 'RTN28NO3  ' .OR. (speci(nldepd(n)) == 'RTN25NO3  ' .OR. (speci(nldepd(n)) == 'RTX28NO3  ' .OR. (speci(nldepd(n)) == 'RTX22NO3  ' .OR. (speci(nldepd(n)) == 'RA22NO3   ' .OR. (speci(nldepd(n)) == 'RA25NO3   ' .OR. (speci(nldepd(n)) == 'RTN23NO3  ' .OR. (speci(nldepd(n)) == 'RU12NO3   ' .OR. (speci(nldepd(n)) == 'RU10NO3   ' .OR. (speci(nldepd(n)) == 'RA13NO3   ' .OR. (speci(nldepd(n)) == 'RA16NO3   ' .OR. speci(nldepd(n)) == 'RA19NO3   '))))))))))))))))))))) then
                                    rsurf(7:13,n) = [618.6, 648.6, 784.3, 888.9, 4000.0, 1290.3, 4000.0]
                                  else
                                    if (speci(nldepd(n)) == 'NH3       ') then
                                      rsurf(7:13,n) = [196.1, 185.8, 196.1, 180.7, 148.9, 213.5, 215.1]
                                    else
                                      if (speci(nldepd(n)) == 'CO        ') then
                                        rsurf(7:13,n) = [1960.0, 4550.0, 4550.0, r_null, r_null, 4550.0, r_null]
                                      else
                                        if (speci(nldepd(n)) == 'HCHO      ') then
                                          rsurf(7:13,n) = [211.6, 203.1, 217.9, 205.1, 182.7, 246.5, 261.8]
                                        else
                                          if (speci(nldepd(n)) == 'MeOH      ' .OR. (speci(nldepd(n)) == 'EtOH      ' .OR. (speci(nldepd(n)) == 'i-PrOH    ' .OR. (speci(nldepd(n)) == 'n-PrOH    ' .OR. (speci(nldepd(n)) == 'AROH14    ' .OR. (speci(nldepd(n)) == 'ARNOH14   ' .OR. (speci(nldepd(n)) == 'AROH17    ' .OR. (speci(nldepd(n)) == 'ARNOH17   ' .OR. (speci(nldepd(n)) == 'IEPOX     ' .OR. (speci(nldepd(n)) == 'HMML      ' .OR. (speci(nldepd(n)) == 'HUCARB9   ' .OR. speci(nldepd(n)) == 'DHCARB9   '))))))))))) then
                                            if (l_fix_improve_drydep) then
                                              rsurf(7:13,n) = [295.6, 282.7, 301.7, 282.1, 245.1, 337.2, 352.2]
                                            end if
                                          else
                                            if (speci(nldepd(n)) == 'MeCHO     ' .OR. (speci(nldepd(n)) == 'EtCHO     ' .OR. (speci(nldepd(n)) == 'MACR      ' .OR. (speci(nldepd(n)) == 'NALD      ' .OR. (speci(nldepd(n)) == 'HOCH2CHO  ' .OR. (speci(nldepd(n)) == 'CARB14    ' .OR. (speci(nldepd(n)) == 'CARB17    ' .OR. (speci(nldepd(n)) == 'CARB11A   ' .OR. (speci(nldepd(n)) == 'UCARB10   ' .OR. (speci(nldepd(n)) == 'UCARB12   ' .OR. (speci(nldepd(n)) == 'UDCARB8   ' .OR. (speci(nldepd(n)) == 'UDCARB11  ' .OR. (speci(nldepd(n)) == 'UDCARB14  ' .OR. (speci(nldepd(n)) == 'TNCARB26  ' .OR. (speci(nldepd(n)) == 'TNCARB10  ' .OR. (speci(nldepd(n)) == 'TNCARB12  ' .OR. (speci(nldepd(n)) == 'TNCARB11  ' .OR. (speci(nldepd(n)) == 'CCARB12   ' .OR. (speci(nldepd(n)) == 'TNCARB15  ' .OR. (speci(nldepd(n)) == 'TXCARB24  ' .OR. (speci(nldepd(n)) == 'TXCARB22  ' .OR. (speci(nldepd(n)) == 'UDCARB17  ' .OR. (speci(nldepd(n)) == 'NOA       ' .OR. (speci(nldepd(n)) == 'NUCARB12  ' .OR. speci(nldepd(n)) == 'ANHY      ')))))))))))))))))))))))) then
                                              rsurf(7:13,n) = [6185.6, 6486.5, 7843.1, 8888.9, 40000.0, 12903.2, 40000.0]
                                            else
                                              if (speci(nldepd(n)) == 'MGLY      ' .OR. (speci(nldepd(n)) == 'CARB3     ' .OR. (speci(nldepd(n)) == 'CARB6     ' .OR. (speci(nldepd(n)) == 'CARB9     ' .OR. (speci(nldepd(n)) == 'CARB12    ' .OR. speci(nldepd(n)) == 'CARB15    '))))) then
                                                rsurf(7:13,n) = [19607.8, 18575.9, 19607.8, 18072.3, 14888.3, 21352.3, 21505.4]
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        enddo
      else
        if (ntype == 17 .OR. ntype == 27) then
          do n = 1, ndepd, 1
            if (speci(nldepd(n)) == 'O3        ' .OR. speci(nldepd(n)) == 'O3S       ') then
              rsurf(7:17,n) = [355.0, 355.0, 309.3, 309.3, 309.3, 324.3, 392.2, 444.4, 2000.0, 645.2, 2000.0]
            else
              if (speci(nldepd(n)) == 'SO2       ') then
                if (ukca_config%l_ukca_dry_dep_so2wet) then
                  rsurf(7:17,n) = [30.0, 209.8, 196.1, 30.0, 196.1, 185.8, 196.1, 180.7, 1.0, 213.5, 215.1]
                else
                  if (l_fix_drydep_so2_water) then
                    rsurf(7:17,n) = [209.8, 209.8, 196.1, 196.1, 196.1, 185.8, 196.1, 180.7, 10.0, 213.5, 215.1]
                  else
                    rsurf(7:17,n) = [209.8, 209.8, 196.1, 196.1, 196.1, 185.8, 196.1, 180.7, 148.9, 213.5, 215.1]
                  end if
                end if
              else
                if (speci(nldepd(n)) == 'NO2       ' .OR. speci(nldepd(n)) == 'NO3       ') then
                  rsurf(7:17,n) = [443.8, 443.8, 386.6, 386.6, 386.6, 405.4, 490.2, 555.6, 2500.0, 806.5, 2500.0]
                else
                  if (speci(nldepd(n)) == 'NO        ') then
                    rsurf(7:17,n) = [3662.7, 3662.7, 2319.6, 2319.6, 2319.6, 2432.4, 2941.2, 3333.3, 15000.0, 4838.7, 15000.0]
                  else
                    if (speci(nldepd(n)) == 'HNO3      ' .OR. (speci(nldepd(n)) == 'HONO2     ' .OR. (speci(nldepd(n)) == 'B2ndry    ' .OR. (speci(nldepd(n)) == 'A2ndry    ' .OR. speci(nldepd(n)) == 'N2O5      ')))) then
                      rsurf(7:17,n) = [13.2, 13.2, 12.0, 12.0, 12.0, 59.1, 65.4, 12.8, 13.9, 16.0, 19.4]
                    else
                      if (speci(nldepd(n)) == 'HCl       ' .OR. (speci(nldepd(n)) == 'HOCl      ' .OR. (speci(nldepd(n)) == 'HBr       ' .OR. speci(nldepd(n)) == 'HOBr      '))) then
                        if (l_fix_improve_drydep) then
                          rsurf(7:17,n) = [13.2, 13.2, 12.0, 12.0, 12.0, 59.1, 65.4, 12.8, 13.9, 16.0, 19.4]
                        end if
                      else
                        if (speci(nldepd(n)) == 'H2SO4     ') then
                          if (l_fix_improve_drydep) then
                            rsurf(7:17,n) = [131.9, 131.9, 120.0, 120.0, 120.0, 118.1, 130.7, 128.5, 138.6, 160.4, 194.2]
                          end if
                        else
                          if (speci(nldepd(n)) == 'HNO4      ' .OR. speci(nldepd(n)) == 'HO2NO2    ') then
                            rsurf(7:17,n) = [26.4, 26.4, 24.0, 24.0, 24.0, 23.6, 26.1, 25.7, 27.7, 32.1, 38.8]
                          else
                            if (speci(nldepd(n)) == 'HONO      ') then
                              rsurf(7:17,n) = [65.9, 65.9, 60.0, 60.0, 60.0, 59.1, 65.4, 64.2, 69.3, 80.2, 97.1]
                            else
                              if (speci(nldepd(n)) == 'H2O2      ' .OR. speci(nldepd(n)) == 'HOOH      ') then
                                rsurf(7:17,n) = [131.9, 131.9, 120.0, 120.0, 120.0, 118.1, 130.7, 128.5, 138.6, 160.4, 194.2]
                              else
                                if (speci(nldepd(n)) == 'PAN       ' .OR. (speci(nldepd(n)) == 'PPAN      ' .OR. (speci(nldepd(n)) == 'PHAN      ' .OR. (speci(nldepd(n)) == 'RU12PAN   ' .OR. speci(nldepd(n)) == 'RTN26PAN  ')))) then
                                  rsurf(7:17,n) = [591.7, 591.7, 515.5, 515.5, 515.5, 540.5, 653.6, 740.7, 3333.3, 1075.3, 3333.3]
                                else
                                  if (speci(nldepd(n)) == 'MPAN      ') then
                                    rsurf(7:17,n) = [1183.4, 1183.4, 1030.9, 1030.9, 1030.9, 1081.1, 1307.2, 1481.5, 6666.7, 2150.5, 6666.7]
                                  else
                                    if (speci(nldepd(n)) == 'OnitU     ' .OR. (speci(nldepd(n)) == 'ISON      ' .OR. (speci(nldepd(n)) == 'ORGNIT    ' .OR. (speci(nldepd(n)) == 'HOC2H4NO3 ' .OR. (speci(nldepd(n)) == 'RTX24NO3  ' .OR. (speci(nldepd(n)) == 'RN9NO3    ' .OR. (speci(nldepd(n)) == 'RN12NO3   ' .OR. (speci(nldepd(n)) == 'RN15NO3   ' .OR. (speci(nldepd(n)) == 'RN18NO3   ' .OR. (speci(nldepd(n)) == 'RU14NO3   ' .OR. (speci(nldepd(n)) == 'RTN28NO3  ' .OR. (speci(nldepd(n)) == 'RTN25NO3  ' .OR. (speci(nldepd(n)) == 'RTX28NO3  ' .OR. (speci(nldepd(n)) == 'RTX22NO3  ' .OR. (speci(nldepd(n)) == 'RA22NO3   ' .OR. (speci(nldepd(n)) == 'RA25NO3   ' .OR. (speci(nldepd(n)) == 'RTN23NO3  ' .OR. (speci(nldepd(n)) == 'RU12NO3   ' .OR. (speci(nldepd(n)) == 'RU10NO3   ' .OR. (speci(nldepd(n)) == 'RA13NO3   ' .OR. (speci(nldepd(n)) == 'RA16NO3   ' .OR. speci(nldepd(n)) == 'RA19NO3   '))))))))))))))))))))) then
                                      rsurf(7:17,n) = [710.1, 710.1, 618.6, 618.6, 618.6, 648.6, 784.3, 888.9, 4000.0, 1290.3, 4000.0]
                                    else
                                      if (speci(nldepd(n)) == 'NH3       ') then
                                        rsurf(7:17,n) = [209.8, 209.8, 196.1, 196.1, 196.1, 185.8, 196.1, 180.7, 148.9, 213.5, 215.1]
                                      else
                                        if (speci(nldepd(n)) == 'CO        ') then
                                          rsurf(7:17,n) = [4550.0, 4550.0, 1960.0, 1960.0, 1960.0, 4550.0, 4550.0, r_null, r_null, 4550.0, r_null]
                                        else
                                          if (speci(nldepd(n)) == 'HCHO      ') then
                                            rsurf(7:17,n) = [228.5, 228.5, 211.6, 211.6, 211.6, 203.1, 217.9, 205.1, 182.7, 246.5, 261.8]
                                          else
                                            if (speci(nldepd(n)) == 'MeOH      ' .OR. (speci(nldepd(n)) == 'EtOH      ' .OR. (speci(nldepd(n)) == 'i-PrOH    ' .OR. (speci(nldepd(n)) == 'n-PrOH    ' .OR. (speci(nldepd(n)) == 'AROH14    ' .OR. (speci(nldepd(n)) == 'ARNOH14   ' .OR. (speci(nldepd(n)) == 'AROH17    ' .OR. (speci(nldepd(n)) == 'ARNOH17   ' .OR. (speci(nldepd(n)) == 'IEPOX     ' .OR. (speci(nldepd(n)) == 'HMML      ' .OR. (speci(nldepd(n)) == 'HUCARB9   ' .OR. speci(nldepd(n)) == 'DHCARB9   '))))))))))) then
                                              if (l_fix_improve_drydep) then
                                                rsurf(7:17,n) = [318.3, 318.3, 295.6, 295.6, 295.6, 282.7, 301.7, 282.1, 245.1, 337.2, 352.2]
                                              end if
                                            else
                                              if (speci(nldepd(n)) == 'MeCHO     ' .OR. (speci(nldepd(n)) == 'EtCHO     ' .OR. (speci(nldepd(n)) == 'MACR      ' .OR. (speci(nldepd(n)) == 'NALD      ' .OR. (speci(nldepd(n)) == 'HOCH2CHO  ' .OR. (speci(nldepd(n)) == 'CARB14    ' .OR. (speci(nldepd(n)) == 'CARB17    ' .OR. (speci(nldepd(n)) == 'CARB11A   ' .OR. (speci(nldepd(n)) == 'UCARB10   ' .OR. (speci(nldepd(n)) == 'UCARB12   ' .OR. (speci(nldepd(n)) == 'UDCARB8   ' .OR. (speci(nldepd(n)) == 'UDCARB11  ' .OR. (speci(nldepd(n)) == 'UDCARB14  ' .OR. (speci(nldepd(n)) == 'TNCARB26  ' .OR. (speci(nldepd(n)) == 'TNCARB10  ' .OR. (speci(nldepd(n)) == 'TNCARB12  ' .OR. (speci(nldepd(n)) == 'TNCARB11  ' .OR. (speci(nldepd(n)) == 'CCARB12   ' .OR. (speci(nldepd(n)) == 'TNCARB15  ' .OR. (speci(nldepd(n)) == 'TXCARB24  ' .OR. (speci(nldepd(n)) == 'TXCARB22  ' .OR. (speci(nldepd(n)) == 'UDCARB17  ' .OR. (speci(nldepd(n)) == 'NOA       ' .OR. (speci(nldepd(n)) == 'NUCARB12  ' .OR. speci(nldepd(n)) == 'ANHY      ')))))))))))))))))))))))) then
                                                rsurf(7:17,n) = [7100.6, 7100.6, 6185.6, 6185.6, 6185.6, 6486.5, 7843.1, 8888.9, 40000.0, 12903.2, 40000.0]
                                              else
                                                if (speci(nldepd(n)) == 'MGLY      ' .OR. (speci(nldepd(n)) == 'CARB3     ' .OR. (speci(nldepd(n)) == 'CARB6     ' .OR. (speci(nldepd(n)) == 'CARB9     ' .OR. (speci(nldepd(n)) == 'CARB12    ' .OR. speci(nldepd(n)) == 'CARB15    '))))) then
                                                  rsurf(7:17,n) = [20979.0, 20979.0, 19607.8, 19607.8, 19607.8, 18575.9, 19607.8, 18072.3, 14888.3, 21352.3, 21505.4]
                                                end if
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          enddo
        end if
      end if
      if (ntype == 27) then
        do n = 1, ndepd, 1
          if (speci(nldepd(n)) == 'O3        ' .OR. speci(nldepd(n)) == 'O3S       ') then
            rsurf(18:27,n) = [2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0]
          else
            if (speci(nldepd(n)) == 'SO2       ') then
              rsurf(18:27,n) = [215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1]
            else
              if (speci(nldepd(n)) == 'NO2       ' .OR. speci(nldepd(n)) == 'NO3       ') then
                rsurf(18:27,n) = [2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0]
              else
                if (speci(nldepd(n)) == 'NO        ') then
                  rsurf(18:27,n) = [15000.0, 15000.0, 15000.0, 15000.0, 15000.0, 15000.0, 15000.0, 15000.0, 15000.0, 15000.0]
                else
                  if (speci(nldepd(n)) == 'HNO3      ' .OR. (speci(nldepd(n)) == 'HONO2     ' .OR. (speci(nldepd(n)) == 'B2ndry    ' .OR. (speci(nldepd(n)) == 'A2ndry    ' .OR. speci(nldepd(n)) == 'N2O5      ')))) then
                    rsurf(18:27,n) = [19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4]
                  else
                    if (speci(nldepd(n)) == 'HCl       ' .OR. (speci(nldepd(n)) == 'HOCl      ' .OR. (speci(nldepd(n)) == 'HBr       ' .OR. speci(nldepd(n)) == 'HOBr      '))) then
                      if (l_fix_improve_drydep) then
                        rsurf(18:27,n) = [19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4]
                      end if
                    else
                      if (speci(nldepd(n)) == 'H2SO4     ') then
                        if (l_fix_improve_drydep) then
                          rsurf(18:27,n) = [194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2]
                        end if
                      else
                        if (speci(nldepd(n)) == 'HNO4      ' .OR. speci(nldepd(n)) == 'HO2NO2    ') then
                          rsurf(18:27,n) = [38.8, 38.8, 38.8, 38.8, 38.8, 38.8, 38.8, 38.8, 38.8, 38.8]
                        else
                          if (speci(nldepd(n)) == 'HONO      ') then
                            rsurf(18:27,n) = [97.1, 97.1, 97.1, 97.1, 97.1, 97.1, 97.1, 97.1, 97.1, 97.1]
                          else
                            if (speci(nldepd(n)) == 'H2O2      ' .OR. speci(nldepd(n)) == 'HOOH      ') then
                              rsurf(18:27,n) = [194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2, 194.2]
                            else
                              if (speci(nldepd(n)) == 'PAN       ' .OR. (speci(nldepd(n)) == 'PPAN      ' .OR. (speci(nldepd(n)) == 'PHAN      ' .OR. (speci(nldepd(n)) == 'RU12PAN   ' .OR. speci(nldepd(n)) == 'RTN26PAN  ')))) then
                                rsurf(18:27,n) = [3333.3, 3333.3, 3333.3, 3333.3, 3333.3, 3333.3, 3333.3, 3333.3, 3333.3, 3333.3]
                              else
                                if (speci(nldepd(n)) == 'MPAN      ') then
                                  rsurf(18:27,n) = [6666.7, 6666.7, 6666.7, 6666.7, 6666.7, 6666.7, 6666.7, 6666.7, 6666.7, 6666.7]
                                else
                                  if (speci(nldepd(n)) == 'OnitU     ' .OR. (speci(nldepd(n)) == 'ISON      ' .OR. (speci(nldepd(n)) == 'ORGNIT    ' .OR. (speci(nldepd(n)) == 'HOC2H4NO3 ' .OR. (speci(nldepd(n)) == 'RTX24NO3  ' .OR. (speci(nldepd(n)) == 'RN9NO3    ' .OR. (speci(nldepd(n)) == 'RN12NO3   ' .OR. (speci(nldepd(n)) == 'RN15NO3   ' .OR. (speci(nldepd(n)) == 'RN18NO3   ' .OR. (speci(nldepd(n)) == 'RU14NO3   ' .OR. (speci(nldepd(n)) == 'RTN28NO3  ' .OR. (speci(nldepd(n)) == 'RTN25NO3  ' .OR. (speci(nldepd(n)) == 'RTX28NO3  ' .OR. (speci(nldepd(n)) == 'RTX22NO3  ' .OR. (speci(nldepd(n)) == 'RA22NO3   ' .OR. (speci(nldepd(n)) == 'RA25NO3   ' .OR. (speci(nldepd(n)) == 'RTN23NO3  ' .OR. (speci(nldepd(n)) == 'RU12NO3   ' .OR. (speci(nldepd(n)) == 'RU10NO3   ' .OR. (speci(nldepd(n)) == 'RA13NO3   ' .OR. (speci(nldepd(n)) == 'RA16NO3   ' .OR. speci(nldepd(n)) == 'RA19NO3   '))))))))))))))))))))) then
                                    rsurf(18:27,n) = [4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0]
                                  else
                                    if (speci(nldepd(n)) == 'NH3       ') then
                                      rsurf(18:27,n) = [215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1, 215.1]
                                    else
                                      if (speci(nldepd(n)) == 'CO        ') then
                                        rsurf(18:27,n) = [r_null, r_null, r_null, r_null, r_null, r_null, r_null, r_null, r_null, r_null]
                                      else
                                        if (speci(nldepd(n)) == 'HCHO      ') then
                                          rsurf(18:27,n) = [261.8, 261.8, 261.8, 261.8, 261.8, 261.8, 261.8, 261.8, 261.8, 261.8]
                                        else
                                          if (speci(nldepd(n)) == 'MeOH      ' .OR. (speci(nldepd(n)) == 'EtOH      ' .OR. (speci(nldepd(n)) == 'i-PrOH    ' .OR. (speci(nldepd(n)) == 'n-PrOH    ' .OR. (speci(nldepd(n)) == 'AROH14    ' .OR. (speci(nldepd(n)) == 'ARNOH14   ' .OR. (speci(nldepd(n)) == 'AROH17    ' .OR. (speci(nldepd(n)) == 'ARNOH17   ' .OR. (speci(nldepd(n)) == 'IEPOX     ' .OR. (speci(nldepd(n)) == 'HMML      ' .OR. (speci(nldepd(n)) == 'HUCARB9   ' .OR. speci(nldepd(n)) == 'DHCARB9   '))))))))))) then
                                            if (l_fix_improve_drydep) then
                                              rsurf(18:27,n) = [352.2, 352.2, 352.2, 352.2, 352.2, 352.2, 352.2, 352.2, 352.2, 352.2]
                                            end if
                                          else
                                            if (speci(nldepd(n)) == 'MeCHO     ' .OR. (speci(nldepd(n)) == 'EtCHO     ' .OR. (speci(nldepd(n)) == 'MACR      ' .OR. (speci(nldepd(n)) == 'NALD      ' .OR. (speci(nldepd(n)) == 'HOCH2CHO  ' .OR. (speci(nldepd(n)) == 'CARB14    ' .OR. (speci(nldepd(n)) == 'CARB17    ' .OR. (speci(nldepd(n)) == 'CARB11A   ' .OR. (speci(nldepd(n)) == 'UCARB10   ' .OR. (speci(nldepd(n)) == 'UCARB12   ' .OR. (speci(nldepd(n)) == 'UDCARB8   ' .OR. (speci(nldepd(n)) == 'UDCARB11  ' .OR. (speci(nldepd(n)) == 'UDCARB14  ' .OR. (speci(nldepd(n)) == 'TNCARB26  ' .OR. (speci(nldepd(n)) == 'TNCARB10  ' .OR. (speci(nldepd(n)) == 'TNCARB12  ' .OR. (speci(nldepd(n)) == 'TNCARB11  ' .OR. (speci(nldepd(n)) == 'CCARB12   ' .OR. (speci(nldepd(n)) == 'TNCARB15  ' .OR. (speci(nldepd(n)) == 'TXCARB24  ' .OR. (speci(nldepd(n)) == 'TXCARB22  ' .OR. (speci(nldepd(n)) == 'UDCARB17  ' .OR. (speci(nldepd(n)) == 'NOA       ' .OR. (speci(nldepd(n)) == 'NUCARB12  ' .OR. speci(nldepd(n)) == 'ANHY      ')))))))))))))))))))))))) then
                                              rsurf(18:27,n) = [40000.0, 40000.0, 40000.0, 40000.0, 40000.0, 40000.0, 40000.0, 40000.0, 40000.0, 40000.0]
                                            else
                                              if (speci(nldepd(n)) == 'MGLY      ' .OR. (speci(nldepd(n)) == 'CARB3     ' .OR. (speci(nldepd(n)) == 'CARB6     ' .OR. (speci(nldepd(n)) == 'CARB9     ' .OR. (speci(nldepd(n)) == 'CARB12    ' .OR. speci(nldepd(n)) == 'CARB15    '))))) then
                                                rsurf(18:27,n) = [21505.4, 21505.4, 21505.4, 21505.4, 21505.4, 21505.4, 21505.4, 21505.4, 21505.4, 21505.4]
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        enddo
      end if
      first = .false.
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_33 = 1, rows, 1
      do idx_34 = 1, row_length, 1
        o3_stom_frac(idx_34,idx_33) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do n = 1, ntype, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          todo(i,k,n) = gsf(i,k,n) > min_tile_frac
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_35 = LBOUND(microb, dim=2), UBOUND(microb, dim=2), 1
      do idx_36 = LBOUND(microb, dim=1), UBOUND(microb, dim=1), 1
        microb(idx_36,idx_35) = rh(idx_36,idx_35) > 0.4 .AND. t0(idx_36,idx_35) > 278.0
      enddo
    enddo
    !$acc end kernels
    do n = 1, ntype, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          p_surf1 = p0(i,k)
          t_surf1 = t0(i,k)
          ustar = u_s(i,k)
          sst = t0(i,k)
          if (todo(i,k,n)) then
            do j = 1, ndepd, 1
              rc(i,k,n,j) = rsurf(n,j)
              if (ukca_config%l_ukca_ddepo3_ocean) then
                if (ntype == 9 .AND. n == 7 .OR. ntype == 13 .AND. n == 11 .OR. ntype == 17 .AND. n == 15 .OR. ntype == 27 .AND. n == 15) then
                  if ((speci(nldepd(j)) == 'O3        ' .OR. speci(nldepd(j)) == 'O3S       ') .AND. gsf(i,k,n) > 0.75) then
                    call ukca_ddepo3_ocean(p_surf1, t_surf1, sst, ustar, rc_ocean)
                    rc(i,k,n,j) = rc_ocean
                  end if
                else
                  if (ntype /= 9 .AND. ntype /= 13 .AND. ntype /= 17 .AND. ntype /= 27) then
                    ! PSyclone CodeBlock (unsupported code) reason:
                    !  - Unsupported statement: Write_Stmt
                    WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                    errcode = 1
                    call ereport(routinename, errcode, cmessage)
                  end if
                end if
              end if
            enddo
          else
            !$acc kernels
            !$acc loop gang vector independent
            do j = 1, ndepd, 1
              rc(i,k,n,j) = r_null
            enddo
            !$acc end kernels
          end if
        enddo
      enddo
    enddo
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do n = 1, npft, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          if (todo(i,k,n) .AND. stcon(i,k,n) > glmin) then
            r_stom(i,k,n,r_stom_no2) = 1.5 * dif(r_stom_no2) / stcon(i,k,n)
            r_stom(i,k,n,r_stom_o3) = dif(r_stom_o3) / stcon(i,k,n)
            r_stom(i,k,n,r_stom_pan) = dif(r_stom_pan) / stcon(i,k,n)
            r_stom(i,k,n,r_stom_so2) = dif(r_stom_so2) / stcon(i,k,n)
            r_stom(i,k,n,r_stom_nh3) = dif(r_stom_nh3) / stcon(i,k,n)
          else
            r_stom(i,k,n,r_stom_no2) = r_null
            r_stom(i,k,n,r_stom_o3) = r_null
            r_stom(i,k,n,r_stom_pan) = r_null
            r_stom(i,k,n,r_stom_so2) = r_null
            r_stom(i,k,n,r_stom_nh3) = r_null
          end if
        enddo
      enddo
    enddo
    !$acc end kernels
    do j = 1, ndepd, 1
      if (speci(nldepd(j)) == 'O3        ') then
        do k = 1, rows, 1
          !$acc kernels
          !$acc loop seq
          do i = 1, row_length, 1
            if (smr(i,k) > soil_moist_thresh) then
              !$acc loop gang vector independent
              do n = 1, npft, 1
                if (todo(i,k,n)) then
                  rc(i,k,n,j) = r_wet_o3
                end if
              enddo
              if (todo(i,k,soil)) then
                rc(i,k,soil,j) = r_wet_o3
              end if
            end if
          enddo
          !$acc end kernels
          if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
            if (l_fix_ukca_h2dd_x) then
              n = shrub
            else
              n = npft
            end if
            !$acc kernels
            !$acc loop gang vector independent
            do i = 1, row_length, 1
              if (sinlat(i,k) > tundra_s_limit) then
                if (todo(i,k,n)) then
                  rc(i,k,n,j) = r_tundra(r_tundra_o3)
                end if
              end if
            enddo
            !$acc end kernels
          else
            if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
              !$acc kernels
              !$acc loop gang vector independent collapse(2)
              do n = shrub_dec, shrub_eg, 1
                do i = 1, row_length, 1
                  if (sinlat(i,k) > tundra_s_limit) then
                    if (todo(i,k,n)) then
                      rc(i,k,n,j) = r_tundra(r_tundra_o3)
                    end if
                  end if
                enddo
              enddo
              !$acc end kernels
            else
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
              errcode = 1
              call ereport(routinename, errcode, cmessage)
            end if
          end if
          !$acc kernels
          !$acc loop gang vector independent
          do i = 1, row_length, 1
            if (sinlat(i,k) > tundra_s_limit) then
              if (todo(i,k,soil)) then
                rc(i,k,soil,j) = r_tundra(r_tundra_o3)
              end if
            end if
          enddo
          !$acc loop gang vector independent collapse(2)
          do n = 1, npft, 1
            do i = 1, row_length, 1
              if (todo(i,k,n) .AND. lai_ft(i,k,n) > 0.0) then
                r_cut_o3(i,n) = cuticle_o3 / lai_ft(i,k,n)
              else
                r_cut_o3(i,n) = r_null
              end if
            enddo
          enddo
          !$acc loop seq
          do n = 1, npft, 1
            !$acc loop gang vector independent
            do i = 1, row_length, 1
              if (todo(i,k,n)) then
                rr = 1.0 / r_stom(i,k,n,r_stom_o3) + 1.0 / r_cut_o3(i,n) + 1.0 / rc(i,k,n,j)
                rc(i,k,n,j) = 1.0 / rr
                o3_stom_frac(i,k) = o3_stom_frac(i,k) + gsf(i,k,n) / (rr * r_stom(i,k,n,r_stom_o3))
              end if
            enddo
          enddo
          !$acc end kernels
        enddo
      else
        if (speci(nldepd(j)) == 'SO2       ') then
          if (ukca_config%l_ukca_dry_dep_so2wet) then
            if (soil > 0) then
              !$acc kernels
              !$acc loop gang vector independent collapse(2)
              do k = 1, rows, 1
                do i = 1, row_length, 1
                  if (smr(i,k) > soil_moist_thresh .OR. surf_wetness(i,k) > tol) then
                    if (t0tile(i,k,soil) > minus1degc) then
                      if (todo(i,k,soil)) then
                        rc(i,k,soil,j) = r_wet_so2
                      end if
                    end if
                  end if
                enddo
              enddo
              !$acc end kernels
            end if
            !$acc kernels
            !$acc loop gang vector independent collapse(3)
            do n = 1, npft, 1
              do k = 1, rows, 1
                do i = 1, row_length, 1
                  if (todo(i,k,n)) then
                    if (smr(i,k) > soil_moist_thresh .OR. surf_wetness(i,k) > tol) then
                      if (t0tile(i,k,n) > minus1degc) then
                        rc(i,k,n,j) = r_wet_so2
                      end if
                    end if
                    if (lai_ft(i,k,n) > tol) then
                      if (t0tile(i,k,n) >= minus1degc) then
                        if (smr(i,k) > soil_moist_thresh .OR. surf_wetness(i,k) > tol) then
                          r_cut_so2 = r_wet_so2
                        else
                          if (n == c3_crop .OR. n == c4_crop) then
                            r_cut_so2 = rc(i,k,n,j)
                          else
                            if (rh(i,k) < rh_thresh) then
                              r_cut_so2 = so2con1 * EXP(so2con3 * rh(i,k))
                            else
                              r_cut_so2 = so2con2 * EXP(so2con4 * rh(i,k))
                            end if
                          end if
                          r_cut_so2 = MAX(r_wet_so2, MIN(r_dry_so2, r_cut_so2))
                        end if
                      else
                        if (t0tile(i,k,n) >= minus5degc) then
                          r_cut_so2 = r_cut_so2_5to_1
                        else
                          r_cut_so2 = r_cut_so2_5degc
                        end if
                      end if
                    else
                      r_cut_so2 = r_null
                    end if
                    rr = 1.0 / r_stom(i,k,n,r_stom_so2) + 1.0 / r_cut_so2 + 1.0 / rc(i,k,n,j)
                    rc(i,k,n,j) = 1.0 / rr
                  end if
                enddo
              enddo
            enddo
            !$acc end kernels
          end if
        else
          if (speci(nldepd(j)) == 'NO2       ') then
            do k = 1, rows, 1
              if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                if (l_fix_ukca_h2dd_x) then
                  n = shrub
                else
                  n = npft
                end if
                !$acc kernels
                !$acc loop gang vector independent
                do i = 1, row_length, 1
                  if (sinlat(i,k) > tundra_s_limit) then
                    if (todo(i,k,n)) then
                      rc(i,k,n,j) = r_tundra(r_tundra_no2)
                    end if
                  end if
                enddo
                !$acc end kernels
              else
                if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                  !$acc kernels
                  !$acc loop gang vector independent collapse(2)
                  do n = shrub_dec, shrub_eg, 1
                    do i = 1, row_length, 1
                      if (sinlat(i,k) > tundra_s_limit) then
                        if (todo(i,k,n)) then
                          rc(i,k,n,j) = r_tundra(r_tundra_no2)
                        end if
                      end if
                    enddo
                  enddo
                  !$acc end kernels
                else
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                  errcode = 1
                  call ereport(routinename, errcode, cmessage)
                end if
              end if
              !$acc kernels
              !$acc loop gang vector independent
              do i = 1, row_length, 1
                if (sinlat(i,k) > tundra_s_limit) then
                  if (todo(i,k,soil)) then
                    rc(i,k,soil,j) = r_tundra(r_tundra_no2)
                  end if
                end if
              enddo
              !$acc loop gang vector independent collapse(2)
              do n = 1, npft, 1
                do i = 1, row_length, 1
                  if (todo(i,k,n)) then
                    rr = 1.0 / r_stom(i,k,n,r_stom_no2) + 1.0 / rc(i,k,n,j)
                    rc(i,k,n,j) = 1.0 / rr
                  end if
                enddo
              enddo
              !$acc end kernels
            enddo
          else
            if (speci(nldepd(j)) == 'PAN       ') then
              do k = 1, rows, 1
                if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                  if (l_fix_ukca_h2dd_x) then
                    n = shrub
                  else
                    n = npft
                  end if
                  !$acc kernels
                  !$acc loop gang vector independent
                  do i = 1, row_length, 1
                    if (sinlat(i,k) > tundra_s_limit) then
                      if (todo(i,k,n)) then
                        rc(i,k,n,j) = r_tundra(r_tundra_pan)
                      end if
                    end if
                  enddo
                  !$acc end kernels
                else
                  if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                    !$acc kernels
                    !$acc loop gang vector independent collapse(2)
                    do n = shrub_dec, shrub_eg, 1
                      do i = 1, row_length, 1
                        if (sinlat(i,k) > tundra_s_limit) then
                          if (todo(i,k,n)) then
                            rc(i,k,n,j) = r_tundra(r_tundra_pan)
                          end if
                        end if
                      enddo
                    enddo
                    !$acc end kernels
                  else
                    ! PSyclone CodeBlock (unsupported code) reason:
                    !  - Unsupported statement: Write_Stmt
                    WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                    errcode = 1
                    call ereport(routinename, errcode, cmessage)
                  end if
                end if
                !$acc kernels
                !$acc loop gang vector independent
                do i = 1, row_length, 1
                  if (sinlat(i,k) > tundra_s_limit) then
                    if (todo(i,k,soil)) then
                      rc(i,k,soil,j) = r_tundra(r_tundra_pan)
                    end if
                  end if
                enddo
                !$acc loop gang vector independent collapse(2)
                do n = 1, npft, 1
                  do i = 1, row_length, 1
                    if (todo(i,k,n)) then
                      rr = 1.0 / r_stom(i,k,n,r_stom_pan) + 1.0 / rc(i,k,n,j)
                      rc(i,k,n,j) = 1.0 / rr
                    end if
                  enddo
                enddo
                !$acc end kernels
              enddo
            else
              if (speci(nldepd(j)) == 'PPAN      ') then
                do k = 1, rows, 1
                  if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                    if (l_fix_ukca_h2dd_x) then
                      n = shrub
                    else
                      n = npft
                    end if
                    !$acc kernels
                    !$acc loop gang vector independent
                    do i = 1, row_length, 1
                      if (sinlat(i,k) > tundra_s_limit) then
                        if (todo(i,k,n)) then
                          rc(i,k,n,j) = r_tundra(r_tundra_pan)
                        end if
                      end if
                    enddo
                    !$acc end kernels
                  else
                    if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                      !$acc kernels
                      !$acc loop gang vector independent collapse(2)
                      do n = shrub_dec, shrub_eg, 1
                        do i = 1, row_length, 1
                          if (sinlat(i,k) > tundra_s_limit) then
                            if (todo(i,k,n)) then
                              rc(i,k,n,j) = r_tundra(r_tundra_pan)
                            end if
                          end if
                        enddo
                      enddo
                      !$acc end kernels
                    else
                      ! PSyclone CodeBlock (unsupported code) reason:
                      !  - Unsupported statement: Write_Stmt
                      WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                      errcode = 1
                      call ereport(routinename, errcode, cmessage)
                    end if
                  end if
                  !$acc kernels
                  !$acc loop gang vector independent
                  do i = 1, row_length, 1
                    if (sinlat(i,k) > tundra_s_limit) then
                      if (todo(i,k,soil)) then
                        rc(i,k,soil,j) = r_tundra(r_tundra_pan)
                      end if
                    end if
                  enddo
                  !$acc loop gang vector independent collapse(2)
                  do n = 1, npft, 1
                    do i = 1, row_length, 1
                      if (todo(i,k,n)) then
                        rr = 1.0 / r_stom(i,k,n,r_stom_pan) + 1.0 / rc(i,k,n,j)
                        rc(i,k,n,j) = 1.0 / rr
                      end if
                    enddo
                  enddo
                  !$acc end kernels
                enddo
              else
                if (speci(nldepd(j)) == 'MPAN      ') then
                  do k = 1, rows, 1
                    if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                      if (l_fix_ukca_h2dd_x) then
                        n = shrub
                      else
                        n = npft
                      end if
                      !$acc kernels
                      !$acc loop gang vector independent
                      do i = 1, row_length, 1
                        if (sinlat(i,k) > tundra_s_limit) then
                          if (todo(i,k,n)) then
                            rc(i,k,n,j) = r_tundra(r_tundra_pan)
                          end if
                        end if
                      enddo
                      !$acc end kernels
                    else
                      if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                        !$acc kernels
                        !$acc loop gang vector independent collapse(2)
                        do n = shrub_dec, shrub_eg, 1
                          do i = 1, row_length, 1
                            if (sinlat(i,k) > tundra_s_limit) then
                              if (todo(i,k,n)) then
                                rc(i,k,n,j) = r_tundra(r_tundra_pan)
                              end if
                            end if
                          enddo
                        enddo
                        !$acc end kernels
                      else
                        ! PSyclone CodeBlock (unsupported code) reason:
                        !  - Unsupported statement: Write_Stmt
                        WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                        errcode = 1
                        call ereport(routinename, errcode, cmessage)
                      end if
                    end if
                    !$acc kernels
                    !$acc loop gang vector independent
                    do i = 1, row_length, 1
                      if (sinlat(i,k) > tundra_s_limit) then
                        if (todo(i,k,soil)) then
                          rc(i,k,soil,j) = r_tundra(r_tundra_pan)
                        end if
                      end if
                    enddo
                    !$acc loop gang vector independent collapse(2)
                    do n = 1, npft, 1
                      do i = 1, row_length, 1
                        if (todo(i,k,n)) then
                          rr = 1.0 / r_stom(i,k,n,r_stom_pan) + 1.0 / rc(i,k,n,j)
                          rc(i,k,n,j) = 1.0 / rr
                        end if
                      enddo
                    enddo
                    !$acc end kernels
                  enddo
                else
                  if (speci(nldepd(j)) == 'ONITU     ') then
                    do k = 1, rows, 1
                      if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                        if (l_fix_ukca_h2dd_x) then
                          n = shrub
                        else
                          n = npft
                        end if
                        !$acc kernels
                        !$acc loop gang vector independent
                        do i = 1, row_length, 1
                          if (sinlat(i,k) > tundra_s_limit) then
                            if (todo(i,k,n)) then
                              rc(i,k,n,j) = r_tundra(r_tundra_pan)
                            end if
                          end if
                        enddo
                        !$acc end kernels
                      else
                        if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                          !$acc kernels
                          !$acc loop gang vector independent collapse(2)
                          do n = shrub_dec, shrub_eg, 1
                            do i = 1, row_length, 1
                              if (sinlat(i,k) > tundra_s_limit) then
                                if (todo(i,k,n)) then
                                  rc(i,k,n,j) = r_tundra(r_tundra_pan)
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        else
                          ! PSyclone CodeBlock (unsupported code) reason:
                          !  - Unsupported statement: Write_Stmt
                          WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                          errcode = 1
                          call ereport(routinename, errcode, cmessage)
                        end if
                      end if
                      !$acc kernels
                      !$acc loop gang vector independent
                      do i = 1, row_length, 1
                        if (sinlat(i,k) > tundra_s_limit) then
                          if (todo(i,k,soil)) then
                            rc(i,k,soil,j) = r_tundra(r_tundra_pan)
                          end if
                        end if
                      enddo
                      !$acc loop gang vector independent collapse(2)
                      do n = 1, npft, 1
                        do i = 1, row_length, 1
                          if (todo(i,k,n)) then
                            rr = 1.0 / r_stom(i,k,n,r_stom_pan) + 1.0 / rc(i,k,n,j)
                            rc(i,k,n,j) = 1.0 / rr
                          end if
                        enddo
                      enddo
                      !$acc end kernels
                    enddo
                  else
                    if (speci(nldepd(j)) == 'H2        ') then
                      if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                        !$acc kernels
                        !$acc loop seq
                        do k = 1, rows, 1
                          if (l_fix_ukca_h2dd_x) then
                            !$acc loop gang vector independent collapse(2)
                            do n = brd_leaf, c3_grass, 1
                              do i = 1, row_length, 1
                                if (todo(i,k,n) .AND. microb(i,k)) then
                                  sm = MAX(smr(i,k), 0.1)
                                  rr = (h2dd_m(n) * sm + h2dd_c(n)) * 1.0e-4
                                  rc(i,k,n,j) = 1.0 / rr
                                end if
                              enddo
                            enddo
                          else
                            !$acc loop gang vector independent collapse(2)
                            do n = 1, npft - 2, 1
                              do i = 1, row_length, 1
                                if (todo(i,k,n) .AND. microb(i,k)) then
                                  sm = MAX(smr(i,k), 0.1)
                                  rc(i,k,n,j) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
                                end if
                              enddo
                            enddo
                          end if
                          if (l_fix_ukca_h2dd_x) then
                            n = c4_grass
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                sm = LOG(MAX(smr(i,k), 0.1))
                                rr = (h2dd_c(c4_grass) + sm * (h2dd_m(c4_grass) + sm * h2dd_q(c4_grass))) * 1.0e-4
                                if (rr > 0.00131) then
                                  rr = 0.00131
                                end if
                                rc(i,k,n,j) = 1.0 / rr
                              end if
                            enddo
                          else
                            n = npft - 1
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                sm = LOG(MAX(smr(i,k), 0.1))
                                rr = (h2dd_c(c4_grass) + sm * (h2dd_m(c4_grass) + sm * h2dd_q(c4_grass))) * 1.0e-4
                                if (rr > 0.00131) then
                                  rr = 0.00131
                                end if
                                rc(i,k,n,j) = 1.0 / rr
                              end if
                            enddo
                          end if
                          if (l_fix_ukca_h2dd_x) then
                            n = shrub
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  rr = h2dd_c(shrub) * 1.0e-4
                                  rc(i,k,n,j) = 1.0 / rr
                                else
                                  sm = MAX(smr(i,k), 0.1)
                                  rr = (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass)) * 1.0e-4
                                  rc(i,k,n,j) = 1.0 / rr
                                end if
                              end if
                            enddo
                          else
                            n = npft
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  rr = 1.0 / h2dd_c(n)
                                  rc(i,k,n,j) = rr
                                else
                                  sm = MAX(smr(i,k), 0.1)
                                  rc(i,k,n,j) = 1.0 / (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass))
                                end if
                              end if
                            enddo
                          end if
                          if (l_fix_ukca_h2dd_x) then
                            n = soil
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  rr = 1.0 / h2dd_c(shrub) * 1.0e-4
                                  rc(i,k,n,j) = 1.0 / rr
                                else
                                  sm = MAX(smr(i,k), 0.1)
                                  rr = (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass)) * 1.0e-4
                                  rc(i,k,n,j) = 1.0 / rr
                                end if
                              end if
                            enddo
                          else
                            n = npft
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,npft) .AND. microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  rr = 1.0 / h2dd_c(npft)
                                  rc(i,k,soil,j) = rr
                                else
                                  sm = MAX(smr(i,k), 0.1)
                                  rc(i,k,soil,j) = 1.0 / (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass))
                                end if
                              end if
                            enddo
                          end if
                        enddo
                        !$acc end kernels
                      else
                        if (ntype == 13) then
                          !$acc kernels
                          !$acc loop seq
                          do k = 1, rows, 1
                            !$acc loop gang vector independent collapse(2)
                            do n = brd_leaf_dec, c3_grass, 1
                              do i = 1, row_length, 1
                                if (todo(i,k,n) .AND. microb(i,k)) then
                                  sm = MAX(smr(i,k), 0.1)
                                  rr = (h2dd_m(n) * sm + h2dd_c(n)) * 1.0e-4
                                  rc(i,k,n,j) = 1.0 / rr
                                end if
                              enddo
                            enddo
                            n = c4_grass
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                sm = LOG(MAX(smr(i,k), 0.1))
                                rr = (h2dd_c(c4_grass) + sm * (h2dd_m(c4_grass) + sm * h2dd_q(c4_grass))) * 1.0e-4
                                if (rr > 0.00131) then
                                  rr = 0.00131
                                end if
                                rc(i,k,n,j) = 1.0 / rr
                              end if
                            enddo
                            !$acc loop gang vector independent collapse(2)
                            do n = shrub_dec, shrub_eg, 1
                              do i = 1, row_length, 1
                                if (todo(i,k,n) .AND. microb(i,k)) then
                                  if (sinlat(i,k) > tundra_s_limit) then
                                    rr = h2dd_c(shrub_dec) * 1.0e-4
                                    rc(i,k,n,j) = 1.0 / rr
                                  else
                                    if (sinlat(i,k) <= tundra_s_limit) then
                                      sm = MAX(smr(i,k), 0.1)
                                      rr = (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass)) * 1.0e-4
                                      rc(i,k,n,j) = 1.0 / rr
                                    end if
                                  end if
                                end if
                              enddo
                            enddo
                            n = soil
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  rr = h2dd_c(shrub_dec) * 1.0e-4
                                  rc(i,k,soil,j) = 1.0 / rr
                                else
                                  if (sinlat(i,k) <= tundra_s_limit) then
                                    sm = MAX(smr(i,k), 0.1)
                                    rr = (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass)) * 1.0e-4
                                    rc(i,k,n,j) = 1.0 / rr
                                  end if
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        else
                          if (ntype == 17 .OR. ntype == 27) then
                            !$acc kernels
                            !$acc loop seq
                            do k = 1, rows, 1
                              !$acc loop gang vector independent collapse(2)
                              do n = brd_leaf_dec, c3_pasture, 1
                                do i = 1, row_length, 1
                                  if (todo(i,k,n) .AND. microb(i,k)) then
                                    sm = MAX(smr(i,k), 0.1)
                                    rr = (h2dd_m(n) * sm + h2dd_c(n)) * 1.0e-4
                                    rc(i,k,n,j) = 1.0 / rr
                                  end if
                                enddo
                              enddo
                              !$acc loop gang vector independent collapse(2)
                              do n = c4_grass, c4_pasture, 1
                                do i = 1, row_length, 1
                                  if (todo(i,k,n) .AND. microb(i,k)) then
                                    sm = LOG(MAX(smr(i,k), 0.1))
                                    rr = (h2dd_c(n) + sm * (h2dd_m(n) + sm * h2dd_q(n))) * 1.0e-4
                                    if (rr > 0.00131) then
                                      rr = 0.00131
                                    end if
                                    rc(i,k,n,j) = 1.0 / rr
                                  end if
                                enddo
                              enddo
                              !$acc loop gang vector independent collapse(2)
                              do n = shrub_dec, shrub_eg, 1
                                do i = 1, row_length, 1
                                  if (todo(i,k,n) .AND. microb(i,k)) then
                                    if (sinlat(i,k) > tundra_s_limit) then
                                      rr = h2dd_c(shrub_dec) * 1.0e-4
                                      rc(i,k,n,j) = 1.0 / rr
                                    else
                                      if (sinlat(i,k) <= tundra_s_limit) then
                                        sm = MAX(smr(i,k), 0.1)
                                        rr = (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass)) * 1.0e-4
                                        rc(i,k,n,j) = 1.0 / rr
                                      end if
                                    end if
                                  end if
                                enddo
                              enddo
                              n = soil
                              !$acc loop gang vector independent
                              do i = 1, row_length, 1
                                if (todo(i,k,n) .AND. microb(i,k)) then
                                  if (sinlat(i,k) > tundra_s_limit) then
                                    rr = h2dd_c(shrub_dec) * 1.0e-4
                                    rc(i,k,soil,j) = 1.0 / rr
                                  else
                                    if (sinlat(i,k) <= tundra_s_limit) then
                                      sm = MAX(smr(i,k), 0.1)
                                      rr = (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass)) * 1.0e-4
                                      rc(i,k,n,j) = 1.0 / rr
                                    end if
                                  end if
                                end if
                              enddo
                            enddo
                            !$acc end kernels
                          else
                            ! PSyclone CodeBlock (unsupported code) reason:
                            !  - Unsupported statement: Write_Stmt
                            WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                            errcode = 1
                            call ereport(routinename, errcode, cmessage)
                          end if
                        end if
                      end if
                    else
                      if (speci(nldepd(j)) == 'CH4       ') then
                        do k = 1, rows, 1
                          do i = 1, row_length, 1
                            if (ukca_config%l_ukca_emsdrvn_ch4) then
                              !$acc kernels
                              !$acc loop gang vector independent
                              do n = 1, ntype, 1
                                rc(i,k,n,j) = 1.0 / r_null
                              enddo
                              !$acc end kernels
                            else
                              !$acc kernels
                              !$acc loop gang vector independent
                              do n = 1, ntype, 1
                                rc(i,k,n,j) = r_null
                              enddo
                              !$acc end kernels
                            end if
                            if (microb(i,k)) then
                              sm = smr(i,k)
                              if (sm < smfrac1) then
                                f_ch4_uptk = sm / smfrac1
                              else
                                if (sm > smfrac2) then
                                  f_ch4_uptk = (smfrac3 - sm) / smfrac4
                                else
                                  f_ch4_uptk = 1.0
                                end if
                              end if
                              if (ukca_config%l_ukca_emsdrvn_ch4) then
                                f_ch4_uptk = MIN(MAX(f_ch4_uptk, 0.0), 1.0)
                              else
                                f_ch4_uptk = MAX(f_ch4_uptk, 0.0)
                              end if
                              if (.NOT.l_fix_ukca_h2dd_x) then
                                if (todo(i,k,1)) then
                                  rc(i,k,1,j) = ch4_up_flux(1) * f_ch4_uptk
                                end if
                                if (todo(i,k,2)) then
                                  rc(i,k,2,j) = ch4_up_flux(2) * f_ch4_uptk
                                end if
                                if (todo(i,k,3)) then
                                  rc(i,k,3,j) = ch4_up_flux(3) * f_ch4_uptk
                                end if
                                if (todo(i,k,4)) then
                                  rc(i,k,4,j) = ch4_up_flux(4) * f_ch4_uptk
                                end if
                              else
                                if (ntype == 9) then
                                  !$acc kernels
                                  !$acc loop gang vector independent
                                  do n = brd_leaf, c4_grass, 1
                                    if (todo(i,k,n)) then
                                      rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
                                    end if
                                  enddo
                                  !$acc end kernels
                                else
                                  if (ntype == 13) then
                                    !$acc kernels
                                    !$acc loop gang vector independent
                                    do n = brd_leaf_dec, c4_grass, 1
                                      if (todo(i,k,n)) then
                                        rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
                                      end if
                                    enddo
                                    !$acc end kernels
                                  else
                                    if (ntype == 17 .OR. ntype == 27) then
                                      !$acc kernels
                                      !$acc loop gang vector independent
                                      do n = brd_leaf_dec, c4_pasture, 1
                                        if (todo(i,k,n)) then
                                          rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
                                        end if
                                      enddo
                                      !$acc end kernels
                                    else
                                      ! PSyclone CodeBlock (unsupported code) reason:
                                      !  - Unsupported statement: Write_Stmt
                                      WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                                      errcode = 1
                                      call ereport(routinename, errcode, cmessage)
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          enddo
                        enddo
                        if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                          if (l_fix_ukca_h2dd_x) then
                            n = shrub
                          else
                            n = npft
                          end if
                          !$acc kernels
                          !$acc loop gang vector independent collapse(2)
                          do k = 1, rows, 1
                            do i = 1, row_length, 1
                              if (microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  ts = t0tile(i,k,n)
                                  rr = ch4dd_tun(fourth_term) + ts * (ch4dd_tun(third_term) + ts * (ch4dd_tun(second_term) + ts * ch4dd_tun(first_term)))
                                  rr = rr * 3600.0
                                  if (todo(i,k,n)) then
                                    rc(i,k,n,j) = MAX(rr, 0.0)
                                  end if
                                else
                                  sm = smr(i,k)
                                  if (sm < smfrac1) then
                                    f_ch4_uptk = sm / smfrac1
                                  else
                                    if (sm > smfrac2) then
                                      f_ch4_uptk = (smfrac3 - sm) / smfrac4
                                    else
                                      f_ch4_uptk = 1.0
                                    end if
                                  end if
                                  if (ukca_config%l_ukca_emsdrvn_ch4) then
                                    f_ch4_uptk = MIN(MAX(f_ch4_uptk, 0.0), 1.0)
                                  else
                                    f_ch4_uptk = MAX(f_ch4_uptk, 0.0)
                                  end if
                                  if (todo(i,k,n)) then
                                    rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
                                  end if
                                end if
                              else
                                if (ukca_config%l_ukca_emsdrvn_ch4) then
                                  rc(i,k,n,j) = 1.0 / r_null
                                else
                                  rc(i,k,n,j) = r_null
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        else
                          if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                            !$acc kernels
                            !$acc loop gang vector independent collapse(3)
                            do n = shrub_dec, shrub_eg, 1
                              do k = 1, rows, 1
                                do i = 1, row_length, 1
                                  if (microb(i,k)) then
                                    if (sinlat(i,k) > tundra_s_limit) then
                                      ts = t0tile(i,k,n)
                                      rr = ch4dd_tun(fourth_term) + ts * (ch4dd_tun(third_term) + ts * (ch4dd_tun(second_term) + ts * ch4dd_tun(first_term)))
                                      rr = rr * 3600.0
                                      if (todo(i,k,n)) then
                                        rc(i,k,n,j) = MAX(rr, 0.0)
                                      end if
                                    else
                                      sm = smr(i,k)
                                      if (sm < smfrac1) then
                                        f_ch4_uptk = sm / smfrac1
                                      else
                                        if (sm > smfrac2) then
                                          f_ch4_uptk = (smfrac3 - sm) / smfrac4
                                        else
                                          f_ch4_uptk = 1.0
                                        end if
                                      end if
                                      if (ukca_config%l_ukca_emsdrvn_ch4) then
                                        f_ch4_uptk = MIN(MAX(f_ch4_uptk, 0.0), 1.0)
                                      else
                                        f_ch4_uptk = MAX(f_ch4_uptk, 0.0)
                                      end if
                                      if (todo(i,k,n)) then
                                        rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
                                      end if
                                    end if
                                  else
                                    if (ukca_config%l_ukca_emsdrvn_ch4) then
                                      rc(i,k,n,j) = 1.0 / r_null
                                    else
                                      rc(i,k,n,j) = r_null
                                    end if
                                  end if
                                enddo
                              enddo
                            enddo
                            !$acc end kernels
                          else
                            ! PSyclone CodeBlock (unsupported code) reason:
                            !  - Unsupported statement: Write_Stmt
                            WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                            errcode = 1
                            call ereport(routinename, errcode, cmessage)
                          end if
                        end if
                        if (l_fix_ukca_h2dd_x) then
                          !$acc kernels
                          n = soil
                          !$acc loop gang vector independent collapse(2)
                          do k = 1, rows, 1
                            do i = 1, row_length, 1
                              if (microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  ts = t0tile(i,k,soil)
                                  rr = ch4dd_tun(fourth_term) + ts * (ch4dd_tun(third_term) + ts * (ch4dd_tun(second_term) + ts * ch4dd_tun(first_term)))
                                  rr = rr * 3600.0
                                  if (todo(i,k,soil)) then
                                    rc(i,k,soil,j) = MAX(rr, 0.0)
                                  end if
                                else
                                  sm = smr(i,k)
                                  if (sm < smfrac1) then
                                    f_ch4_uptk = sm / smfrac1
                                  else
                                    if (sm > smfrac2) then
                                      f_ch4_uptk = (smfrac3 - sm) / smfrac4
                                    else
                                      f_ch4_uptk = 1.0
                                    end if
                                  end if
                                  if (ukca_config%l_ukca_emsdrvn_ch4) then
                                    f_ch4_uptk = MIN(MAX(f_ch4_uptk, 0.0), 1.0)
                                  else
                                    f_ch4_uptk = MAX(f_ch4_uptk, 0.0)
                                  end if
                                  if (todo(i,k,soil)) then
                                    rc(i,k,soil,j) = ch4_up_flux(soil) * f_ch4_uptk
                                  end if
                                end if
                              else
                                if (ukca_config%l_ukca_emsdrvn_ch4) then
                                  rc(i,k,n,j) = 1.0 / r_null
                                else
                                  rc(i,k,n,j) = r_null
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        else
                          !$acc kernels
                          n = npft
                          !$acc loop gang vector independent collapse(2)
                          do k = 1, rows, 1
                            do i = 1, row_length, 1
                              if (microb(i,k)) then
                                if (sinlat(i,k) > tundra_s_limit) then
                                  ts = t0tile(i,k,npft)
                                  rr = ch4dd_tun(fourth_term) + ts * (ch4dd_tun(third_term) + ts * (ch4dd_tun(second_term) + ts * ch4dd_tun(first_term)))
                                  rr = rr * 3600.0
                                  if (todo(i,k,soil)) then
                                    rc(i,k,soil,j) = MAX(rr, 0.0)
                                  end if
                                else
                                  sm = smr(i,k)
                                  if (sm < smfrac1) then
                                    f_ch4_uptk = sm / smfrac1
                                  else
                                    if (sm > smfrac2) then
                                      f_ch4_uptk = (smfrac3 - sm) / smfrac4
                                    else
                                      f_ch4_uptk = 1.0
                                    end if
                                  end if
                                  if (ukca_config%l_ukca_emsdrvn_ch4) then
                                    f_ch4_uptk = MIN(MAX(f_ch4_uptk, 0.0), 1.0)
                                  else
                                    f_ch4_uptk = MAX(f_ch4_uptk, 0.0)
                                  end if
                                  if (todo(i,k,soil)) then
                                    rc(i,k,soil,j) = ch4_up_flux(soil) * f_ch4_uptk
                                  end if
                                end if
                              else
                                if (ukca_config%l_ukca_emsdrvn_ch4) then
                                  rc(i,k,n,j) = 1.0 / r_null
                                else
                                  rc(i,k,n,j) = r_null
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        end if
                        if (ukca_config%l_ukca_emsdrvn_ch4) then
                          !$acc kernels
                          !$acc loop gang vector independent collapse(3)
                          do n = 1, ntype, 1
                            do k = 1, rows, 1
                              do i = 1, row_length, 1
                                if (todo(i,k,n)) then
                                  if (microb(i,k)) then
                                    rr = rc(i,k,n,j)
                                    if (rr > 0.0) then
                                      if (n == soil) then
                                        rc(i,k,n,j) = (p0(i,k) * mml) / (rmol * t0tile(i,k,n) * rr * tar_scaling)
                                      else
                                        if (n == urban .OR. n == lake .OR. n >= ice) then
                                          rc(i,k,n,j) = r_null
                                        else
                                          rc(i,k,n,j) = (p0(i,k) * mml) / (rmol * t0tile(i,k,n) * rr)
                                        end if
                                      end if
                                    else
                                      rc(i,k,n,j) = r_null
                                    end if
                                  else
                                    rc(i,k,n,j) = r_null
                                  end if
                                  if (n == urban) then
                                    rc(i,k,n,j) = r_null
                                  end if
                                  if (n == lake) then
                                    rc(i,k,n,j) = r_null
                                  end if
                                  if (n >= ice) then
                                    rc(i,k,n,j) = r_null
                                  end if
                                else
                                  rc(i,k,n,j) = r_null
                                end if
                              enddo
                            enddo
                          enddo
                          !$acc end kernels
                        else
                          !$acc kernels
                          !$acc loop seq
                          do k = 1, rows, 1
                            !$acc loop gang vector independent collapse(2)
                            do n = 1, npft, 1
                              do i = 1, row_length, 1
                                if (todo(i,k,n) .AND. microb(i,k)) then
                                  rr = rc(i,k,n,j)
                                  if (rr > 0.0) then
                                    rc(i,k,n,j) = p0(i,k) * mml / (rmol * t0tile(i,k,n) * rr)
                                  else
                                    rc(i,k,n,j) = r_null
                                  end if
                                end if
                              enddo
                            enddo
                            n = soil
                            !$acc loop gang vector independent
                            do i = 1, row_length, 1
                              if (todo(i,k,n) .AND. microb(i,k)) then
                                rr = rc(i,k,n,j)
                                if (rr > 0.0) then
                                  rc(i,k,n,j) = p0(i,k) * mml / (rmol * t0tile(i,k,n) * rr * tar_scaling)
                                else
                                  rc(i,k,n,j) = r_null
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        end if
                      else
                        if (speci(nldepd(j)) == 'CO        ') then
                          if (.NOT.l_fix_ukca_h2dd_x .OR. ntype == 9) then
                            if (l_fix_ukca_h2dd_x) then
                              n = shrub
                            else
                              n = npft
                            end if
                            !$acc kernels
                            !$acc loop gang vector independent collapse(2)
                            do k = 1, rows, 1
                              do i = 1, row_length, 1
                                if (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) then
                                  if (todo(i,k,n)) then
                                    rc(i,k,n,j) = r_tundra(r_tundra_co)
                                  end if
                                end if
                              enddo
                            enddo
                            !$acc end kernels
                          else
                            if (ntype == 13 .OR. ntype == 17 .OR. ntype == 27) then
                              !$acc kernels
                              !$acc loop gang vector independent collapse(3)
                              do n = shrub_dec, shrub_eg, 1
                                do k = 1, rows, 1
                                  do i = 1, row_length, 1
                                    if (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) then
                                      if (todo(i,k,n)) then
                                        rc(i,k,n,j) = r_tundra(r_tundra_co)
                                      end if
                                    end if
                                  enddo
                                enddo
                              enddo
                              !$acc end kernels
                            else
                              ! PSyclone CodeBlock (unsupported code) reason:
                              !  - Unsupported statement: Write_Stmt
                              WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                              errcode = 1
                              call ereport(routinename, errcode, cmessage)
                            end if
                          end if
                          !$acc kernels
                          !$acc loop gang vector independent collapse(2)
                          do k = 1, rows, 1
                            do i = 1, row_length, 1
                              if (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) then
                                if (todo(i,k,soil)) then
                                  rc(i,k,soil,j) = r_tundra(r_tundra_co)
                                end if
                              end if
                            enddo
                          enddo
                          !$acc end kernels
                        else
                          if ((speci(nldepd(j)) == 'HNO3      ' .OR. speci(nldepd(j)) == 'HONO2     ' .OR. speci(nldepd(j)) == 'ISON      ') .OR. (speci(nldepd(j)) == 'HCl       ' .OR. speci(nldepd(j)) == 'HOCl      ' .OR. speci(nldepd(j)) == 'HBr       ' .OR. speci(nldepd(j)) == 'HOBr      ') .AND. l_fix_improve_drydep) then
                            if (.NOT.l_fix_ukca_h2dd_x .OR. (ntype == 9 .OR. ntype == 13 .OR. ntype == 17)) then
                              if (l_fix_ukca_h2dd_x) then
                                n = ice
                              else
                                n = ntype
                              end if
                              !$acc kernels
                              !$acc loop gang vector independent collapse(2)
                              do k = 1, rows, 1
                                do i = 1, row_length, 1
                                  if (todo(i,k,n)) then
                                    ts = MAX(t0tile(i,k,n), 252.0)
                                    rr = hno3dd_ice(third_term) + ts * (hno3dd_ice(second_term) + ts * hno3dd_ice(first_term))
                                    rc(i,k,n,j) = MAX(rr, 10.0)
                                  end if
                                enddo
                              enddo
                              !$acc end kernels
                            else
                              if (ntype == 27) then
                                !$acc kernels
                                !$acc loop gang vector independent collapse(3)
                                do n = ice, ice + n_elev_ice, 1
                                  do k = 1, rows, 1
                                    do i = 1, row_length, 1
                                      if (todo(i,k,n)) then
                                        ts = MAX(t0tile(i,k,n), 252.0)
                                        rr = hno3dd_ice(third_term) + ts * (hno3dd_ice(second_term) + ts * hno3dd_ice(first_term))
                                        rc(i,k,n,j) = MAX(rr, 10.0)
                                      end if
                                    enddo
                                  enddo
                                enddo
                                !$acc end kernels
                              else
                                ! PSyclone CodeBlock (unsupported code) reason:
                                !  - Unsupported statement: Write_Stmt
                                WRITE(cmessage, '(A)') 'UKCA does not handle flexible tiles yet: ' // 'Dry deposition needs extending. ' // 'Please use standard tile configuration'
                                errcode = 1
                                call ereport(routinename, errcode, cmessage)
                              end if
                            end if
                          else
                            if ((speci(nldepd(j)) == 'ORGNIT    ' .OR. speci(nldepd(j)) == 'BSOA      ' .OR. speci(nldepd(j)) == 'ASOA      ' .OR. speci(nldepd(j)) == 'ISOSOA    ') .OR. speci(nldepd(j)) == 'Sec_Org   ' .AND. l_fix_improve_drydep .OR. speci(nldepd(j)) == 'SEC_ORG_I ' .AND. l_fix_improve_drydep) then
                              !$acc kernels
                              !$acc loop gang vector independent collapse(3)
                              do n = 1, urban, 1
                                do k = 1, rows, 1
                                  do i = 1, row_length, 1
                                    if (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) then
                                      rr = 1.0 / so4_vd(i,k)
                                      rc(i,k,n,j) = rr
                                    end if
                                  enddo
                                enddo
                              enddo
                              n = soil
                              !$acc loop gang vector independent collapse(2)
                              do k = 1, rows, 1
                                do i = 1, row_length, 1
                                  if (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) then
                                    rr = 1.0 / so4_vd(i,k)
                                    rc(i,k,n,j) = rr
                                  end if
                                enddo
                              enddo
                              !$acc end kernels
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_surfddr

end module ukca_surfddr_mod
