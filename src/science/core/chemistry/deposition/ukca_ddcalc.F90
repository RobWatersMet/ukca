module ukca_ddcalc_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DDCALC_MOD'
  public

  contains
  subroutine ukca_ddcalc(row_length, rows, bl_levels, ntype, npft, timestep, dzl, zbl, gsf, ra, rb, rc, nlev_with_ddep, zdryrt, len_stashwork, stashwork)
    use asad_mod, only : jpdd, ndepd, nldepd, speci
    use parkind1, only : jpim, jprb
    use ukca_config_specification_mod, only : ukca_config
    use ukca_um_legacy_mod, only : copydiag, len_stlist, num_stash_pseudo, pdims, set_pseudo_list, sf, si, si_last, stash_pseudo_levels, stindex, stlist
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: bl_levels
    integer, intent(in) :: ntype
    integer, intent(in) :: npft
    real, intent(in) :: timestep
    real, dimension(row_length,rows), intent(in) :: zbl
    real, dimension(row_length,rows,bl_levels), intent(in) :: dzl
    real, dimension(row_length,rows,ntype), intent(in) :: gsf
    real, dimension(row_length,rows,ntype), intent(in) :: ra
    real, dimension(row_length,rows,jpdd), intent(in) :: rb
    real, dimension(row_length,rows,ntype,jpdd), intent(in) :: rc
    integer, intent(in) :: len_stashwork
    real, dimension(len_stashwork), intent(inout) :: stashwork
    integer, dimension(row_length,rows), intent(out) :: nlev_with_ddep
    real, dimension(row_length,rows,jpdd), intent(out) :: zdryrt
    integer :: i
    integer :: j
    integer :: loop
    integer :: n
    integer :: pslevel
    integer :: pslevel_out
    integer :: si_start
    integer :: si_stop
    logical, dimension(ntype) :: plltile
    real :: dd
    real, dimension(row_length,rows) :: layer_depth
    real, dimension(row_length,rows,ntype,jpdd) :: vd
    real, save :: r_nodep = 1.0e40
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DDCALC'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do loop = 1, jpdd, 1
      do j = 1, rows, 1
        do i = 1, row_length, 1
          zdryrt(i,j,loop) = 0.0
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(4)
    do loop = 1, jpdd, 1
      do n = 1, ntype, 1
        do j = 1, rows, 1
          do i = 1, row_length, 1
            vd(i,j,n,loop) = 0.0
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_ddep_lev1) then
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx = LBOUND(nlev_with_ddep, dim=2), UBOUND(nlev_with_ddep, dim=2), 1
        do idx_1 = LBOUND(nlev_with_ddep, dim=1), UBOUND(nlev_with_ddep, dim=1), 1
          nlev_with_ddep(idx_1,idx) = 1
        enddo
      enddo
      !$acc loop gang vector independent collapse(2)
      do idx_2 = LBOUND(layer_depth, dim=2), UBOUND(layer_depth, dim=2), 1
        do idx_3 = LBOUND(layer_depth, dim=1), UBOUND(layer_depth, dim=1), 1
          layer_depth(idx_3,idx_2) = dzl(idx_3,idx_2,1)
        enddo
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do j = 1, rows, 1
        do i = 1, row_length, 1
          layer_depth(i,j) = dzl(i,j,1)
          nlev_with_ddep(i,j) = 1
        enddo
      enddo
      !$acc loop seq
      do loop = 2, bl_levels, 1
        !$acc loop gang vector independent collapse(2)
        do j = 1, rows, 1
          do i = 1, row_length, 1
            dd = layer_depth(i,j) + dzl(i,j,loop)
            if (dd < zbl(i,j)) then
              layer_depth(i,j) = dd
              nlev_with_ddep(i,j) = loop
            end if
          enddo
        enddo
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(4)
    do loop = 1, ndepd, 1
      do n = 1, npft, 1
        do j = 1, rows, 1
          do i = 1, row_length, 1
            if (rc(i,j,n,loop) < r_nodep .AND. gsf(i,j,n) > 0.0) then
              vd(i,j,n,loop) = 1.0 / (ra(i,j,n) + 2.0 * rb(i,j,loop) + rc(i,j,n,loop))
            end if
          enddo
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(4)
    do loop = 1, ndepd, 1
      do n = npft + 1, ntype, 1
        do j = 1, rows, 1
          do i = 1, row_length, 1
            if (rc(i,j,n,loop) < r_nodep .AND. gsf(i,j,n) > 0.0) then
              vd(i,j,n,loop) = 1.0 / (ra(i,j,n) + rb(i,j,loop) + rc(i,j,n,loop))
            end if
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_scale_ppe) then
      !$acc kernels
      !$acc loop gang vector independent
      do j = 1, ndepd, 1
        if (speci(nldepd(j)) == 'SO2       ') then
          vd(:,:,:,j) = vd(:,:,:,j) * ukca_config%dry_depvel_so2_scaling
        end if
      enddo
      !$acc end kernels
    end if
    if (ukca_config%l_enable_diag_um) then
      do loop = 1, ndepd, 1
        if (speci(nldepd(loop)) == 'CH4       ') then
          if (sf(432,50)) then
            call set_pseudo_list(ntype, len_stlist, stlist(:,stindex(1,432,50,1)), plltile, stash_pseudo_levels, num_stash_pseudo)
            pslevel_out = 0
            do pslevel = 1, ntype, 1
              if (plltile(pslevel)) then
                pslevel_out = pslevel_out + 1
                si_start = si(432,50,1) + (pslevel_out - 1) * pdims%i_end * pdims%j_end
                si_stop = si(432,50,1) + pslevel_out * pdims%i_end * pdims%j_end - 1
                call copydiag(stashwork(si_start:si_stop), vd(:,:,pslevel,loop), pdims%i_end, pdims%j_end)
              end if
            enddo
          end if
        else
          if (speci(nldepd(loop)) == 'O3        ') then
            if (sf(433,50)) then
              call set_pseudo_list(ntype, len_stlist, stlist(:,stindex(1,433,50,1)), plltile, stash_pseudo_levels, num_stash_pseudo)
              pslevel_out = 0
              do pslevel = 1, ntype, 1
                if (plltile(pslevel)) then
                  pslevel_out = pslevel_out + 1
                  si_start = si(433,50,1) + (pslevel_out - 1) * pdims%i_end * pdims%j_end
                  si_stop = si(433,50,1) + pslevel_out * pdims%i_end * pdims%j_end - 1
                  call copydiag(stashwork(si_start:si_stop), vd(:,:,pslevel,loop), pdims%i_end, pdims%j_end)
                end if
              enddo
            end if
          else
            if (speci(nldepd(loop)) == 'HONO2     ') then
              if (sf(434,50)) then
                call set_pseudo_list(ntype, len_stlist, stlist(:,stindex(1,434,50,1)), plltile, stash_pseudo_levels, num_stash_pseudo)
                pslevel_out = 0
                do pslevel = 1, ntype, 1
                  if (plltile(pslevel)) then
                    pslevel_out = pslevel_out + 1
                    si_start = si(434,50,1) + (pslevel_out - 1) * pdims%i_end * pdims%j_end
                    si_stop = si(434,50,1) + pslevel_out * pdims%i_end * pdims%j_end - 1
                    call copydiag(stashwork(si_start:si_stop), vd(:,:,pslevel,loop), pdims%i_end, pdims%j_end)
                  end if
                enddo
              end if
            end if
          end if
        end if
      enddo
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(4)
    do loop = 1, ndepd, 1
      do n = 1, ntype, 1
        do j = 1, rows, 1
          do i = 1, row_length, 1
            if (vd(i,j,n,loop) > 0.0) then
              zdryrt(i,j,loop) = zdryrt(i,j,loop) + gsf(i,j,n) * (1.0 - EXP(-vd(i,j,n,loop) * timestep / layer_depth(i,j)))
            end if
          enddo
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do loop = 1, ndepd, 1
      do j = 1, rows, 1
        do i = 1, row_length, 1
          zdryrt(i,j,loop) = -LOG(1.0 - zdryrt(i,j,loop)) / timestep
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ukca_config%l_enable_diag_um) then
      do loop = 1, ndepd, 1
        if (speci(nldepd(loop)) == 'CH4       ') then
          if (sf(435,50)) then
            call copydiag(stashwork(si(435,50,1):si_last(435,50,1)), zdryrt(:,:,loop), row_length, rows)
          end if
        else
          if (speci(nldepd(loop)) == 'O3        ') then
            if (sf(436,50)) then
              call copydiag(stashwork(si(436,50,1):si_last(436,50,1)), zdryrt(:,:,loop), row_length, rows)
            end if
          else
            if (speci(nldepd(loop)) == 'HONO2     ') then
              if (sf(437,50)) then
                call copydiag(stashwork(si(437,50,1):si_last(437,50,1)), zdryrt(:,:,loop), row_length, rows)
              end if
            end if
          end if
        end if
      enddo
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_ddcalc

end module ukca_ddcalc_mod
