module ukca_sediment_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_SEDIMENT_MOD'
  public

  contains
  subroutine ukca_sediment(rows, row_length, model_levels, shno3, qcf, r_theta_levels, mass, secs_per_step, l_stratosphere)
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    real, parameter :: fv1 = 0.000463
    real, parameter :: fv2 = 0.01736
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: rows
    integer, intent(in) :: row_length
    integer, intent(in) :: model_levels
    real, intent(in) :: secs_per_step
    real, dimension(row_length,rows,0:model_levels), intent(in) :: r_theta_levels
    real, dimension(row_length,rows,model_levels), intent(in) :: mass
    real, dimension(row_length,rows,model_levels), intent(in) :: qcf
    logical, dimension(row_length,rows,model_levels), intent(in) :: l_stratosphere
    real, dimension(row_length,rows,model_levels), intent(inout) :: shno3
    integer :: i
    integer :: j
    integer :: k
    real, dimension(row_length,rows) :: thick
    real, dimension(row_length,rows) :: sed_mmr
    real, dimension(row_length,rows) :: frac
    real, dimension(row_length,rows) :: psc1
    real, dimension(row_length,rows) :: psc2
    real, dimension(row_length,rows,model_levels) :: sh2o
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SEDIMENT'
    integer :: widx3
    integer :: widx2
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
    integer :: widx2_1
    integer :: widx1_1
    integer :: widx2_2
    integer :: widx1_2
    integer :: widx2_3
    integer :: widx1_3

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do k = 1, model_levels, 1
      do j = 1, rows, 1
        do i = 1, row_length, 1
          sh2o(i,j,k) = qcf(i,j,k)
        enddo
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do widx3 = 1, model_levels, 1
      do widx2 = 1, rows, 1
        do widx1 = 1, row_length, 1
          if (.NOT.l_stratosphere(widx1,widx2,widx3)) then
            sh2o(widx1,widx2,widx3) = 0.0
          end if
        enddo
      enddo
    enddo
    !$acc end kernels
    do k = 2, model_levels, 1
      if (MAXVAL(shno3(:,:,k)) > 0.0) then
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx = 1, rows, 1
          do idx_1 = 1, row_length, 1
            psc1(idx_1,idx) = 0.0
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_2 = 1, rows, 1
          do idx_3 = 1, row_length, 1
            psc2(idx_3,idx_2) = 0.0
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do widx2_1 = 1, SIZE(shno3, dim=2), 1
          do widx1_1 = 1, SIZE(shno3, dim=1), 1
            if (shno3(widx1_1,widx2_1,k) > 0.0) then
              psc1(widx1_1,widx2_1) = 1.0
            end if
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do widx2_2 = 1, SIZE(sh2o, dim=2), 1
          do widx1_2 = 1, SIZE(sh2o, dim=1), 1
            if (sh2o(widx1_2,widx2_2,k) > 0.0) then
              psc2(widx1_2,widx2_2) = 1.0
            end if
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_4 = 1, rows, 1
          do idx_5 = 1, row_length, 1
            thick(idx_5,idx_4) = r_theta_levels(idx_5,idx_4,k) - r_theta_levels(idx_5,idx_4,k - 1)
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_6 = 1, rows, 1
          do idx_7 = 1, row_length, 1
            frac(idx_7,idx_6) = (psc1(idx_7,idx_6) * fv1 + psc2(idx_7,idx_6) * fv2) * secs_per_step / thick(idx_7,idx_6)
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do widx2_3 = 1, rows, 1
          do widx1_3 = 1, row_length, 1
            if (frac(widx1_3,widx2_3) > 1.0) then
              frac(widx1_3,widx2_3) = 1.0
            end if
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_8 = 1, rows, 1
          do idx_9 = 1, row_length, 1
            sed_mmr(idx_9,idx_8) = shno3(idx_9,idx_8,k) * frac(idx_9,idx_8)
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_10 = LBOUND(shno3, dim=2), UBOUND(shno3, dim=2), 1
          do idx_11 = LBOUND(shno3, dim=1), UBOUND(shno3, dim=1), 1
            shno3(idx_11,idx_10,k) = shno3(idx_11,idx_10,k) - sed_mmr(idx_11,idx_10)
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_12 = LBOUND(shno3, dim=2), UBOUND(shno3, dim=2), 1
          do idx_13 = LBOUND(shno3, dim=1), UBOUND(shno3, dim=1), 1
            shno3(idx_13,idx_12,k - 1) = shno3(idx_13,idx_12,k - 1) + sed_mmr(idx_13,idx_12) * mass(idx_13,idx_12,k) / mass(idx_13,idx_12,k - 1)
          enddo
        enddo
        !$acc end kernels
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_sediment

end module ukca_sediment_mod
