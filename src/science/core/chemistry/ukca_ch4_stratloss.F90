module ukca_ch4_stratloss_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CH4_STRATLOSS_MOD'
  public

  contains
  subroutine ukca_ch4_stratloss(n_be_calls, n_pnts, vol, dts, y, stratloss)
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    use ukca_config_constants_mod, only : avogadro
    use nvtx
    real, parameter :: loss_rate = 2.0e-7
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_be_calls
    integer, intent(in) :: n_pnts
    real, intent(in) :: dts
    real, dimension(n_pnts), intent(in) :: vol
    real, dimension(n_pnts), intent(inout) :: y
    real, dimension(n_pnts), intent(inout) :: stratloss
    integer :: n
    real, dimension(n_pnts) :: p
    real, dimension(n_pnts) :: l
    real, dimension(n_pnts) :: yp
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CH4_STRATLOSS'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop seq
    do n = 1, n_be_calls, 1
      !$acc loop gang vector independent
      do idx = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx) = y(idx)
      enddo
      !$acc loop gang vector independent
      do idx_1 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
        p(idx_1) = 0.0
      enddo
      !$acc loop gang vector independent
      do idx_2 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
        l(idx_2) = loss_rate
      enddo
      !$acc loop gang vector independent
      do idx_3 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_3) = (yp(idx_3) + dts * p(idx_3)) / (1.0 + dts * l(idx_3))
      enddo
      !$acc loop gang vector independent
      do idx_4 = LBOUND(stratloss, dim=1), UBOUND(stratloss, dim=1), 1
        stratloss(idx_4) = stratloss(idx_4) + loss_rate * y(idx_4) * dts * vol(idx_4) / avogadro
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_ch4_stratloss

end module ukca_ch4_stratloss_mod
