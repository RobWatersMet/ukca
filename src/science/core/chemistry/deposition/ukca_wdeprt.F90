module ukca_wdeprt_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_WDEPRT_MOD'
  public

  contains
  subroutine ukca_wdeprt(p_fieldda, model_levels, drain, crain, t, lat, tstep, wetrt, h_plus_2d_arr)
    use asad_mod, only : ddhr, dhr, jpdw, k298, kd298
    use ukca_um_legacy_mod, only : exp_v, log_v
    use ukca_config_constants_mod, only : rmol
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    real, parameter :: fac = 0.1
    real, parameter :: fc = 0.3
    real, parameter :: csca = 4.7
    real, parameter :: dsca = 2.4
    real, parameter :: clw = 1.0e-6
    real, parameter :: nlat = 65.0
    real, parameter :: slat = -65.0
    real, parameter :: tmax = 273.15
    real, parameter :: tmin = 253.15
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: p_fieldda
    integer, intent(in) :: model_levels
    real, dimension(p_fieldda,model_levels), intent(in) :: drain
    real, dimension(p_fieldda,model_levels), intent(in) :: crain
    real, dimension(p_fieldda,model_levels), intent(in) :: t
    real, dimension(p_fieldda), intent(in) :: lat
    real, intent(in) :: tstep
    real, dimension(p_fieldda,model_levels), intent(in) :: h_plus_2d_arr
    real, dimension(p_fieldda,model_levels,jpdw), intent(out) :: wetrt
    integer :: first_point
    integer :: last_point
    integer :: i
    integer :: k
    integer :: ns
    real :: rgas
    real, dimension(p_fieldda,model_levels) :: tmp1
    real :: tmp2
    real :: kaq
    real :: rdw
    real :: flim
    real, dimension(p_fieldda,model_levels,jpdw) :: rwcon
    real, dimension(p_fieldda,model_levels,jpdw) :: rwdyn
    real, dimension(p_fieldda,model_levels,jpdw) :: hcoef
    real, dimension(p_fieldda,model_levels,jpdw) :: faq
    real, dimension(p_fieldda,model_levels) :: ddrain
    real, dimension(p_fieldda,model_levels) :: ccrain
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_WDEPRT'
    real, dimension(p_fieldda,model_levels,jpdw) :: rcw
    real, dimension(p_fieldda,model_levels,jpdw) :: tmp4
    real, dimension(p_fieldda,model_levels,jpdw) :: tmp5
    real, dimension(p_fieldda,model_levels,jpdw) :: rcw_out
    real, dimension(p_fieldda,model_levels,jpdw) :: tmp4_out
    real, dimension(p_fieldda,model_levels,jpdw) :: tmp5_out
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    first_point = 1
    last_point = p_fieldda
    rgas = rmol / 100.0
    !NO!$acc loop gang vector independent collapse(3)
    do idx = LBOUND(wetrt, dim=3), UBOUND(wetrt, dim=3), 1
      do idx_1 = LBOUND(wetrt, dim=2), UBOUND(wetrt, dim=2), 1
        do idx_2 = LBOUND(wetrt, dim=1), UBOUND(wetrt, dim=1), 1
          wetrt(idx_2,idx_1,idx) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_3 = LBOUND(rwdyn, dim=3), UBOUND(rwdyn, dim=3), 1
      do idx_4 = LBOUND(rwdyn, dim=2), UBOUND(rwdyn, dim=2), 1
        do idx_5 = LBOUND(rwdyn, dim=1), UBOUND(rwdyn, dim=1), 1
          rwdyn(idx_5,idx_4,idx_3) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_6 = LBOUND(rwcon, dim=3), UBOUND(rwcon, dim=3), 1
      do idx_7 = LBOUND(rwcon, dim=2), UBOUND(rwcon, dim=2), 1
        do idx_8 = LBOUND(rwcon, dim=1), UBOUND(rwcon, dim=1), 1
          rwcon(idx_8,idx_7,idx_6) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_9 = 1, model_levels, 1
      do idx_10 = 1, p_fieldda, 1
        ddrain(idx_10,idx_9) = drain(idx_10,idx_9)
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_11 = 1, model_levels, 1
      do idx_12 = 1, p_fieldda, 1
        ccrain(idx_12,idx_11) = crain(idx_12,idx_11)
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do k = 1, model_levels, 1
      do i = first_point, last_point, 1
        tmp1(i,k) = (298.0 - t(i,k)) / (t(i,k) * 298.0)
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do ns = 1, jpdw, 1
      do k = 1, model_levels, 1
        do i = first_point, last_point, 1
          tmp5(i,k,ns) = ddhr(ns,1) * tmp1(i,k)
          tmp4(i,k,ns) = dhr(ns) * tmp1(i,k)
        enddo
      enddo
    enddo
    !NO!$acc end kernels
    call exp_v(jpdw * model_levels * p_fieldda, tmp5, tmp5_out)
    call exp_v(jpdw * model_levels * p_fieldda, tmp4, tmp4_out)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(3)
    do ns = 1, jpdw, 1
      do k = 1, model_levels, 1
        do i = first_point, last_point, 1
          kaq = kd298(ns,1) * tmp5_out(i,k,ns)
          hcoef(i,k,ns) = k298(ns) * tmp4_out(i,k,ns) * (1.0 + kaq / h_plus_2d_arr(i,k))
          tmp2 = clw * hcoef(i,k,ns) * rgas * t(i,k)
          faq(i,k,ns) = tmp2 / (1.0 + tmp2)
          rcw(i,k,ns) = ccrain(i,k) * fac * csca * faq(i,k,ns)
          rcw(i,k,ns) = rcw(i,k,ns) / fc
          rcw(i,k,ns) = -tstep * rcw(i,k,ns)
        enddo
      enddo
    enddo
    !NO!$acc end kernels
    call exp_v(jpdw * model_levels * p_fieldda, rcw, rcw_out)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(3)
    do idx_13 = 1, jpdw, 1
      do idx_14 = 1, model_levels, 1
        do idx_15 = 1, p_fieldda, 1
          rcw(idx_15,idx_14,idx_13) = 1 - fc + fc * rcw_out(idx_15,idx_14,idx_13)
        enddo
      enddo
    enddo
    !NO!$acc end kernels
    call log_v(jpdw * model_levels * p_fieldda, rcw, rcw_out)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(3)
    do ns = 1, jpdw, 1
      do k = 1, model_levels, 1
        do i = first_point, last_point, 1
          rwcon(i,k,ns) = -rcw_out(i,k,ns) / tstep
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do ns = 1, jpdw, 1
      do k = 1, model_levels, 1
        do i = first_point, last_point, 1
          rdw = ddrain(i,k) * fac * dsca * faq(i,k,ns)
          if (lat(i) >= nlat .OR. lat(i) <= slat) then
            if (t(i,k) >= tmax) then
              flim = 1.0
            else
              if (t(i,k) < tmin) then
                flim = 0.0
              else
                flim = 1.0 + 0.05 * (t(i,k) - tmax)
              end if
            end if
            rdw = rdw * flim
          end if
          rwdyn(i,k,ns) = rdw
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do ns = 1, jpdw, 1
      do k = 1, model_levels, 1
        do i = first_point, last_point, 1
          wetrt(i,k,ns) = rwcon(i,k,ns) + rwdyn(i,k,ns)
        enddo
      enddo
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_wdeprt

end module ukca_wdeprt_mod
