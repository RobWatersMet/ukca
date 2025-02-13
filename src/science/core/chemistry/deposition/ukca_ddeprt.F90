module ukca_ddeprt_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DDEPRT_MOD'
  public

  contains
  subroutine ukca_ddeprt(p_fieldda, bl_levels, i_day_number, i_mon, i_hour, r_minute, secs_per_step, lon, lat, tanlat, dzl, z0m, u_s, temp, dryrt)
    use ukca_um_legacy_mod, only : vkman
    use ukca_constants, only : fxb, fxc, pi_over_180
    use asad_mod, only : depvel, jpdd
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    real, parameter :: midday = 12.0
    real, parameter :: temp_summer = 275.0
    real, parameter :: z0m_min = 1.0e-3
    real, parameter :: z0m_max = 1.0e-1
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: p_fieldda
    integer, intent(in) :: bl_levels
    integer, intent(in) :: i_day_number
    integer, intent(in) :: i_mon
    integer, intent(in) :: i_hour
    real, intent(in) :: r_minute
    real, intent(in) :: secs_per_step
    real, dimension(p_fieldda), intent(in) :: lon
    real, dimension(p_fieldda), intent(in) :: lat
    real, dimension(p_fieldda), intent(in) :: tanlat
    real, dimension(p_fieldda,bl_levels), intent(in) :: dzl
    real, dimension(p_fieldda), intent(in) :: z0m
    real, dimension(p_fieldda), intent(in) :: u_s
    real, dimension(p_fieldda), intent(in) :: temp
    real, dimension(p_fieldda,jpdd), intent(out) :: dryrt
    integer :: i
    integer :: ns
    integer :: iday
    integer :: isum
    integer :: icat
    real :: dawn
    real :: dusk
    real :: timel
    real :: vdep1
    real :: vdeph
    real :: tgmt
    real :: tan_declin
    real :: cs_hour_ang
    real, dimension(p_fieldda) :: tloc
    real, dimension(p_fieldda) :: daylen
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DDEPRT'
    integer :: idx
    integer :: idx_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(dryrt, dim=2), UBOUND(dryrt, dim=2), 1
      do idx_1 = LBOUND(dryrt, dim=1), UBOUND(dryrt, dim=1), 1
        dryrt(idx_1,idx) = 0.0
      enddo
    enddo
    !$acc end kernels
    do i = 1, p_fieldda, 1
      tgmt = REAL(i_hour) + r_minute / 60.0 + secs_per_step * 0.5 / 3600.0
      tloc(i) = tgmt + 24.0 * lon(i) / 360.0
      if (tloc(i) > 24.0) then
        tloc(i) = tloc(i) - 24.0
      end if
      tan_declin = TAN(fxb * SIN(pi_over_180 * (266.0 + i_day_number)))
      if (ABS(tan_declin) < EPSILON(0.0)) then
        cs_hour_ang = 0.0
      else
        cs_hour_ang = MAX(-1.0, MIN(1.0, -tanlat(i) * tan_declin))
      end if
      daylen(i) = fxc * ACOS(cs_hour_ang)
    enddo
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do ns = 1, jpdd, 1
      do i = 1, p_fieldda, 1
        dawn = midday - daylen(i) / 2.0
        dusk = midday + daylen(i) / 2.0
        timel = tloc(i)
        if (timel < dawn .OR. timel > dusk) then
          iday = 0
        else
          iday = 1
        end if
        if (temp(i) > temp_summer) then
          isum = 1
        else
          if (temp(i) <= temp_summer) then
            isum = 0
          end if
        end if
        if (z0m(i) < z0m_min) then
          icat = 1
        else
          if (z0m(i) > z0m_max) then
            icat = 2
          else
            if (z0m(i) <= z0m_max .OR. z0m(i) >= z0m_min) then
              icat = 3
            end if
          end if
        end if
        if (lat(i) >= 0.0) then
          if ((i_mon == 12 .OR. i_mon <= 3) .AND. lat(i) >= 45.0) then
            icat = 5
          end if
          if ((i_mon == 11 .OR. i_mon == 4) .AND. lat(i) >= 55.0) then
            icat = 5
          end if
          if ((i_mon == 10 .OR. i_mon == 5 .OR. i_mon == 6) .AND. lat(i) >= 65.0) then
            icat = 5
          end if
          if ((i_mon >= 7 .AND. i_mon <= 9) .AND. lat(i) >= 70.0) then
            icat = 5
          end if
        else
          if (lat(i) < 0.0) then
            if ((i_mon == 12 .OR. i_mon <= 4) .AND. lat(i) <= (-75.0)) then
              icat = 5
            end if
            if ((i_mon >= 5 .AND. i_mon <= 11) .AND. lat(i) <= (-70.0)) then
              icat = 5
            end if
          end if
        end if
        if (iday == 1 .AND. isum == 1) then
          vdep1 = depvel(1,icat,ns) / 100.0
        else
          if (iday == 0 .AND. isum == 1) then
            vdep1 = depvel(2,icat,ns) / 100.0
          else
            if (iday == 1 .AND. isum == 0) then
              vdep1 = depvel(4,icat,ns) / 100.0
            else
              if (iday == 0 .AND. isum == 0) then
                vdep1 = depvel(5,icat,ns) / 100.0
              end if
            end if
          end if
        end if
        vdeph = vdep1 / (1.0 + vdep1 * LOG(dzl(i,1) / 2.0) / (vkman * u_s(i)))
        dryrt(i,ns) = vdeph / dzl(i,1)
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_ddeprt

end module ukca_ddeprt_mod
