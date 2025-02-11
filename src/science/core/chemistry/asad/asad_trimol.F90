module asad_trimol_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_TRIMOL_MOD'
  public

  contains
  subroutine asad_trimol(n_points)
    use asad_mod, only : at, f, jpcspf, jptk, ntrkx, peps, rk, specf, spt, t, t300, tnd, wp
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ukca_um_legacy_mod, only : exp_v, oneover_v, powr_v
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, save :: ih2o
    integer :: iho2
    integer :: in2o5
    integer :: ino2no3
    integer :: idmsoho2_a
    integer :: idmsoho2_b
    integer :: j
    integer :: jl
    integer :: jtr
    integer :: jr
    real :: zo
    real :: zi
    real :: zfc
    real :: zr
    real, allocatable, dimension(:), save :: nf
    real :: nf2
    logical, save :: first = .true.
    logical, save :: first_pass = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_TRIMOL'
    real, dimension(n_points) :: tmp
    real, dimension(n_points) :: tmp_out
    real, dimension(n_points) :: inv_t_power4
    real, dimension(n_points) :: inv_t_power7
    real, dimension(n_points) :: inv_t_power4_out
    real, dimension(n_points) :: inv_t_power7_out
    real, dimension(n_points) :: t_power3
    real, dimension(n_points) :: t_power6
    real, dimension(n_points) :: inv_t
    integer :: idx
    integer :: idx_1
    integer :: idx_2

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    iho2 = 0
    in2o5 = 0
    ino2no3 = 0
    idmsoho2_a = 0
    idmsoho2_b = 0
    call oneover_v(n_points, t, inv_t)
    if (first_pass) then
      if (first) then
        ALLOCATE(nf(1:jptk))
        !$acc kernels
        !$acc loop gang vector independent
        do idx = LBOUND(nf, dim=1), UBOUND(nf, dim=1), 1
          nf(idx) = 1.0
        enddo
        !$acc end kernels
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
        WHERE (at(:, 1) > 1E-3) nf(:) = 0.75 - 1.27 * LOG10(at(:, 1))
        !$acc kernels
        ih2o = 0
        !$acc loop seq
        do jtr = 1, jpcspf, 1
          if (specf(jtr) == 'H2O       ') then
            ih2o = jtr
          end if
        enddo
        first = .false.
        !$acc end kernels
      end if
      first_pass = .false.
    end if
    do j = 1, jptk, 1
      jr = ntrkx(j)
      if (spt(j,1) == 'HO2    ' .AND. spt(j,2) == 'HO2    ') then
        iho2 = jr
      end if
      in2o5 = 0
      if (spt(j,1) == 'N2O5      ') then
        in2o5 = jr
      end if
      ino2no3 = 0
      if (spt(j,1) == 'NO2       ' .AND. spt(j,2) == 'NO3       ' .OR. spt(j,2) == 'NO2       ' .AND. spt(j,1) == 'NO3       ') then
        ino2no3 = jr
      end if
      idmsoho2_a = 0
      if (spt(j,1) == 'DMS       ' .AND. spt(j,2) == 'OH        ' .AND. spt(j,3) == 'DMSO      ' .AND. spt(j,4) == 'HO2       ') then
        idmsoho2_a = jr
      end if
      idmsoho2_b = 0
      if (spt(j,1) == 'DMS       ' .AND. spt(j,2) == 'OH        ' .AND. spt(j,3) == 'SO2       ' .AND. spt(j,4) == 'MeOO      ') then
        idmsoho2_b = jr
      end if
      call powr_v(n_points, t300, at(j,3), t_power3)
      call powr_v(n_points, t300, at(j,6), t_power6)
      inv_t_power4(:) = -at(j,4) * inv_t(:)
      inv_t_power7(:) = -at(j,7) * inv_t(:)
      call exp_v(n_points, inv_t_power4, inv_t_power4_out)
      call exp_v(n_points, inv_t_power7, inv_t_power7_out)
      !$acc kernels
      !$acc loop gang vector independent
      do jl = 1, n_points, 1
        zo = at(j,2) * t_power3(jl) * inv_t_power4_out(jl) * tnd(jl)
        zi = at(j,5) * t_power6(jl) * inv_t_power7_out(jl)
        if (zo < peps) then
          rk(jl,jr) = zi
        else
          if (zi < peps) then
            rk(jl,jr) = zo
          else
            if (idmsoho2_a /= 0 .OR. idmsoho2_b /= 0) then
              zo = zo * 0.2095
              zi = zi * 0.2095 * tnd(jl)
              rk(jl,jr) = zo / (1 + zi)
            else
              nf2 = nf(j)
              if (at(j,1) <= 1.0) then
                zfc = at(j,1)
              else
                if (in2o5 /= 0 .OR. ino2no3 /= 0) then
                  zfc = 2.5 * EXP(-1950.0 / t(jl)) + 0.9 * EXP(-t(jl) / at(j,1))
                else
                  zfc = EXP(-t(jl) / at(j,1))
                  nf2 = 0.75 - 1.27 * LOG10(zfc)
                end if
              end if
              zr = zo / zi
              rk(jl,jr) = (zo / (1.0 + zr)) * zfc ** (1.0 / (1.0 + (LOG10(zr) / nf2) ** 2))
            end if
          end if
        end if
      enddo
      !$acc end kernels
    enddo
    if (ih2o /= 0 .AND. iho2 /= 0) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_1 = 1, n_points, 1
        tmp(idx_1) = 2200.0 * inv_t(idx_1)
      enddo
      !$acc end kernels
      call exp_v(n_points, tmp, tmp_out)
      !$acc kernels
      !$acc loop gang vector independent
      do jl = 1, n_points, 1
        rk(jl,iho2) = rk(jl,iho2) * (1.0 + 1.4e-21 * f(jl,ih2o) * tmp_out(jl))
      enddo
      !$acc end kernels
    else
      if (ih2o == 0 .AND. iho2 /= 0) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_2 = 1, n_points, 1
          tmp(idx_2) = 2200.0 * inv_t(idx_2)
        enddo
        !$acc end kernels
        call exp_v(n_points, tmp, tmp_out)
        !$acc kernels
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          rk(jl,iho2) = rk(jl,iho2) * (1.0 + 1.4e-21 * wp(jl) * tnd(jl) * tmp_out(jl))
        enddo
        !$acc end kernels
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_trimol

end module asad_trimol_mod
