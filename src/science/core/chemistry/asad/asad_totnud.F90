module asad_totnud_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_TOTNUD_MOD'
  public

  contains
  subroutine asad_totnud(n_points)
    use asad_mod, only : p, pmin, pmintnd, t, tnd
    use ukca_config_constants_mod, only : boltzmann
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer :: jl
    real :: zb
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_TOTNUD'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    zb = boltzmann * 1.0e6
    !$acc loop gang vector independent
    do jl = 1, n_points, 1
      tnd(jl) = p(jl) / (zb * t(jl))
      pmintnd(jl) = pmin * tnd(jl)
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_totnud

end module asad_totnud_mod
