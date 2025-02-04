module ukca_photol_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_PHOTOL_MOD'
  public

  contains
  subroutine ukca_photol(prt, n_points)
    use asad_mod, only : jppj, nprkx, rk
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, dimension(n_points,jppj), intent(in) :: prt
    integer :: jl
    integer :: jr
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_PHOTOL'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do jl = 1, n_points, 1
      do jr = 1, jppj, 1
        rk(jl,nprkx(jr)) = prt(jl,jr)
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_photol

end module ukca_photol_mod
