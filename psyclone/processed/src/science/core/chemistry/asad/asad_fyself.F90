module asad_fyself_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_FYSELF_MOD'
  public

  contains
  subroutine asad_fyself(n_points)
    use asad_mod, only : jpnr, nlstst, nspi, nstst, qa, rk
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer :: j
    integer :: jl
    integer :: jr
    integer :: js
    integer :: isp
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_FYSELF'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    do j = 1, nstst, 1
      js = nlstst(j)
      do jl = 1, n_points, 1
        qa(jl,js) = 0.0
      enddo
    enddo
    do jr = 1, jpnr, 1
      isp = nspi(jr,1)
      if (isp == nspi(jr,2)) then
        do jl = 1, n_points, 1
          qa(jl,isp) = qa(jl,isp) + 2.0 * rk(jl,jr)
        enddo
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine asad_fyself

end module asad_fyself_mod
