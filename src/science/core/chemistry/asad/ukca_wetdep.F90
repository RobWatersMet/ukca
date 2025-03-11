module ukca_wetdep_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_WETDEP_MOD'
  public

  contains
  subroutine ukca_wetdep(wetrt, n_points)
    use asad_mod, only : dpw, jpdw, ndepw, nldepw
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, dimension(n_points,jpdw), intent(in) :: wetrt
    integer :: i
    integer :: js
    integer :: nspec
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_WETDEP'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop seq
    do js = 1, ndepw, 1
      nspec = nldepw(js)
      !NO!$acc loop gang vector independent
      do i = 1, n_points, 1
        dpw(i,nspec) = wetrt(i,js)
      enddo
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_wetdep

end module ukca_wetdep_mod
