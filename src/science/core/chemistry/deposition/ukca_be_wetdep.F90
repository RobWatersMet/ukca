module ukca_be_wetdep_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_BE_WETDEP_MOD'
  public

  contains
  subroutine ukca_be_wetdep(n_pnts, wetrt, owetrt)
    use asad_mod, only : jpdw, jpspec, ndepw, nldepw
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_pnts
    real, dimension(n_pnts,jpdw), intent(in) :: wetrt
    real, dimension(n_pnts,jpspec), intent(out) :: owetrt
    integer :: i
    integer :: js
    integer :: nspec
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_BE_WETDEP'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do js = 1, jpspec, 1
      do i = 1, n_pnts, 1
        owetrt(i,js) = 0.0
      enddo
    enddo
    !$acc loop seq
    do js = 1, ndepw, 1
      nspec = nldepw(js)
      !$acc loop gang vector independent
      do i = 1, n_pnts, 1
        owetrt(i,nspec) = wetrt(i,js)
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_be_wetdep

end module ukca_be_wetdep_mod
