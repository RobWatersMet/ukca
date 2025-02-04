module ukca_drydep_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DRYDEP_MOD'
  public

  contains
  subroutine ukca_drydep(nlev, dryrt, n_points)
    use asad_mod, only : dpd, jpdd, ndepd, nldepd
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nlev
    integer, intent(in) :: n_points
    real, dimension(n_points,jpdd), intent(in) :: dryrt
    integer :: i
    integer :: js
    integer :: nspec
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DRYDEP'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (nlev == 1 .OR. ukca_config%l_ukca_intdd) then
      !$acc kernels
      !$acc loop seq
      do js = 1, ndepd, 1
        nspec = nldepd(js)
        !$acc loop gang vector independent
        do i = 1, n_points, 1
          dpd(i,nspec) = dryrt(i,js)
        enddo
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop seq
      do js = 1, ndepd, 1
        nspec = nldepd(js)
        !$acc loop gang vector independent
        do i = 1, n_points, 1
          dpd(i,nspec) = 0.0
        enddo
      enddo
      !$acc end kernels
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_drydep

end module ukca_drydep_mod
