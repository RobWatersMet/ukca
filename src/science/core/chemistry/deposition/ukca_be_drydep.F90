module ukca_be_drydep_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_BE_DRYDEP_MOD'
  public

  contains
  subroutine ukca_be_drydep(k, n_pnts, nlev_with_ddep2, dryrt, odryrt)
    use asad_mod, only : jpdd, jpspec, ndepd, nldepd
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer, parameter :: bottom_level = 1
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: k
    integer, intent(in) :: n_pnts
    integer, dimension(n_pnts), intent(in) :: nlev_with_ddep2
    real, dimension(n_pnts,jpdd), intent(in) :: dryrt
    real, dimension(n_pnts,jpspec), intent(out) :: odryrt
    integer :: i
    integer :: js
    integer :: nspec
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_BE_DRYDEP'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do js = 1, jpspec, 1
      do i = 1, n_pnts, 1
        odryrt(i,js) = 0.0
      enddo
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_intdd) then
      !$acc kernels
      !$acc loop seq
      do i = 1, n_pnts, 1
        if (k <= nlev_with_ddep2(i)) then
          !$acc loop seq
          do js = 1, ndepd, 1
            nspec = nldepd(js)
            odryrt(i,nspec) = dryrt(i,js)
          enddo
        end if
      enddo
      !$acc end kernels
    else
      if (k == bottom_level) then
        !$acc kernels
        !$acc loop seq
        do js = 1, ndepd, 1
          nspec = nldepd(js)
          !$acc loop gang vector independent
          do i = 1, n_pnts, 1
            odryrt(i,nspec) = dryrt(i,js)
          enddo
        enddo
        !$acc end kernels
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_be_drydep

end module ukca_be_drydep_mod
