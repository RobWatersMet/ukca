module ukca_inddep_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_INDDEP_MOD'
  public

  contains
  subroutine ukca_inddep()
    use asad_mod, only : depvel, jddepc, jddept, jpdd
    use ukca_chem_defs_mod, only : depvel_defs
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: errcode
    integer :: ns
    integer :: nt
    integer :: nc
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_INDDEP'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (SIZE(depvel_defs) /= SIZE(depvel)) then
      cmessage(:) = 'sizes of depvel_defs and depvel are inconsistent'
      errcode = 1
      call ereport('UKCA_INDDEP', errcode, cmessage)
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do ns = 1, jpdd, 1
      do nc = 1, jddepc, 1
        do nt = 1, jddept, 1
          depvel(nt,nc,ns) = depvel_defs(nt,nc,ns)
        enddo
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_inddep

end module ukca_inddep_mod
