module asad_posthet_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_POSTHET_MOD'
  public

  contains
  subroutine asad_posthet()
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: icode
    CHARACTER(LEN = errormessagelength) :: cmessage
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_POSTHET'
    real(kind=jprb) :: zhook_handle

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    cmessage(:) = 'Routine should not be callable'
    icode = 1
    call ereport(routinename, icode, cmessage)
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_posthet

end module asad_posthet_mod
