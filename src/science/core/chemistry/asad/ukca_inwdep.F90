module ukca_inwdep_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_INWDEP_MOD'
  public

  contains
  subroutine ukca_inwdep()
    use asad_mod, only : ct_ddhr, ct_dhr, ct_k298, ct_kd298, ddhr, dhr, jpdw, k298, kd298
    use ukca_chem_defs_mod, only : henry_defs
    use ukca_chem_offline, only : henry_defs_const, nwet_constant
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: errcode
    integer :: ns
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_INWDEP'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (SIZE(henry_defs) /= jpdw * 6) then
      cmessage(:) = 'jpdw and henry_defs are inconsistent'
      errcode = 1
      call ereport('UKCA_INWDEP', errcode, cmessage)
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do ns = 1, jpdw, 1
      k298(ns) = henry_defs(1,ns)
      dhr(ns) = henry_defs(2,ns)
      kd298(ns,1) = henry_defs(3,ns)
      ddhr(ns,1) = henry_defs(4,ns)
      kd298(ns,2) = henry_defs(5,ns)
      ddhr(ns,2) = henry_defs(6,ns)
    enddo
    !NO!$acc end kernels
    if (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) then
      if (nwet_constant > 0) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do ns = 1, nwet_constant, 1
          ct_k298(ns) = henry_defs_const(1,ns)
          ct_dhr(ns) = henry_defs_const(2,ns)
          ct_kd298(ns,1) = henry_defs_const(3,ns)
          ct_ddhr(ns,1) = henry_defs_const(4,ns)
          ct_kd298(ns,2) = henry_defs_const(5,ns)
          ct_ddhr(ns,2) = henry_defs_const(6,ns)
        enddo
        !NO!$acc end kernels
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_inwdep

end module ukca_inwdep_mod
