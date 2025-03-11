module asad_findreaction_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_FINDREACTION_MOD'
  public

  contains
  function asad_findreaction(react1, react2, products, np, reactions, idx_in, nr, ns)
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr
    integer, dimension(nr), intent(in) :: idx_in
    integer, intent(in) :: ns
    integer, intent(in) :: np
    CHARACTER(LEN = 10), INTENT(IN) :: react1
    CHARACTER(LEN = 10), INTENT(IN) :: react2
    CHARACTER(LEN = 10), INTENT(IN) :: products(np)
    CHARACTER(LEN = 10), INTENT(IN) :: reactions(nr, ns)
    integer :: j
    integer :: jp
    integer :: jp2
    integer :: jr
    integer :: nreacts
    integer :: nprods
    integer :: iprods
    integer :: ipargs
    CHARACTER(LEN = 10) :: prods(ns - 2)
    CHARACTER(LEN = 10) :: inprods(ns - 2)
    CHARACTER(LEN = 10) :: empty_ten
    logical, dimension(ns - 2) :: found
    logical :: thelast
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_FINDREACTION'
    integer :: asad_findreaction
    integer :: idx

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    asad_findreaction = 0
    empty_ten(:) = '         '
    nreacts = SIZE(reactions, 1)
    nprods = SIZE(reactions, 2) - 2
    ipargs = 0
    inprods(:) = empty_ten(:)
    do j = 1, np, 1
      if (products(j) /= empty_ten(:)) then
        ipargs = ipargs + 1
        inprods(ipargs) = products(j)
      end if
    enddo
    j = 1
    do while (.true.)
      jr = idx_in(j)
      if (reactions(j,1) == react1(:) .AND. reactions(j,2) == react2(:) .OR. reactions(j,1) == react2(:) .AND. reactions(j,2) == react1(:)) then
        prods(:) = empty_ten(:)
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx = 1, ns - 2, 1
          found(idx) = .false.
        enddo
        !NO!$acc end kernels
        iprods = 0
        do jp = 1, nprods, 1
          if (reactions(j,2 + jp) /= empty_ten(:)) then
            iprods = iprods + 1
            prods(iprods) = reactions(j,2 + jp)
          end if
        enddo
        if (iprods /= ipargs) then
          j = j + 1
          if (j > nreacts) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Cycle_Stmt
          CYCLE
        end if
        do jp = 1, iprods, 1
          do jp2 = 1, iprods, 1
            if (.NOT.found(jp) .AND. inprods(jp) == prods(jp2)) then
              found(jp) = .true.
              prods(jp2) = empty_ten(:)
            end if
          enddo
        enddo
        thelast = .true.
        !NO!$acc kernels
        !NO!$acc loop seq
        do jp = 1, iprods, 1
          thelast = thelast .AND. found(jp)
        enddo
        !NO!$acc end kernels
        if (thelast) then
          asad_findreaction = jr
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Goto_Stmt
          GO TO 9999
        end if
      end if
      j = j + 1
      if (j > nreacts) then
        if (lhook) then
          call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
        end if
        call nvtxEndRange()
        return
      end if
    end do
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - Unsupported labelled statement
    9999 CONTINUE
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end function asad_findreaction

end module asad_findreaction_mod
