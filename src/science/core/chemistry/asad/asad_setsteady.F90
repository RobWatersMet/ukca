module asad_setsteady_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_SETSTEADY_MOD'
  public

  contains
  subroutine asad_setsteady()
    use asad_mod, only : ctype, h_in_ss, jpmsp, jpnr, jpspec, jpss, jpssr, n_in_ss, nspi, nssi, nsspi, nsspt, nssri, nssrt, nssrx, nsst, o1d_in_ss, o3p_in_ss, speci
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use ukca_um_legacy_mod, only : mype
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: js
    integer :: jr
    integer :: jp
    integer :: ix
    integer :: i
    integer :: j
    integer :: k
    integer :: issr
    integer, dimension(jpss) :: indss
    integer, dimension(jpssr) :: issx
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_SETSTEADY'
    integer :: idx

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    nsst = 0
    nssi = 0
    !$acc loop gang vector independent
    do idx = 1, jpss, 1
      indss(idx) = 0
    enddo
    nssrt = 0
    nssrx = 0
    nsspt = 0
    nssri = 0
    nsspi = 0
    ix = 0
    !$acc loop seq
    do js = 1, jpspec, 1
      if (speci(js) == 'O(1D)     ' .AND. o1d_in_ss) then
        indss(1) = js
        ix = ix + 1
      else
        if (speci(js) == 'O(3P)     ' .AND. o3p_in_ss) then
          indss(2) = js
          ix = ix + 1
        else
          if (speci(js) == 'N         ' .AND. n_in_ss) then
            indss(3) = js
            ix = ix + 1
          else
            if (speci(js) == 'H         ' .AND. h_in_ss) then
              indss(4) = js
              ix = ix + 1
            else
              if (ctype(js) == 'SS') then
                ix = ix + 1
                indss(ix) = js
              end if
            end if
          end if
        end if
      end if
    enddo
    !$acc loop seq
    do i = 1, jpss, 1
      if (indss(i) /= 0) then
        if (ctype(indss(i)) == 'SS') then
          nsst = nsst + 1
          nssi(nsst) = indss(i)
        end if
      end if
    enddo
    !$acc end kernels
    if (nsst == 0) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Goto_Stmt
      GO TO 9999
    end if
    if (nsst > jpss) then
      errcode = 1
      cmessage(:) = 'Too many steady state species - increase jpss in asad_mod'
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I0)') 'INCREASE JPSS TO: ', nsst
      call umprint(ummessage, src='asad_setsteady')
      call ereport('ASAD_SETSTEADY', errcode, cmessage)
    end if
    !$acc kernels
    !$acc loop seq
    do jr = 1, jpnr, 1
      !$acc loop seq
      do jp = 1, jpmsp, 1
        !$acc loop gang vector independent
        do ix = 1, nsst, 1
          if (nspi(jr,jp) == nssi(ix)) then
            if (jp <= 2) then
              nssrt(ix) = nssrt(ix) + 1
              nssri(ix,nssrt(ix)) = jr
              nssrx(ix,nssrt(ix)) = 3 - jp
            else
              nsspt(ix) = nsspt(ix) + 1
              nsspi(ix,nsspt(ix)) = jr
            end if
          end if
        enddo
      enddo
    enddo
    !$acc end kernels
    do ix = 1, nsst, 1
      i = MAX(nssrt(ix), nsspt(ix))
      if (i > jpssr) then
        errcode = 2
        cmessage(:) = ' Too many reactions for steady state ' // speci(nssi(ix)) // '- increase jpssr in asad_mod'
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I0)') 'INCREASE JPSSR TO: ', i
        call umprint(ummessage, src='asad_setsteady')
        call ereport('ASAD_SETSTEADY', errcode, cmessage)
      end if
    enddo
    do ix = 1, nsst, 1
      !$acc kernels
      issr = 0
      !$acc loop gang vector independent
      do i = 1, jpssr, 1
        issx(i) = 0
      enddo
      !$acc end kernels
      do i = 1, nssrt(ix), 1
        if (nspi(nssri(ix,i),nssrx(ix,i)) == nssi(ix)) then
          if (issr == 0) then
            cmessage(:) = ' Steady state ' // speci(nssi(ix)) // ' has self reaction - disallowed'
            errcode = -1
            call ereport('ASAD_SETSTEADY', errcode, cmessage)
          end if
          issr = issr + 1
          issx(issr) = i
        end if
      enddo
      if (issr /= 0) then
        !$acc kernels
        nssrt(ix) = nssrt(ix) - issr
        !$acc loop seq
        do i = 1, nssrt(ix), 1
          k = i
          !$acc loop seq
          do j = 1, issr, 1
            if (k >= issx(j)) then
              k = k + 1
            end if
          enddo
          nssri(ix,i) = nssri(ix,k)
          nssrx(ix,i) = nssrx(ix,k)
        enddo
        !$acc end kernels
      end if
    enddo
    do ix = 1, nsst, 1
      if (mype == 0 .AND. printstatus >= prstatus_oper) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(3A)') ' Putting ', speci(nssi(ix)), 'in steady state:'
        call umprint(ummessage, src='asad_setsteady')
      end if
    enddo
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - Unsupported labelled statement
    9999 CONTINUE
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_setsteady

end module asad_setsteady_mod
