module asad_spmjpdriv_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_SPMJPDRIV_MOD'
  public

  contains
  subroutine asad_spmjpdriv(ix, jy, nlev, n_points, num_iter)
    use asad_mod, only : cdt, f, jpcspf, jpspec, ltrig, ncsteps, nitfg, speci, y
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use ukca_um_legacy_mod, only : mype
    use errormessagelength_mod, only : errormessagelength
    use asad_diffun_mod, only : asad_diffun
    use asad_spimpmjp_mod, only : asad_spimpmjp
    use asad_ftoy_mod, only : asad_ftoy
    integer, parameter :: max_redo = 128
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    integer, intent(in) :: nlev
    integer, intent(out) :: num_iter
    integer :: exit_code
    integer :: ncst
    integer :: iredo
    integer :: jl
    integer :: iter
    integer :: location
    integer :: solver_iter
    integer :: jit
    integer :: errcode
    logical, save :: not_first_call = .false.
    CHARACTER(LEN = errormessagelength) :: cmessage
    real :: ctrd
    real, dimension(n_points,jpcspf) :: f_initial
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_SPMJPDRIV'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    ncst = ncsteps
    ctrd = cdt
    ltrig = .false.
    call asad_diffun(n_points)
    !NO!$acc kernels
    iredo = 1
    !NO!$acc loop gang vector independent collapse(2)
    do idx = LBOUND(f_initial, dim=2), UBOUND(f_initial, dim=2), 1
      do idx_1 = 1, n_points, 1
        f_initial(idx_1,idx) = f(idx_1,idx + (LBOUND(f, dim=2) - LBOUND(f_initial, dim=2)))
      enddo
    enddo
    !NO!$acc end kernels
    if (ukca_config%l_ukca_asad_columns) then
      location = ix + (jy - 1) * ukca_config%row_length
    else
      location = nlev
    end if
    num_iter = 0
    iter = 1
    do while (iter <= iredo)
      call asad_spimpmjp(exit_code, ix, jy, nlev, n_points, location, solver_iter)
      num_iter = num_iter + solver_iter
      if (exit_code == 4) then
        if (ltrig) then
          errcode = 1
          cmessage(:) = 'Slow-converging system, Set printstatus for Jacobian debug'
          do jl = 1, jpspec, 1
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(a4,i6,a12,2e14.5,i12)') 'y: ', jl, speci(jl), MAXVAL(y(:, jl)), MINVAL(y(:, jl)), SIZE(y(:, jl))
            call umprint(ummessage, src='asad_spmjpdriv')
          enddo
          call ereport('ASAD_SPMJPDRIV', errcode, cmessage)
        end if
        !NO!$acc kernels
        ltrig = .true.
        !NO!$acc loop gang vector independent collapse(2)
        do idx_2 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
          do idx_3 = 1, n_points, 1
            f(idx_3,idx_2) = f_initial(idx_3,idx_2 + (LBOUND(f_initial, dim=2) - LBOUND(f, dim=2)))
          enddo
        enddo
        !NO!$acc end kernels
        jit = 0
        call asad_ftoy(not_first_call, nitfg, jit, n_points, ix, jy, nlev)
        call asad_diffun(n_points)
        iter = 1
      else
        if (exit_code > 1) then
          ncsteps = ncsteps * 2
          cdt = cdt / 2.0
          iredo = iredo * 2
          if (ukca_config%l_ukca_debug_asad) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, "('ASAD: failed to converge at location  = ',I0)") location
            call umprint(ummessage, src='asad_spmjpdriv')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A,I0,A,I0,A,E18.8)') 'ASAD: halving timestep: ncsteps = ', ncsteps, ' iredo = ', iredo, ' cdt = ', cdt
            call umprint(ummessage, src='asad_spmjpdriv')
          end if
          if (cdt < 1.0e-05) then
            errcode = 2
            cmessage(:) = ' Time step now too short'
            call ereport('ASAD_SPMJPDRIV', errcode, cmessage)
          end if
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do idx_4 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
            do idx_5 = 1, n_points, 1
              f(idx_5,idx_4) = f_initial(idx_5,idx_4 + (LBOUND(f_initial, dim=2) - LBOUND(f, dim=2)))
            enddo
          enddo
          !NO!$acc end kernels
          if (iredo >= max_redo) then
            if (printstatus >= prstatus_oper) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, "(' Resetting array after',i4,' iterations')") iredo
              call umprint(ummessage, src='asad_spmjpdriv')
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, "('NO CONVERGENCE location: ',i4,' pe: ',i4)") location, mype
              call umprint(ummessage, src='asad_spmjpdriv')
            end if
            ncsteps = ncst
            cdt = ctrd
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          jit = 0
          call asad_ftoy(not_first_call, nitfg, jit, n_points, ix, jy, nlev)
          call asad_diffun(n_points)
          iter = 1
        else
          iter = iter + 1
        end if
      end if
    end do
    if (iredo > 1) then
      if (iredo > 2) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, "('   No. iterations =',i2)") iredo
        call umprint(ummessage, src='asad_spmjpdriv')
      end if
      ncsteps = ncst
      cdt = ctrd
    end if
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - Unsupported labelled statement
    9999 CONTINUE
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_spmjpdriv

end module asad_spmjpdriv_mod
