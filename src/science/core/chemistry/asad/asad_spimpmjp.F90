module asad_spimpmjp_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_SPIMPMJP_MOD'
  public

  contains
  subroutine forward_euler(n_points, f, f_initial, f_min, nonzero_map, spfj)
    use asad_mod, only : cdt, fdot, jpcspf, spfjsize_max
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, dimension(n_points,jpcspf), intent(out) :: f
    real, dimension(n_points,jpcspf), intent(in) :: f_initial
    real, intent(in) :: f_min
    integer, dimension(jpcspf,jpcspf), intent(in) :: nonzero_map
    real, dimension(n_points,spfjsize_max), intent(in) :: spfj
    integer :: jl
    integer :: jtr
    integer :: ip
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'FORWARD_EULER'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop seq
    do jtr = 1, jpcspf, 1
      ip = nonzero_map(jtr,jtr)
      !$acc loop gang vector independent
      do jl = 1, n_points, 1
        if (spfj(jl,ip) > 0.0) then
          f(jl,jtr) = f_initial(jl,jtr) + cdt * fdot(jl,jtr)
        else
          f(jl,jtr) = f_initial(jl,jtr) + (cdt * fdot(jl,jtr)) / (1.0 - cdt * spfj(jl,ip))
        end if
        if (f(jl,jtr) < f_min) then
          f(jl,jtr) = f_min
        end if
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine forward_euler
  subroutine calc_residual_error(n_points, residual_error, g_f, f_min)
    use asad_mod, only : jpcspf, nlf, prod, slos
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, intent(out) :: residual_error
    real, dimension(n_points,jpcspf), intent(in) :: g_f
    real, intent(in) :: f_min
    real, dimension(n_points,jpcspf) :: tmprc
    integer :: jl
    integer :: jtr
    integer :: j
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'CALC_RESIDUAL_ERROR'
    integer :: idx

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop seq
    do jtr = 1, jpcspf, 1
      j = nlf(jtr)
      !NO!$acc loop gang vector independent
      do idx = 1, n_points, 1
        tmprc(idx,jtr) = prod(idx,j) + slos(idx,j)
      enddo
    enddo
    residual_error = 0.0
    !NO!$acc loop seq
    do jl = 1, n_points, 1
      !NO!$acc loop seq
      do jtr = 1, jpcspf, 1
        if (ABS(tmprc(jl,jtr)) > f_min) then
          residual_error = MAX(residual_error, ABS(g_f(jl,jtr) / tmprc(jl,jtr)))
        end if
      enddo
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine calc_residual_error
  subroutine calc_error_norm(n_points, error_norm, f, f_incr, f_min)
    use asad_mod, only : jpcspf
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, intent(out) :: error_norm
    real, dimension(n_points,jpcspf), intent(in) :: f
    real, dimension(n_points,jpcspf), intent(in) :: f_incr
    real, intent(in) :: f_min
    integer :: jl
    integer :: jtr
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'CALC_ERROR_NORM'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    error_norm = 0.0
    !$acc loop seq
    do jtr = 1, jpcspf, 1
      !$acc loop seq
      do jl = 1, n_points, 1
        if (ABS(f_incr(jl,jtr)) > 1.0e-16) then
          error_norm = MAX(error_norm, ABS(f_incr(jl,jtr) / MAX(f(jl,jtr), f_min)))
        end if
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine calc_error_norm
  subroutine asad_spimpmjp(exit_code, ix, jy, nlev, n_points, location, solver_iter)
    use asad_mod, only : cdt, f, fdot, fj, jpcspf, ltrig, modified_map, nitnr, nonzero_map, nonzero_map_unordered, nstst, peps, ptol, spfj, y
    use asad_sparse_vars, only : setup_spfuljac, spfuljac, splinslv2, spresolv2
    use ukca_config_specification_mod, only : ukca_config
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    use umPrintMgr, only : printstatus, prstatus_diag, prstatus_oper, ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    use ukca_um_legacy_mod, only : mype
    use asad_diffun_mod, only : asad_diffun
    use asad_fuljac_mod, only : asad_fuljac
    use asad_steady_mod, only : asad_steady
    use asad_ftoy_mod, only : asad_ftoy
    integer, parameter :: maxneg = 2000
    real, parameter :: f_max = 1.0e30
    integer, parameter :: ltrig_iter = 51
    integer, parameter :: incr_limiter = 7
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    integer, intent(in) :: nlev
    integer, intent(in) :: location
    integer, intent(out) :: exit_code
    integer, intent(out) :: solver_iter
    integer :: iter
    integer :: jl
    integer :: jtr
    integer :: ifi
    integer :: i
    integer :: itr
    integer :: count_negatives
    real :: ztmp
    real :: f_min
    real :: reltol_residual_error
    real :: reltol_error
    real :: rafmin
    real :: rafmax
    real :: deltt
    real :: damp1
    real :: error_norm
    real :: residual_error
    logical, save :: not_first_call = .false.
    real, dimension(n_points,jpcspf) :: f_initial
    real, dimension(n_points,jpcspf) :: f_incr
    real, dimension(n_points,jpcspf) :: g_f
    real, dimension(n_points,jpcspf) :: g_f_old
    real, dimension(n_points,jpcspf) :: g_ftmp
    real, dimension(n_points,jpcspf) :: delta_g
    real :: coeff
    real :: denom
    real, dimension(jpcspf) :: zsum
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_SPIMPMJP'
    CHARACTER(LEN = errormessagelength) :: cmessage1 = "(1x,i2,20(1x,1pG12.4))"
    CHARACTER(LEN = errormessagelength) :: cmessage2 = "(1x,a3,20(1x,1pG12.4))"
    logical, save :: first = .true.
    logical, save :: first_pass = .true.
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19
    integer :: idx_20
    integer :: idx_21

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    f_min = SQRT(peps)
    reltol_residual_error = 1.0e-10
    reltol_error = 10 * ptol
    rafmin = 1.0e-01
    rafmax = 1.0e+04
    exit_code = 0
    error_norm = 1.0
    deltt = 1.0 / cdt
    damp1 = 0.5
    !NO!$acc loop gang vector independent collapse(2)
    do idx = 1, jpcspf, 1
      do idx_1 = 1, n_points, 1
        f_initial(idx_1,idx) = f(idx_1 + (LBOUND(f, dim=1) - 1),idx + (LBOUND(f, dim=2) - 1))
      enddo
    enddo
    !NO!$acc end kernels
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (f(:, :) < f_min) f = f_min
    if (nstst /= 0) then
      call asad_steady(n_points)
    end if
    print *, "RW - part 1 - min value of nonzero_map = ", minval(nonzero_map)
    if (first_pass) then
      if (first) then
        call setup_spfuljac()
        first = .false.
        print *, "RW - part 2 - min value of nonzero_map = ", minval(nonzero_map)
      end if
      first_pass = .false.
    end if
    print *, "RW - part 3 - min value of nonzero_map = ", minval(nonzero_map)
    call spfuljac(n_points, cdt, f_min, nonzero_map, spfj)
    call forward_euler(n_points, f, f_initial, f_min, nonzero_map, spfj)
    do iter = 1, ukca_config%nrsteps, 1
      ifi = 0
      if (iter == 1) then
        ifi = nitnr
      end if
      call asad_ftoy(not_first_call, ifi, iter, n_points, ix, jy, nlev)
      if (nstst /= 0 .AND. ifi == 0) then
        call asad_steady(n_points)
      end if
      if (ltrig .AND. printstatus >= prstatus_oper) then
        do jl = 1, n_points, 1
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, "('Point: ',i4)") jl
          call umprint(ummessage, src='asad_spimpmjp')
          if (iter == 1) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, cmessage1) iter - 1, (f_initial(jl, jtr), jtr = 1, jpcspf)
            call umprint(ummessage, src='asad_spimpmjp')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, cmessage1) iter - 1, ((f_initial(jl, jtr) + cdt * fdot(jl, jtr)), jtr = 1, jpcspf)
            call umprint(ummessage, src='asad_spimpmjp')
          end if
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage1) iter - 1, (f(jl, jtr), jtr = 1, jpcspf), (y(jl, i), i = 1, 2)
          call umprint(ummessage, src='asad_spimpmjp')
        enddo
      end if
      call asad_diffun(n_points)
      if (error_norm < reltol_error) then
        exit_code = 0
        solver_iter = iter - 1
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Goto_Stmt
        GO TO 9999
      end if
      if (iter == ukca_config%nrsteps) then
        exit_code = 2
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Goto_Stmt
        GO TO 9999
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx_2 = 1, jpcspf, 1
        do idx_3 = 1, n_points, 1
          g_f(idx_3,idx_2) = (f(idx_3 + (LBOUND(f, dim=1) - 1),idx_2 + (LBOUND(f, dim=2) - 1)) - f_initial(idx_3,idx_2)) * deltt - fdot(idx_3 + (LBOUND(fdot, dim=1) - 1),idx_2 + (LBOUND(fdot, dim=2) - 1))
        enddo
      enddo
      !NO!$acc end kernels
      call calc_residual_error(n_points, residual_error, g_f, f_min)
      if (residual_error < reltol_residual_error) then
        exit_code = 0
        solver_iter = iter
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Goto_Stmt
        GO TO 9999
      end if
      do jl = 1, n_points, 1
        do jtr = 1, jpcspf, 1
          if (f(jl,jtr) > f_max) then
            exit_code = 3
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
        enddo
      enddo
      call spfuljac(n_points, cdt, f_min, nonzero_map, spfj)
      if (ltrig .AND. printstatus == prstatus_diag) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, "('Iteration ',i4)") iter
        call umprint(ummessage, src='asad_spimpmjp')
        do jl = 1, n_points, 1
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, "('Point: ',i4)") jl
          call umprint(ummessage, src='asad_spimpmjp')
          !NO!$acc kernels
          !NO!$acc loop seq
          do idx_4 = LBOUND(fj, dim=3), UBOUND(fj, dim=3), 1
            !NO!$acc loop seq
            do idx_5 = LBOUND(fj, dim=2), UBOUND(fj, dim=2), 1
              fj(jl,idx_5,idx_4) = 0.0
            enddo
          enddo
          !NO!$acc loop seq
          do jtr = 1, jpcspf, 1
            !NO!$acc loop seq
            do itr = 1, jpcspf, 1
              if (nonzero_map(jtr,itr) > 0) then
                fj(jl,jtr,itr) = spfj(jl,nonzero_map(jtr,itr))
              end if
            enddo
          enddo
          !NO!$acc end kernels
          do jtr = 1, jpcspf, 1
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, cmessage1) jtr, (fj(jl, jtr, itr), itr = 1, jpcspf)
            call umprint(ummessage, src='asad_spimpmjp')
          enddo
        enddo
      end if
      call splinslv2(n_points, g_f, f_incr, f_min, f_max, nonzero_map_unordered, modified_map, spfj)
      if (ltrig .AND. printstatus == prstatus_diag) then
        do jl = 1, n_points, 1
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, "('Point: ',i4)") jl
          call umprint(ummessage, src='asad_spimpmjp')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage2) 'G_f', (G_f(jl, jtr), jtr = 1, jpcspf)
          call umprint(ummessage, src='asad_spimpmjp')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage2) 'fdt', (fdot(jl, jtr), jtr = 1, jpcspf)
          call umprint(ummessage, src='asad_spimpmjp')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage2) 'del', ((f(jl, jtr) - f_initial(jl, jtr)) * deltt, jtr = 1, jpcspf)
          call umprint(ummessage, src='asad_spimpmjp')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage2) 'f  ', (f(jl, jtr), jtr = 1, jpcspf)
          call umprint(ummessage, src='asad_spimpmjp')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage2) 'f_incr', (f_incr(jl, jtr), jtr = 1, jpcspf)
          call umprint(ummessage, src='asad_spimpmjp')
        enddo
        call asad_fuljac(n_points)
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, "('Iteration ',i4)") iter
        call umprint(ummessage, src='asad_spimpmjp')
        do jl = 1, n_points, 1
          do jtr = 1, jpcspf, 1
            !NO!$acc kernels
            zsum(jtr) = 0.0
            !NO!$acc loop seq
            do itr = 1, jpcspf, 1
              zsum(jtr) = zsum(jtr) + fj(jl,jtr,itr) * f_incr(jl,itr)
            enddo
            !NO!$acc end kernels
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, cmessage1) jtr, (fj(jl, jtr, itr) * f_incr(jl, itr), itr = 1, jpcspf)
            call umprint(ummessage, src='asad_spimpmjp')
          enddo
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, cmessage2) 'sum', (zsum(jtr), jtr = 1, jpcspf)
          call umprint(ummessage, src='asad_spimpmjp')
        enddo
      end if
      if (iter == 1) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(2)
        do idx_6 = 1, jpcspf, 1
          do idx_7 = 1, n_points, 1
            f_incr(idx_7,idx_6) = damp1 * f_incr(idx_7,idx_6)
          enddo
        enddo
        !NO!$acc end kernels
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx_8 = 1, jpcspf, 1
        do idx_9 = 1, n_points, 1
          f_incr(idx_9,idx_8) = MIN(MAX(f_incr(idx_9,idx_8), -f_max), f_max)
        enddo
      enddo
      !NO!$acc end kernels
      call calc_error_norm(n_points, error_norm, f, f_incr, f_min)
      count_negatives = 0
      do jtr = 1, jpcspf, 1
        do jl = 1, n_points, 1
          ztmp = f(jl,jtr) + f_incr(jl,jtr)
          if (iter < incr_limiter) then
            ztmp = MAX(rafmin * f(jl,jtr), MIN(rafmax * f(jl,jtr), ztmp))
          end if
          if (ztmp == 0.0) then
            ztmp = f_min
          end if
          if (ztmp < 0.0) then
            ztmp = f_min
            count_negatives = count_negatives + 1
          end if
          f(jl,jtr) = ztmp
        enddo
      enddo
      if (count_negatives > maxneg) then
        exit_code = 2
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Goto_Stmt
        GO TO 9999
      end if
      if (ukca_config%l_ukca_quasinewton .AND. error_norm >= reltol_error) then
        if (iter >= ukca_config%i_ukca_quasinewton_start .AND. iter <= ukca_config%i_ukca_quasinewton_end) then
          call asad_ftoy(not_first_call, ifi, iter, n_points, ix, jy, nlev)
          call asad_steady(n_points)
          call asad_diffun(n_points)
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do idx_10 = 1, jpcspf, 1
            do idx_11 = 1, n_points, 1
              g_f_old(idx_11,idx_10) = g_f(idx_11,idx_10)
            enddo
          enddo
          !NO!$acc loop gang vector independent collapse(2)
          do idx_12 = 1, jpcspf, 1
            do idx_13 = 1, n_points, 1
              g_f(idx_13,idx_12) = (f(idx_13 + (LBOUND(f, dim=1) - 1),idx_12 + (LBOUND(f, dim=2) - 1)) - f_initial(idx_13,idx_12)) * deltt - fdot(idx_13 + (LBOUND(fdot, dim=1) - 1),idx_12 + (LBOUND(fdot, dim=2) - 1))
            enddo
          enddo
          !NO!$acc loop gang vector independent collapse(2)
          do idx_14 = 1, jpcspf, 1
            do idx_15 = 1, n_points, 1
              delta_g(idx_15,idx_14) = g_f(idx_15,idx_14) - g_f_old(idx_15,idx_14)
            enddo
          enddo
          coeff = 0.0
          denom = 0.0
          !NO!$acc loop seq
          do jl = 1, n_points, 1
            !NO!$acc loop seq
            do jtr = 1, jpcspf, 1
              coeff = coeff + g_f(jl,jtr) * delta_g(jl,jtr)
              denom = denom + delta_g(jl,jtr) * delta_g(jl,jtr)
            enddo
          enddo
          coeff = 1.0 - coeff / denom
          !NO!$acc loop gang vector independent collapse(2)
          do jl = 1, n_points, 1
            do jtr = 1, jpcspf, 1
              g_ftmp(jl,jtr) = g_f(jl,jtr) * coeff
            enddo
          enddo
          !NO!$acc end kernels
          call spresolv2(n_points, g_ftmp, f_incr, f_min, modified_map, spfj)
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(2)
          do idx_16 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
            do idx_17 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
              f(idx_17,idx_16) = f(idx_17,idx_16) + f_incr(idx_17 + (1 - LBOUND(f, dim=1)),idx_16 + (1 - LBOUND(f, dim=2)))
            enddo
          enddo
          !NO!$acc loop gang vector independent collapse(2)
          do idx_18 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
            do idx_19 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
              f(idx_19,idx_18) = ABS(f(idx_19,idx_18))
            enddo
          enddo
          !NO!$acc end kernels
        end if
      end if
    enddo
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - Unsupported labelled statement
    9999 CONTINUE
    if (exit_code /= 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx_20 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
        do idx_21 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
          f(idx_21,idx_20) = f_initial(idx_21 + (1 - LBOUND(f, dim=1)),idx_20 + (1 - LBOUND(f, dim=2)))
        enddo
      enddo
      !NO!$acc end kernels
      solver_iter = iter
      if (count_negatives > maxneg) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, "('Negatives - exceeds maxneg')")
        call umprint(ummessage, src='asad_spimpmjp')
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, "(1x,'Too many negatives (>',i4,') in spimpmjp (iter',i3," // "')  lon=',i3,'  lat=',i3,'; halving step')") maxneg, iter, location, mype
        call umprint(ummessage, src='asad_spimpmjp')
      end if
      if (iter >= ltrig_iter) then
        exit_code = 4
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, "('Convergence problems (',i3,1x,'iter) at location=',i3,' pe=',i3)") iter, location, mype
        call umprint(ummessage, src='asad_spimpmjp')
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_spimpmjp

end module asad_spimpmjp_mod
