module asad_impact_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_IMPACT_MOD'
  public

  contains
  subroutine asad_impact(n_points, ix, jy, nlev)
    use asad_mod, only : cdt, ctype, ej, f, fdot, jpcspf, jpif, linfam, ljacx, madvtr, majors, moffam, nitnr, nlpdv, nltrim, nprkx, nspi, peps, pmintnd, ptol, rk, y
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use asad_diffun_mod, only : asad_diffun
    use asad_inimpct_mod, only : asad_inimpct
    use asad_jac_mod, only : asad_jac
    use asad_ftoy_mod, only : asad_ftoy
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    integer, intent(in) :: nlev
    integer, dimension(n_points) :: lphot
    integer :: inl
    integer :: nl
    integer :: j
    integer :: jit
    integer :: jl
    integer :: jr
    integer :: jtr
    integer :: jt1
    integer :: isp
    integer :: ifi
    integer :: ntr1
    integer :: ntr2
    integer :: ipos0
    integer :: ipos
    integer :: ir
    integer :: irk
    integer :: njr
    integer :: ireac
    integer :: iprod
    real, dimension(n_points,jpcspf) :: zf
    real, dimension(n_points,jpcspf) :: zprf
    real, dimension(n_points,jpcspf) :: binv
    real, dimension(n_points,jpcspf) :: dely
    real, dimension(n_points) :: dd
    real, dimension(n_points) :: corrn
    logical :: gconv
    logical :: gfam
    logical, save :: gfirst = .true.
    logical, save :: first_pass = .true.
    logical, save :: not_first_call = .false.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_IMPACT'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    inl = 0
    !$acc loop seq
    do jl = 1, n_points, 1
      if (rk(jl,nprkx(1)) > peps) then
        inl = inl + 1
        lphot(inl) = jl
      end if
    enddo
    !$acc end kernels
    if (first_pass) then
      if (gfirst) then
        gfirst = .false.
        call asad_inimpct()
      end if
      first_pass = .false.
    end if
    !$acc kernels
    nl = n_points
    !$acc loop seq
    do jtr = 1, jpcspf, 1
      isp = majors(jtr)
      !$acc loop gang vector independent
      do jl = 1, n_points, 1
        zf(jl,jtr) = f(jl,jtr)
        if (ej(jl,jtr) > 0.0) then
          f(jl,jtr) = zf(jl,jtr) + cdt * fdot(jl,jtr)
        else
          f(jl,jtr) = zf(jl,jtr) + cdt * fdot(jl,jtr) / (1.0 - cdt * ej(jl,jtr))
        end if
        if (linfam(jl,jtr)) then
          f(jl,jtr) = y(jl,isp)
        end if
      enddo
    enddo
    !$acc end kernels
    do jit = 1, ukca_config%nrsteps, 1
      if (jit < 4) then
        ifi = 0
      else
        ifi = nitnr
      end if
      call asad_ftoy(not_first_call, ifi, jit, n_points, ix, jy, nlev)
      call asad_diffun(nl)
      call asad_jac(n_points)
      !$acc kernels
      !$acc loop seq
      do jtr = 1, jpcspf, 1
        isp = majors(jtr)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          zprf(jl,jtr) = f(jl,jtr)
          binv(jl,jtr) = 1.0 / (1.0 - cdt * ej(jl,jtr))
          dely(jl,jtr) = binv(jl,jtr) * (zf(jl,jtr) - f(jl,jtr) + cdt * fdot(jl,jtr))
          f(jl,jtr) = f(jl,jtr) + dely(jl,jtr)
          if (linfam(jl,jtr)) then
            f(jl,jtr) = y(jl,isp)
          end if
        enddo
      enddo
      !$acc end kernels
      if (ljacx .AND. inl > 0) then
        ntr1 = nltrim(0,1)
        do jt1 = 1, ntr1, 1
          !$acc kernels
          ntr2 = nltrim(jt1,2)
          ipos0 = nltrim(jt1,3)
          ipos = ipos0
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            corrn(jl) = 0.0
          enddo
          !$acc end kernels
          do while (.true.)
            !$acc kernels
            ir = nlpdv(ipos,1)
            irk = nprkx(ir)
            isp = nspi(irk,1)
            njr = nlpdv(ipos,2)
            ireac = madvtr(isp)
            gfam = ireac == 0
            !$acc loop gang vector independent
            do j = 1, inl, 1
              dd(j) = 0.0
            enddo
            !$acc end kernels
            if (gfam) then
              !$acc kernels
              ireac = moffam(isp)
              !$acc loop seq
              do jr = 1, njr, 1
                ir = nlpdv(ipos,1)
                irk = nprkx(ir)
                isp = nspi(irk,1)
                !$acc loop gang vector independent
                do j = 1, inl, 1
                  jl = lphot(j)
                  dd(j) = dd(j) + rk(jl,irk) * y(jl,isp)
                enddo
                ipos = ipos + 1
              enddo
              !$acc loop gang vector independent
              do j = 1, inl, 1
                jl = lphot(j)
                dd(j) = dd(j) / zprf(jl,ireac)
              enddo
              !$acc end kernels
            else
              !$acc kernels
              !$acc loop seq
              do jr = 1, njr, 1
                ir = nlpdv(ipos,1)
                irk = nprkx(ir)
                isp = nspi(irk,1)
                !$acc loop gang vector independent
                do j = 1, inl, 1
                  jl = lphot(j)
                  if (linfam(jl,ireac)) then
                    dd(j) = dd(j) + rk(jl,irk) * y(jl,isp) / zprf(jl,ireac)
                  else
                    dd(j) = dd(j) + rk(jl,irk)
                  end if
                enddo
                ipos = ipos + 1
              enddo
              !$acc end kernels
            end if
            !$acc kernels
            !$acc loop gang vector independent
            do j = 1, inl, 1
              jl = lphot(j)
              corrn(j) = corrn(j) + dd(j) * dely(jl,ireac)
            enddo
            !$acc end kernels
            if (ipos - ipos0 >= ntr2) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Exit_Stmt
              EXIT
            end if
          end do
          iprod = nltrim(jt1,1)
          isp = majors(iprod)
          if (ctype(isp) == jpif) then
            !$acc kernels
            !$acc loop seq
            do j = 1, inl, 1
              jl = lphot(j)
              iprod = nltrim(jt1,1)
              if (linfam(jl,iprod)) then
                iprod = moffam(isp)
              end if
              f(jl,iprod) = f(jl,iprod) + cdt * binv(jl,iprod) * corrn(j)
            enddo
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop seq
            do j = 1, inl, 1
              jl = lphot(j)
              f(jl,iprod) = f(jl,iprod) + cdt * binv(jl,iprod) * corrn(j)
            enddo
            !$acc end kernels
          end if
        enddo
      end if
      gconv = .true.
      do jtr = 1, jpcspf, 1
        !$acc kernels
        !$acc loop seq
        do jl = 1, n_points, 1
          if (ABS(f(jl,jtr) - zprf(jl,jtr)) > ptol * f(jl,jtr) .AND. f(jl,jtr) > pmintnd(jl)) then
            gconv = .false.
          end if
        enddo
        !$acc end kernels
        if (.NOT.gconv) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Exit_Stmt
          EXIT
        end if
      enddo
      if (gconv) then
        if (lhook) then
          call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
        end if
        call nvtxEndRange()
        return
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_impact

end module asad_impact_mod
