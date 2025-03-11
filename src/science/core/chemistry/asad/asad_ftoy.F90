module asad_ftoy_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_FTOY_MOD'
  public

  contains
  subroutine asad_ftoy(ofirst, iter, num_iter, n_points, ix, jy, nlev)
    use asad_mod, only : cdt, ctype, f, ftol, ilft, ilstmin, jpfm, jpif, jpna, jpro2, jpspec, jsro2, linfam, madvtr, majors, moffam, nlfro2, nlmajmin, nlstst, nodd, nstst, peps, pmintnd, prod, qa, ratio, slos, y
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : umprint
    use errormessagelength_mod, only : errormessagelength
    use asad_fyfixr_mod, only : asad_fyfixr
    use asad_fyinit_mod, only : asad_fyinit
    use asad_prls_mod, only : asad_prls
    use asad_fyself_mod, only : asad_fyself
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    integer, intent(in) :: nlev
    integer, intent(inout) :: iter
    integer, intent(in) :: num_iter
    logical, intent(in) :: ofirst
    real, dimension(n_points) :: fro2
    real, dimension(n_points) :: fro2_old
    real :: damp_ro2
    integer :: j
    integer :: jit
    integer :: jl
    integer :: js
    integer :: ifam
    integer :: imaj
    integer :: nl
    integer :: istart
    integer :: iend
    integer :: iodd
    integer :: itr
    integer :: icode
    integer :: iro2
    integer :: errcode
    integer, dimension(nstst) :: errcodes
    integer, save :: istmin
    integer, save :: ift
    CHARACTER(LEN = errormessagelength) :: cmessage
    real :: zthresh
    real :: sl
    real, dimension(n_points,jpspec) :: zy
    real, dimension(n_points) :: zb
    real, dimension(n_points) :: zc
    real, dimension(n_points) :: zd
    logical :: gconv
    logical, save :: gonce = .true.
    logical, save :: first_pass = .true.
    logical, save :: gdepem = .false.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_FTOY'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    errcode = 0
    !NO!$acc loop gang vector independent
    do idx = LBOUND(errcodes, dim=1), UBOUND(errcodes, dim=1), 1
      errcodes(idx) = 0
    enddo
    !NO!$acc end kernels
    if (first_pass) then
      if (gonce) then
        gonce = .false.
        istmin = 0
        ift = 0
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do j = 1, jpspec, 1
          ilstmin(j) = 0
          ilft(j) = 0
        enddo
        !NO!$acc end kernels
        !NO!$acc kernels
        !NO!$acc loop seq
        do j = 1, nstst, 1
          js = nlstst(j)
          if (ctype(js) == jpfm(:)) then
            ifam = moffam(js)
            imaj = majors(ifam)
            if (imaj /= js) then
              istmin = istmin + 1
              ilstmin(istmin) = js
            end if
          else
            if (ctype(js) == jpif(:)) then
              istmin = istmin + 1
              ilstmin(istmin) = js
              ift = ift + 1
              ilft(ift) = js
            else
              if (ctype(js) == jpna(:)) then
                istmin = istmin + 1
                ilstmin(istmin) = js
              else
                errcodes(j) = 1
              end if
            end if
          end if
        enddo
        !NO!$acc end kernels
        errcode = MAXVAL(errcodes)
      end if
      first_pass = .false.
    end if
    if (errcode /= 0) then
      call umprint('**** ASAD ERROR in ASAD_FTOY!! ', src='asad_ftoy')
      call umprint('ASAD_FTOY found an unexpected species type', src='asad_ftoy')
      call umprint('in the species list nlstst ', src='asad_ftoy')
      cmessage(:) = 'Found unexpected species type'
      call ereport('ASAD_FTOY', js, cmessage)
    end if
    call asad_fyinit(ofirst, n_points, ix, jy, nlev)
    if (ukca_config%l_ukca_ro2_perm) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_1 = 1, n_points, 1
        fro2_old(idx_1) = y(idx_1,jsro2)
      enddo
      !NO!$acc loop gang vector independent
      do idx_2 = 1, n_points, 1
        fro2(idx_2) = 0.0
      enddo
      !NO!$acc loop seq
      do j = 1, jpro2, 1
        iro2 = nlfro2(j)
        !NO!$acc loop gang vector independent
        do idx_3 = 1, n_points, 1
          fro2(idx_3) = fro2(idx_3) + f(idx_3,iro2)
        enddo
      enddo
      !NO!$acc end kernels
      if (num_iter <= 0) then
        damp_ro2 = 1.0
      else
        if (num_iter == 1) then
          damp_ro2 = 0.5
        else
          if (num_iter == 2) then
            damp_ro2 = 0.9
          else
            damp_ro2 = 0.99
          end if
        end if
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_4 = 1, n_points, 1
        y(idx_4,jsro2) = fro2_old(idx_4) + damp_ro2 * (fro2(idx_4) - fro2_old(idx_4))
      enddo
      !NO!$acc end kernels
    end if
    if (nstst == 0) then
      if (lhook) then
        call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
      end if
      call nvtxEndRange()
      return
    end if
    !NO!$acc kernels
    nl = n_points
    !NO!$acc loop seq
    do j = 1, nstst, 1
      js = nlstst(j)
      !NO!$acc loop gang vector independent
      do jl = 1, n_points, 1
        zy(jl,js) = 0.0
      enddo
    enddo
    !NO!$acc end kernels
    if (ofirst .AND. iter < 5) then
      iter = 5
      icode = -1
      cmessage(:) = 'iter too low on first call, resetting to 5'
      call ereport('ASAD_FTOY', icode, cmessage)
    end if
    zthresh = 2.0 / cdt
    if (ofirst) then
      call asad_fyself(n_points)
    end if
    if (iter == 0) then
      call asad_fyfixr(n_points)
      if (lhook) then
        call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
      end if
      call nvtxEndRange()
      return
    end if
    do jit = 1, iter, 1
      if (jit > 1) then
        !NO!$acc kernels
        !NO!$acc loop seq
        do j = 1, nstst, 1
          js = nlstst(j)
          !NO!$acc loop gang vector independent
          do jl = 1, n_points, 1
            if (ABS(y(jl,js)) < peps .AND. prod(jl,js) > peps) then
              y(jl,js) = peps
            end if
          enddo
        enddo
        !NO!$acc end kernels
      end if
      call asad_prls(nl, istmin, ilstmin, gdepem)
      !NO!$acc kernels
      !NO!$acc loop seq
      do j = 1, nstst, 1
        js = nlstst(j)
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          ratio(jl,js) = 0.0
        enddo
      enddo
      !NO!$acc loop seq
      do j = 1, istmin, 1
        js = ilstmin(j)
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          if (y(jl,js) < peps) then
            sl = 0.0
          else
            sl = slos(jl,js) / y(jl,js)
          end if
          if (qa(jl,js) > peps) then
            zb(jl) = sl - qa(jl,js) * y(jl,js)
            zc(jl) = prod(jl,js)
            zd(jl) = zb(jl) * zb(jl) + 4.0 * qa(jl,js) * zc(jl)
            if (zd(jl) > 0.0) then
              y(jl,js) = (zb(jl) - SQRT(zd(jl))) / (-2.0 * qa(jl,js))
              if (y(jl,js) < 0.0) then
                y(jl,js) = 10.0 * peps
              end if
            else
              y(jl,js) = zb(jl) / (-2.0 * qa(jl,js))
            end if
          else
            if (qa(jl,js) <= peps .AND. sl > peps) then
              y(jl,js) = prod(jl,js) / sl
            else
              y(jl,js) = 0.0
            end if
          end if
        enddo
      enddo
      istart = nlmajmin(3)
      iend = nlmajmin(4)
      !NO!$acc loop seq
      do j = nlmajmin(3), nlmajmin(4), 1
        js = nlmajmin(j)
        iodd = nodd(js)
        ifam = moffam(js)
        imaj = 0
        if (ifam /= 0) then
          imaj = majors(ifam)
        end if
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          if (y(jl,imaj) > peps) then
            ratio(jl,js) = y(jl,js) / y(jl,imaj)
            ratio(jl,imaj) = ratio(jl,imaj) + iodd * ratio(jl,js)
          else
            ratio(jl,js) = 0.0
          end if
        enddo
      enddo
      !NO!$acc end kernels
      do j = 1, ift, 1
        js = ilft(j)
        iodd = nodd(js)
        ifam = moffam(js)
        itr = madvtr(js)
        imaj = majors(ifam)
        if (ofirst .AND. jit == 1) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do jl = 1, n_points, 1
            if (y(jl,js) > peps) then
              sl = slos(jl,js) / y(jl,js)
            else
              sl = 0.0
            end if
            linfam(jl,itr) = sl > zthresh
            if (linfam(jl,itr)) then
              f(jl,ifam) = f(jl,ifam) + iodd * f(jl,itr)
            end if
          enddo
          !NO!$acc end kernels
        end if
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          if (linfam(jl,itr)) then
            if (y(jl,imaj) > peps) then
              ratio(jl,js) = y(jl,js) / y(jl,imaj)
              ratio(jl,imaj) = ratio(jl,imaj) + iodd * ratio(jl,js)
            else
              ratio(jl,js) = 0.0
            end if
          else
            y(jl,js) = f(jl,itr)
          end if
        enddo
        !NO!$acc end kernels
      enddo
      !NO!$acc kernels
      istart = nlmajmin(1)
      iend = nlmajmin(2)
      !NO!$acc loop seq
      do j = istart, iend, 1
        js = nlmajmin(j)
        iodd = nodd(js)
        ifam = moffam(js)
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          ratio(jl,js) = 1.0 / (iodd + ratio(jl,js))
          y(jl,js) = f(jl,ifam) * ratio(jl,js)
        enddo
      enddo
      istart = nlmajmin(3)
      iend = nlmajmin(4)
      !NO!$acc loop seq
      do j = istart, iend, 1
        js = nlmajmin(j)
        ifam = moffam(js)
        imaj = majors(ifam)
        if (ctype(js) /= jpif(:)) then
          !NO!$acc loop gang vector independent
          do jl = 1, n_points, 1
            y(jl,js) = y(jl,imaj) * ratio(jl,js)
          enddo
        else
          itr = madvtr(js)
          !NO!$acc loop gang vector independent
          do jl = 1, n_points, 1
            if (linfam(jl,itr)) then
              y(jl,js) = y(jl,imaj) * ratio(jl,js)
            end if
          enddo
        end if
      enddo
      gconv = .true.
      !NO!$acc end kernels
      do j = 1, nstst, 1
        !NO!$acc kernels
        js = nlstst(j)
        !NO!$acc loop seq
        do jl = 1, n_points, 1
          if (ABS(y(jl,js) - zy(jl,js)) > ftol * y(jl,js) .AND. y(jl,js) > pmintnd(jl)) then
            gconv = .false.
          end if
        enddo
        !NO!$acc end kernels
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
      !NO!$acc kernels
      !NO!$acc loop seq
      do j = 1, nstst, 1
        js = nlstst(j)
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          zy(jl,js) = y(jl,js)
        enddo
      enddo
      !NO!$acc end kernels
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_ftoy

end module asad_ftoy_mod
