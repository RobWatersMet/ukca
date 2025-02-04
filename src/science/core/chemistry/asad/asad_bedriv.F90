module asad_bedriv_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_BEDRIV_MOD'
  public

  contains
  subroutine asad_bedriv(nslon, nslat, n_points, nlev)
    use asad_mod, only : cdt, ctype, dpd, dpw, f, frpx, jpcspf, jpif, jpmsp, jpna, jpnr, jpoo, jppj, jpsp, jpspec, jptk, ldepd, ldepw, nfrpx, nnfrp, nprkx, nspi, nstst, ntabfp, ntrkx, nuni, p, rk, speci, spj, spt, y, ydot
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use ukca_um_legacy_mod, only : mype
    use umPrintMgr, only : printstatus, prstatus_diag, ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    use asad_steady_mod, only : asad_steady
    use asad_ftoy_mod, only : asad_ftoy
    integer, parameter :: prodmax = 50
    integer, parameter :: besteps = 8
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: nlev
    integer, intent(in) :: nslon
    integer, intent(in) :: nslat
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer, save :: knspec
    integer :: jit
    integer :: j
    integer :: jp
    integer :: jr
    integer :: js
    integer :: k
    integer :: k1
    integer :: i
    integer :: ix
    integer :: jr1
    integer :: jr2
    integer :: jr3
    integer :: jr4
    integer :: jr5
    integer, save :: in2o5 = 0
    integer, save :: ino3 = 0
    integer, save :: ino2 = 0
    integer, save :: pn2o5 = 0
    integer, save :: tn2o5 = 0
    integer, save :: tno2no3 = 0
    logical, save :: first = .true.
    logical, save :: first_pass = .true.
    logical, save :: not_first_call = .false.
    integer, dimension(jpspec) :: kspec
    real, dimension(n_points) :: pdn
    real, dimension(n_points) :: l
    real, dimension(n_points) :: r1
    real, dimension(n_points) :: r2
    real, dimension(n_points) :: l1
    real, dimension(n_points) :: l2
    real, dimension(n_points) :: l3
    integer, dimension(jpspec,prodmax) :: uprodreac
    integer, dimension(jpspec,prodmax) :: ulossreac
    integer, dimension(jpspec,prodmax) :: bprodreac
    integer, dimension(jpspec,prodmax) :: blossreac
    integer, dimension(jpspec,prodmax) :: blosspartner
    integer, dimension(jpspec) :: nuintprod
    integer, dimension(jpspec) :: nbintprod
    integer, dimension(jpspec) :: nfracprod
    integer, dimension(jpspec) :: nuloss
    integer, dimension(jpspec) :: nbloss
    real, dimension(jpspec,prodmax) :: frac_prod
    integer :: iftr
    integer, dimension(jpcspf) :: ilftra
    integer :: iter
    integer, dimension(0:prodmax), save :: mod5
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_BEDRIV'
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
    integer :: idx_22
    integer :: idx_23
    integer :: idx_24
    integer :: idx_25
    integer :: idx_26
    integer :: idx_27
    integer :: idx_28
    integer :: idx_29
    integer :: idx_30
    integer :: idx_31
    integer :: idx_32
    integer :: idx_33
    integer :: idx_34
    integer :: idx_35
    integer :: idx_36
    integer :: idx_37
    integer :: idx_38
    integer :: idx_39
    integer :: idx_40
    integer :: idx_41
    integer :: idx_42
    integer :: idx_43
    integer :: idx_44
    integer :: idx_45
    integer :: idx_46
    integer :: idx_47
    integer :: idx_48
    integer :: idx_49
    integer :: idx_50

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first_pass) then
      if (first) then
        !$acc kernels
        first = .false.
        !$acc loop seq
        do js = 1, jpspec, 1
          if (speci(js) == 'N2O5      ') then
            in2o5 = js
          else
            if (speci(js) == 'NO3       ') then
              ino3 = js
            else
              if (speci(js) == 'NO2       ') then
                ino2 = js
              end if
            end if
          end if
        enddo
        !$acc loop seq
        do js = 1, jptk, 1
          if (spt(js,1) == 'N2O5      ') then
            tn2o5 = ntrkx(js)
          else
            if (spt(js,1) == 'NO2       ') then
              if (spt(js,2) == 'NO3       ') then
                tno2no3 = ntrkx(js)
              end if
            else
              if (spt(js,1) == 'NO3       ') then
                if (spt(js,2) == 'NO2       ') then
                  tno2no3 = ntrkx(js)
                end if
              end if
            end if
          end if
        enddo
        !$acc loop seq
        do js = 1, jppj, 1
          if (spj(js,1) == 'N2O5      ') then
            pn2o5 = nprkx(js)
          end if
        enddo
        knspec = 0
        !$acc loop seq
        do js = 1, jpspec, 1
          if (ctype(js) == jpna) then
            knspec = knspec + 1
            kspec(knspec) = js
          end if
        enddo
        !$acc loop seq
        do js = 1, jpspec, 1
          if (ctype(js) == jpsp) then
            knspec = knspec + 1
            kspec(knspec) = js
          end if
        enddo
        !$acc end kernels
        if (mype == 0 .AND. printstatus >= prstatus_diag) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') 'ASAD_BEDRIV:'
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A5,I4,A5,I4,A5,I4)') 'N2O5=', in2o5, ' NO3 ', ino3, ' NO2 ', ino2
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A6,I4,1X,A8,I4,1X,A5,I4)') 'tn2o5=', tn2o5, 'tno2no3=', tno2no3, 'pn2o5 ', pn2o5
          call umprint(ummessage, src='asad_bedriv')
        end if
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx = 1, prodmax, 1
          do idx_1 = 1, jpspec, 1
            uprodreac(idx_1,idx) = 0
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_2 = 1, prodmax, 1
          do idx_3 = 1, jpspec, 1
            ulossreac(idx_3,idx_2) = 0
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_4 = 1, prodmax, 1
          do idx_5 = 1, jpspec, 1
            bprodreac(idx_5,idx_4) = 0
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_6 = 1, prodmax, 1
          do idx_7 = 1, jpspec, 1
            blossreac(idx_7,idx_6) = 0
          enddo
        enddo
        !$acc loop gang vector independent collapse(2)
        do idx_8 = 1, prodmax, 1
          do idx_9 = 1, jpspec, 1
            frac_prod(idx_9,idx_8) = 1.0
          enddo
        enddo
        !$acc loop gang vector independent
        do idx_10 = 1, jpspec, 1
          nuintprod(idx_10) = 0
        enddo
        !$acc loop gang vector independent
        do idx_11 = 1, jpspec, 1
          nbintprod(idx_11) = 0
        enddo
        !$acc loop gang vector independent
        do idx_12 = 1, jpspec, 1
          nfracprod(idx_12) = 0
        enddo
        !$acc loop gang vector independent
        do idx_13 = 1, jpspec, 1
          nuloss(idx_13) = 0
        enddo
        !$acc loop gang vector independent
        do idx_14 = 1, jpspec, 1
          nbloss(idx_14) = 0
        enddo
        !$acc end kernels
        do jr = 1, jpnr, 1
          if (nfrpx(jr) == 0) then
            do jp = 3, jpmsp, 1
              js = nspi(jr,jp)
              if (js > 0) then
                if (jr <= nuni) then
                  nuintprod(js) = nuintprod(js) + 1
                  if (nuintprod(js) > prodmax) then
                    errcode = 1
                    cmessage(:) = 'Increase prodmax.'
                    call ereport('ASAD_BEDRIV', errcode, cmessage)
                  end if
                  uprodreac(js,nuintprod(js)) = jr
                else
                  nbintprod(js) = nbintprod(js) + 1
                  if (nbintprod(js) > prodmax) then
                    errcode = 1
                    cmessage(:) = 'Increase prodmax.'
                    call ereport('ASAD_BEDRIV', errcode, cmessage)
                  end if
                  bprodreac(js,nbintprod(js)) = jr
                end if
              end if
            enddo
          end if
          do jp = 1, 2, 1
            js = nspi(jr,jp)
            if (js > 0) then
              if (jr <= nuni) then
                nuloss(js) = nuloss(js) + 1
                if (nuloss(js) > prodmax) then
                  errcode = 1
                  cmessage(:) = 'Increase prodmax.'
                  call ereport('ASAD_BEDRIV', errcode, cmessage)
                end if
                ulossreac(js,nuloss(js)) = jr
              else
                nbloss(js) = nbloss(js) + 1
                if (nbloss(js) > prodmax) then
                  errcode = 1
                  cmessage(:) = 'Increase prodmax.'
                  call ereport('ASAD_BEDRIV', errcode, cmessage)
                end if
                blossreac(js,nbloss(js)) = jr
                if (jp == 1) then
                  blosspartner(js,nbloss(js)) = nspi(jr,2)
                else
                  blosspartner(js,nbloss(js)) = nspi(jr,1)
                end if
              end if
            end if
          enddo
        enddo
        !$acc kernels
        !$acc loop gang vector independent
        do idx_15 = 1, jpspec, 1
          nfracprod(idx_15) = nbintprod(idx_15)
        enddo
        !$acc end kernels
        do i = 1, nnfrp, 1
          js = ntabfp(i,1)
          jr = ntabfp(i,2)
          ix = ntabfp(i,3)
          nfracprod(js) = nfracprod(js) + 1
          if (nfracprod(js) > prodmax) then
            errcode = 1
            cmessage(:) = 'Increase prodmax.'
            call ereport('ASAD_BEDRIV', errcode, cmessage)
          end if
          bprodreac(js,nfracprod(js)) = jr
          frac_prod(js,nfracprod(js)) = frpx(ix)
        enddo
        if (mype == 0) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'NUINTPROD = ', nuintprod
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'NBINTPROD = ', nbintprod
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'NFRACPROD = ', nfracprod
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A9,I4)') 'NULOSS = ', nuloss
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A9,I4)') 'NBLOSS = ', nbloss
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'UPRODREAC = ', uprodreac(1, :)
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'BPRODREAC = ', bprodreac(1, :)
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'ULOSSREAC = ', ulossreac(1, :)
          call umprint(ummessage, src='asad_bedriv')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A12,I4)') 'BLOSSREAC = ', blossreac(1, :)
          call umprint(ummessage, src='asad_bedriv')
        end if
        !$acc kernels
        iftr = 0
        !$acc loop seq
        do js = 1, jpspec, 1
          if (ctype(js) == jpif .OR. ctype(js) == jpsp .OR. ctype(js) == jpoo) then
            iftr = iftr + 1
            ilftra(iftr) = js
          end if
        enddo
        !$acc loop gang vector independent
        do k = 0, prodmax, 1
          mod5(k) = MOD(k, 5)
        enddo
        !$acc end kernels
      end if
      first_pass = .false.
    end if
    iter = 1
    jit = 0
    call asad_ftoy(not_first_call, iter, jit, n_points, nslon, nslat, nlev)
    if (nstst /= 0) then
      call asad_steady(n_points)
    end if
    ydot = y
    do jit = 1, besteps, 1
      do j = 1, knspec, 1
        js = kspec(j)
        if (js /= in2o5) then
          k1 = nuintprod(js)
          if (mod5(k1) == 4) then
            !$acc kernels
            jr1 = uprodreac(js,k1)
            jr2 = uprodreac(js,k1 - 1)
            jr3 = uprodreac(js,k1 - 2)
            jr4 = uprodreac(js,k1 - 3)
            !$acc loop gang vector independent
            do idx_16 = 1, n_points, 1
              pdn(idx_16) = rk(idx_16 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_16 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) + rk(idx_16 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_16 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) + rk(idx_16 + (LBOUND(rk, dim=1) - 1),jr3) * y(idx_16 + (LBOUND(y, dim=1) - 1),nspi(jr3,1)) + rk(idx_16 + (LBOUND(rk, dim=1) - 1),jr4) * y(idx_16 + (LBOUND(y, dim=1) - 1),nspi(jr4,1))
            enddo
            !$acc end kernels
          else
            if (mod5(k1) == 3) then
              !$acc kernels
              jr1 = uprodreac(js,k1)
              jr2 = uprodreac(js,k1 - 1)
              jr3 = uprodreac(js,k1 - 3)
              !$acc loop gang vector independent
              do idx_17 = 1, n_points, 1
                pdn(idx_17) = rk(idx_17 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_17 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) + rk(idx_17 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_17 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) + rk(idx_17 + (LBOUND(rk, dim=1) - 1),jr3) * y(idx_17 + (LBOUND(y, dim=1) - 1),nspi(jr3,1))
              enddo
              !$acc end kernels
            else
              if (mod5(k1) == 2) then
                !$acc kernels
                jr1 = uprodreac(js,k1)
                jr2 = uprodreac(js,k1 - 1)
                !$acc loop gang vector independent
                do idx_18 = 1, n_points, 1
                  pdn(idx_18) = rk(idx_18 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_18 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) + rk(idx_18 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_18 + (LBOUND(y, dim=1) - 1),nspi(jr2,1))
                enddo
                !$acc end kernels
              else
                if (mod5(k1) == 1) then
                  !$acc kernels
                  jr1 = uprodreac(js,k1)
                  !$acc loop gang vector independent
                  do idx_19 = 1, n_points, 1
                    pdn(idx_19) = rk(idx_19 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_19 + (LBOUND(y, dim=1) - 1),nspi(jr1,1))
                  enddo
                  !$acc end kernels
                else
                  if (mod5(k1) == 0) then
                    !$acc kernels
                    !$acc loop gang vector independent
                    do idx_20 = 1, n_points, 1
                      pdn(idx_20) = 0.0
                    enddo
                    !$acc end kernels
                  end if
                end if
              end if
            end if
          end if
          !$acc kernels
          !$acc loop seq
          do k = 1, k1 - 4, 5
            jr1 = uprodreac(js,k)
            jr2 = uprodreac(js,k + 1)
            jr3 = uprodreac(js,k + 2)
            jr4 = uprodreac(js,k + 3)
            jr5 = uprodreac(js,k + 4)
            !$acc loop gang vector independent
            do idx_21 = 1, n_points, 1
              pdn(idx_21) = pdn(idx_21) + rk(idx_21 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_21 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) + rk(idx_21 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_21 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) + rk(idx_21 + (LBOUND(rk, dim=1) - 1),jr3) * y(idx_21 + (LBOUND(y, dim=1) - 1),nspi(jr3,1)) + rk(idx_21 + (LBOUND(rk, dim=1) - 1),jr4) * y(idx_21 + (LBOUND(y, dim=1) - 1),nspi(jr4,1)) + rk(idx_21 + (LBOUND(rk, dim=1) - 1),jr5) * y(idx_21 + (LBOUND(y, dim=1) - 1),nspi(jr5,1))
            enddo
          enddo
          !$acc end kernels
          k = nbintprod(js)
          do while (k >= 5)
            !$acc kernels
            jr1 = bprodreac(js,k)
            jr2 = bprodreac(js,k - 1)
            jr3 = bprodreac(js,k - 2)
            jr4 = bprodreac(js,k - 3)
            jr5 = bprodreac(js,k - 4)
            !$acc loop gang vector independent
            do idx_22 = 1, n_points, 1
              pdn(idx_22) = pdn(idx_22) + rk(idx_22 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr1,2)) + rk(idx_22 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr2,2)) + rk(idx_22 + (LBOUND(rk, dim=1) - 1),jr3) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr3,1)) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr3,2)) + rk(idx_22 + (LBOUND(rk, dim=1) - 1),jr4) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr4,1)) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr4,2)) + rk(idx_22 + (LBOUND(rk, dim=1) - 1),jr5) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr5,1)) * y(idx_22 + (LBOUND(y, dim=1) - 1),nspi(jr5,2))
            enddo
            k = k - 5
            !$acc end kernels
          end do
          if (k == 4) then
            !$acc kernels
            jr1 = bprodreac(js,1)
            jr2 = bprodreac(js,2)
            jr3 = bprodreac(js,3)
            jr4 = bprodreac(js,4)
            !$acc loop gang vector independent
            do idx_23 = 1, n_points, 1
              pdn(idx_23) = pdn(idx_23) + rk(idx_23 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr1,2)) + rk(idx_23 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr2,2)) + rk(idx_23 + (LBOUND(rk, dim=1) - 1),jr3) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr3,1)) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr3,2)) + rk(idx_23 + (LBOUND(rk, dim=1) - 1),jr4) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr4,1)) * y(idx_23 + (LBOUND(y, dim=1) - 1),nspi(jr4,2))
            enddo
            !$acc end kernels
          else
            if (k == 3) then
              !$acc kernels
              jr1 = bprodreac(js,1)
              jr2 = bprodreac(js,2)
              jr3 = bprodreac(js,3)
              !$acc loop gang vector independent
              do idx_24 = 1, n_points, 1
                pdn(idx_24) = pdn(idx_24) + rk(idx_24 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_24 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) * y(idx_24 + (LBOUND(y, dim=1) - 1),nspi(jr1,2)) + rk(idx_24 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_24 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) * y(idx_24 + (LBOUND(y, dim=1) - 1),nspi(jr2,2)) + rk(idx_24 + (LBOUND(rk, dim=1) - 1),jr3) * y(idx_24 + (LBOUND(y, dim=1) - 1),nspi(jr3,1)) * y(idx_24 + (LBOUND(y, dim=1) - 1),nspi(jr3,2))
              enddo
              !$acc end kernels
            else
              if (k == 2) then
                !$acc kernels
                jr1 = bprodreac(js,1)
                jr2 = bprodreac(js,2)
                !$acc loop gang vector independent
                do idx_25 = 1, n_points, 1
                  pdn(idx_25) = pdn(idx_25) + rk(idx_25 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_25 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) * y(idx_25 + (LBOUND(y, dim=1) - 1),nspi(jr1,2)) + rk(idx_25 + (LBOUND(rk, dim=1) - 1),jr2) * y(idx_25 + (LBOUND(y, dim=1) - 1),nspi(jr2,1)) * y(idx_25 + (LBOUND(y, dim=1) - 1),nspi(jr2,2))
                enddo
                !$acc end kernels
              else
                if (k == 1) then
                  !$acc kernels
                  jr1 = bprodreac(js,1)
                  !$acc loop gang vector independent
                  do idx_26 = 1, n_points, 1
                    pdn(idx_26) = pdn(idx_26) + rk(idx_26 + (LBOUND(rk, dim=1) - 1),jr1) * y(idx_26 + (LBOUND(y, dim=1) - 1),nspi(jr1,1)) * y(idx_26 + (LBOUND(y, dim=1) - 1),nspi(jr1,2))
                  enddo
                  !$acc end kernels
                end if
              end if
            end if
          end if
          !$acc kernels
          !$acc loop seq
          do k = nbintprod(js) + 1, nfracprod(js), 1
            jr = bprodreac(js,k)
            !$acc loop gang vector independent
            do idx_27 = 1, n_points, 1
              pdn(idx_27) = pdn(idx_27) + rk(idx_27 + (LBOUND(rk, dim=1) - 1),jr) * y(idx_27 + (LBOUND(y, dim=1) - 1),nspi(jr,1)) * y(idx_27 + (LBOUND(y, dim=1) - 1),nspi(jr,2)) * frac_prod(js,k)
            enddo
          enddo
          k1 = nuloss(js)
          !$acc end kernels
          if (mod5(k1) == 4) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_28 = 1, n_points, 1
              l(idx_28) = rk(idx_28 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1)) + rk(idx_28 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1 - 1)) + rk(idx_28 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1 - 2)) + rk(idx_28 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1 - 3))
            enddo
            !$acc end kernels
          else
            if (mod5(k1) == 3) then
              !$acc kernels
              !$acc loop gang vector independent
              do idx_29 = 1, n_points, 1
                l(idx_29) = rk(idx_29 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1)) + rk(idx_29 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1 - 1)) + rk(idx_29 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1 - 2))
              enddo
              !$acc end kernels
            else
              if (mod5(k1) == 2) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_30 = 1, n_points, 1
                  l(idx_30) = rk(idx_30 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1)) + rk(idx_30 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1 - 1))
                enddo
                !$acc end kernels
              else
                if (mod5(k1) == 1) then
                  !$acc kernels
                  !$acc loop gang vector independent
                  do idx_31 = 1, n_points, 1
                    l(idx_31) = rk(idx_31 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k1))
                  enddo
                  !$acc end kernels
                else
                  if (mod5(k1) == 0) then
                    !$acc kernels
                    !$acc loop gang vector independent
                    do idx_32 = 1, n_points, 1
                      l(idx_32) = 0.0
                    enddo
                    !$acc end kernels
                  end if
                end if
              end if
            end if
          end if
          !$acc kernels
          !$acc loop seq
          do k = 1, k1 - 4, 5
            !$acc loop gang vector independent
            do idx_33 = 1, n_points, 1
              l(idx_33) = l(idx_33) + rk(idx_33 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k)) + rk(idx_33 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k + 1)) + rk(idx_33 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k + 2)) + rk(idx_33 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k + 3)) + rk(idx_33 + (LBOUND(rk, dim=1) - 1),ulossreac(js,k + 4))
            enddo
          enddo
          !$acc end kernels
          k = nbloss(js)
          do while (k >= 5)
            !$acc kernels
            !$acc loop gang vector independent
            do idx_34 = 1, n_points, 1
              l(idx_34) = l(idx_34) + rk(idx_34 + (LBOUND(rk, dim=1) - 1),blossreac(js,k)) * y(idx_34 + (LBOUND(y, dim=1) - 1),blosspartner(js,k)) + rk(idx_34 + (LBOUND(rk, dim=1) - 1),blossreac(js,k - 1)) * y(idx_34 + (LBOUND(y, dim=1) - 1),blosspartner(js,k - 1)) + rk(idx_34 + (LBOUND(rk, dim=1) - 1),blossreac(js,k - 2)) * y(idx_34 + (LBOUND(y, dim=1) - 1),blosspartner(js,k - 2)) + rk(idx_34 + (LBOUND(rk, dim=1) - 1),blossreac(js,k - 3)) * y(idx_34 + (LBOUND(y, dim=1) - 1),blosspartner(js,k - 3)) + rk(idx_34 + (LBOUND(rk, dim=1) - 1),blossreac(js,k - 4)) * y(idx_34 + (LBOUND(y, dim=1) - 1),blosspartner(js,k - 4))
            enddo
            k = k - 5
            !$acc end kernels
          end do
          if (k == 4) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_35 = 1, n_points, 1
              l(idx_35) = l(idx_35) + rk(idx_35 + (LBOUND(rk, dim=1) - 1),blossreac(js,4)) * y(idx_35 + (LBOUND(y, dim=1) - 1),blosspartner(js,4)) + rk(idx_35 + (LBOUND(rk, dim=1) - 1),blossreac(js,3)) * y(idx_35 + (LBOUND(y, dim=1) - 1),blosspartner(js,3)) + rk(idx_35 + (LBOUND(rk, dim=1) - 1),blossreac(js,2)) * y(idx_35 + (LBOUND(y, dim=1) - 1),blosspartner(js,2)) + rk(idx_35 + (LBOUND(rk, dim=1) - 1),blossreac(js,1)) * y(idx_35 + (LBOUND(y, dim=1) - 1),blosspartner(js,1))
            enddo
            !$acc end kernels
          else
            if (k == 3) then
              !$acc kernels
              !$acc loop gang vector independent
              do idx_36 = 1, n_points, 1
                l(idx_36) = l(idx_36) + rk(idx_36 + (LBOUND(rk, dim=1) - 1),blossreac(js,3)) * y(idx_36 + (LBOUND(y, dim=1) - 1),blosspartner(js,3)) + rk(idx_36 + (LBOUND(rk, dim=1) - 1),blossreac(js,2)) * y(idx_36 + (LBOUND(y, dim=1) - 1),blosspartner(js,2)) + rk(idx_36 + (LBOUND(rk, dim=1) - 1),blossreac(js,1)) * y(idx_36 + (LBOUND(y, dim=1) - 1),blosspartner(js,1))
              enddo
              !$acc end kernels
            else
              if (k == 2) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_37 = 1, n_points, 1
                  l(idx_37) = l(idx_37) + rk(idx_37 + (LBOUND(rk, dim=1) - 1),blossreac(js,2)) * y(idx_37 + (LBOUND(y, dim=1) - 1),blosspartner(js,2)) + rk(idx_37 + (LBOUND(rk, dim=1) - 1),blossreac(js,1)) * y(idx_37 + (LBOUND(y, dim=1) - 1),blosspartner(js,1))
                enddo
                !$acc end kernels
              else
                if (k == 1) then
                  !$acc kernels
                  !$acc loop gang vector independent
                  do idx_38 = 1, n_points, 1
                    l(idx_38) = l(idx_38) + rk(idx_38 + (LBOUND(rk, dim=1) - 1),blossreac(js,1)) * y(idx_38 + (LBOUND(y, dim=1) - 1),blosspartner(js,1))
                  enddo
                  !$acc end kernels
                end if
              end if
            end if
          end if
          if (ldepd(js)) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_39 = 1, n_points, 1
              l(idx_39) = l(idx_39) + dpd(idx_39 + (LBOUND(dpd, dim=1) - 1),js)
            enddo
            !$acc end kernels
          end if
          if (ldepw(js)) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_40 = 1, n_points, 1
              l(idx_40) = l(idx_40) + dpw(idx_40 + (LBOUND(dpw, dim=1) - 1),js)
            enddo
            !$acc end kernels
          end if
          if (ctype(js) == jpna) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_41 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
              y(idx_41,js) = pdn(idx_41 + (1 - LBOUND(y, dim=1))) / l(idx_41 + (1 - LBOUND(y, dim=1)))
            enddo
            !$acc end kernels
          else
            if (js /= ino3) then
              y(:,js) = (ydot(:,js) + cdt * pdn(:)) / (1.0 + cdt * l(:))
            else
              !$acc kernels
              !$acc loop gang vector independent
              do idx_42 = 1, n_points, 1
                r1(idx_42) = rk(idx_42 + (LBOUND(rk, dim=1) - 1),pn2o5) + rk(idx_42 + (LBOUND(rk, dim=1) - 1),tn2o5)
              enddo
              !$acc loop gang vector independent
              do idx_43 = 1, n_points, 1
                pdn(idx_43) = pdn(idx_43) - r1(idx_43) * y(idx_43 + (LBOUND(y, dim=1) - 1),in2o5)
              enddo
              !$acc loop gang vector independent
              do idx_44 = 1, n_points, 1
                r2(idx_44) = rk(idx_44 + (LBOUND(rk, dim=1) - 1),tno2no3) * y(idx_44 + (LBOUND(y, dim=1) - 1),ino2)
              enddo
              !$acc loop gang vector independent
              do idx_45 = 1, n_points, 1
                l1(idx_45) = 0.0
              enddo
              !$acc loop seq
              do k = 1, nuloss(in2o5), 1
                !$acc loop gang vector independent
                do idx_46 = 1, n_points, 1
                  l1(idx_46) = l1(idx_46) + rk(idx_46 + (LBOUND(rk, dim=1) - 1),ulossreac(in2o5,k))
                enddo
              enddo
              !$acc loop seq
              do k = 1, nbloss(in2o5), 1
                !$acc loop gang vector independent
                do idx_47 = 1, n_points, 1
                  l1(idx_47) = l1(idx_47) + rk(idx_47 + (LBOUND(rk, dim=1) - 1),blossreac(in2o5,k)) * y(idx_47 + (LBOUND(y, dim=1) - 1),blosspartner(in2o5,k))
                enddo
              enddo
              !$acc end kernels
              if (ldepd(in2o5)) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_48 = 1, n_points, 1
                  l1(idx_48) = l1(idx_48) + dpd(idx_48 + (LBOUND(dpd, dim=1) - 1),in2o5)
                enddo
                !$acc end kernels
              end if
              if (ldepw(in2o5)) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_49 = 1, n_points, 1
                  l1(idx_49) = l1(idx_49) + dpw(idx_49 + (LBOUND(dpw, dim=1) - 1),in2o5)
                enddo
                !$acc end kernels
              end if
              !$acc kernels
              l2(:) = 1.0 + cdt * l(:)
              l3(:) = 1.0 + cdt * l1(:)
              y(:,ino3) = (l3(:) * (ydot(:,js) + p * cdt) + r1(:) * cdt * ydot(:,in2o5)) / (l3(:) * l2(:) - r1(:) * r2(:) * cdt * cdt)
              y(:,in2o5) = (ydot(:,in2o5) + r2(:) * cdt * y(:,js)) / l3(:)
              !$acc end kernels
            end if
          end if
        end if
      enddo
    enddo
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do j = 1, jpcspf, 1
      do idx_50 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_50,j) = y(idx_50 + (LBOUND(y, dim=1) - LBOUND(f, dim=1)),ilftra(j))
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_bedriv

end module asad_bedriv_mod
