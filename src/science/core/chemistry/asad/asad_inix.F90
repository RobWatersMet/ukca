module asad_inix_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_INIX_MOD'
  public

  contains
  subroutine asad_inix()
    use asad_mod, only : ctype, jpfm, jpfrpx, jpif, jpmsp, jpna, jpnr, jpoo, jpsp, jpspb, jpspec, ldepd, ldepw, ndepd, ndepw, nf, nfrpx, ngrp, nlall, nldepd, nldepw, nldepx, nlf, nlstst, nnfrp, nprdx1, nprdx2, nprdx3, nspi, nstst, ntabfp
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: j
    integer :: jr
    integer :: js
    integer :: jp
    integer :: iprd
    integer :: igrp3
    integer :: iloss
    integer :: ngmax
    integer :: iss
    integer :: ispec
    integer :: jdry
    integer :: jwet
    integer :: jpnpx3
    integer, dimension(3) :: igrp
    integer, dimension(jpspec) :: idepd
    integer, dimension(jpspec) :: idepw
    CHARACTER(LEN = 2) :: itype
    CHARACTER(LEN = errormessagelength) :: cmessage
    logical :: gdry
    logical :: gwet
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INIX'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    jpnpx3 = jpnr / (3 * 3) + 3 * 3
    nstst = 0
    nf = 0
    do js = 1, jpspec, 1
      itype(:) = ctype(js)
      nlall(js) = js
      if (itype(:) == jpfm .OR. itype(:) == jpif .OR. itype(:) == jpna) then
        nstst = nstst + 1
        nlstst(nstst) = js
      end if
      if (itype(:) == jpfm .OR. itype(:) == jpif .OR. itype(:) == jpsp .OR. itype(:) == jpoo) then
        nf = nf + 1
        nlf(nf) = js
      end if
    enddo
    !$acc kernels
    !$acc loop seq
    do js = 1, jpspec, 1
      iprd = 0
      igrp3 = 0
      !$acc loop seq
      do jr = 1, jpnr, 1
        if (nfrpx(jr) == 0) then
          !$acc loop seq
          do jp = 3, jpmsp, 1
            if (nspi(jr,jp) == js) then
              iprd = iprd + 1
              igrp(iprd) = jr
              if (iprd == 3) then
                iprd = 0
                igrp3 = igrp3 + 1
                !$acc loop gang vector independent
                do j = 1, 3, 1
                  nprdx3(j,igrp3,js) = igrp(j)
                enddo
              end if
            end if
          enddo
        end if
      enddo
      ngrp(js,3) = igrp3
      if (iprd == 2) then
        ngrp(js,2) = 1
        nprdx2(1,js) = igrp(1)
        nprdx2(2,js) = igrp(2)
      else
        if (iprd == 1) then
          ngrp(js,1) = 1
          nprdx1(js) = igrp(1)
        end if
      end if
      iloss = 0
      igrp3 = 0
      !$acc loop seq
      do jr = 1, jpnr, 1
        !$acc loop seq
        do jp = 1, 2, 1
          if (nspi(jr,jp) == js) then
            iloss = iloss + 1
            igrp(iloss) = jr
            if (iloss == 3) then
              iloss = 0
              igrp3 = igrp3 + 1
              !$acc loop gang vector independent
              do j = 1, 3, 1
                nprdx3(j,igrp3,jpspec + js) = igrp(j)
              enddo
            end if
          end if
        enddo
      enddo
      ngrp(jpspec + js,3) = igrp3
      if (iloss == 2) then
        ngrp(js + jpspec,2) = 1
        nprdx2(1,jpspec + js) = igrp(1)
        nprdx2(2,jpspec + js) = igrp(2)
      else
        if (iloss == 1) then
          ngrp(js + jpspec,1) = 1
          nprdx1(jpspec + js) = igrp(1)
        end if
      end if
    enddo
    !$acc end kernels
    ngmax = ngrp(1,3)
    !$acc kernels
    !$acc loop seq
    do js = 2, jpspec, 1
      ngmax = MAX(ngmax, ngrp(js,3))
    enddo
    !$acc end kernels
    if (ngmax > jpnpx3) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(2A,I0)') '** ASAD ERROR: The parameter jpnpx3 is too', ' small. Please change it not less than ', ngmax
      call umprint(ummessage, src='asad_inix')
      cmessage(:) = ' Parameter jpnpx3 is too small'
      call ereport('ASAD_INIX', ngmax, cmessage)
    end if
    nnfrp = 0
    do jr = 1, jpnr, 1
      if (nfrpx(jr) /= 0) then
        do jp = 3, jpspb, 1
          iss = nspi(jr,jp)
          if (iss /= 0) then
            nnfrp = nnfrp + 1
            if (nnfrp > jpfrpx) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2A)') '*ASAD ERROR: The parameter jpfrpx is', ' too small. Please increase it.'
              call umprint(ummessage, src='asad_inix')
              cmessage(:) = ' Parameter jpfrpx is too small for no. ' // 'of fractional products'
              call ereport('ASAD_INIX', nnfrp, cmessage)
            end if
            ntabfp(nnfrp,1) = iss
            ntabfp(nnfrp,2) = jr
            ntabfp(nnfrp,3) = nfrpx(jr) + (jp - 3)
          end if
        enddo
      end if
    enddo
    !$acc kernels
    ndepd = 0
    ndepw = 0
    !$acc loop seq
    do js = 1, jpspec, 1
      if (ldepd(js)) then
        ndepd = ndepd + 1
        nldepd(ndepd) = js
      end if
      if (ldepw(js)) then
        ndepw = ndepw + 1
        nldepw(ndepw) = js
      end if
    enddo
    !$acc loop gang vector independent
    do js = 1, jpspec, 1
      idepd(js) = nldepd(js)
      idepw(js) = nldepw(js)
    enddo
    nldepx(1) = 7
    nldepx(2) = nldepx(1) - 1
    do ispec = 1, jpspec, 1
      gdry = .false.
      gwet = .false.
      !$acc loop seq
      do j = 1, jpspec, 1
        if (idepd(j) == ispec) then
          gdry = .true.
          jdry = j
        end if
        if (idepw(j) == ispec) then
          gwet = .true.
          jwet = j
        end if
      enddo
      if (gdry .AND. gwet) then
        nldepx(2) = nldepx(2) + 1
        nldepx(nldepx(2)) = ispec
        idepd(jdry) = 0
        idepw(jwet) = 0
      end if
    enddo
    nldepx(3) = nldepx(2) + 1
    nldepx(4) = nldepx(3) - 1
    do j = 1, jpspec, 1
      if (idepd(j) /= 0) then
        nldepx(4) = nldepx(4) + 1
        nldepx(nldepx(4)) = idepd(j)
      end if
    enddo
    nldepx(5) = nldepx(4) + 1
    nldepx(6) = nldepx(5) - 1
    do j = 1, jpspec, 1
      if (idepw(j) /= 0) then
        nldepx(6) = nldepx(6) + 1
        nldepx(nldepx(6)) = idepw(j)
      end if
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_inix

end module asad_inix_mod
