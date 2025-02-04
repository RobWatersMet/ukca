module asad_inijac_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_INIJAC_MOD'
  public

  contains
  subroutine asad_inijac()
    use asad_mod, only : ctype, jpcspf, jpfm, jpmsp, jpna, jpnr, jppjac, jpspec, madvtr, moffam, n_in_ss, njacx1, njacx2, njacx3, njcgrp, njcoss, njcoth, nmpjac, nmsjac, nmzjac, nodd, npjac1, nsjac1, nsph, nspho2, nspi, nspn, nspno, nspo1d, nspo3, nspo3p, nspoh, nstst, ntrho2, ntrno, ntro3, ntroh, nzjac1, o3p_in_ss, specf, speci
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, dimension(jpmsp) :: ifam
    integer, dimension(jpmsp) :: itr
    integer, dimension(2 * jpnr,jpcspf) :: ifneg
    integer, dimension(jppjac,jpcspf) :: ifpos
    integer, dimension(jpcspf) :: ifn
    integer, dimension(jpcspf) :: ifp
    integer, dimension(jpcspf) :: ifz
    integer, dimension(jpnr,jpcspf) :: ifzer
    integer, dimension(jpcspf) :: ifs
    integer, dimension(jpnr,jpcspf) :: ifsss
    integer :: ir1
    integer :: ir2
    integer :: ip1
    integer :: ip2
    integer :: ip3
    integer :: iss
    integer :: njx
    integer :: j
    integer :: jc
    integer :: je
    integer :: jg
    integer :: jr
    integer :: jtr
    integer :: i
    integer :: irem
    integer :: jpnjcx3
    CHARACTER(LEN = 2) :: itype(jpmsp)
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INIJAC'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    jpnjcx3 = jpnr / (3 * 3) + 3 * 3
    !$acc loop gang vector independent
    do idx = LBOUND(ifn, dim=1), UBOUND(ifn, dim=1), 1
      ifn(idx) = 0
    enddo
    !$acc loop gang vector independent
    do idx_1 = LBOUND(ifp, dim=1), UBOUND(ifp, dim=1), 1
      ifp(idx_1) = 0
    enddo
    !$acc loop gang vector independent
    do idx_2 = LBOUND(ifz, dim=1), UBOUND(ifz, dim=1), 1
      ifz(idx_2) = 0
    enddo
    !$acc loop gang vector independent
    do idx_3 = LBOUND(ifs, dim=1), UBOUND(ifs, dim=1), 1
      ifs(idx_3) = 0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_4 = LBOUND(ifzer, dim=2), UBOUND(ifzer, dim=2), 1
      do idx_5 = LBOUND(ifzer, dim=1), UBOUND(ifzer, dim=1), 1
        ifzer(idx_5,idx_4) = 0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_6 = LBOUND(ifsss, dim=2), UBOUND(ifsss, dim=2), 1
      do idx_7 = LBOUND(ifsss, dim=1), UBOUND(ifsss, dim=1), 1
        ifsss(idx_7,idx_6) = 0
      enddo
    enddo
    !$acc end kernels
    do jr = 1, jpnr, 1
      !$acc kernels
      ir1 = nspi(jr,1)
      ir2 = nspi(jr,2)
      ip1 = nspi(jr,3)
      ip2 = nspi(jr,4)
      ip3 = nspi(jr,5)
      !$acc loop gang vector independent
      do j = 1, jpmsp, 1
        ifam(j) = 0
        itype(j) = '  '
        itr(j) = 0
      enddo
      !$acc end kernels
      if (ir1 /= 0) then
        ifam(1) = moffam(ir1)
        itype(1) = ctype(ir1)
        itr(1) = madvtr(ir1)
      end if
      if (ir2 /= 0) then
        ifam(2) = moffam(ir2)
        itype(2) = ctype(ir2)
        itr(2) = madvtr(ir2)
      end if
      if (ip1 /= 0) then
        ifam(3) = moffam(ip1)
        itype(3) = ctype(ip1)
        itr(3) = madvtr(ip1)
      end if
      if (ip2 /= 0) then
        ifam(4) = moffam(ip2)
        itype(4) = ctype(ip2)
        itr(4) = madvtr(ip2)
      end if
      if (ip3 /= 0) then
        ifam(5) = moffam(ip3)
        itype(5) = ctype(ip3)
        itr(5) = madvtr(ip3)
      end if
      if (ifam(1) /= 0) then
        iss = -nodd(ir1)
        if (ifam(2) == ifam(1) .AND. itype(2) == jpfm) then
          iss = iss - nodd(ir2)
        end if
        if (ifam(3) == ifam(1) .AND. itype(3) == jpfm) then
          iss = iss + nodd(ip1)
        end if
        if (ifam(4) == ifam(1) .AND. itype(4) == jpfm) then
          iss = iss + nodd(ip2)
        end if
        if (ifam(5) == ifam(1) .AND. itype(5) == jpfm) then
          iss = iss + nodd(ip3)
        end if
        if (iss < 0) then
          do je = 1, ABS(iss), 1
            ifn(ifam(1)) = ifn(ifam(1)) + 1
            if (ifn(ifam(1)) > 2 * jpnr) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifneg dimension.'
              call umprint(ummessage, src='asad_inijac')
              cmessage(:) = 'ARRAY TOO SMALL: Increase ifneg dimension'
              errcode = je
              call ereport('ASAD_INIJAC', errcode, cmessage)
            end if
            ifneg(ifn(ifam(1)),ifam(1)) = jr
          enddo
        else
          if (iss > 0) then
            do je = 1, iss, 1
              ifp(ifam(1)) = ifp(ifam(1)) + 1
              if (ifp(ifam(1)) > jppjac) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifpos dimension.'
                call umprint(ummessage, src='asad_inijac')
                cmessage(:) = 'ARRAY TOO SMALL: Increase ifpos dimension'
                errcode = je
                call ereport('ASAD_INIJAC', errcode, cmessage)
              end if
              ifpos(ifp(ifam(1)),ifam(1)) = jr
            enddo
          end if
        end if
      end if
      if (ifam(2) /= 0) then
        iss = -nodd(ir2)
        if (ifam(1) == ifam(2) .AND. itype(2) == jpfm) then
          iss = iss - nodd(ir1)
        end if
        if (ifam(3) == ifam(2) .AND. itype(3) == jpfm) then
          iss = iss + nodd(ip1)
        end if
        if (ifam(4) == ifam(2) .AND. itype(4) == jpfm) then
          iss = iss + nodd(ip2)
        end if
        if (ifam(5) == ifam(2) .AND. itype(5) == jpfm) then
          iss = iss + nodd(ip3)
        end if
        if (iss < 0) then
          do je = 1, ABS(iss), 1
            ifn(ifam(2)) = ifn(ifam(2)) + 1
            if (ifn(ifam(2)) > 2 * jpnr) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifneg dimension.'
              call umprint(ummessage, src='asad_inijac')
              cmessage(:) = 'ARRAY TOO SMALL: Increase ifneg'
              errcode = je
              call ereport('ASAD_INIJAC', errcode, cmessage)
            end if
            ifneg(ifn(ifam(2)),ifam(2)) = jr
          enddo
        else
          if (iss > 0) then
            do je = 1, iss, 1
              ifp(ifam(2)) = ifp(ifam(2)) + 1
              if (ifp(ifam(2)) > jppjac) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifpos dimension.'
                call umprint(ummessage, src='asad_inijac')
                cmessage(:) = 'ARRAY TOO SMALL: Increase ifpos dimension'
                errcode = je
                call ereport('ASAD_INIJAC', errcode, cmessage)
              end if
              ifpos(ifp(ifam(2)),ifam(2)) = jr
            enddo
          end if
        end if
      end if
      if (itr(1) /= 0) then
        iss = -1
        if (itr(2) == itr(1) .AND. ir2 /= 0) then
          iss = iss - 1
        end if
        if (itr(3) == itr(1) .AND. ip1 /= 0) then
          iss = iss + 1
        end if
        if (itr(4) == itr(1) .AND. ip2 /= 0) then
          iss = iss + 1
        end if
        if (itr(5) == itr(1) .AND. ip3 /= 0) then
          iss = iss + 1
        end if
        if (iss < 0) then
          do je = 1, ABS(iss), 1
            ifn(itr(1)) = ifn(itr(1)) + 1
            if (ifn(itr(1)) > 2 * jpnr) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifneg dimension.'
              call umprint(ummessage, src='asad_inijac')
              cmessage(:) = 'ARRAY TOO SMALL: Increase ifneg'
              errcode = je
              call ereport('ASAD_INIJAC', errcode, cmessage)
            end if
            ifneg(ifn(itr(1)),itr(1)) = jr
          enddo
        else
          if (iss > 0) then
            do je = 1, iss, 1
              ifp(itr(1)) = ifp(itr(1)) + 1
              if (ifp(itr(1)) > jppjac) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifpos dimension.'
                call umprint(ummessage, src='asad_inijac')
                cmessage(:) = 'ARRAY TOO SMALL: Increase ifpos dimension'
                errcode = je
                call ereport('ASAD_INIJAC', errcode, cmessage)
              end if
              ifpos(ifp(itr(1)),itr(1)) = jr
            enddo
          end if
        end if
      end if
      if (itr(2) /= 0) then
        iss = -1
        if (itr(2) == itr(1) .AND. ir1 /= 0) then
          iss = iss - 1
        end if
        if (itr(3) == itr(2) .AND. ip1 /= 0) then
          iss = iss + 1
        end if
        if (itr(4) == itr(2) .AND. ip2 /= 0) then
          iss = iss + 1
        end if
        if (itr(5) == itr(2) .AND. ip3 /= 0) then
          iss = iss + 1
        end if
        if (iss < 0) then
          do je = 1, ABS(iss), 1
            ifn(itr(2)) = ifn(itr(2)) + 1
            if (ifn(itr(2)) > 2 * jpnr) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifneg dimension.'
              call umprint(ummessage, src='asad_inijac')
              cmessage(:) = 'ARRAY TOO SMALL: Increase igneg dimension'
              errcode = je
              call ereport('ASAD_INIJAC', errcode, cmessage)
            end if
            ifneg(ifn(itr(2)),itr(2)) = jr
          enddo
        else
          if (iss > 0) then
            do je = 1, iss, 1
              ifp(itr(2)) = ifp(itr(2)) + 1
              if (ifp(itr(2)) > jppjac) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(2A)') 'INTERNAL ASAD ERROR: array too small ', ' in inijac.f. Increase ifpos dimension.'
                call umprint(ummessage, src='asad_inijac')
                cmessage(:) = 'ARRAY TOO SMALL: Increase ifpos dimension'
                errcode = je
                call ereport('ASAD_INIJAC', errcode, cmessage)
              end if
              ifpos(ifp(itr(2)),itr(2)) = jr
            enddo
          end if
        end if
      end if
      !$acc kernels
      !$acc loop gang vector independent
      do j = 1, jpmsp, 1
        njcoth(jr,j) = itr(j)
      enddo
      !$acc loop seq
      do j = 1, 2, 1
        if (itr(j) /= 0) then
          ifz(itr(j)) = ifz(itr(j)) + 1
          ifzer(ifz(itr(j)),itr(j)) = jr
        end if
      enddo
      !$acc loop seq
      do j = 1, jpmsp, 1
        njcoss(jr,j) = 0
        if (itype(j) == jpna) then
          if (speci(nspi(jr,j)) == 'O(1D)') then
            if (j <= 2) then
              ifs(1) = ifs(1) + 1
              ifsss(ifs(1),1) = jr
            end if
            njcoss(jr,j) = 1
          else
            if (o3p_in_ss .AND. speci(nspi(jr,j)) == 'O(3P)     ') then
              if (j <= 2) then
                ifs(2) = ifs(2) + 1
                ifsss(ifs(2),2) = jr
              end if
              njcoss(jr,j) = 2
            else
              if (n_in_ss .AND. speci(nspi(jr,j)) == 'N         ') then
                if (j <= 2) then
                  ifs(3) = ifs(3) + 1
                  ifsss(ifs(3),3) = jr
                end if
                njcoss(jr,j) = 3
              end if
            end if
          end if
        end if
      enddo
      !$acc end kernels
    enddo
    do jc = 1, jpcspf, 1
      njx = ifn(jc) / 3
      if (njx > jpnjcx3) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(3A,I0)') 'FATAL ASAD ERROR: array njacx3 too ', 'small in inijac.f. Increase jpnjcx3 to at ', 'least ', njx
        call umprint(ummessage, src='asad_inijac')
        cmessage(:) = 'ARRAY NJACX3 TOO SMALL: Increase jpnjcx3'
        errcode = jc
        call ereport('ASAD_INIJAC', errcode, cmessage)
      end if
      !$acc kernels
      i = 1
      njcgrp(jc,3) = njx
      !$acc loop seq
      do jg = 1, njcgrp(jc,3), 1
        njacx3(1,jg,jc) = ifneg(i,jc)
        njacx3(2,jg,jc) = ifneg(i + 1,jc)
        njacx3(3,jg,jc) = ifneg(i + 2,jc)
        i = i + 3
      enddo
      irem = MOD(ifn(jc), 3)
      !$acc end kernels
      if (irem == 2) then
        njcgrp(jc,2) = 1
        njacx2(1,jc) = ifneg(i,jc)
        njacx2(2,jc) = ifneg(i + 1,jc)
        i = i + 2
      else
        if (irem == 1) then
          njcgrp(jc,1) = 1
          njacx1(jc) = ifneg(i,jc)
        end if
      end if
    enddo
    do jc = 1, jpcspf, 1
      nmpjac(jc) = ifp(jc)
      if (ifp(jc) > jppjac) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(3A,I0)') 'FATAL ASAD ERROR: array npjac1 too ', ' small in injac.f. Increase parameter jppjac to at ', ' least ', ifp(jc)
        call umprint(ummessage, src='asad_inijac')
        cmessage(:) = 'ARRAY NPJAC1 TOO SMALL: Increase jppjac'
        errcode = jc
        call ereport('ASAD_INIJAC', errcode, cmessage)
      end if
      !$acc kernels
      !$acc loop gang vector independent
      do jr = 1, ifp(jc), 1
        npjac1(jr,jc) = ifpos(jr,jc)
      enddo
      !$acc end kernels
    enddo
    !$acc kernels
    !$acc loop seq
    do jc = 1, jpcspf, 1
      nmzjac(jc) = ifz(jc)
      !$acc loop gang vector independent
      do jr = 1, ifz(jc), 1
        nzjac1(jr,jc) = ifzer(jr,jc)
      enddo
    enddo
    !$acc loop seq
    do jc = 1, nstst, 1
      nmsjac(jc) = ifs(jc)
      !$acc loop gang vector independent
      do jr = 1, ifs(jc), 1
        nsjac1(jr,jc) = ifsss(jr,jc)
      enddo
    enddo
    ntro3 = 0
    ntroh = 0
    ntrho2 = 0
    nspo3 = 0
    nspo1d = 0
    nspo3p = 0
    nspoh = 0
    nspho2 = 0
    nspno = 0
    nspn = 0
    nsph = 0
    !$acc loop seq
    do jtr = 1, jpcspf, 1
      if (specf(jtr) == 'O3        ') then
        ntro3 = jtr
      end if
      if (specf(jtr) == 'OH        ') then
        ntroh = jtr
      end if
      if (specf(jtr) == 'HO2       ') then
        ntrho2 = jtr
      end if
      if (specf(jtr) == 'NO        ') then
        ntrno = jtr
      end if
    enddo
    !$acc loop seq
    do jtr = 1, jpspec, 1
      if (speci(jtr) == 'O3        ') then
        nspo3 = jtr
      end if
      if (speci(jtr) == 'O(1D)     ') then
        nspo1d = jtr
      end if
      if (speci(jtr) == 'O(3P)     ') then
        nspo3p = jtr
      end if
      if (speci(jtr) == 'OH        ') then
        nspoh = jtr
      end if
      if (speci(jtr) == 'HO2       ') then
        nspho2 = jtr
      end if
      if (speci(jtr) == 'NO        ') then
        nspno = jtr
      end if
      if (speci(jtr) == 'N         ') then
        nspn = jtr
      end if
      if (speci(jtr) == 'H         ') then
        nsph = jtr
      end if
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_inijac

end module asad_inijac_mod
