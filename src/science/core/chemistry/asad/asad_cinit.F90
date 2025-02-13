module asad_cinit_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_CINIT_MOD'
  public

  contains
  subroutine asad_cinit()
    use asad_mod, only : cdt, ctype, dpd, dpw, ej, emr, f, fdot, frpb, frph, frpj, frpt, frpx, jpcspf, jpfm, jpfrpx, jpif, jpnr, jppjac, jpspec, kfphot, ldepd, ldepw, linfam, ljacx, method, nbrkx, ncsteps, ndepd, ndepw, nfphot, nfrpx, ngrp, nhrkx, nit0, nitfg, nitnr, njacx1, njacx2, njacx3, njcgrp, nlall, nldepd, nldepw, nldepx, nlf, nlmajmin, nlstst, nmpjac, nmsjac, nmzjac, npdfr, npjac1, nprdx1, nprdx2, nprdx3, nprkx, nrsteps_max, nsjac1, nspi, ntabfp, ntabpd, ntrkx, nzjac1, peps, prk, prod, rk, slos, spb, sph, spj, spt, y, ydot, ztabpd
    use ukca_config_specification_mod, only : int_method_be_explicit, int_method_nr, ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    use asad_inijac_mod, only : asad_inijac
    use asad_inix_mod, only : asad_inix
    use asad_inrats_mod, only : asad_inrats
    use asad_setsteady_mod, only : asad_setsteady
    use ukca_inddep_mod, only : ukca_inddep
    use ukca_inwdep_mod, only : ukca_inwdep
    integer, parameter :: nit0_max = 50
    integer, parameter :: nitfg_max = 50
    integer, parameter :: nitnr_max = 50
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: j
    integer :: jc
    integer :: jp
    integer :: js
    integer :: jx
    integer :: jpnpx3
    integer :: errcode
    real :: sfmin
    CHARACTER(LEN = 10), PARAMETER :: nullx = '          '
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_CINIT'
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
    integer :: idx_51
    integer :: idx_52
    integer :: idx_53
    integer :: idx_54
    integer :: idx_55
    integer :: idx_56
    integer :: idx_57
    integer :: idx_58
    integer :: idx_59
    integer :: idx_60
    integer :: idx_61
    integer :: idx_62
    integer :: idx_63
    integer :: idx_64

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    ljacx = .true.
    jpnpx3 = jpnr / (3 * 3) + 3 * 3
    if (method /= int_method_be_explicit) then
      if (ukca_config%nrsteps < 0 .OR. ukca_config%nrsteps > nrsteps_max) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(cmessage, '(2(A,1x,I0))') 'NRSTEPS is not in range - 0:', nrsteps_max, '. NRSTEPS=', ukca_config % nrsteps
        errcode = 1
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      if (nit0 < 0 .OR. nit0 > nit0_max) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(cmessage, '(2(A,1x,I0))') 'NIT0 is not in range - 0:', nit0_max, '. NIT0=', nit0
        errcode = 2
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      if (nitfg < 0 .OR. nitfg > nitfg_max) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(cmessage, '(2(A,1x,I0))') 'NITFG is not in range - 0:', nitfg_max, '. NITFG=', nitfg
        errcode = 3
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      if (nitnr < 0 .OR. nitnr > nitnr_max) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(cmessage, '(2(A,1x,I0))') 'NITNR is not in range - 0:', nitnr_max, '. NITNR=', nitnr
        errcode = 4
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    end if
    if (kfphot < 0 .AND. ABS(kfphot) > INT(ukca_config%timestep)) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I0,A,A)') '**CINIT WARNING: VALUE OF KFPHOT ', kfphot, 'EXCEEDS THE MODEL TIMESTEP. ROUTINE PHOTOL WILL ONLY', ' BE CALLED ONCE.'
      call umprint(ummessage, src=routinename)
      nfphot = 0
    else
      if (kfphot > 0 .AND. kfphot > ncsteps) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I0,A,A)') '**CINIT WARNING: FREQUENCY KFPHOT ', kfphot, ' EXCEEDS THE TOTAL NUMBER OF CHEMICAL SUBSTEPS. ROUTINE ', ' PHOTOL WILL BE CALLED ONCE ONLY.'
        call umprint(ummessage, src=routinename)
        nfphot = 0
      else
        if (kfphot < 0) then
          nfphot = INT(ABS(kfphot) / cdt)
        else
          nfphot = kfphot
        end if
      end if
    end if
    sfmin = TINY(1.0)
    sfmin = 10.0 ** (INT(LOG10(sfmin)) + 1)
    peps = 1.0e19 * sfmin
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
      do idx_1 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_1,idx) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_2 = LBOUND(fdot, dim=2), UBOUND(fdot, dim=2), 1
      do idx_3 = LBOUND(fdot, dim=1), UBOUND(fdot, dim=1), 1
        fdot(idx_3,idx_2) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_4 = LBOUND(ej, dim=2), UBOUND(ej, dim=2), 1
      do idx_5 = LBOUND(ej, dim=1), UBOUND(ej, dim=1), 1
        ej(idx_5,idx_4) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_6 = LBOUND(linfam, dim=2), UBOUND(linfam, dim=2), 1
      do idx_7 = LBOUND(linfam, dim=1), UBOUND(linfam, dim=1), 1
        linfam(idx_7,idx_6) = .false.
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_8 = LBOUND(y, dim=2), UBOUND(y, dim=2), 1
      do idx_9 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_9,idx_8) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_10 = LBOUND(ydot, dim=2), UBOUND(ydot, dim=2), 1
      do idx_11 = LBOUND(ydot, dim=1), UBOUND(ydot, dim=1), 1
        ydot(idx_11,idx_10) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_12 = LBOUND(prod, dim=2), UBOUND(prod, dim=2), 1
      do idx_13 = LBOUND(prod, dim=1), UBOUND(prod, dim=1), 1
        prod(idx_13,idx_12) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_14 = LBOUND(slos, dim=2), UBOUND(slos, dim=2), 1
      do idx_15 = LBOUND(slos, dim=1), UBOUND(slos, dim=1), 1
        slos(idx_15,idx_14) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_16 = LBOUND(dpd, dim=2), UBOUND(dpd, dim=2), 1
      do idx_17 = LBOUND(dpd, dim=1), UBOUND(dpd, dim=1), 1
        dpd(idx_17,idx_16) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_18 = LBOUND(dpw, dim=2), UBOUND(dpw, dim=2), 1
      do idx_19 = LBOUND(dpw, dim=1), UBOUND(dpw, dim=1), 1
        dpw(idx_19,idx_18) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_20 = LBOUND(emr, dim=2), UBOUND(emr, dim=2), 1
      do idx_21 = LBOUND(emr, dim=1), UBOUND(emr, dim=1), 1
        emr(idx_21,idx_20) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_22 = LBOUND(rk, dim=2), UBOUND(rk, dim=2), 1
      do idx_23 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
        rk(idx_23,idx_22) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_24 = LBOUND(prk, dim=2), UBOUND(prk, dim=2), 1
      do idx_25 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
        prk(idx_25,idx_24) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_26 = LBOUND(nspi, dim=2), UBOUND(nspi, dim=2), 1
      do idx_27 = LBOUND(nspi, dim=1), UBOUND(nspi, dim=1), 1
        nspi(idx_27,idx_26) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop seq
    do js = 1, jpspec, 1
      ngrp(js,1) = 0
      ngrp(js,2) = 0
      ngrp(js,3) = 0
      nprdx2(1,js) = 0
      nprdx2(2,js) = 0
      nprdx1(js) = 0
      ngrp(js + jpspec,1) = 0
      ngrp(js + jpspec,2) = 0
      ngrp(js + jpspec,3) = 0
      nprdx2(1,js + jpspec) = 0
      nprdx2(2,js + jpspec) = 0
      nprdx1(js + jpspec) = 0
      nlall(js) = 0
      nlstst(js) = 0
      nlf(js) = 0
      nlmajmin(js) = 0
      nldepd(js) = 0
      nldepw(js) = 0
      nldepx(js) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do js = 1, 2 * jpspec, 1
      do jx = 1, jpnpx3, 1
        nprdx3(1,jx,js) = 0
        nprdx3(2,jx,js) = 0
        nprdx3(3,jx,js) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_28 = LBOUND(nbrkx, dim=1), UBOUND(nbrkx, dim=1), 1
      nbrkx(idx_28) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_29 = LBOUND(ntrkx, dim=1), UBOUND(ntrkx, dim=1), 1
      ntrkx(idx_29) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_30 = LBOUND(nprkx, dim=1), UBOUND(nprkx, dim=1), 1
      nprkx(idx_30) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_31 = LBOUND(nhrkx, dim=1), UBOUND(nhrkx, dim=1), 1
      nhrkx(idx_31) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_32 = LBOUND(njacx3, dim=3), UBOUND(njacx3, dim=3), 1
      do idx_33 = LBOUND(njacx3, dim=2), UBOUND(njacx3, dim=2), 1
        njacx3(1,idx_33,idx_32) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_34 = LBOUND(njacx3, dim=3), UBOUND(njacx3, dim=3), 1
      do idx_35 = LBOUND(njacx3, dim=2), UBOUND(njacx3, dim=2), 1
        njacx3(2,idx_35,idx_34) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_36 = LBOUND(njacx3, dim=3), UBOUND(njacx3, dim=3), 1
      do idx_37 = LBOUND(njacx3, dim=2), UBOUND(njacx3, dim=2), 1
        njacx3(3,idx_37,idx_36) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_38 = LBOUND(njcgrp, dim=1), UBOUND(njcgrp, dim=1), 1
      njcgrp(idx_38,1) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_39 = LBOUND(njcgrp, dim=1), UBOUND(njcgrp, dim=1), 1
      njcgrp(idx_39,2) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_40 = LBOUND(njcgrp, dim=1), UBOUND(njcgrp, dim=1), 1
      njcgrp(idx_40,3) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_41 = LBOUND(njacx2, dim=2), UBOUND(njacx2, dim=2), 1
      njacx2(1,idx_41) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_42 = LBOUND(njacx2, dim=2), UBOUND(njacx2, dim=2), 1
      njacx2(2,idx_42) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_43 = LBOUND(njacx1, dim=1), UBOUND(njacx1, dim=1), 1
      njacx1(idx_43) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop seq
    do jc = 1, jpcspf, 1
      nmpjac(jc) = 0
      !$acc loop gang vector independent
      do jp = 1, jppjac, 1
        npjac1(jp,jc) = 0
      enddo
    enddo
    !$acc end kernels
    spb(:,:) = nullx(:)
    spt(:,:) = nullx(:)
    spj(:,:) = nullx(:)
    sph(:,:) = nullx(:)
    !$acc kernels
    !$acc loop gang vector independent
    do idx_44 = LBOUND(frpb, dim=1), UBOUND(frpb, dim=1), 1
      frpb(idx_44) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_45 = LBOUND(frpt, dim=1), UBOUND(frpt, dim=1), 1
      frpt(idx_45) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_46 = LBOUND(frpj, dim=1), UBOUND(frpj, dim=1), 1
      frpj(idx_46) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_47 = LBOUND(frph, dim=1), UBOUND(frph, dim=1), 1
      frph(idx_47) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_48 = LBOUND(frpx, dim=1), UBOUND(frpx, dim=1), 1
      frpx(idx_48) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_49 = LBOUND(nfrpx, dim=1), UBOUND(nfrpx, dim=1), 1
      nfrpx(idx_49) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_50 = 1, jpfrpx, 1
      ntabfp(idx_50,1) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_51 = 1, jpfrpx, 1
      ntabfp(idx_51,2) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_52 = 1, jpfrpx, 1
      ntabfp(idx_52,3) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_53 = LBOUND(nmzjac, dim=1), UBOUND(nmzjac, dim=1), 1
      nmzjac(idx_53) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_54 = LBOUND(nmsjac, dim=1), UBOUND(nmsjac, dim=1), 1
      nmsjac(idx_54) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_55 = LBOUND(nzjac1, dim=2), UBOUND(nzjac1, dim=2), 1
      do idx_56 = LBOUND(nzjac1, dim=1), UBOUND(nzjac1, dim=1), 1
        nzjac1(idx_56,idx_55) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_57 = LBOUND(nsjac1, dim=2), UBOUND(nsjac1, dim=2), 1
      do idx_58 = LBOUND(nsjac1, dim=1), UBOUND(nsjac1, dim=1), 1
        nsjac1(idx_58,idx_57) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_59 = LBOUND(ntabpd, dim=2), UBOUND(ntabpd, dim=2), 1
      do idx_60 = LBOUND(ntabpd, dim=1), UBOUND(ntabpd, dim=1), 1
        ntabpd(idx_60,idx_59) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_61 = LBOUND(ztabpd, dim=2), UBOUND(ztabpd, dim=2), 1
      do idx_62 = LBOUND(ztabpd, dim=1), UBOUND(ztabpd, dim=1), 1
        ztabpd(idx_62,idx_61) = 0.0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_63 = LBOUND(npdfr, dim=2), UBOUND(npdfr, dim=2), 1
      do idx_64 = LBOUND(npdfr, dim=1), UBOUND(npdfr, dim=1), 1
        npdfr(idx_64,idx_63) = 0
      enddo
    enddo
    !$acc end kernels
    call asad_inrats()
    do js = 1, jpspec, 1
      if (ldepd(js) .AND. ctype(js)(1 : 1) == 'C') then
        cmessage(:) = 'Dry deposition turned on for constant species'
        errcode = js
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
      if (ldepw(js) .AND. ctype(js)(1 : 1) == 'C') then
        cmessage(:) = 'Wet deposition turned on for constant species'
        errcode = js
        call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
      end if
    enddo
    if (method == int_method_nr) then
      call asad_setsteady()
    end if
    if (method >= 10) then
      do j = 1, jpspec, 1
        if (ctype(j) == jpfm(:) .OR. ctype(j) == jpif(:)) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,A,A,I0,A,I0)') '*** ASAD ERROR: You cannot use families ', ' with one of the stiff integrators. If method  >=  10 ', ' you cannot have species specified as ', jpfm, ' or ', jpif
          call umprint(ummessage, src=routinename)
          cmessage(:) = 'ASAD ABORTED'
          errcode = j
          call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
        end if
      enddo
    end if
    call asad_inix()
    call asad_inijac()
    if (ndepw /= 0) then
      call ukca_inwdep()
    end if
    if (ndepd /= 0) then
      call ukca_inddep()
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_cinit

end module asad_cinit_mod
