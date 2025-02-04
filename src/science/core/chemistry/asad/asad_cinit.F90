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
    f = 0.0
    fdot = 0.0
    ej = 0.0
    linfam = .false.
    y = 0.0
    ydot = 0.0
    prod = 0.0
    slos = 0.0
    dpd = 0.0
    dpw = 0.0
    emr = 0.0
    rk = 0.0
    prk = 0.0
    nspi = 0
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
    nbrkx = 0
    ntrkx = 0
    nprkx = 0
    nhrkx = 0
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(njacx3, dim=3), UBOUND(njacx3, dim=3), 1
      do idx_1 = LBOUND(njacx3, dim=2), UBOUND(njacx3, dim=2), 1
        njacx3(1,idx_1,idx) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_2 = LBOUND(njacx3, dim=3), UBOUND(njacx3, dim=3), 1
      do idx_3 = LBOUND(njacx3, dim=2), UBOUND(njacx3, dim=2), 1
        njacx3(2,idx_3,idx_2) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_4 = LBOUND(njacx3, dim=3), UBOUND(njacx3, dim=3), 1
      do idx_5 = LBOUND(njacx3, dim=2), UBOUND(njacx3, dim=2), 1
        njacx3(3,idx_5,idx_4) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_6 = LBOUND(njcgrp, dim=1), UBOUND(njcgrp, dim=1), 1
      njcgrp(idx_6,1) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_7 = LBOUND(njcgrp, dim=1), UBOUND(njcgrp, dim=1), 1
      njcgrp(idx_7,2) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_8 = LBOUND(njcgrp, dim=1), UBOUND(njcgrp, dim=1), 1
      njcgrp(idx_8,3) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_9 = LBOUND(njacx2, dim=2), UBOUND(njacx2, dim=2), 1
      njacx2(1,idx_9) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_10 = LBOUND(njacx2, dim=2), UBOUND(njacx2, dim=2), 1
      njacx2(2,idx_10) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_11 = LBOUND(njacx1, dim=1), UBOUND(njacx1, dim=1), 1
      njacx1(idx_11) = 0
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
    do idx_12 = LBOUND(frpb, dim=1), UBOUND(frpb, dim=1), 1
      frpb(idx_12) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_13 = LBOUND(frpt, dim=1), UBOUND(frpt, dim=1), 1
      frpt(idx_13) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_14 = LBOUND(frpj, dim=1), UBOUND(frpj, dim=1), 1
      frpj(idx_14) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_15 = LBOUND(frph, dim=1), UBOUND(frph, dim=1), 1
      frph(idx_15) = 0.0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_16 = LBOUND(frpx, dim=1), UBOUND(frpx, dim=1), 1
      frpx(idx_16) = 0.0
    enddo
    !$acc end kernels
    nfrpx = 0
    !$acc kernels
    !$acc loop gang vector independent
    do idx_17 = 1, jpfrpx, 1
      ntabfp(idx_17,1) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_18 = 1, jpfrpx, 1
      ntabfp(idx_18,2) = 0
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent
    do idx_19 = 1, jpfrpx, 1
      ntabfp(idx_19,3) = 0
    enddo
    !$acc end kernels
    nmzjac = 0
    nmsjac = 0
    nzjac1 = 0
    nsjac1 = 0
    ntabpd = 0
    ztabpd = 0.0
    npdfr = 0
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
        if (ctype(j) == jpfm .OR. ctype(j) == jpif) then
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
