module asad_inrats_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_INRATS_MOD'
  integer, save, public :: iadv
  integer, save, public :: ispf
  public

  contains
  subroutine asad_inrats_set_sp_lists()
    use asad_mod, only : advt, ctype, family, h_in_ss, jpctr, jpfm, jpif, jpna, jpoo, jpro2, jpsp, jpspec, n_in_ss, nadvt, nlfro2, nlnaro2, nltr3, nltrf, nnaf, nodd, nro2, ntr3, ntrf, o1d_in_ss, o3p_in_ss, specf, speci, spro2
    use ukca_chem_defs_mod, only : chch_defs
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: errcode
    integer :: inadv
    integer :: iro2
    integer :: jadv
    integer :: js
    integer :: k
    CHARACTER(LEN = errormessagelength) :: cmessage
    logical :: l_fa
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INRATS_SET_SP_LISTS'
    integer :: idx
    integer :: idx_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    o1d_in_ss = .false.
    o3p_in_ss = .false.
    n_in_ss = .false.
    h_in_ss = .false.
    !$acc loop seq
    do k = 1, jpspec, 1
      if (chch_defs(k)%speci == 'O(1D)     ' .AND. chch_defs(k)%ctype(1:2) == jpna) then
        o1d_in_ss = .true.
      end if
      if (chch_defs(k)%speci == 'O(3P)     ' .AND. chch_defs(k)%ctype(1:2) == jpna) then
        o3p_in_ss = .true.
      end if
      if (chch_defs(k)%speci == 'N         ' .AND. chch_defs(k)%ctype(1:2) == jpna) then
        n_in_ss = .true.
      end if
      if (chch_defs(k)%speci == 'H         ' .AND. chch_defs(k)%ctype(1:2) == jpna) then
        h_in_ss = .true.
      end if
    enddo
    !$acc end kernels
    if (.NOT.o1d_in_ss .AND. (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR. ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat)) then
      cmessage(:) = ' O(1D) is not a Steady-State species'
      errcode = 1
      call ereport(routinename, errcode, cmessage)
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(nltrf, dim=1), UBOUND(nltrf, dim=1), 1
      nltrf(idx) = 0
    enddo
    !$acc loop gang vector independent
    do idx_1 = LBOUND(nltr3, dim=1), UBOUND(nltr3, dim=1), 1
      nltr3(idx_1) = 0
    enddo
    !$acc end kernels
    if (SIZE(chch_defs) /= jpspec) then
      errcode = 1
      cmessage(:) = ' jpspec and chch_defs are inconsistent'
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') cmessage
      call umprint(ummessage, src=routinename)
      call ereport(routinename, errcode, cmessage)
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do k = 1, jpspec, 1
      speci(k) = chch_defs(k)%speci
      nodd(k) = chch_defs(k)%nodd
      ctype(k) = chch_defs(k)%ctype(1:2)
      family(k) = chch_defs(k)%family
    enddo
    iadv = 0
    inadv = 0
    ntrf = 0
    nnaf = 0
    ntr3 = 0
    iro2 = 0
    nro2 = 0
    ispf = 0
    !$acc end kernels
    do js = 1, jpspec, 1
      if (ctype(js) /= jpfm .AND. ctype(js) /= jpif) then
        family(js) = '          '
      end if
      if (ctype(js) == jpsp) then
        iadv = iadv + 1
        ntrf = ntrf + 1
        ispf = ispf + 1
        if (iadv > jpctr) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') '** ASAD ERROR in subroutine ' // RoutineName
          call umprint(ummessage, src=routinename)
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,I3)') '** Parameter jpctr is too low; found', iadv
          call umprint(ummessage, src=routinename)
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,I3)') '** tracers so far with ', jpspec - js
          call umprint(ummessage, src=routinename)
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') '** species to check.'
          call umprint(ummessage, src=routinename)
          cmessage(:) = 'ASAD ERROR: jpctr is too low'
          call ereport(routinename, iadv, cmessage)
        end if
        advt(iadv) = speci(js)
        nltrf(ntrf) = ispf
        specf(ntrf) = speci(js)
      else
        if (ctype(js) == jpoo) then
          iro2 = iro2 + 1
          nro2 = nro2 + 1
          spro2(iro2) = speci(js)
          if (iro2 > jpro2) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') '** ASAD ERROR in subroutine ' // RoutineName
            call umprint(ummessage, src=routinename)
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A,I3)') '** Parameter jpro2 is too low; found', iro2
            call umprint(ummessage, src=routinename)
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A,I3)') '** tracers so far with ', jpspec - js
            call umprint(ummessage, src=routinename)
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') '** species to check.'
            call umprint(ummessage, src=routinename)
            cmessage(:) = 'ASAD ERROR: jpro2 is too low'
            call ereport(routinename, iadv, cmessage)
          end if
          ntrf = ntrf + 1
          ispf = ispf + 1
          nltrf(ntrf) = ispf
          specf(ntrf) = speci(js)
          nlfro2(iro2) = ispf
          if (ukca_config%l_ukca_ro2_ntp) then
            inadv = inadv + 1
            nnaf = nnaf + 1
            nadvt(inadv) = speci(js)
            nlnaro2(iro2) = inadv
          else
            iadv = iadv + 1
            advt(iadv) = speci(js)
          end if
        else
          if (ctype(js) == jpna) then
            inadv = inadv + 1
            nnaf = nnaf + 1
            nadvt(inadv) = speci(js)
          else
            if (ctype(js) == jpfm .OR. ctype(js) == jpif) then
              cmessage(:) = ' Family chemistry not available in this version'
              errcode = js
              call ereport('ASAD_INRATS', errcode, cmessage)
              if (ctype(js) == jpif) then
                iadv = iadv + 1
                ntr3 = ntr3 + 1
                if (iadv > jpctr) then
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') '** ASAD ERROR in subroutine ' // RoutineName
                  call umprint(ummessage, src=routinename)
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,I3)') '** Parameter jpctr is too low; found', iadv
                  call umprint(ummessage, src=routinename)
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,I3)') '** tracers so far with ', jpspec - js
                  call umprint(ummessage, src=routinename)
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') '** species to check.'
                  call umprint(ummessage, src=routinename)
                  cmessage(:) = 'ERROR in jpctr'
                  errcode = js
                  call ereport(routinename, errcode, cmessage)
                end if
                advt(iadv) = speci(js)
                nltr3(ntr3) = iadv
              end if
              l_fa = .true.
              do jadv = 1, iadv, 1
                if (family(js) == advt(jadv)) then
                  l_fa = .false.
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Exit_Stmt
                  EXIT
                end if
              enddo
              if (l_fa) then
                iadv = iadv + 1
                ntrf = ntrf + 1
                ispf = ispf + 1
                if (iadv > jpctr) then
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') '***** ASAD ERROR in subroutine ' // RoutineName
                  call umprint(ummessage, src=routinename)
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,I3)') '***** Param jpctr is too low; found', iadv
                  call umprint(ummessage, src=routinename)
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,I3)') '***** tracers so far with ', jpspec - js
                  call umprint(ummessage, src=routinename)
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') '***** species to check.'
                  call umprint(ummessage, src=routinename)
                  cmessage(:) = 'INRATS ERROR : jpctr is too low'
                  call ereport(routinename, iadv, cmessage)
                end if
                advt(iadv) = family(js)
                nltrf(ntrf) = ispf
              end if
            end if
          end if
        end if
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_inrats_set_sp_lists
  subroutine asad_inrats()
    use asad_mod, only : ab, advt, at, ctype, family, frpb, frph, frpj, frpt, frpx, jpbk, jpcspf, jpctr, jpfrpb, jpfrph, jpfrpj, jpfrpt, jphk, jpif, jpnr, jppj, jpspb, jpspec, jpsph, jpspj, jpspt, jptk, madvtr, majors, moffam, nadvt, nbrkx, nfrpx, nhrkx, nlmajmin, nnaf, nprkx, nspi, ntrf, ntrkx, nuni, spb, specf, speci, sph, spj, spt
    use ukca_chem_defs_mod, only : ratb_defs, rath_defs, ratj_defs, ratt_defs
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use ukca_um_legacy_mod, only : mype
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: errcode
    integer, dimension(jpbk + 1,jpspb) :: ispb
    integer, dimension(jptk + 1,jpspt) :: ispt
    integer, dimension(jppj + 1,jpspj) :: ispj
    integer, dimension(jphk + 1,jpsph) :: isph
    integer, dimension(jpbk + 1) :: ifrpbx
    integer, dimension(jppj + 1) :: ifrpjx
    integer, dimension(jptk + 1) :: ifrptx
    integer, dimension(jphk + 1) :: ifrphx
    integer, dimension(jpspec) :: ilmin
    integer, dimension(jpspec) :: ilmaj
    integer :: imajor
    integer :: iminor
    integer :: ix
    integer :: icount
    integer :: ifam
    integer :: imaj
    integer :: iflag
    integer :: j
    integer :: jf
    integer :: jb
    integer :: jctr
    integer :: jcspf
    integer :: jh
    integer :: jj
    integer :: jp
    integer :: jr
    integer :: js
    integer :: jspb
    integer :: jsph
    integer :: jspj
    integer :: jspt
    integer :: jt
    integer :: k
    integer :: ind
    CHARACTER(LEN = 10), PARAMETER :: nullx = '          '
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INRATS'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(ispb, dim=2), UBOUND(ispb, dim=2), 1
      do idx_1 = LBOUND(ispb, dim=1), UBOUND(ispb, dim=1), 1
        ispb(idx_1,idx) = 0
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_2 = LBOUND(ifrpbx, dim=1), UBOUND(ifrpbx, dim=1), 1
      ifrpbx(idx_2) = 0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_3 = LBOUND(ispt, dim=2), UBOUND(ispt, dim=2), 1
      do idx_4 = LBOUND(ispt, dim=1), UBOUND(ispt, dim=1), 1
        ispt(idx_4,idx_3) = 0
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_5 = LBOUND(ifrptx, dim=1), UBOUND(ifrptx, dim=1), 1
      ifrptx(idx_5) = 0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_6 = LBOUND(ispj, dim=2), UBOUND(ispj, dim=2), 1
      do idx_7 = LBOUND(ispj, dim=1), UBOUND(ispj, dim=1), 1
        ispj(idx_7,idx_6) = 0
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_8 = LBOUND(ifrpjx, dim=1), UBOUND(ifrpjx, dim=1), 1
      ifrpjx(idx_8) = 0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_9 = LBOUND(isph, dim=2), UBOUND(isph, dim=2), 1
      do idx_10 = LBOUND(isph, dim=1), UBOUND(isph, dim=1), 1
        isph(idx_10,idx_9) = 0
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_11 = LBOUND(ifrphx, dim=1), UBOUND(ifrphx, dim=1), 1
      ifrphx(idx_11) = 0
    enddo
    !$acc loop gang vector independent collapse(2)
    do jcspf = 1, ntrf, 1
      do js = 1, jpspec, 1
        if (family(js) == specf(jcspf) .AND. ctype(js) /= jpif) then
          majors(jcspf) = js
        end if
        if (speci(js) == specf(jcspf)) then
          majors(jcspf) = js
        end if
      enddo
    enddo
    !$acc end kernels
    do js = 1, jpspec, 1
      moffam(js) = 0
      madvtr(js) = 0
      do jcspf = 1, ntrf, 1
        if (family(js) == specf(jcspf)) then
          moffam(js) = jcspf
        end if
        if (speci(js) == specf(jcspf)) then
          madvtr(js) = jcspf
        end if
        if (family(js) == specf(jcspf) .AND. js > majors(jcspf)) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') '** ASAD ERROR: '
          call umprint(ummessage, src='asad_inrats')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') 'RE-ORDER SPECIES FILE SO THAT THE MAJOR '
          call umprint(ummessage, src='asad_inrats')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A)') 'SPECIES OF A FAMILY OCCURS AFTER THE OTHERS'
          call umprint(ummessage, src='asad_inrats')
          cmessage(:) = 'INRATS ERROR : Order of species is incorrect'
          errcode = jcspf
          call ereport('ASAD_INRATS', errcode, cmessage)
        end if
      enddo
    enddo
    !$acc kernels
    nlmajmin(1) = 5
    imajor = 0
    iminor = 0
    !$acc loop seq
    do js = 1, jpspec, 1
      ifam = moffam(js)
      imaj = 0
      if (ifam /= 0) then
        imaj = majors(ifam)
        if (imaj /= js) then
          iminor = iminor + 1
          ilmin(iminor) = js
        else
          imajor = imajor + 1
          ilmaj(imajor) = js
        end if
      end if
    enddo
    nlmajmin(2) = nlmajmin(1) + imajor - 1
    nlmajmin(3) = nlmajmin(1) + imajor
    nlmajmin(4) = nlmajmin(3) + iminor - 1
    !$acc loop seq
    do j = nlmajmin(1), nlmajmin(2), 1
      nlmajmin(j) = ilmaj(j - nlmajmin(1) + 1)
    enddo
    !$acc loop seq
    do j = nlmajmin(3), nlmajmin(4), 1
      nlmajmin(j) = ilmin(j - nlmajmin(3) + 1)
    enddo
    !$acc end kernels
    if (iadv /= jpctr) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(2A)') '** ASAD ERROR: Number of advected tracers', ' specified in chch_defs does not match jpctr'
      call umprint(ummessage, src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I3,A,I3)') 'Found ', iadv, ' but expected ', jpctr
      call umprint(ummessage, src='asad_inrats')
      cmessage(:) = 'INRATS ERROR : iadv and jpctr do not match'
      call ereport('ASAD_INRATS', iadv, cmessage)
    end if
    if (ispf /= jpcspf) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(2A)') '** ASAD ERROR: Number of species treated by asad', ' specified in chch_defs does not match jpcspf'
      call umprint(ummessage, src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I3,A,I3)') 'Found ', ispf, ' but expected ', jpcspf
      call umprint(ummessage, src='asad_inrats')
      cmessage(:) = 'INRATS ERROR : ispf and jpcspf do not match'
      call ereport('ASAD_INRATS', ispf, cmessage)
    end if
    if (mype == 0 .AND. printstatus >= prstatus_oper) then
      call umprint('', src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') '  ***  CHEMISTRY INFORMATION  ***'
      call umprint(ummessage, src='asad_inrats')
      call umprint('', src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'ASAD IS TREATING ADVECTED TRACERS IN THE ORDER:'
      call umprint(ummessage, src='asad_inrats')
      call umprint('', src='asad_inrats')
      do jctr = 1, jpctr, 1
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I3,A,A10)') '  ', jctr, ' ', advt(jctr)
        call umprint(ummessage, src='asad_inrats')
      enddo
      call umprint('', src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'ASAD IS TREATING NON-ADVECTED TRACERS IN THE ORDER:'
      call umprint(ummessage, src='asad_inrats')
      call umprint('', src='asad_inrats')
      do jctr = 1, nnaf, 1
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I2,A,A10)') '  ', jctr, ' ', nadvt(jctr)
        call umprint(ummessage, src='asad_inrats')
      enddo
      call umprint('', src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'ASAD IS TREATING CHEMICALLY ACTIVE SPECIES IN THE ORDER:'
      call umprint(ummessage, src='asad_inrats')
      call umprint('', src='asad_inrats')
      do jcspf = 1, ntrf, 1
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I2,A,A10)') '  ', jcspf, ' ', specf(jcspf)
        call umprint(ummessage, src='asad_inrats')
      enddo
      call umprint('', src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'IF THE TRACERS WERE NOT INITIALISED IN THIS '
      call umprint(ummessage, src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'ORDER THEN THE MODEL RESULTS ARE WORTHLESS '
      call umprint(ummessage, src='asad_inrats')
      call umprint('', src='asad_inrats')
      iflag = 0
      do jcspf = 1, jpcspf, 1
        if (specf(jcspf) /= speci(majors(jcspf))) then
          if (iflag == 0) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') 'THE MAJOR MEMBER OF EACH OF THE FAMILIES'
            call umprint(ummessage, src='asad_inrats')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') 'IS GIVEN BELOW. IF THIS IS NOT ACCEPTABLE,'
            call umprint(ummessage, src='asad_inrats')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') 'THEN YOU MUST REORDER THE SPECIES IN' // ' chch_defs'
            call umprint(ummessage, src='asad_inrats')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)') 'SO THE MAJOR SPECIES FOLLOWS THE OTHERS.'
            call umprint(ummessage, src='asad_inrats')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A)')
            call umprint(ummessage, src='asad_inrats')
            iflag = 1
          end if
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A10,1X,A10)') specf(jcspf), speci(majors(jcspf))
          call umprint(ummessage, src='asad_inrats')
        end if
      enddo
    end if
    if (SIZE(ratb_defs) /= jpbk) then
      errcode = 1
      cmessage(:) = 'size of ratb_defs is inconsistent with jpbk'
      call ereport('ASAD_INRATS', errcode, cmessage)
    end if
    !$acc kernels
    icount = 1
    !$acc loop seq
    do k = 1, jpbk, 1
      spb(k,1) = ratb_defs(k)%react1
      spb(k,2) = ratb_defs(k)%react2
      spb(k,3) = ratb_defs(k)%prod1
      spb(k,4) = ratb_defs(k)%prod2
      spb(k,5) = ratb_defs(k)%prod3
      spb(k,6) = ratb_defs(k)%prod4
      ab(k,1) = ratb_defs(k)%k0
      ab(k,2) = ratb_defs(k)%alpha
      ab(k,3) = ratb_defs(k)%beta
      if (ratb_defs(k)%pyield1 > 1e-18) then
        ifrpbx(k) = icount
        frpb(icount) = ratb_defs(k)%pyield1
        if (spb(k,4) /= nullx(:)) then
          frpb(icount + 1) = ratb_defs(k)%pyield2
        end if
        if (spb(k,5) /= nullx(:)) then
          frpb(icount + 2) = ratb_defs(k)%pyield3
        end if
        if (spb(k,6) /= nullx(:)) then
          frpb(icount + 3) = ratb_defs(k)%pyield4
        end if
        icount = icount + 4
      end if
    enddo
    !$acc loop gang vector independent collapse(3)
    do jb = 1, jpbk, 1
      do js = 1, jpspec, 1
        do jspb = 1, jpspb, 1
          if (speci(js) == spb(jb,jspb)) then
            ispb(jb,jspb) = js
          end if
        enddo
      enddo
    enddo
    !$acc loop gang vector independent
    do jf = 1, jpfrpb, 1
      frpx(jf) = frpb(jf)
    enddo
    !$acc end kernels
    if (SIZE(ratt_defs) /= jptk) then
      errcode = 1
      cmessage(:) = 'size of ratt_defs is inconsistent with jptk'
      call ereport('ASAD_INRATS', errcode, cmessage)
    end if
    !$acc kernels
    icount = 1
    !$acc loop seq
    do k = 1, jptk, 1
      spt(k,1) = ratt_defs(k)%react1
      spt(k,2) = ratt_defs(k)%react2
      spt(k,3) = ratt_defs(k)%prod1
      spt(k,4) = ratt_defs(k)%prod2
      at(k,1) = ratt_defs(k)%f
      at(k,2) = ratt_defs(k)%k1
      at(k,3) = ratt_defs(k)%alpha1
      at(k,4) = ratt_defs(k)%beta1
      at(k,5) = ratt_defs(k)%k2
      at(k,6) = ratt_defs(k)%alpha2
      at(k,7) = ratt_defs(k)%beta2
      if (ratt_defs(k)%pyield1 > 1e-18) then
        ifrptx(k) = icount
        frpt(icount) = ratt_defs(k)%pyield1
        if (spt(k,4) /= nullx(:)) then
          frpt(icount + 1) = ratt_defs(k)%pyield2
        end if
        icount = icount + 2
      end if
    enddo
    !$acc loop gang vector independent collapse(3)
    do jt = 1, jptk, 1
      do js = 1, jpspec, 1
        do jspt = 1, jpspt, 1
          if (speci(js) == spt(jt,jspt)) then
            ispt(jt,jspt) = js
          end if
        enddo
      enddo
    enddo
    !$acc loop seq
    do jf = 1, jpfrpt, 1
      ind = jf + jpfrpb
      frpx(ind) = frpt(jf)
    enddo
    !$acc end kernels
    if (SIZE(ratj_defs) /= jppj) then
      errcode = 1
      cmessage(:) = 'size of ratj_defs is not equal to jppj'
      call ereport('ASAD_INRATS', errcode, cmessage)
    end if
    !$acc kernels
    icount = 1
    !$acc loop seq
    do k = 1, jppj, 1
      spj(k,1) = ratj_defs(k)%react1
      spj(k,2) = ratj_defs(k)%react2
      spj(k,3) = ratj_defs(k)%prod1
      spj(k,4) = ratj_defs(k)%prod2
      spj(k,5) = ratj_defs(k)%prod3
      spj(k,6) = ratj_defs(k)%prod4
      if (ratj_defs(k)%pyield1 > 1e-18) then
        ifrpjx(k) = icount
        frpj(icount) = ratj_defs(k)%pyield1
        if (spj(k,4) /= nullx(:)) then
          frpj(icount + 1) = ratj_defs(k)%pyield2
        end if
        if (spj(k,5) /= nullx(:)) then
          frpj(icount + 2) = ratj_defs(k)%pyield3
        end if
        if (spj(k,6) /= nullx(:)) then
          frpj(icount + 3) = ratj_defs(k)%pyield4
        end if
        icount = icount + 4
      end if
    enddo
    !$acc loop gang vector independent collapse(3)
    do jj = 1, jppj, 1
      do js = 1, jpspec, 1
        do jspj = 1, jpspj, 1
          if (speci(js) == spj(jj,jspj)) then
            ispj(jj,jspj) = js
          end if
        enddo
      enddo
    enddo
    !$acc loop seq
    do jf = 1, jpfrpj, 1
      ind = jf + jpfrpb + jpfrpt
      frpx(ind) = frpj(jf)
    enddo
    !$acc end kernels
    if (jphk > 0) then
      if (SIZE(rath_defs) /= jphk) then
        errcode = 1
        cmessage(:) = 'size of rath_defs is not equal to jphk'
        call ereport('ASAD_INRATS', errcode, cmessage)
      end if
      !$acc kernels
      icount = 1
      !$acc loop seq
      do k = 1, jphk, 1
        sph(k,1) = rath_defs(k)%react1
        sph(k,2) = rath_defs(k)%react2
        sph(k,3) = rath_defs(k)%prod1
        sph(k,4) = rath_defs(k)%prod2
        sph(k,5) = rath_defs(k)%prod3
        sph(k,6) = rath_defs(k)%prod4
        if (rath_defs(k)%pyield1 > 1e-18) then
          ifrphx(k) = icount
          frph(icount) = rath_defs(k)%pyield1
          if (sph(k,4) /= nullx(:)) then
            frph(icount + 1) = rath_defs(k)%pyield2
          end if
          if (sph(k,5) /= nullx(:)) then
            frph(icount + 2) = rath_defs(k)%pyield3
          end if
          if (sph(k,6) /= nullx(:)) then
            frph(icount + 3) = rath_defs(k)%pyield4
          end if
          icount = icount + 4
        end if
      enddo
      !$acc loop gang vector independent collapse(3)
      do jh = 1, jphk, 1
        do js = 1, jpspec, 1
          do jsph = 1, jpsph, 1
            if (speci(js) == sph(jh,jsph)) then
              isph(jh,jsph) = js
            end if
          enddo
        enddo
      enddo
      !$acc loop seq
      do jf = 1, jpfrph, 1
        ind = jf + jpfrpb + jpfrpt + jpfrpj
        frpx(ind) = frph(jf)
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    nuni = 0
    !$acc loop seq
    do jr = 1, jpbk, 1
      if (ispb(jr,2) == 0) then
        nuni = nuni + 1
        nbrkx(jr) = nuni
        !$acc loop gang vector independent
        do jp = 1, jpspb, 1
          nspi(nuni,jp) = ispb(jr,jp)
        enddo
        if (ifrpbx(jr) /= 0) then
          nfrpx(nuni) = ifrpbx(jr)
        end if
      end if
    enddo
    !$acc loop seq
    do jr = 1, jptk, 1
      if (ispt(jr,2) == 0) then
        nuni = nuni + 1
        ntrkx(jr) = nuni
        !$acc loop gang vector independent
        do jp = 1, jpspt, 1
          nspi(nuni,jp) = ispt(jr,jp)
        enddo
        if (ifrptx(jr) /= 0) then
          nfrpx(nuni) = ifrptx(jr) + jpfrpb
        end if
      end if
    enddo
    !$acc loop seq
    do jr = 1, jppj, 1
      if (ispj(jr,2) == 0) then
        nuni = nuni + 1
        nprkx(jr) = nuni
        !$acc loop gang vector independent
        do jp = 1, jpspj, 1
          nspi(nuni,jp) = ispj(jr,jp)
        enddo
        if (ifrpjx(jr) /= 0) then
          nfrpx(nuni) = ifrpjx(jr) + jpfrpb + jpfrpt
        end if
      end if
    enddo
    !$acc end kernels
    if (jphk > 0) then
      !$acc kernels
      !$acc loop seq
      do jr = 1, jphk, 1
        if (isph(jr,2) == 0) then
          nuni = nuni + 1
          nhrkx(jr) = nuni
          !$acc loop gang vector independent
          do jp = 1, jpsph, 1
            nspi(nuni,jp) = isph(jr,jp)
          enddo
          if (ifrphx(jr) /= 0) then
            nfrpx(nuni) = ifrphx(jr) + jpfrpb + jpfrpt + jpfrpj
          end if
        end if
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    ix = nuni
    !$acc loop seq
    do jr = 1, jpbk, 1
      if (ispb(jr,2) /= 0) then
        ix = ix + 1
        nbrkx(jr) = ix
        !$acc loop gang vector independent
        do jp = 1, jpspb, 1
          nspi(ix,jp) = ispb(jr,jp)
        enddo
        if (ifrpbx(jr) /= 0) then
          nfrpx(ix) = ifrpbx(jr)
        end if
      end if
    enddo
    !$acc loop seq
    do jr = 1, jptk, 1
      if (ispt(jr,2) /= 0) then
        ix = ix + 1
        ntrkx(jr) = ix
        !$acc loop gang vector independent
        do jp = 1, jpspt, 1
          nspi(ix,jp) = ispt(jr,jp)
        enddo
        if (ifrptx(jr) /= 0) then
          nfrpx(ix) = ifrptx(jr) + jpfrpb
        end if
      end if
    enddo
    !$acc loop seq
    do jr = 1, jppj, 1
      if (ispj(jr,2) /= 0) then
        ix = ix + 1
        nprkx(jr) = ix
        !$acc loop gang vector independent
        do jp = 1, jpspj, 1
          nspi(ix,jp) = ispj(jr,jp)
        enddo
        if (ifrpjx(jr) /= 0) then
          nfrpx(ix) = ifrpjx(jr) + jpfrpb + jpfrpt
        end if
      end if
    enddo
    !$acc end kernels
    if (jphk > 0) then
      !$acc kernels
      !$acc loop seq
      do jr = 1, jphk, 1
        if (isph(jr,2) /= 0) then
          ix = ix + 1
          nhrkx(jr) = ix
          !$acc loop gang vector independent
          do jp = 1, jpsph, 1
            nspi(ix,jp) = isph(jr,jp)
          enddo
          if (ifrphx(jr) /= 0) then
            nfrpx(ix) = ifrphx(jr) + jpfrpb + jpfrpt + jpfrpj
          end if
        end if
      enddo
      !$acc end kernels
    end if
    if (ix /= jpnr) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(2A)') '*** INTERNAL ASAD ERROR: Number of reactions', ' placed in nspi array does not equal jpnr. '
      call umprint(ummessage, src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(2A)') '                         Check that reaction', ' files and value of jpnr in UKCA namelist are'
      call umprint(ummessage, src='asad_inrats')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(2A)') '                         consistent. Found: ', ix, ' jpnr: ', jpnr
      call umprint(ummessage, src='asad_inrats')
      cmessage(:) = 'No of rxns in nspi array is not equal to jpnr'
      call ereport('ASAD_INRATS', ix, cmessage)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_inrats

end module asad_inrats_mod
