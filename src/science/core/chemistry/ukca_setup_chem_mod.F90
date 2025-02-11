module ukca_setup_chem_mod
  implicit none
  interface psyclone_internal_cmp
    procedure :: psyclone_cmp_int, psyclone_cmp_logical, psyclone_cmp_char
  end interface psyclone_internal_cmp
  character(len = *), parameter, private :: ModuleName = 'UKCA_SETUP_CHEM_MOD'
  public

  private :: psyclone_cmp_int, psyclone_cmp_logical, psyclone_cmp_char, psyclone_internal_cmp

  contains
  pure function psyclone_cmp_int(op1, op2)
    integer, intent(in) :: op1
    integer, intent(in) :: op2
    logical :: psyclone_cmp_int

    psyclone_cmp_int = op1 == op2

  end function psyclone_cmp_int
  pure function psyclone_cmp_logical(op1, op2)
    logical, intent(in) :: op1
    logical, intent(in) :: op2
    logical :: psyclone_cmp_logical

    psyclone_cmp_logical = op1 .EQV. op2

  end function psyclone_cmp_logical
  pure function psyclone_cmp_char(op1, op2)
    CHARACTER(LEN = *), INTENT(IN) :: op1
    CHARACTER(LEN = *), INTENT(IN) :: op2
    logical :: psyclone_cmp_char

    psyclone_cmp_char = op1(:) == op2(:)

  end function psyclone_cmp_char
  subroutine ukca_setup_chem()
    use ukca_config_specification_mod, only : glomap_config, i_blomquist, i_du_2mode, i_liss_merlivat, i_nightingale, i_ukca_chem_cristrat, i_ukca_chem_off, i_ukca_chem_offline, i_ukca_chem_offline_be, i_ukca_chem_raq, i_ukca_chem_strat, i_ukca_chem_strattrop, i_ukca_chem_trop, i_ukca_chem_tropisop, i_wanninkhof, int_method_be_explicit, int_method_none, int_method_nr, ukca_config
    use asad_mod, only : jpbk, jpcspf, jpctr, jpdd, jpdw, jphk, jpnr, jppj, jpspec, jptk
    use ukca_missing_data_mod, only : imdi
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_diag, ummessage, umprint
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use ukca_chem_offline, only : ndry_offline, nwet_offline
    use ukca_chem_aer, only : ndry_aer_aer, nwet_aer_aer
    use ukca_chem1_dat, only : ndry_trop, nwet_trop
    use ukca_chem_raq, only : ndry_raq, nwet_raq
    use ukca_chem_raqaero_mod, only : ndry_raqaero, nwet_raqaero
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SETUP_CHEM'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_off)) then
      ukca_config%l_ukca_chem = .false.
      ukca_config%ukca_int_method = int_method_none
      jpctr = 0
      jpspec = 0
      jpbk = 0
      jptk = 0
      jppj = 0
      jphk = 0
      jpdd = 0
      jpdw = 0
      jpcspf = 0
    else
      if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_trop)) then
        ukca_config%l_ukca_chem = .true.
        ukca_config%l_ukca_trop = .true.
        ukca_config%ukca_int_method = int_method_be_explicit
        if (ukca_config%l_ukca_chem_aero) then
          ukca_config%l_ukca_aerchem = .true.
          jpctr = 33
          jpspec = 53
          jpbk = 88
          jptk = 14
          jppj = 20
          jphk = 0
          jpdd = ndry_aer_aer
          jpdw = nwet_aer_aer
        else
          jpctr = 26
          jpspec = 46
          jpbk = 88
          jptk = 14
          jppj = 20
          jphk = 0
          jpdd = ndry_trop
          jpdw = nwet_trop
        end if
        jpcspf = jpctr
      else
        if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_raq)) then
          ukca_config%l_ukca_chem = .true.
          ukca_config%l_ukca_raq = .true.
          ukca_config%l_ukca_raqaero = .false.
          ukca_config%ukca_int_method = int_method_be_explicit
          jpctr = 40
          jpspec = 58
          jpbk = 113
          jptk = 12
          jppj = 23
          jphk = 0
          jpdd = ndry_raq
          jpdw = nwet_raq
          if (ukca_config%l_ukca_classic_hetchem) then
            jphk = jphk + 2
          end if
          if (ukca_config%l_ukca_chem_aero) then
            if (ukca_config%l_ukca_classic_hetchem) then
              cmessage(:) = 'RAQ-AERO is not compatibile with the CLASSIC ' // 'heterogeneous chemistry'
              errcode = 113
              call ereport('ukca_setup_chem', errcode, cmessage)
            end if
            ukca_config%l_ukca_raq = .false.
            ukca_config%l_ukca_raqaero = .true.
            jpctr = jpctr + 8
            jpspec = jpspec + 8
            jpbk = jpbk + 10
            jptk = jptk + 2
            jppj = jppj + 0
            jphk = jphk + 2
            jpdd = ndry_raqaero
            jpdw = nwet_raqaero
          end if
          jpcspf = jpctr
        else
          if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_offline_be)) then
            ukca_config%l_ukca_chem = .true.
            ukca_config%l_ukca_offline_be = .true.
            ukca_config%ukca_int_method = int_method_be_explicit
            jpctr = 7
            jpspec = 11
            jpbk = 9
            jptk = 1
            jppj = 0
            jphk = 3
            jpdd = ndry_offline
            jpdw = nwet_offline
            jpcspf = jpctr
          else
            if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_tropisop)) then
              ukca_config%l_ukca_chem = .true.
              ukca_config%l_ukca_tropisop = .true.
              ukca_config%ukca_int_method = int_method_nr
              if (ukca_config%l_ukca_chem_aero) then
                ukca_config%l_ukca_achem = .true.
                ukca_config%l_ukca_nr_aqchem = .true.
              else
                if (ukca_config%l_ukca_trophet) then
                  call umprint('trop het chem requires aerosol chemistry on', src='ukca_setup_chem_mod')
                  cmessage(:) = 'Unsupported option choice'
                  errcode = ABS(ukca_config%i_ukca_chem)
                  call ereport('ukca_setup_chem', errcode, cmessage)
                end if
              end if
            else
              if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_strattrop)) then
                ukca_config%l_ukca_chem = .true.
                ukca_config%l_ukca_strattrop = .true.
                ukca_config%l_ukca_advh2o = .true.
                ukca_config%ukca_int_method = int_method_nr
                if (ukca_config%l_ukca_chem_aero) then
                  ukca_config%l_ukca_achem = .true.
                  ukca_config%l_ukca_nr_aqchem = .true.
                else
                  if (ukca_config%l_ukca_trophet) then
                    call umprint('trop het chem requires aerosol chemistry on', src='ukca_setup_chem_mod')
                    cmessage(:) = 'Unsupported option choice'
                    errcode = ABS(ukca_config%i_ukca_chem)
                    call ereport('ukca_setup_chem', errcode, cmessage)
                  end if
                end if
              else
                if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_strat)) then
                  ukca_config%l_ukca_chem = .true.
                  ukca_config%l_ukca_strat = .true.
                  ukca_config%l_ukca_advh2o = .true.
                  ukca_config%ukca_int_method = int_method_nr
                  if (ukca_config%l_ukca_chem_aero) then
                    ukca_config%l_ukca_achem = .true.
                    ukca_config%l_ukca_nr_aqchem = .true.
                  end if
                  if (ukca_config%l_ukca_trophet) then
                    call umprint('Strat chem does not support trop het chem', src='ukca_setup_chem_mod')
                    cmessage(:) = 'Unsupported option choice'
                    errcode = ABS(ukca_config%i_ukca_chem)
                    call ereport('ukca_setup_chem', errcode, cmessage)
                  end if
                else
                  if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_offline)) then
                    ukca_config%l_ukca_chem = .true.
                    ukca_config%l_ukca_offline = .true.
                    ukca_config%l_ukca_nr_aqchem = .true.
                    ukca_config%l_ukca_advh2o = .false.
                    ukca_config%ukca_int_method = int_method_nr
                  else
                    if (psyclone_internal_cmp(ukca_config%i_ukca_chem, i_ukca_chem_cristrat)) then
                      ukca_config%l_ukca_chem = .true.
                      ukca_config%l_ukca_cristrat = .true.
                      ukca_config%l_ukca_advh2o = .true.
                      ukca_config%ukca_int_method = int_method_nr
                      if (ukca_config%l_ukca_chem_aero) then
                        ukca_config%l_ukca_achem = .true.
                        ukca_config%l_ukca_nr_aqchem = .true.
                      else
                        if (ukca_config%l_ukca_trophet) then
                          call umprint('trop het chem requires aerosol chemistry on', src='ukca_setup_chem_mod')
                          cmessage(:) = 'Unsupported option choice'
                          errcode = ABS(ukca_config%i_ukca_chem)
                          call ereport('ukca_setup_chem', errcode, cmessage)
                        end if
                      end if
                    else
                      ! PSyclone CodeBlock (unsupported code) reason:
                      !  - Unsupported statement: Write_Stmt
                      WRITE(umMessage, '(A,1X,I6)') 'Unknown chemistry scheme: ', ukca_config % i_ukca_chem
                      call umprint(ummessage, src='ukca_setup_chem_mod')
                      cmessage(:) = 'Unknown chemistry scheme'
                      errcode = ABS(ukca_config%i_ukca_chem)
                      call ereport('ukca_setup_chem', errcode, cmessage)
                    end if
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if
    jpnr = jpbk + jptk + jppj + jphk
    if (ukca_config%ukca_int_method == int_method_nr .OR. ukca_config%i_ukca_chem == i_ukca_chem_offline_be) then
      if (ukca_config%chem_timestep == imdi) then
        cmessage(:) = ' The chemical timestep has not been set'
        errcode = 1
        call ereport('ukca_setup_chem', errcode, cmessage)
      end if
    end if
    if (ukca_config%ukca_int_method == int_method_be_explicit .AND. ukca_config%i_ukca_chem /= i_ukca_chem_offline_be .AND. ukca_config%i_ukca_chem /= i_ukca_chem_raq) then
      if (ukca_config%l_ukca_chem_aero .OR. ukca_config%l_ukca_trophet) then
        call umprint('This BE scheme does not support add-on aerosol ' // 'chemistry or tropospheric heterogeneous chemistry', src='ukca_setup_chem_mod')
        cmessage(:) = 'Unsupported option combination'
        errcode = ABS(ukca_config%i_ukca_chem)
        call ereport('ukca_setup_chem', errcode, cmessage)
      end if
    end if
    if (ukca_config%l_ukca_trophet .AND. (.NOT.ukca_config%l_ukca_mode)) then
      call umprint('Need Glomap MODE on for ' // 'tropospheric heterogeneous chemistry', src='ukca_setup_chem_mod')
      cmessage(:) = 'Unsupported option combination'
      errcode = 1
      call ereport('ukca_setup_chem', errcode, cmessage)
    end if
    if (ukca_config%l_ukca_mode) then
      if (ukca_config%ukca_int_method == int_method_nr) then
        if (.NOT.ukca_config%l_ukca_nr_aqchem) then
          call umprint('Need aerosol chemistry on for ' // 'Glomap MODE', src='ukca_setup_chem_mod')
          cmessage(:) = 'Unsupported option combination'
          errcode = 2
          call ereport('ukca_setup_chem', errcode, cmessage)
        end if
      else
        if (ukca_config%ukca_int_method == int_method_none) then
          if (.NOT.glomap_config%i_mode_setup == i_du_2mode) then
            call umprint('Need chemistry scheme for all aerosol ' // 'setups except 2-mode dust', src='ukca_setup_chem_mod')
            cmessage(:) = 'Unsupported option combination'
            errcode = 3
            call ereport('ukca_setup_chem', errcode, cmessage)
          end if
        else
          if (.NOT.(ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_offline_be .OR. ukca_config%l_ukca_raqaero)) then
            call umprint('Need aerosol chemistry on for ' // 'Glomap MODE', src='ukca_setup_chem_mod')
            cmessage(:) = 'Unsupported option combination'
            errcode = 4
            call ereport('ukca_setup_chem', errcode, cmessage)
          end if
        end if
      end if
    end if
    if (ukca_config%l_ukca_chem_aero .AND. ukca_config%l_seawater_dms) then
      if (psyclone_internal_cmp(ukca_config%i_ukca_dms_flux, i_liss_merlivat)) then
      else
        if (psyclone_internal_cmp(ukca_config%i_ukca_dms_flux, i_wanninkhof)) then
        else
          if (psyclone_internal_cmp(ukca_config%i_ukca_dms_flux, i_nightingale)) then
          else
            if (psyclone_internal_cmp(ukca_config%i_ukca_dms_flux, i_blomquist)) then
            else
              call umprint('Marine DMS emissions expected but no UKCA DMS scheme is selected ' // 'i_ukca_dms_flux should be 1,2,3 or 4')
              errcode = 5
              cmessage(:) = 'RUN_UKCA: DMS flux scheme not specified'
              call ereport('UKCA_SETUP_CHEM_MOD', errcode, cmessage)
            end if
          end if
        end if
      end if
    end if
    if (printstatus >= prstatus_diag) then
      call umprint('Logicals and parameters set by ukca_setup_chem:', src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_chem      = ', ukca_config % l_ukca_chem
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_trop      = ', ukca_config % l_ukca_trop
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_achem   = ', ukca_config % l_ukca_achem
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_nr_aqchem = ', ukca_config % l_ukca_nr_aqchem
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_aerchem   = ', ukca_config % l_ukca_aerchem
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_tropisop  = ', ukca_config % l_ukca_tropisop
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_strattrop = ', ukca_config % l_ukca_strattrop
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_strat     = ', ukca_config % l_ukca_strat
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_offline   = ', ukca_config % l_ukca_offline
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,L7)') 'l_ukca_cristrat  = ', ukca_config % l_ukca_cristrat
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'ukca_int_method  = ', ukca_config % ukca_int_method
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpctr            = ', jpctr
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpcspf           = ', jpcspf
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpspec           = ', jpspec
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpbk             = ', jpbk
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jptk             = ', jptk
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jppj             = ', jppj
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jphk             = ', jphk
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpnr             = ', jpnr
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpdd             = ', jpdd
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I5)') 'jpdw             = ', jpdw
      call umprint(ummessage, src='ukca_setup_chem_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I2)') 'i_ukca_dms_flux  = ', ukca_config % i_ukca_dms_flux
      call umprint(ummessage, src='ukca_setup_chem_mod')
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if

  end subroutine ukca_setup_chem

end module ukca_setup_chem_mod
