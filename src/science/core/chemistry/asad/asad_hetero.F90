module asad_hetero_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_HETERO_MOD'
  public

  contains
  subroutine asad_hetero(n_points, cld_f, cld_l, rc_het, h_plus_1d_arr)
    use asad_findreaction_mod, only : asad_findreaction
    use asad_mod, only : ih2o2_oh, ih_h2o2, ih_hno3, ih_o3, ih_o3_const, ih_so2, ihno3_oh, iho2_h, ihso3_h2o2, ihso3_o3, in2o5_h, iso3_o3, jpbk, jpdw, jpeq, jphk, jpspb, jpsph, nbrkx, nhrkx, p, rk, spb, sph, t, tnd
    use ukca_config_specification_mod, only : ukca_config
    use ukca_chem_offline, only : nwet_constant
    use ukca_fdiss_constant_mod, only : ukca_fdiss_constant
    use ukca_config_constants_mod, only : avogadro, rho_water
    use ukca_constants, only : h_plus, m_air
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    use ukca_fdiss_mod, only : ukca_fdiss
    real, parameter :: qcl_min = 1.0e-12
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, dimension(n_points), intent(in) :: cld_f
    real, dimension(n_points), intent(in) :: cld_l
    real, dimension(n_points,2), intent(in) :: rc_het
    real, dimension(n_points), intent(in) :: h_plus_1d_arr
    real, dimension(n_points) :: vr
    real, dimension(n_points,jpdw,jpeq + 1) :: fdiss
    real, allocatable, dimension(:,:,:) :: fdiss_constant
    real, dimension(n_points) :: fdiss_o3
    integer, save :: icode = 0
    CHARACTER(LEN = errormessagelength) :: cmessage
    CHARACTER(LEN = 10) :: prods(2)
    logical, save :: first = .true.
    logical, save :: first_pass = .true.
    logical, dimension(n_points) :: todo
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_HETERO'
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
    integer :: widx1
    integer :: widx1_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first_pass) then
      if (first) then
        if (ukca_config%l_ukca_achem) then
          if (ih_o3 == 0 .OR. ih_h2o2 == 0 .OR. ih_so2 == 0 .OR. ih_hno3 == 0) then
            cmessage(:) = ' Indicies for Aqueous chemistry uninitialised' // ' - O3, H2O2, SO2, and HNO3 must be made ' // ' soluble species in chch_defs array'
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A9,I5)') 'ih_o3:   ', ih_o3
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A9,I5)') 'ih_h2o2: ', ih_h2o2
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A9,I5)') 'ih_hno3: ', ih_hno3
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A9,I5)') 'ih_so2:  ', ih_so2
            call umprint(ummessage, src='asad_hetero')
            icode = 1
            call ereport('ASAD_HETERO', icode, cmessage)
          end if
        end if
        if (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) then
          if (ih_h2o2 == 0 .OR. ih_so2 == 0) then
            cmessage(:) = ' Indices for Aqueous chemistry uninitialised' // ' - H2O2, and SO2, must be made ' // ' soluble species in chch_defs array'
            icode = 1
          end if
          if (icode /= 0) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A10,I5)') 'ih_h2o2: ', ih_h2o2
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A10,I5)') 'ih_so2:  ', ih_so2
            call umprint(ummessage, src='asad_hetero')
            call ereport('ASAD_HETERO', icode, cmessage)
          end if
        end if
        if (ukca_config%l_ukca_trophet .AND. (.NOT.ukca_config%l_ukca_mode)) then
          cmessage(:) = ' Tropospheric heterogeneous chemistry is flagged' // ' but MODE aerosol scheme is not in use'
          icode = 1
          call ereport('ASAD_HETERO', icode, cmessage)
        end if
        ihso3_h2o2 = 0
        iso3_o3 = 0
        ihso3_o3 = 0
        ih2o2_oh = 0
        ihno3_oh = 0
        in2o5_h = 0
        iho2_h = 0
        if (ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_offline_be) then
          prods(:) = ['NULL0     ', '          ']
          ihso3_h2o2 = asad_findreaction('SO2       ', 'H2O2      ', prods, 2, sph, nhrkx, jphk + 1, jpsph)
          prods(:) = ['NULL1     ', '          ']
          ihso3_o3 = asad_findreaction('SO2       ', 'O3        ', prods, 2, sph, nhrkx, jphk + 1, jpsph)
          prods(:) = ['NULL2     ', '          ']
          iso3_o3 = asad_findreaction('SO2       ', 'O3        ', prods, 2, sph, nhrkx, jphk + 1, jpsph)
          if (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) then
            prods(:) = ['H2O       ', '          ']
          else
            prods(:) = ['H2O       ', 'HO2       ']
          end if
          ih2o2_oh = asad_findreaction('H2O2      ', 'OH        ', prods, 2, spb, nbrkx, jpbk + 1, jpspb)
          prods(:) = ['H2O       ', 'NO3       ']
          ihno3_oh = asad_findreaction('HONO2     ', 'OH        ', prods, 2, spb, nbrkx, jpbk + 1, jpspb)
          icode = 0
          if (ihso3_h2o2 == 0 .OR. iso3_o3 == 0 .OR. ih2o2_oh == 0 .OR. ihso3_o3 == 0) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A12,I5)') 'ihso3_h2o2: ', ihso3_h2o2
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A12,I5)') 'ihso3_o3: ', ihso3_o3
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A12,I5)') 'iso3_o3: ', iso3_o3
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A12,I5)') 'ih2o2_oh: ', ih2o2_oh
            call umprint(ummessage, src='asad_hetero')
            icode = 1
          end if
          if (ukca_config%l_ukca_achem .AND. ihno3_oh == 0) then
            icode = 1
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A12,I5)') 'ihno3_oh: ', ihno3_oh
            call umprint(ummessage, src='asad_hetero')
          end if
          if (icode > 0) then
            cmessage(:) = ' Heterogeneous chemistry called, but eqns' // ' not found - see output'
            call ereport('ASAD_HETERO', icode, cmessage)
          end if
        end if
        if (ukca_config%l_ukca_trophet) then
          prods(:) = ['HONO2     ', '          ']
          in2o5_h = asad_findreaction('N2O5      ', '          ', prods, 2, sph, nhrkx, jphk + 1, jpsph)
          prods(:) = ['H2O2      ', '          ']
          iho2_h = asad_findreaction('HO2       ', '          ', prods, 2, sph, nhrkx, jphk + 1, jpsph)
          if (iho2_h == 0 .OR. in2o5_h == 0) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A9,I5)') 'in2o5_h: ', in2o5_h
            call umprint(ummessage, src='asad_hetero')
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(A9,I5)') 'iho2_h: ', iho2_h
            call umprint(ummessage, src='asad_hetero')
            cmessage(:) = ' Tropospheric heterogeneous chemistry is flagged,' // ' but equations not found - see output'
            icode = 1
            call ereport('ASAD_HETERO', icode, cmessage)
          end if
        end if
        first = .false.
      end if
    end if
    if ((ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_offline_be) .AND. ANY(cld_l > qcl_min)) then
      call ukca_fdiss(n_points, qcl_min, t, p, cld_l, fdiss, h_plus_1d_arr)
      !$acc kernels
      !$acc loop gang vector independent
      do idx = LBOUND(todo, dim=1), UBOUND(todo, dim=1), 1
        todo(idx) = cld_l(idx) > qcl_min
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop gang vector independent collapse(3)
      do idx_1 = LBOUND(fdiss, dim=3), UBOUND(fdiss, dim=3), 1
        do idx_2 = LBOUND(fdiss, dim=2), UBOUND(fdiss, dim=2), 1
          do idx_3 = LBOUND(fdiss, dim=1), UBOUND(fdiss, dim=1), 1
            fdiss(idx_3,idx_2,idx_1) = 0.0
          enddo
        enddo
      enddo
      !$acc loop gang vector independent
      do idx_4 = LBOUND(todo, dim=1), UBOUND(todo, dim=1), 1
        todo(idx_4) = .false.
      enddo
      !$acc end kernels
    end if
    if (ANY(cld_l > qcl_min)) then
      if ((ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) .AND. nwet_constant > 0) then
        ALLOCATE(fdiss_constant(1:n_points,1:nwet_constant,1:jpeq + 1))
        call ukca_fdiss_constant(n_points, qcl_min, t, p, cld_l, fdiss_constant, h_plus_1d_arr)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_5 = LBOUND(fdiss_o3, dim=1), UBOUND(fdiss_o3, dim=1), 1
          fdiss_o3(idx_5) = fdiss_constant(idx_5 + (LBOUND(fdiss_constant, dim=1) - LBOUND(fdiss_o3, dim=1)),ih_o3_const,1)
        enddo
        !$acc end kernels
        DEALLOCATE(fdiss_constant)
      else
        !$acc kernels
        !$acc loop gang vector independent
        do idx_6 = LBOUND(fdiss_o3, dim=1), UBOUND(fdiss_o3, dim=1), 1
          fdiss_o3(idx_6) = fdiss(idx_6,ih_o3,1)
        enddo
        !$acc end kernels
      end if
    end if
    if (ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_offline_be) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (todo(:))
        vr(:) = cld_l(:) * tnd(:) * m_air * (1E6 / avogadro) / rho_water
        rk(:, ihso3_h2o2) = 2.1295E+14 * EXP(- 4430.0 / t(:)) * (H_plus_1d_arr(:) / (1.0 + 13.0 * H_plus_1d_arr(:))) * cld_f(:) * fdiss(:, ih_so2, 2) * fdiss(:, ih_h2o2, 1) * 1000.0 / (avogadro * vr(:))
        rk(:, ihso3_o3) = 4.0113E+13 * EXP(- 5530.0 / t(:)) * cld_f(:) * fdiss(:, ih_so2, 2) * fdiss_o3(:) * 1000.0 / (avogadro * vr(:))
        rk(:, iso3_o3) = 7.43E+16 * EXP(- 5280.0 / t(:)) * cld_f(:) * fdiss(:, ih_so2, 3) * fdiss_o3(:) * 1000.0 / (avogadro * vr(:))
        rk(:, ih2o2_oh) = rk(:, ih2o2_oh) * (1.0 - (fdiss(:, ih_h2o2, 1) + fdiss(:, ih_h2o2, 2)) * cld_f(:))
      ELSEWHERE
        rk(:, ihso3_h2o2) = 0.0
        rk(:, ihso3_o3) = 0.0
        rk(:, iso3_o3) = 0.0
      END WHERE
    end if
    if (ukca_config%l_ukca_achem) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (todo(:))
        rk(:, ihno3_oh) = rk(:, ihno3_oh) * (1.0 - (fdiss(:, ih_hno3, 1) + fdiss(:, ih_hno3, 2)) * cld_f(:))
      END WHERE
    end if
    if (ukca_config%l_ukca_trophet) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_7 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
        rk(idx_7,in2o5_h) = rc_het(idx_7 + (LBOUND(rc_het, dim=1) - LBOUND(rk, dim=1)),1)
      enddo
      !$acc loop gang vector independent
      do idx_8 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
        rk(idx_8,iho2_h) = rc_het(idx_8 + (LBOUND(rc_het, dim=1) - LBOUND(rk, dim=1)),2)
      enddo
      !$acc end kernels
    else
      if (in2o5_h > 0) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_9 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_9,in2o5_h) = 0.0
        enddo
        !$acc end kernels
      end if
      if (iho2_h > 0) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_10 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_10,iho2_h) = 0.0
        enddo
        !$acc end kernels
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_hetero

end module asad_hetero_mod
