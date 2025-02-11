module ukca_trop_hetchem_mod
  use ukca_config_constants_mod, only : avogadro, rmol
  use ukca_constants, only : m_air, m_ho2, m_n2o5, pi
  use ukca_mode_setup, only : cp_bc, cp_cl, cp_du, cp_nh4, cp_nn, cp_no3, cp_oc, cp_so, cp_su, nmodes
  use ukca_config_specification_mod, only : glomap_variables, ukca_config
  use yomhook, only : dr_hook, lhook
  use parkind1, only : jpim, jprb
  use ereport_mod, only : ereport
  use umPrintMgr, only : printstatus, prstatus_diag, ummessage, umprint
  use errormessagelength_mod, only : errormessagelength
  implicit none
  integer, parameter, public :: ihet_n2o5 = 1
  integer, parameter, public :: ihet_ho2_ho2 = 2
  integer, private :: errcode
  character(len = *), parameter, private :: ModuleName = 'UKCA_TROP_HETCHEM_MOD'
  private

  public :: cal_n2o5, cal_loss1k, ukca_trop_hetchem

  contains
  subroutine ukca_trop_hetchem(nbox, nhet, temp, rh, aird, pvol, wetdp, sarea, het_rates)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nbox
    integer, intent(in) :: nhet
    real, dimension(nbox), intent(in) :: temp
    real, dimension(nbox), intent(in) :: rh
    real, dimension(nbox), intent(in) :: aird
    real, dimension(nbox,nmodes,glomap_variables%ncp), intent(in) :: pvol
    real, dimension(nbox,nmodes), intent(in) :: wetdp
    real, dimension(nbox,nmodes), intent(in) :: sarea
    real, dimension(nbox,nhet), intent(out) :: het_rates
    LOGICAL, POINTER :: component(:, :)
    LOGICAL, POINTER :: mode(:)
    INTEGER, POINTER :: ncp
    integer :: i
    integer :: imode
    integer :: icp
    real, dimension(nhet) :: molec_weight
    real, dimension(nbox,nhet,glomap_variables%ncp) :: sticking_cf
    real, dimension(nbox,nmodes) :: total_vol
    real, dimension(nbox,nmodes,glomap_variables%ncp) :: frac_vol
    real, dimension(nbox,nhet,nmodes) :: mode_gamma
    real, dimension(nbox,nhet,nmodes) :: mode_rate
    real, dimension(nbox) :: total_sa
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_TROP_HETCHEM'
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
    integer :: widx1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    component => glomap_variables%component
    mode => glomap_variables%mode
    ncp => glomap_variables%ncp
    if (printstatus >= prstatus_diag) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'UKCA_TROP_HETCHEM Routine:'
      call umprint(ummessage, src='ukca_trop_hetchem')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I5)') 'nbox: ', nbox
      call umprint(ummessage, src='ukca_trop_hetchem')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,3E15.6)') 'temp: ', MAXVAL(temp), MINVAL(temp), SUM(temp) / SIZE(temp)
      call umprint(ummessage, src='ukca_trop_hetchem')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,3E15.6)') 'rh  : ', MAXVAL(rh), MINVAL(rh), SUM(rh) / SIZE(rh)
      call umprint(ummessage, src='ukca_trop_hetchem')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,3E15.6)') 'aird: ', MAXVAL(aird), MINVAL(aird), SUM(aird) / SIZE(aird)
      call umprint(ummessage, src='ukca_trop_hetchem')
      do imode = 1, nmodes, 1
        if (mode(imode)) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,I0,3E15.6)') 'sarea: ', imode, MAXVAL(sarea(:, imode)), MINVAL(sarea(:, imode)), SUM(sarea(:, imode)) / SIZE(sarea(:, 1))
          call umprint(ummessage, src='ukca_trop_hetchem')
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,I0,3E15.6)') 'wetdp: ', imode, MAXVAL(wetdp(:, imode)), MINVAL(wetdp(:, imode)), SUM(wetdp(:, imode)) / SIZE(wetdp(:, 1))
          call umprint(ummessage, src='ukca_trop_hetchem')
        end if
      enddo
      do imode = 1, nmodes, 1
        if (mode(imode)) then
          do icp = 1, ncp, 1
            if (component(imode,icp) .AND. icp /= cp_nh4) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(A,2I0,3E15.6)') 'pvol: ', imode, icp, MAXVAL(pvol(:, imode, icp)), MINVAL(pvol(:, imode, icp)), SUM(pvol(:, imode, icp)) / SIZE(pvol(:, 1, 1))
              call umprint(ummessage, src='ukca_trop_hetchem')
            end if
          enddo
        end if
      enddo
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(molec_weight, dim=1), UBOUND(molec_weight, dim=1), 1
      molec_weight(idx) = 9.0e5
    enddo
    molec_weight(ihet_n2o5) = m_n2o5
    molec_weight(ihet_ho2_ho2) = m_ho2
    !$acc end kernels
    if (ANY(molec_weight > 1.0e5)) then
      cmessage(:) = ' Molecular weights not defined for all cases'
      errcode = 1
      call ereport('UKCA_TROP_HETCHEM', errcode, cmessage)
    end if
    if (printstatus >= prstatus_diag) then
      call umprint('Sticking_cf:')
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx_1 = LBOUND(sticking_cf, dim=3), UBOUND(sticking_cf, dim=3), 1
      do idx_2 = LBOUND(sticking_cf, dim=2), UBOUND(sticking_cf, dim=2), 1
        do idx_3 = LBOUND(sticking_cf, dim=1), UBOUND(sticking_cf, dim=1), 1
          sticking_cf(idx_3,idx_2,idx_1) = 0.0e0
        enddo
      enddo
    enddo
    !$acc end kernels
    do i = 1, nhet, 1
      if (i == ihet_n2o5) then
        do icp = 1, ncp, 1
          if (icp /= cp_nh4) then
            sticking_cf(:,i,icp) = cal_n2o5(icp, nbox, temp(:), rh(:))
            if (printstatus >= prstatus_diag) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2I0,3E15.6)') i, icp, MAXVAL(sticking_cf(:, i, icp)), MINVAL(sticking_cf(:, i, icp)), SUM(sticking_cf(:, i, icp)) / SIZE(sticking_cf(:, i, icp))
              call umprint(ummessage, src='ukca_trop_hetchem')
            end if
          end if
        enddo
      else
        if (i == ihet_ho2_ho2) then
          do icp = 1, ncp, 1
            if (icp /= cp_nh4) then
              sticking_cf(:,i,icp) = cal_ho2(icp, nbox, temp(:), rh(:))
              if (printstatus >= prstatus_diag) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(2I0,E15.6)') i, icp, MAXVAL(sticking_cf(:, i, icp))
                call umprint(ummessage, src='ukca_trop_hetchem')
              end if
            end if
          enddo
        else
          cmessage(:) = ' Error in calculating sticking_cf'
          errcode = 1
          call ereport('UKCA_TROP_HETCHEM', errcode, cmessage)
        end if
      end if
    enddo
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_4 = LBOUND(total_vol, dim=2), UBOUND(total_vol, dim=2), 1
      do idx_5 = LBOUND(total_vol, dim=1), UBOUND(total_vol, dim=1), 1
        total_vol(idx_5,idx_4) = 0.0e0
      enddo
    enddo
    !$acc loop seq
    do imode = 1, nmodes, 1
      if (mode(imode)) then
        !$acc loop seq
        do icp = 1, ncp, 1
          if (component(imode,icp) .AND. icp /= cp_nh4) then
            !$acc loop gang vector independent
            do idx_6 = LBOUND(total_vol, dim=1), UBOUND(total_vol, dim=1), 1
              total_vol(idx_6,imode) = total_vol(idx_6,imode) + pvol(idx_6,imode,icp)
            enddo
          end if
        enddo
      end if
    enddo
    !$acc end kernels
    if (printstatus >= prstatus_diag) then
      call umprint('frac_vol:')
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx_7 = LBOUND(frac_vol, dim=3), UBOUND(frac_vol, dim=3), 1
      do idx_8 = LBOUND(frac_vol, dim=2), UBOUND(frac_vol, dim=2), 1
        do idx_9 = LBOUND(frac_vol, dim=1), UBOUND(frac_vol, dim=1), 1
          frac_vol(idx_9,idx_8,idx_7) = 0.0e0
        enddo
      enddo
    enddo
    !$acc end kernels
    do imode = 1, nmodes, 1
      if (mode(imode)) then
        do icp = 1, ncp, 1
          if (component(imode,icp) .AND. icp /= cp_nh4) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_10 = LBOUND(frac_vol, dim=1), UBOUND(frac_vol, dim=1), 1
              frac_vol(idx_10,imode,icp) = pvol(idx_10,imode,icp) / total_vol(idx_10,imode)
            enddo
            !$acc end kernels
            if (printstatus >= prstatus_diag) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Write_Stmt
              WRITE(umMessage, '(2I0,3E15.6)') imode, icp, MAXVAL(frac_vol(:, imode, icp)), MINVAL(frac_vol(:, imode, icp)), SUM(frac_vol(:, imode, icp)) / SIZE(frac_vol(:, imode, icp))
              call umprint(ummessage, src='ukca_trop_hetchem')
            end if
          end if
        enddo
      end if
    enddo
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx_11 = LBOUND(mode_gamma, dim=3), UBOUND(mode_gamma, dim=3), 1
      do idx_12 = LBOUND(mode_gamma, dim=2), UBOUND(mode_gamma, dim=2), 1
        do idx_13 = LBOUND(mode_gamma, dim=1), UBOUND(mode_gamma, dim=1), 1
          mode_gamma(idx_13,idx_12,idx_11) = 0.0
        enddo
      enddo
    enddo
    !$acc loop seq
    do i = 1, nhet, 1
      !$acc loop seq
      do imode = 1, nmodes, 1
        if (mode(imode)) then
          !$acc loop seq
          do icp = 1, ncp, 1
            if (component(imode,icp) .AND. icp /= cp_nh4) then
              !$acc loop gang vector independent
              do idx_14 = LBOUND(mode_gamma, dim=1), UBOUND(mode_gamma, dim=1), 1
                mode_gamma(idx_14,i,imode) = mode_gamma(idx_14,i,imode) + frac_vol(idx_14,imode,icp) * sticking_cf(idx_14,i,icp)
              enddo
            end if
          enddo
        end if
      enddo
    enddo
    !$acc end kernels
    if (printstatus >= prstatus_diag) then
      call umprint('Mode_rate:')
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx_15 = LBOUND(mode_rate, dim=3), UBOUND(mode_rate, dim=3), 1
      do idx_16 = LBOUND(mode_rate, dim=2), UBOUND(mode_rate, dim=2), 1
        do idx_17 = LBOUND(mode_rate, dim=1), UBOUND(mode_rate, dim=1), 1
          mode_rate(idx_17,idx_16,idx_15) = 0.0e0
        enddo
      enddo
    enddo
    !$acc end kernels
    do i = 1, nhet, 1
      do imode = 1, nmodes, 1
        if (mode(imode)) then
          mode_rate(:,i,imode) = cal_loss1k(nbox, sarea(:,imode), wetdp(:,imode), aird(:), mode_gamma(:,i,imode), temp(:), molec_weight(i))
          if (printstatus >= prstatus_diag) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Write_Stmt
            WRITE(umMessage, '(I0,3E15.6)') imode, MAXVAL(mode_rate(:, i, imode)), MINVAL(mode_rate(:, i, imode)), SUM(mode_rate(:, i, imode)) / SIZE(mode_rate(:, i, imode))
            call umprint(ummessage, src='ukca_trop_hetchem')
          end if
        end if
      enddo
    enddo
    !$acc kernels
    !$acc loop gang vector independent
    do idx_18 = LBOUND(total_sa, dim=1), UBOUND(total_sa, dim=1), 1
      total_sa(idx_18) = 0.0e0
    enddo
    !$acc loop seq
    do imode = 1, nmodes, 1
      if (mode(imode)) then
        !$acc loop gang vector independent
        do idx_19 = LBOUND(total_sa, dim=1), UBOUND(total_sa, dim=1), 1
          total_sa(idx_19) = total_sa(idx_19) + sarea(idx_19,imode)
        enddo
      end if
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_20 = LBOUND(het_rates, dim=2), UBOUND(het_rates, dim=2), 1
      do idx_21 = LBOUND(het_rates, dim=1), UBOUND(het_rates, dim=1), 1
        het_rates(idx_21,idx_20) = 0.0e0
      enddo
    enddo
    !$acc loop seq
    do i = 1, nhet, 1
      !$acc loop seq
      do imode = 1, nmodes, 1
        if (mode(imode)) then
          !$acc loop gang vector independent
          do widx1 = 1, nbox, 1
            if (total_sa(widx1) > 1.0e-30) then
              het_rates(widx1,i) = het_rates(widx1,i) + mode_rate(widx1,i,imode)
            end if
          enddo
        end if
      enddo
    enddo
    !$acc end kernels
    if (printstatus >= prstatus_diag) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,3E15.6)') 'rc_n2o5: ', MAXVAL(het_rates(:, 1)), MINVAL(het_rates(:, 1)), SUM(het_rates(:, 1)) / SIZE(het_rates(:, 1))
      call umprint(ummessage, src='ukca_trop_hetchem')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,3E15.6)') 'rc_ho2:  ', MAXVAL(het_rates(:, 2)), MINVAL(het_rates(:, 2)), SUM(het_rates(:, 2)) / SIZE(het_rates(:, 2))
      call umprint(ummessage, src='ukca_trop_hetchem')
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_trop_hetchem
  function cal_n2o5(aerotype, nbox, temp, rh)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: aerotype
    integer, intent(in) :: nbox
    real, dimension(nbox), intent(in) :: temp
    real, dimension(nbox), intent(in) :: rh
    real, dimension(nbox) :: cal_n2o5
    real, dimension(nbox) :: ttemp
    real, dimension(nbox) :: rh_p
    real :: maxtemp
    real, dimension(nbox) :: lambda
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'CAL_N2O5'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: widx1
    integer :: widx1_1
    integer :: widx1_2
    integer :: widx1_3
    integer :: widx1_4

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(rh_p, dim=1), UBOUND(rh_p, dim=1), 1
      rh_p(idx) = rh(idx) * 100.0
    enddo
    !$acc end kernels
    if (aerotype == cp_du) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_1 = LBOUND(cal_n2o5, dim=1), UBOUND(cal_n2o5, dim=1), 1
        cal_n2o5(idx_1) = 0.1
      enddo
      !$acc end kernels
    else
      if (aerotype == cp_su) then
        !$acc kernels
        !$acc loop gang vector independent
        do idx_2 = LBOUND(cal_n2o5, dim=1), UBOUND(cal_n2o5, dim=1), 1
          cal_n2o5(idx_2) = 2.79e-4 + rh_p(idx_2) * (1.30e-4 + rh_p(idx_2) * (-3.43e-6 + rh_p(idx_2) * 7.52e-8))
        enddo
        maxtemp = 282.0 - 294.0
        !$acc loop gang vector independent
        do widx1 = 1, nbox, 1
          if (temp(widx1) < 282.0e0) then
            ttemp(widx1) = maxtemp
          else
            ttemp(widx1) = temp(widx1) - 294.0
          end if
        enddo
        !$acc loop gang vector independent
        do idx_3 = LBOUND(cal_n2o5, dim=1), UBOUND(cal_n2o5, dim=1), 1
          cal_n2o5(idx_3) = cal_n2o5(idx_3) * 10.0 ** (-2.0 - 4.0e-2 * ttemp(idx_3)) / 1.0e-2
        enddo
        !$acc end kernels
      else
        if (aerotype == cp_bc) then
          !$acc kernels
          !$acc loop gang vector independent
          do idx_4 = LBOUND(cal_n2o5, dim=1), UBOUND(cal_n2o5, dim=1), 1
            cal_n2o5(idx_4) = 0.005e0
          enddo
          !$acc end kernels
        else
          if (aerotype == cp_oc .OR. aerotype == cp_so) then
            !$acc kernels
            !$acc loop gang vector independent
            do widx1_1 = 1, nbox, 1
              if (rh_p(widx1_1) >= 57.0e0) then
                cal_n2o5(widx1_1) = 0.03e0
              else
                cal_n2o5(widx1_1) = rh_p(widx1_1) * 5.2e-4
              end if
            enddo
            !$acc end kernels
          else
            if (aerotype == cp_cl) then
              !$acc kernels
              !$acc loop gang vector independent
              do widx1_2 = 1, nbox, 1
                if (rh_p(widx1_2) >= 62.0) then
                  cal_n2o5(widx1_2) = 0.03e0
                else
                  cal_n2o5(widx1_2) = rh_p(widx1_2) * 0.0005e0
                end if
              enddo
              !$acc end kernels
            else
              if (aerotype == cp_no3) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_5 = LBOUND(lambda, dim=1), UBOUND(lambda, dim=1), 1
                  lambda(idx_5) = -8.107744 + 0.04902 * rh_p(idx_5)
                enddo
                !$acc loop gang vector independent
                do idx_6 = LBOUND(cal_n2o5, dim=1), UBOUND(cal_n2o5, dim=1), 1
                  cal_n2o5(idx_6) = 1.0 / (1.0 + EXP(-1.0 * lambda(idx_6)))
                enddo
                !$acc loop gang vector independent
                do widx1_3 = 1, nbox, 1
                  if (cal_n2o5(widx1_3) > 0.0154) then
                    cal_n2o5(widx1_3) = 0.0154
                  end if
                enddo
                !$acc end kernels
              else
                if (aerotype == cp_nn) then
                  !$acc kernels
                  !$acc loop gang vector independent
                  do widx1_4 = 1, nbox, 1
                    if (rh_p(widx1_4) < 50.0) then
                      cal_n2o5(widx1_4) = 0.0018
                    else
                      if (rh_p(widx1_4) >= 50.0 .AND. rh_p(widx1_4) < 60.0) then
                        cal_n2o5(widx1_4) = 0.0018 + (0.0032 - 0.0018) * (rh_p(widx1_4) - 50.0) / 10.0
                      else
                        if (rh_p(widx1_4) >= 60.0 .AND. rh_p(widx1_4) < 90.0) then
                          cal_n2o5(widx1_4) = 0.0032 + (0.0230 - 0.0032) * (rh_p(widx1_4) - 60.0) / 30.0
                        else
                          if (rh_p(widx1_4) >= 90.0) then
                            cal_n2o5(widx1_4) = 0.023
                          end if
                        end if
                      end if
                    end if
                  enddo
                  !$acc end kernels
                else
                  cmessage(:) = ' Unknown aerosol surface for N2O5 hydrolysis'
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') cmessage
                  call umprint(ummessage, src='ukca_trop_hetchem')
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,I0)') 'AEROSOL TYPE =', aerotype
                  call umprint(ummessage, src='ukca_trop_hetchem')
                  errcode = aerotype
                  call ereport('UKCA_TROP_HETCHEM_MOD:CAL_N2O5', errcode, cmessage)
                end if
              end if
            end if
          end if
        end if
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end function cal_n2o5
  function cal_ho2(aerotype, nbox, temp, rh)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: aerotype
    integer, intent(in) :: nbox
    real, dimension(nbox), intent(in) :: temp
    real, dimension(nbox), intent(in) :: rh
    real, dimension(nbox) :: cal_ho2
    real, dimension(nbox) :: rh_p
    real, dimension(nbox) :: factt
    real, dimension(nbox) :: facth
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'CAL_HO2'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: widx1
    integer :: widx1_1
    integer :: widx1_2

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(rh_p, dim=1), UBOUND(rh_p, dim=1), 1
      rh_p(idx) = rh(idx) * 100.0
    enddo
    !$acc end kernels
    if (aerotype == cp_du) then
      !$acc kernels
      !$acc loop gang vector independent
      do widx1 = 1, nbox, 1
        if (rh_p(widx1) >= 50.0) then
          cal_ho2(widx1) = 0.1e0
        else
          cal_ho2(widx1) = 0.05e0
        end if
      enddo
      !$acc end kernels
    else
      if (aerotype == cp_su) then
        !$acc kernels
        !$acc loop gang vector independent
        do widx1_1 = 1, nbox, 1
          if (rh_p(widx1_1) <= 35.0) then
            cal_ho2(widx1_1) = 0.01e0
          else
            factt(widx1_1) = 5.66e-5 * EXP(1560.0 / temp(widx1_1)) / 0.11
            facth(widx1_1) = (-2.86 * EXP(-0.078 * rh_p(widx1_1)) + 0.192) / 0.11
            cal_ho2(widx1_1) = 0.11 * factt(widx1_1) * facth(widx1_1)
          end if
        enddo
        !$acc end kernels
      else
        if (aerotype == cp_bc) then
          !$acc kernels
          !$acc loop gang vector independent
          do idx_1 = LBOUND(cal_ho2, dim=1), UBOUND(cal_ho2, dim=1), 1
            cal_ho2(idx_1) = 0.01e0
          enddo
          !$acc end kernels
        else
          if (aerotype == cp_oc .OR. aerotype == cp_so) then
            !$acc kernels
            !$acc loop gang vector independent
            do idx_2 = LBOUND(cal_ho2, dim=1), UBOUND(cal_ho2, dim=1), 1
              cal_ho2(idx_2) = 2.5e-2
            enddo
            !$acc end kernels
          else
            if (aerotype == cp_cl) then
              !$acc kernels
              !$acc loop gang vector independent
              do widx1_2 = 1, nbox, 1
                if (rh_p(widx1_2) >= 62.0) then
                  cal_ho2(widx1_2) = 0.05e0
                else
                  cal_ho2(widx1_2) = 5.66e-5 * EXP(1560.0 / temp(widx1_2))
                end if
              enddo
              !$acc end kernels
            else
              if (aerotype == cp_no3) then
                !$acc kernels
                !$acc loop gang vector independent
                do idx_3 = LBOUND(cal_ho2, dim=1), UBOUND(cal_ho2, dim=1), 1
                  cal_ho2(idx_3) = 0.005
                enddo
                !$acc end kernels
              else
                if (aerotype == cp_nn) then
                  !$acc kernels
                  !$acc loop gang vector independent
                  do idx_4 = LBOUND(cal_ho2, dim=1), UBOUND(cal_ho2, dim=1), 1
                    cal_ho2(idx_4) = 0.015
                  enddo
                  !$acc end kernels
                else
                  cmessage(:) = 'Unknown aerosol surface for HO2 self reaction'
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A)') cmessage, ' aerotype=', aerotype
                  call umprint(ummessage, src='ukca_trop_hetchem')
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,I0)') 'AEROSOL TYPE =', aerotype
                  call umprint(ummessage, src='ukca_trop_hetchem')
                  errcode = aerotype
                  call ereport('UKCA_TROP_HETCHEM_MOD:CAL_HO2', errcode, cmessage)
                end if
              end if
            end if
          end if
        end if
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end function cal_ho2
  function cal_loss1k(nbox, sarea, wetdp, aird, mode_gamma, temp, molec_wt_s)
    real, parameter :: cf = 1.0e-2
    real, parameter :: dq = 4.5e-10
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nbox
    real, dimension(nbox), intent(in) :: sarea
    real, dimension(nbox), intent(in) :: wetdp
    real, dimension(nbox), intent(in) :: aird
    real, dimension(nbox), intent(in) :: mode_gamma
    real, dimension(nbox), intent(in) :: temp
    real, intent(in) :: molec_wt_s
    real, dimension(nbox) :: cal_loss1k
    real :: cal_loss1k_min
    real :: sarea_min
    real, dimension(nbox) :: dfkg
    real, dimension(nbox) :: rho_a
    real :: mw_s
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'CAL_LOSS1K'
    integer :: widx1
    integer :: idx

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    mw_s = molec_wt_s / 1.0e3
    rho_a(:) = aird(:) * m_air * 1.0e6 / avogadro
    if (printstatus >= prstatus_diag) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,2E15.6)') 'RHO_A: ', MAXVAL(rho_a), MINVAL(rho_a)
      call umprint(ummessage, src='ukca_trop_hetchem')
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(dfkg, dim=1), UBOUND(dfkg, dim=1), 1
      dfkg(idx) = 0.0
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_classic_hetchem) then
      sarea_min = 1.0e-30
      cal_loss1k_min = 0.0
    else
      sarea_min = 0.0
      cal_loss1k_min = 1.0e-3
    end if
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (sarea < sarea_min .OR. wetdp < 1.0E-30 .OR. mode_gamma < 1.0E-30)
      cal_loss1k = cal_loss1k_min
    ELSEWHERE
      dfkg(:) = (3.0E4 / (8.0 * avogadro * dq ** 2 * rho_a(:))) * SQRT((rmol * temp(:) * m_air / (2.0 * pi)) * (mw_s + m_air) / mw_s)
      cal_loss1k(:) = sarea(:) / (wetdp(:) / dfkg(:) + (cf / mode_gamma(:)) * SQRT(2.0 * pi * mw_s / (rmol * temp(:))))
    END WHERE
    if (printstatus >= prstatus_diag) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,2E15.6)') 'dfkg: ', MAXVAL(dfkg), MINVAL(dfkg)
      call umprint(ummessage, src='ukca_trop_hetchem')
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end function cal_loss1k

end module ukca_trop_hetchem_mod
