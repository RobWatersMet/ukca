module ukca_stratf_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_STRATF_MOD'
  public

  contains
  subroutine ukca_stratf(row_length, rows, model_levels, ntracer, o33d, tracer)
    use ukca_config_specification_mod, only : ukca_config
    use ukca_tropopause, only : tropopause_level
    use ukca_cspecies, only : n_hono2, n_o3, n_o3s
    use ukca_constants, only : c_hno3, c_n
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_diag, ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    real, parameter :: o3_hno3_ratio = 1.0 / 1000.0
    integer, parameter :: no_above_trop1 = 3
    integer, parameter :: no_above_trop2 = 10
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: ntracer
    real, dimension(row_length,rows,model_levels), intent(in) :: o33d
    real, dimension(row_length,rows,model_levels,ntracer), intent(inout) :: tracer
    integer :: l
    logical, dimension(row_length,rows,model_levels) :: mask
    real, dimension(row_length,rows,model_levels) :: hno33d
    integer :: no_above_trop
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_STRATF'
    integer :: widx3
    integer :: widx2
    integer :: widx1
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: widx3_1
    integer :: widx2_1
    integer :: widx1_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    if (ukca_config%nlev_above_trop_o3_env > 0) then
      no_above_trop = ukca_config%nlev_above_trop_o3_env
    else
      if (model_levels == 38 .OR. model_levels == 60) then
        no_above_trop = no_above_trop1
      else
        if (model_levels == 63 .OR. model_levels == 70 .OR. model_levels == 85) then
          no_above_trop = no_above_trop2
        else
          errcode = 1
          cmessage(:) = 'Levels above tropopause must be set explicitly at this resolution'
          call ereport(routinename, errcode, cmessage)
        end if
      end if
    end if
    if (printstatus == prstatus_diag) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'UKCA_STRATF:'
      call umprint(ummessage, src='ukca_stratf')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I0)') 'no_above_trop: ', no_above_trop
      call umprint(ummessage, src='ukca_stratf')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') ' '
      call umprint(ummessage, src='ukca_stratf')
    end if
    if (no_above_trop >= model_levels) then
      if (lhook) then
        call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
      end if
      return
    end if
    !$acc kernels
    hno33d(:,:,:) = o33d(:,:,:) * o3_hno3_ratio * c_hno3 / c_n
    !$acc loop gang vector independent collapse(3)
    do idx = LBOUND(mask, dim=3), UBOUND(mask, dim=3), 1
      do idx_1 = LBOUND(mask, dim=2), UBOUND(mask, dim=2), 1
        do idx_2 = LBOUND(mask, dim=1), UBOUND(mask, dim=1), 1
          mask(idx_2,idx_1,idx) = .false.
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do l = model_levels, 1, -1
      do idx_3 = LBOUND(mask, dim=2), UBOUND(mask, dim=2), 1
        do idx_4 = LBOUND(mask, dim=1), UBOUND(mask, dim=1), 1
          mask(idx_4,idx_3,l) = tropopause_level(idx_4 + (LBOUND(tropopause_level, dim=1) - LBOUND(mask, dim=1)),idx_3 + (LBOUND(tropopause_level, dim=2) - LBOUND(mask, dim=2))) + no_above_trop <= l
        enddo
      enddo
    enddo
    !$acc end kernels
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (mask(:, :, :))
      tracer(:, :, :, n_o3) = o33d(:, :, :)
      tracer(:, :, :, n_hono2) = hno33d(:, :, :)
    END WHERE
    if (n_o3s > 0) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (mask(:, :, :))
        tracer(:, :, :, n_o3s) = o33d(:, :, :)
      END WHERE
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_stratf

end module ukca_stratf_mod
