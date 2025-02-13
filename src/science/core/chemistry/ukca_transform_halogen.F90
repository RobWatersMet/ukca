module ukca_transform_halogen_mod
  use ukca_constants, only : c_br, c_brcl, c_bro, c_brono2, c_ccl4, c_cf2cl2, c_cf2clbr, c_cf2clcfcl2, c_cf3br, c_cfcl3, c_ch2br2, c_chf2cl, c_cl, c_cl2o2, c_clo, c_clono2, c_hbr, c_hcl, c_ho2no2, c_hobr, c_hocl, c_hono, c_hono2, c_ison, c_mebr, c_meccl3, c_mecl, c_meono2, c_mpan, c_n, c_n2o5, c_no, c_no2, c_no3, c_oclo, c_pan, c_ppan
  implicit none
  integer, parameter, private :: unrestricted = 0
  integer, parameter, private :: source_gas = 1
  integer, parameter, private :: mixed_compound = 2
  integer, parameter, private :: tr_max = 30
  integer, save, private :: n_hcl = 0
  integer, save, private :: n_bro = 0
  integer, save, private :: n_no2 = 0
  integer, save, private :: n_cl_tracers = 0
  integer, save, private :: n_br_tracers = 0
  integer, save, private :: n_n_tracers = 0
  integer, dimension(tr_max), save, private :: cl_tracers = 0
  integer, dimension(tr_max), save, private :: br_tracers = 0
  integer, dimension(tr_max), save, private :: n_tracers = 0
  integer, dimension(tr_max), save, private :: s_cl_tracers = unrestricted
  integer, dimension(tr_max), save, private :: s_br_tracers = unrestricted
  integer, dimension(tr_max), save, private :: s_n_tracers = unrestricted
  real, dimension(tr_max), save, private :: c_cl_tracers
  real, dimension(tr_max), save, private :: c_br_tracers
  real, dimension(tr_max), save, private :: c_n_tracers
  character(len = *), parameter, private :: ModuleName = 'UKCA_TRANSFORM_HALOGEN_MOD'
  private

  public :: ukca_transform_halogen

  contains
  subroutine ukca_transform_halogen(rows, row_length, model_levels, tr_ukca, forward, tracers)
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: rows
    integer, intent(in) :: row_length
    integer, intent(in) :: model_levels
    integer, intent(in) :: tr_ukca
    logical, intent(in) :: forward
    real, dimension(row_length,rows,model_levels,tr_ukca), intent(inout) :: tracers
    logical, save :: first = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_TRANSFORM_HALOGEN'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first) then
      call ukca_setup_lumping()
      first = .false.
    end if
    if (forward) then
      call transform(row_length, rows, model_levels, tr_ukca, n_n_tracers, n_tracers, c_n_tracers, s_n_tracers, n_no2, c_no2, forward, tracers)
      call transform(row_length, rows, model_levels, tr_ukca, n_cl_tracers, cl_tracers, c_cl_tracers, s_cl_tracers, n_hcl, c_hcl, forward, tracers)
      call transform(row_length, rows, model_levels, tr_ukca, n_br_tracers, br_tracers, c_br_tracers, s_br_tracers, n_bro, c_bro, forward, tracers)
    else
      call transform(row_length, rows, model_levels, tr_ukca, n_br_tracers, br_tracers, c_br_tracers, s_br_tracers, n_bro, c_bro, forward, tracers)
      call transform(row_length, rows, model_levels, tr_ukca, n_cl_tracers, cl_tracers, c_cl_tracers, s_cl_tracers, n_hcl, c_hcl, forward, tracers)
      call transform(row_length, rows, model_levels, tr_ukca, n_n_tracers, n_tracers, c_n_tracers, s_n_tracers, n_no2, c_no2, forward, tracers)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_transform_halogen
  subroutine transform(row_length, rows, model_levels, tr_ukca, n_tracers, iposition, c_tracers, s_tracers, major, c_major, forward, tracers)
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    real, parameter :: rafeps = 1e-150
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: tr_ukca
    integer, intent(in) :: n_tracers
    integer, intent(in) :: major
    real, intent(in) :: c_major
    logical, intent(in) :: forward
    integer, dimension(n_tracers), intent(in) :: iposition
    real, dimension(n_tracers), intent(in) :: c_tracers
    integer, dimension(n_tracers), intent(in) :: s_tracers
    real, dimension(row_length,rows,model_levels,tr_ukca), intent(inout) :: tracers
    integer :: i
    real, allocatable, dimension(:,:,:) :: total_all
    real, allocatable, dimension(:,:,:) :: total_source
    real, allocatable, dimension(:,:,:) :: total_nochange
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'TRANSFORM'
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
    integer :: widx3
    integer :: widx2
    integer :: widx1
    integer :: widx3_1
    integer :: widx2_1
    integer :: widx1_1
    integer :: widx3_2
    integer :: widx2_2
    integer :: widx1_2
    integer :: widx3_1_1
    integer :: widx2_1_1
    integer :: widx1_1_1
    integer :: widx3_1_2
    integer :: widx2_1_2
    integer :: widx1_1_2
    integer :: widx3_2_1
    integer :: widx2_2_1
    integer :: widx1_2_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    ALLOCATE(total_all(1:row_length,1:rows,1:model_levels))
    ALLOCATE(total_source(1:row_length,1:rows,1:model_levels))
    ALLOCATE(total_nochange(1:row_length,1:rows,1:model_levels))
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx = LBOUND(total_all, dim=3), UBOUND(total_all, dim=3), 1
      do idx_1 = LBOUND(total_all, dim=2), UBOUND(total_all, dim=2), 1
        do idx_2 = LBOUND(total_all, dim=1), UBOUND(total_all, dim=1), 1
          total_all(idx_2,idx_1,idx) = 0.0
        enddo
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx_3 = LBOUND(total_source, dim=3), UBOUND(total_source, dim=3), 1
      do idx_4 = LBOUND(total_source, dim=2), UBOUND(total_source, dim=2), 1
        do idx_5 = LBOUND(total_source, dim=1), UBOUND(total_source, dim=1), 1
          total_source(idx_5,idx_4,idx_3) = 0.0
        enddo
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx_6 = LBOUND(total_nochange, dim=3), UBOUND(total_nochange, dim=3), 1
      do idx_7 = LBOUND(total_nochange, dim=2), UBOUND(total_nochange, dim=2), 1
        do idx_8 = LBOUND(total_nochange, dim=1), UBOUND(total_nochange, dim=1), 1
          total_nochange(idx_8,idx_7,idx_6) = 0.0
        enddo
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    !$acc loop seq
    do i = 1, n_tracers, 1
      if (s_tracers(i) == unrestricted) then
        !$acc loop gang vector independent collapse(3)
        do idx_9 = LBOUND(total_all, dim=3), UBOUND(total_all, dim=3), 1
          do idx_10 = LBOUND(total_all, dim=2), UBOUND(total_all, dim=2), 1
            do idx_11 = LBOUND(total_all, dim=1), UBOUND(total_all, dim=1), 1
              total_all(idx_11,idx_10,idx_9) = total_all(idx_11,idx_10,idx_9) + tracers(idx_11 + (LBOUND(tracers, dim=1) - LBOUND(total_all, dim=1)),idx_10 + (LBOUND(tracers, dim=2) - LBOUND(total_all, dim=2)),idx_9 + (LBOUND(tracers, dim=3) - LBOUND(total_all, dim=3)),iposition(i)) / c_tracers(i)
            enddo
          enddo
        enddo
      else
        if (s_tracers(i) == source_gas) then
          !$acc loop gang vector independent collapse(3)
          do idx_12 = LBOUND(total_source, dim=3), UBOUND(total_source, dim=3), 1
            do idx_13 = LBOUND(total_source, dim=2), UBOUND(total_source, dim=2), 1
              do idx_14 = LBOUND(total_source, dim=1), UBOUND(total_source, dim=1), 1
                total_source(idx_14,idx_13,idx_12) = total_source(idx_14,idx_13,idx_12) + tracers(idx_14 + (LBOUND(tracers, dim=1) - LBOUND(total_source, dim=1)),idx_13 + (LBOUND(tracers, dim=2) - LBOUND(total_source, dim=2)),idx_12 + (LBOUND(tracers, dim=3) - LBOUND(total_source, dim=3)),iposition(i)) / c_tracers(i)
              enddo
            enddo
          enddo
        else
          if (s_tracers(i) == mixed_compound) then
            !$acc loop gang vector independent collapse(3)
            do idx_15 = LBOUND(total_nochange, dim=3), UBOUND(total_nochange, dim=3), 1
              do idx_16 = LBOUND(total_nochange, dim=2), UBOUND(total_nochange, dim=2), 1
                do idx_17 = LBOUND(total_nochange, dim=1), UBOUND(total_nochange, dim=1), 1
                  total_nochange(idx_17,idx_16,idx_15) = total_nochange(idx_17,idx_16,idx_15) + tracers(idx_17 + (LBOUND(tracers, dim=1) - LBOUND(total_nochange, dim=1)),idx_16 + (LBOUND(tracers, dim=2) - LBOUND(total_nochange, dim=2)),idx_15 + (LBOUND(tracers, dim=3) - LBOUND(total_nochange, dim=3)),iposition(i)) / c_tracers(i)
                enddo
              enddo
            enddo
          end if
        end if
      end if
    enddo
    !$acc end kernels
    if (forward) then
      !$acc kernels
      !$acc loop gang vector independent collapse(3)
      do idx_18 = LBOUND(tracers, dim=3), UBOUND(tracers, dim=3), 1
        do idx_19 = LBOUND(tracers, dim=2), UBOUND(tracers, dim=2), 1
          do idx_20 = LBOUND(tracers, dim=1), UBOUND(tracers, dim=1), 1
            tracers(idx_20,idx_19,idx_18,major) = tracers(idx_20,idx_19,idx_18,major) - c_major * total_nochange(idx_20 + (LBOUND(total_nochange, dim=1) - LBOUND(tracers, dim=1)),idx_19 + (LBOUND(total_nochange, dim=2) - LBOUND(tracers, dim=2)),idx_18 + (LBOUND(total_nochange, dim=3) - LBOUND(tracers, dim=3)))
          enddo
        enddo
      enddo
      !$acc loop seq
      do i = 1, n_tracers, 1
        if (s_tracers(i) < mixed_compound) then
          !$acc loop gang vector independent collapse(3)
          do widx3_2 = 1, SIZE(tracers, dim=3), 1
            do widx2_2 = 1, SIZE(tracers, dim=2), 1
              do widx1_2 = 1, SIZE(tracers, dim=1), 1
                if (tracers(widx1_2,widx2_2,widx3_2,major) < rafeps) then
                  tracers(widx1_2,widx2_2,widx3_2,iposition(i)) = rafeps
                end if
              enddo
            enddo
          enddo
        end if
      enddo
      !$acc loop gang vector independent collapse(3)
      do widx3 = 1, SIZE(tracers, dim=3), 1
        do widx2 = 1, SIZE(tracers, dim=2), 1
          do widx1 = 1, SIZE(tracers, dim=1), 1
            if (tracers(widx1,widx2,widx3,major) < rafeps) then
              tracers(widx1,widx2,widx3,major) = rafeps
              total_source(LBOUND(total_source, dim=1) + widx1 - 1,LBOUND(total_source, dim=2) + widx2 - 1,LBOUND(total_source, dim=3) + widx3 - 1) = rafeps
              total_all(LBOUND(total_all, dim=1) + widx1 - 1,LBOUND(total_all, dim=2) + widx2 - 1,LBOUND(total_all, dim=3) + widx3 - 1) = rafeps
            end if
          enddo
        enddo
      enddo
      !$acc loop gang vector independent collapse(3)
      do idx_21 = LBOUND(total_nochange, dim=3), UBOUND(total_nochange, dim=3), 1
        do idx_22 = LBOUND(total_nochange, dim=2), UBOUND(total_nochange, dim=2), 1
          do idx_23 = LBOUND(total_nochange, dim=1), UBOUND(total_nochange, dim=1), 1
            total_nochange(idx_23,idx_22,idx_21) = total_source(idx_23 + (LBOUND(total_source, dim=1) - LBOUND(total_nochange, dim=1)),idx_22 + (LBOUND(total_source, dim=2) - LBOUND(total_nochange, dim=2)),idx_21 + (LBOUND(total_source, dim=3) - LBOUND(total_nochange, dim=3))) * c_major / MAX(tracers(idx_23 + (LBOUND(tracers, dim=1) - LBOUND(total_nochange, dim=1)),idx_22 + (LBOUND(tracers, dim=2) - LBOUND(total_nochange, dim=2)),idx_21 + (LBOUND(tracers, dim=3) - LBOUND(total_nochange, dim=3)),major), rafeps)
          enddo
        enddo
      enddo
      !$acc loop seq
      do i = 1, n_tracers, 1
        if (s_tracers(i) == unrestricted) then
          !$acc loop gang vector independent collapse(3)
          do widx3_1_1 = 1, SIZE(total_nochange(:,:,:), dim=3), 1
            do widx2_1_1 = 1, SIZE(total_nochange(:,:,:), dim=2), 1
              do widx1_1_1 = 1, SIZE(total_nochange(:,:,:), dim=1), 1
                if (total_nochange(LBOUND(total_nochange, dim=1) + widx1_1_1 - 1,LBOUND(total_nochange, dim=2) + widx2_1_1 - 1,LBOUND(total_nochange, dim=3) + widx3_1_1 - 1) > 1.0) then
                  tracers(widx1_1_1,widx2_1_1,widx3_1_1,iposition(i)) = rafeps
                end if
              enddo
            enddo
          enddo
        else
          if (s_tracers(i) == source_gas) then
            !$acc loop gang vector independent collapse(3)
            do widx3_1_2 = 1, SIZE(total_nochange(:,:,:), dim=3), 1
              do widx2_1_2 = 1, SIZE(total_nochange(:,:,:), dim=2), 1
                do widx1_1_2 = 1, SIZE(total_nochange(:,:,:), dim=1), 1
                  if (total_nochange(LBOUND(total_nochange, dim=1) + widx1_1_2 - 1,LBOUND(total_nochange, dim=2) + widx2_1_2 - 1,LBOUND(total_nochange, dim=3) + widx3_1_2 - 1) > 1.0) then
                    tracers(widx1_1_2,widx2_1_2,widx3_1_2,iposition(i)) = tracers(widx1_1_2,widx2_1_2,widx3_1_2,iposition(i)) / total_nochange(LBOUND(total_nochange, dim=1) + widx1_1_2 - 1,LBOUND(total_nochange, dim=2) + widx2_1_2 - 1,LBOUND(total_nochange, dim=3) + widx3_1_2 - 1)
                  end if
                enddo
              enddo
            enddo
          end if
        end if
      enddo
      !$acc loop gang vector independent collapse(3)
      do widx3_1 = 1, SIZE(total_nochange(:,:,:), dim=3), 1
        do widx2_1 = 1, SIZE(total_nochange(:,:,:), dim=2), 1
          do widx1_1 = 1, SIZE(total_nochange(:,:,:), dim=1), 1
            if (total_nochange(LBOUND(total_nochange, dim=1) + widx1_1 - 1,LBOUND(total_nochange, dim=2) + widx2_1 - 1,LBOUND(total_nochange, dim=3) + widx3_1 - 1) > 1.0) then
              total_all(LBOUND(total_all, dim=1) + widx1_1 - 1,LBOUND(total_all, dim=2) + widx2_1 - 1,LBOUND(total_all, dim=3) + widx3_1 - 1) = rafeps
            end if
          enddo
        enddo
      enddo
      !$acc loop gang vector independent collapse(3)
      do idx_24 = LBOUND(tracers, dim=3), UBOUND(tracers, dim=3), 1
        do idx_25 = LBOUND(tracers, dim=2), UBOUND(tracers, dim=2), 1
          do idx_26 = LBOUND(tracers, dim=1), UBOUND(tracers, dim=1), 1
            tracers(idx_26,idx_25,idx_24,major) = MAX(tracers(idx_26,idx_25,idx_24,major) - c_major * total_source(idx_26 + (LBOUND(total_source, dim=1) - LBOUND(tracers, dim=1)),idx_25 + (LBOUND(total_source, dim=2) - LBOUND(tracers, dim=2)),idx_24 + (LBOUND(total_source, dim=3) - LBOUND(tracers, dim=3))), rafeps)
          enddo
        enddo
      enddo
      !$acc loop gang vector independent collapse(3)
      do idx_27 = LBOUND(total_nochange, dim=3), UBOUND(total_nochange, dim=3), 1
        do idx_28 = LBOUND(total_nochange, dim=2), UBOUND(total_nochange, dim=2), 1
          do idx_29 = LBOUND(total_nochange, dim=1), UBOUND(total_nochange, dim=1), 1
            total_nochange(idx_29,idx_28,idx_27) = total_all(idx_29 + (LBOUND(total_all, dim=1) - LBOUND(total_nochange, dim=1)),idx_28 + (LBOUND(total_all, dim=2) - LBOUND(total_nochange, dim=2)),idx_27 + (LBOUND(total_all, dim=3) - LBOUND(total_nochange, dim=3))) * c_major / MAX(tracers(idx_29 + (LBOUND(tracers, dim=1) - LBOUND(total_nochange, dim=1)),idx_28 + (LBOUND(tracers, dim=2) - LBOUND(total_nochange, dim=2)),idx_27 + (LBOUND(tracers, dim=3) - LBOUND(total_nochange, dim=3)),major), rafeps)
          enddo
        enddo
      enddo
      !$acc loop seq
      do i = 1, n_tracers, 1
        if (s_tracers(i) == unrestricted) then
          !$acc loop gang vector independent collapse(3)
          do widx3_2_1 = 1, SIZE(total_nochange(:,:,:), dim=3), 1
            do widx2_2_1 = 1, SIZE(total_nochange(:,:,:), dim=2), 1
              do widx1_2_1 = 1, SIZE(total_nochange(:,:,:), dim=1), 1
                if (total_nochange(LBOUND(total_nochange, dim=1) + widx1_2_1 - 1,LBOUND(total_nochange, dim=2) + widx2_2_1 - 1,LBOUND(total_nochange, dim=3) + widx3_2_1 - 1) > 1.0) then
                  tracers(widx1_2_1,widx2_2_1,widx3_2_1,iposition(i)) = tracers(widx1_2_1,widx2_2_1,widx3_2_1,iposition(i)) / total_nochange(LBOUND(total_nochange, dim=1) + widx1_2_1 - 1,LBOUND(total_nochange, dim=2) + widx2_2_1 - 1,LBOUND(total_nochange, dim=3) + widx3_2_1 - 1)
                end if
              enddo
            enddo
          enddo
        end if
      enddo
      !$acc loop gang vector independent collapse(3)
      do idx_30 = LBOUND(tracers, dim=3), UBOUND(tracers, dim=3), 1
        do idx_31 = LBOUND(tracers, dim=2), UBOUND(tracers, dim=2), 1
          do idx_32 = LBOUND(tracers, dim=1), UBOUND(tracers, dim=1), 1
            tracers(idx_32,idx_31,idx_30,major) = MAX(tracers(idx_32,idx_31,idx_30,major) - c_major * total_all(idx_32 + (LBOUND(total_all, dim=1) - LBOUND(tracers, dim=1)),idx_31 + (LBOUND(total_all, dim=2) - LBOUND(tracers, dim=2)),idx_30 + (LBOUND(total_all, dim=3) - LBOUND(tracers, dim=3))), rafeps)
          enddo
        enddo
      enddo
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop gang vector independent collapse(3)
      do idx_33 = LBOUND(tracers, dim=3), UBOUND(tracers, dim=3), 1
        do idx_34 = LBOUND(tracers, dim=2), UBOUND(tracers, dim=2), 1
          do idx_35 = LBOUND(tracers, dim=1), UBOUND(tracers, dim=1), 1
            tracers(idx_35,idx_34,idx_33,major) = tracers(idx_35,idx_34,idx_33,major) + c_major * (total_all(idx_35 + (LBOUND(total_all, dim=1) - LBOUND(tracers, dim=1)),idx_34 + (LBOUND(total_all, dim=2) - LBOUND(tracers, dim=2)),idx_33 + (LBOUND(total_all, dim=3) - LBOUND(tracers, dim=3))) + total_source(idx_35 + (LBOUND(total_source, dim=1) - LBOUND(tracers, dim=1)),idx_34 + (LBOUND(total_source, dim=2) - LBOUND(tracers, dim=2)),idx_33 + (LBOUND(total_source, dim=3) - LBOUND(tracers, dim=3))) + total_nochange(idx_35 + (LBOUND(total_nochange, dim=1) - LBOUND(tracers, dim=1)),idx_34 + (LBOUND(total_nochange, dim=2) - LBOUND(tracers, dim=2)),idx_33 + (LBOUND(total_nochange, dim=3) - LBOUND(tracers, dim=3))))
          enddo
        enddo
      enddo
      !$acc end kernels
    end if
    DEALLOCATE(total_all)
    DEALLOCATE(total_source)
    DEALLOCATE(total_nochange)
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine transform
  subroutine ukca_setup_lumping()
    use asad_mod, only : advt, jpctr
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: i
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SETUP_LUMPING'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    do i = 1, jpctr, 1
      if (advt(i) == 'HCl       ') then
        n_hcl = i
      else
        if (advt(i) == 'Cl        ') then
          n_cl_tracers = n_cl_tracers + 1
          c_cl_tracers(n_cl_tracers) = c_cl
          cl_tracers(n_cl_tracers) = i
        else
          if (advt(i) == 'ClO       ' .OR. advt(i) == 'Clx       ') then
            n_cl_tracers = n_cl_tracers + 1
            c_cl_tracers(n_cl_tracers) = c_clo
            cl_tracers(n_cl_tracers) = i
          else
            if (advt(i) == 'Cl2O2     ') then
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cl2o2 / 2.0
              cl_tracers(n_cl_tracers) = i
            else
              if (advt(i) == 'HOCl      ') then
                n_cl_tracers = n_cl_tracers + 1
                c_cl_tracers(n_cl_tracers) = c_hocl
                cl_tracers(n_cl_tracers) = i
              else
                if (advt(i) == 'OClO      ') then
                  n_cl_tracers = n_cl_tracers + 1
                  c_cl_tracers(n_cl_tracers) = c_oclo
                  cl_tracers(n_cl_tracers) = i
                else
                  if (advt(i) == 'ClONO2    ') then
                    n_cl_tracers = n_cl_tracers + 1
                    c_cl_tracers(n_cl_tracers) = c_clono2
                    cl_tracers(n_cl_tracers) = i
                    n_n_tracers = n_n_tracers + 1
                    c_n_tracers(n_n_tracers) = c_clono2
                    n_tracers(n_n_tracers) = i
                    s_n_tracers(n_n_tracers) = mixed_compound
                  else
                    if (advt(i) == 'CFCl3     ') then
                      n_cl_tracers = n_cl_tracers + 1
                      c_cl_tracers(n_cl_tracers) = c_cfcl3 / 3.0
                      cl_tracers(n_cl_tracers) = i
                      s_cl_tracers(n_cl_tracers) = source_gas
                    else
                      if (advt(i) == 'CF2Cl2    ') then
                        n_cl_tracers = n_cl_tracers + 1
                        c_cl_tracers(n_cl_tracers) = c_cf2cl2 / 2.0
                        cl_tracers(n_cl_tracers) = i
                        s_cl_tracers(n_cl_tracers) = source_gas
                      else
                        if (advt(i) == 'MeCl      ') then
                          n_cl_tracers = n_cl_tracers + 1
                          c_cl_tracers(n_cl_tracers) = c_mecl
                          cl_tracers(n_cl_tracers) = i
                          s_cl_tracers(n_cl_tracers) = source_gas
                        else
                          if (advt(i) == 'CCl4      ') then
                            n_cl_tracers = n_cl_tracers + 1
                            c_cl_tracers(n_cl_tracers) = c_ccl4 / 4.0
                            cl_tracers(n_cl_tracers) = i
                            s_cl_tracers(n_cl_tracers) = source_gas
                          else
                            if (advt(i) == 'CF2ClCFCl2') then
                              n_cl_tracers = n_cl_tracers + 1
                              c_cl_tracers(n_cl_tracers) = c_cf2clcfcl2 / 3.0
                              cl_tracers(n_cl_tracers) = i
                              s_cl_tracers(n_cl_tracers) = source_gas
                            else
                              if (advt(i) == 'CHF2Cl    ') then
                                n_cl_tracers = n_cl_tracers + 1
                                c_cl_tracers(n_cl_tracers) = c_chf2cl
                                cl_tracers(n_cl_tracers) = i
                                s_cl_tracers(n_cl_tracers) = source_gas
                              else
                                if (advt(i) == 'MeCCl3    ') then
                                  n_cl_tracers = n_cl_tracers + 1
                                  c_cl_tracers(n_cl_tracers) = c_meccl3 / 3.0
                                  cl_tracers(n_cl_tracers) = i
                                  s_cl_tracers(n_cl_tracers) = source_gas
                                else
                                  if (advt(i) == 'BrCl      ') then
                                    n_cl_tracers = n_cl_tracers + 1
                                    c_cl_tracers(n_cl_tracers) = c_brcl
                                    cl_tracers(n_cl_tracers) = i
                                    s_cl_tracers(n_cl_tracers) = mixed_compound
                                    n_br_tracers = n_br_tracers + 1
                                    c_br_tracers(n_br_tracers) = c_brcl
                                    br_tracers(n_br_tracers) = i
                                  else
                                    if (advt(i) == 'CF2ClBr   ') then
                                      n_cl_tracers = n_cl_tracers + 1
                                      c_cl_tracers(n_cl_tracers) = c_cf2clbr
                                      cl_tracers(n_cl_tracers) = i
                                      s_cl_tracers(n_cl_tracers) = mixed_compound
                                      n_br_tracers = n_br_tracers + 1
                                      c_br_tracers(n_br_tracers) = c_cf2clbr
                                      br_tracers(n_br_tracers) = i
                                      s_br_tracers(n_br_tracers) = source_gas
                                    else
                                      if (advt(i) == 'BrO       ' .OR. advt(i) == 'Brx       ') then
                                        n_bro = i
                                      else
                                        if (advt(i) == 'Br        ') then
                                          n_br_tracers = n_br_tracers + 1
                                          c_br_tracers(n_br_tracers) = c_br
                                          br_tracers(n_br_tracers) = i
                                        else
                                          if (advt(i) == 'HBr       ') then
                                            n_br_tracers = n_br_tracers + 1
                                            c_br_tracers(n_br_tracers) = c_hbr
                                            br_tracers(n_br_tracers) = i
                                          else
                                            if (advt(i) == 'HOBr      ') then
                                              n_br_tracers = n_br_tracers + 1
                                              c_br_tracers(n_br_tracers) = c_hobr
                                              br_tracers(n_br_tracers) = i
                                            else
                                              if (advt(i) == 'BrONO2    ') then
                                                n_br_tracers = n_br_tracers + 1
                                                c_br_tracers(n_br_tracers) = c_brono2
                                                br_tracers(n_br_tracers) = i
                                                n_n_tracers = n_n_tracers + 1
                                                c_n_tracers(n_n_tracers) = c_brono2
                                                n_tracers(n_n_tracers) = i
                                                s_n_tracers(n_n_tracers) = mixed_compound
                                              else
                                                if (advt(i) == 'MeBr      ') then
                                                  n_br_tracers = n_br_tracers + 1
                                                  c_br_tracers(n_br_tracers) = c_mebr
                                                  br_tracers(n_br_tracers) = i
                                                  s_br_tracers(n_br_tracers) = source_gas
                                                else
                                                  if (advt(i) == 'CF3Br     ') then
                                                    n_br_tracers = n_br_tracers + 1
                                                    c_br_tracers(n_br_tracers) = c_cf3br
                                                    br_tracers(n_br_tracers) = i
                                                    s_br_tracers(n_br_tracers) = source_gas
                                                  else
                                                    if (advt(i) == 'CH2Br2    ') then
                                                      n_br_tracers = n_br_tracers + 1
                                                      c_br_tracers(n_br_tracers) = c_ch2br2 / 2.0
                                                      br_tracers(n_br_tracers) = i
                                                      s_br_tracers(n_br_tracers) = 1
                                                    else
                                                      if (advt(i) == 'N         ') then
                                                        n_n_tracers = n_n_tracers + 1
                                                        c_n_tracers(n_n_tracers) = c_n
                                                        n_tracers(n_n_tracers) = i
                                                      else
                                                        if (advt(i) == 'NO        ') then
                                                          n_n_tracers = n_n_tracers + 1
                                                          c_n_tracers(n_n_tracers) = c_no
                                                          n_tracers(n_n_tracers) = i
                                                        else
                                                          if (advt(i) == 'NO2       ' .OR. advt(i) == 'NOx       ') then
                                                            n_no2 = i
                                                          else
                                                            if (advt(i) == 'NO3       ') then
                                                              n_n_tracers = n_n_tracers + 1
                                                              c_n_tracers(n_n_tracers) = c_no3
                                                              n_tracers(n_n_tracers) = i
                                                            else
                                                              if (advt(i) == 'N2O5      ') then
                                                                n_n_tracers = n_n_tracers + 1
                                                                c_n_tracers(n_n_tracers) = c_n2o5 / 2.0
                                                                n_tracers(n_n_tracers) = i
                                                              else
                                                                if (advt(i) == 'HO2NO2    ') then
                                                                  n_n_tracers = n_n_tracers + 1
                                                                  c_n_tracers(n_n_tracers) = c_ho2no2
                                                                  n_tracers(n_n_tracers) = i
                                                                else
                                                                  if (advt(i) == 'HONO2     ') then
                                                                    n_n_tracers = n_n_tracers + 1
                                                                    c_n_tracers(n_n_tracers) = c_hono2
                                                                    n_tracers(n_n_tracers) = i
                                                                  else
                                                                    if (advt(i) == 'HONO      ') then
                                                                      n_n_tracers = n_n_tracers + 1
                                                                      c_n_tracers(n_n_tracers) = c_hono
                                                                      n_tracers(n_n_tracers) = i
                                                                    else
                                                                      if (advt(i) == 'MeONO2    ') then
                                                                        n_n_tracers = n_n_tracers + 1
                                                                        c_n_tracers(n_n_tracers) = c_meono2
                                                                        n_tracers(n_n_tracers) = i
                                                                      else
                                                                        if (advt(i) == 'PAN       ') then
                                                                          n_n_tracers = n_n_tracers + 1
                                                                          c_n_tracers(n_n_tracers) = c_pan
                                                                          n_tracers(n_n_tracers) = i
                                                                        else
                                                                          if (advt(i) == 'PPAN      ') then
                                                                            n_n_tracers = n_n_tracers + 1
                                                                            c_n_tracers(n_n_tracers) = c_ppan
                                                                            n_tracers(n_n_tracers) = i
                                                                          else
                                                                            if (advt(i) == 'ISON      ') then
                                                                              n_n_tracers = n_n_tracers + 1
                                                                              c_n_tracers(n_n_tracers) = c_ison
                                                                              n_tracers(n_n_tracers) = i
                                                                            else
                                                                              if (advt(i) == 'MPAN      ') then
                                                                                n_n_tracers = n_n_tracers + 1
                                                                                c_n_tracers(n_n_tracers) = c_mpan
                                                                                n_tracers(n_n_tracers) = i
                                                                              else
                                                                                ummessage(:) = 'Species: ' // advt(i) // ' not included in CASE'
                                                                                if (printstatus > prstatus_oper) then
                                                                                  call umprint(ummessage, src=ModuleName // ':' // RoutineName)
                                                                                end if
                                                                              end if
                                                                            end if
                                                                          end if
                                                                        end if
                                                                      end if
                                                                    end if
                                                                  end if
                                                                end if
                                                              end if
                                                            end if
                                                          end if
                                                        end if
                                                      end if
                                                    end if
                                                  end if
                                                end if
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
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

  end subroutine ukca_setup_lumping

end module ukca_transform_halogen_mod
