module ukca_chem_offline
  use ukca_chem_defs_mod, only : chch_t, ratb_t, rath_t, ratj_t, ratt_t
  implicit none
  integer, parameter, private :: naq_offline = 3
  integer, parameter, public :: ndry_offline = 6
  integer, parameter, public :: nwet_offline = 4
  integer, parameter, public :: nwet_constant = 1
  type(ratb_t), dimension(9), parameter, private :: ratb_defs_a = [ratb_t('DMS       ', 'OH        ', 'SO2       ', '          ', '          ', '          ', 9.60E-12, 0.00, 240.00, 0.000, 0.000, 0.000, 0.000), ratb_t('DMS       ', 'OH        ', 'SO2       ', 'DMSO      ', '          ', '          ', 3.04E-12, 0.00, - 350.00, 0.600, 0.400, 0.000, 0.000), ratb_t('DMS       ', 'NO3       ', 'SO2       ', '          ', '          ', '          ', 1.90E-13, 0.00, - 500.00, 0.000, 0.000, 0.000, 0.000), ratb_t('DMSO      ', 'OH        ', 'SO2       ', '          ', '          ', '          ', 5.80E-11, 0.00, 0.00, 0.600, 0.000, 0.000, 0.000), ratb_t('Monoterp  ', 'OH        ', 'Sec_Org   ', '          ', '          ', '          ', 1.20E-11, 0.00, - 444.00, 0.130, 0.000, 0.000, 0.000), ratb_t('Monoterp  ', 'O3        ', 'Sec_Org   ', '          ', '          ', '          ', 1.01E-15, 0.00, 732.00, 0.130, 0.000, 0.000, 0.000), ratb_t('Monoterp  ', 'NO3       ', 'Sec_Org   ', '          ', '          ', '          ', 1.19E-12, 0.00, - 925.00, 0.130, 0.000, 0.000, 0.000), ratb_t('HO2       ', 'HO2       ', 'H2O2      ', '          ', '          ', '          ', 2.20E-13, 0.00, - 600.00, 0.000, 0.000, 0.000, 0.000), ratb_t('OH        ', 'H2O2      ', 'H2O       ', '          ', '          ', '          ', 2.90E-12, 0.00, 160.00, 0.000, 0.000, 0.000, 0.000)]
  type(chch_t), allocatable, dimension(:), public :: chch_defs_offline
  type(ratb_t), allocatable, dimension(:), public :: ratb_defs_offline
  type(ratt_t), allocatable, dimension(:), public :: ratt_defs_offline
  type(ratj_t), allocatable, dimension(:), public :: ratj_defs_offline
  type(rath_t), allocatable, dimension(:), public :: rath_defs_offline
  real, allocatable, dimension(:,:,:), public :: depvel_defs_offline
  real, allocatable, dimension(:,:), public :: henry_defs_offline
  integer, private :: n_tracers
  real, allocatable, dimension(:,:), save, public :: o3_offline_diag
  real, allocatable, dimension(:,:), save, public :: oh_offline_diag
  real, allocatable, dimension(:,:), save, public :: no3_offline_diag
  real, allocatable, dimension(:,:), save, public :: ho2_offline_diag
  type(chch_t), dimension(11), save, private :: chch_defs_a = [chch_t(1, 'H2O2      ', 1, 'TR        ', '          ', 1, 1), chch_t(2, 'DMS       ', 1, 'TR        ', '          ', 0, 0), chch_t(3, 'SO2       ', 1, 'TR        ', '          ', 1, 1), chch_t(4, 'H2SO4     ', 1, 'TR        ', '          ', 1, 0), chch_t(5, 'DMSO      ', 1, 'TR        ', '          ', 1, 1), chch_t(6, 'Monoterp  ', 1, 'TR        ', '          ', 1, 0), chch_t(7, 'Sec_Org   ', 1, 'TR        ', '          ', 1, 1), chch_t(8, 'O3        ', 1, 'CF        ', '          ', 0, 0), chch_t(9, 'OH        ', 1, 'CF        ', '          ', 0, 0), chch_t(10, 'NO3       ', 1, 'CF        ', '          ', 0, 0), chch_t(11, 'HO2       ', 1, 'CF        ', '          ', 0, 0)]
  type(rath_t), dimension(naq_offline), save, private :: rath_defs_a = [rath_t('SO2       ', 'H2O2      ', 'NULL0     ', '          ', '          ', '          ', 0.000, 0.000, 0.000, 0.000), rath_t('SO2       ', 'O3        ', 'NULL1     ', '          ', '          ', '          ', 0.000, 0.000, 0.000, 0.000), rath_t('SO2       ', 'O3        ', 'NULL2     ', '          ', '          ', '          ', 0.000, 0.000, 0.000, 0.000)]
  type(ratt_t), dimension(1), save, private :: ratt_defs_a = [ratt_t('SO2       ', 'OH        ', 'H2SO4     ', '          ', 0.6, 3.00E-31, - 3.30, 0.0, 1.50E-12, 0.00, 0.0, 0.000, 0.000)]
  real, dimension(6,5,ndry_offline), save, private :: depvel_defs_a = RESHAPE([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.25, 0.16, 0.71, 0.28, 0.12, 0.20, 1.25, 0.53, 0.89, 0.83, 0.78, 0.81, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10], [6, 5, ndry_offline])
  real, dimension(6,nwet_offline), save, private :: henry_defs_a = RESHAPE([0.830E+05, 0.740E+04, 0.240E-11, - 0.373E+04, 0.000E+00, 0.000E+00, 0.123E+01, 0.302E+04, 0.123E-01, 0.201E+04, 0.600E-07, 0.112E+04, 0.500E+05, 0.6425E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.100E+06, 0.120E+02, 0.000E-00, 0.000E+00, 0.000E+00, 0.000E+00], [6, nwet_offline])
  real, dimension(6,nwet_constant), save, public :: henry_defs_const = RESHAPE([0.113E-01, 0.230E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00], [6, nwet_constant])
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEM_OFFLINE'
  private

  public :: ukca_init_offline

  contains
  subroutine ukca_init_offline()
    use asad_mod, only : ct_ddhr, ct_dhr, ct_k298, ct_kd298, ih_o3_const, jddepc, jddept, jpbk, jpctr, jpdd, jpdw, jpeq, jphk, jppj, jpsp, jpspec, jptk
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    integer :: i
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_INIT_OFFLINE'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    ih_o3_const = 1
    n_tracers = 0
    !$acc loop seq
    do i = 1, SIZE(chch_defs_a), 1
      if (chch_defs_a(i)%ctype == jpsp) then
        n_tracers = n_tracers + 1
      end if
    enddo
    !$acc end kernels
    if (n_tracers /= jpctr) then
      cmessage(:) = 'No. of tracers found not equal to jpctr'
      errcode = 1
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A40,A12,I5,A9,I5)') cmessage, ' n_tracers:', n_tracers, ' jpctr: ', jpctr
      call umprint(ummessage, src='UKCA_INIT_OFFLINE')
      call ereport('UKCA_INIT_OFFLINE', errcode, cmessage)
    end if
    ALLOCATE(chch_defs_offline(1:jpspec))
    ALLOCATE(ratb_defs_offline(1:jpbk))
    ALLOCATE(ratj_defs_offline(1:jppj))
    ALLOCATE(ratt_defs_offline(1:jptk))
    ALLOCATE(rath_defs_offline(1:jphk))
    ALLOCATE(depvel_defs_offline(1:jddept,1:jddepc,1:jpdd))
    ALLOCATE(henry_defs_offline(1:6,1:jpdw))
    if (nwet_constant > 0) then
      if (.NOT.ALLOCATED(ct_k298)) then
        ALLOCATE(ct_k298(1:nwet_constant))
      end if
      if (.NOT.ALLOCATED(ct_dhr)) then
        ALLOCATE(ct_dhr(1:nwet_constant))
      end if
      if (.NOT.ALLOCATED(ct_kd298)) then
        ALLOCATE(ct_kd298(1:nwet_constant,1:jpeq))
      end if
      if (.NOT.ALLOCATED(ct_ddhr)) then
        ALLOCATE(ct_ddhr(1:nwet_constant,1:jpeq))
      end if
    end if
    chch_defs_offline(:) = [chch_defs_a]
    ratb_defs_offline(:) = [ratb_defs_a]
    ratt_defs_offline(:) = [ratt_defs_a]
    rath_defs_offline(:) = [rath_defs_a]
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do idx = LBOUND(depvel_defs_offline, dim=3), UBOUND(depvel_defs_offline, dim=3), 1
      do idx_1 = LBOUND(depvel_defs_offline, dim=2), UBOUND(depvel_defs_offline, dim=2), 1
        do idx_2 = LBOUND(depvel_defs_offline, dim=1), UBOUND(depvel_defs_offline, dim=1), 1
          depvel_defs_offline(idx_2,idx_1,idx) = depvel_defs_a(idx_2 + (1 - LBOUND(depvel_defs_offline, dim=1)),idx_1 + (1 - LBOUND(depvel_defs_offline, dim=2)),idx + (1 - LBOUND(depvel_defs_offline, dim=3)))
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_3 = LBOUND(henry_defs_offline, dim=2), UBOUND(henry_defs_offline, dim=2), 1
      do idx_4 = LBOUND(henry_defs_offline, dim=1), UBOUND(henry_defs_offline, dim=1), 1
        henry_defs_offline(idx_4,idx_3) = henry_defs_a(idx_4 + (1 - LBOUND(henry_defs_offline, dim=1)),idx_3 + (1 - LBOUND(henry_defs_offline, dim=2)))
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_init_offline

end module ukca_chem_offline
