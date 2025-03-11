module ukca_diurnal_oxidant
  use yomhook, only : dr_hook, lhook
  use nvtx
  use parkind1, only : jpim, jprb
  implicit none
  integer(kind=jpim), parameter, private :: zhook_in = 0
  integer(kind=jpim), parameter, private :: zhook_out = 1
  real, allocatable, dimension(:), save, public :: cosza
  real, allocatable, dimension(:), save, public :: intcosza
  real, allocatable, dimension(:), save, public :: daylength1d
  character(len = *), parameter, private :: ModuleName = 'UKCA_DIURNAL_OXIDANT'
  private

  public :: ukca_set_diurnal_ox, ukca_set_diurnal_ox_col, ukca_int_cosz, dealloc_diurnal_oxidant

  contains
  subroutine ukca_set_diurnal_ox(species, field, n_points, nlev)
    use ukca_chem_offline, only : ho2_offline_diag, no3_offline_diag, o3_offline_diag, oh_offline_diag
    use ukca_constants, only : pi, rsec_per_day
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    real, parameter :: daylen_cutoff = 3600.0
    integer, intent(in) :: n_points
    integer, intent(in) :: nlev
    CHARACTER(LEN = 10), INTENT(IN) :: species
    real, dimension(n_points), intent(inout) :: field
    real, dimension(n_points) :: daylen_secs
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SET_DIURNAL_OX'
    integer :: widx1
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: widx1_1
    integer :: widx1_2

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (species(1:4) == 'OH  ' .OR. species(1:4) == 'HO2 ') then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE ((cosza > 0.0) .AND. (intcosza > 1.0E-1))
        field(:) = field(:) * rsec_per_day * (cosza(:) / intcosza(:))
      ELSEWHERE
        field(:) = 0.0
      END WHERE
    else
      if (species(1:4) == 'NO3 ') then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx = 1, n_points, 1
          daylen_secs(idx) = daylength1d(idx + (LBOUND(daylength1d, dim=1) - 1)) * rsec_per_day / (2.0 * pi)
        enddo
        !NO!$acc end kernels
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
        !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
        WHERE (daylen_secs < daylen_cutoff)
          daylen_secs = rsec_per_day
        END WHERE
        WHERE ((cosza > 0.0) .AND. (intcosza > 1.0E-1) .AND. (daylen_secs < rsec_per_day))
          field(:) = 0.0
        END WHERE
      else
        if (species(1:4) == 'O3  ') then
        else
          cmessage(:) = 'ukca_set_diurnal_ox: Species: ' // species // ' not found'
          errcode = 1
          call ereport('ukca_set_diurnal_ox', errcode, cmessage)
        end if
      end if
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do widx1 = 1, n_points, 1
      if (field(widx1) < 0.0) then
        field(widx1) = 0.0
      end if
    enddo
    !NO!$acc end kernels
    if (species(1:4) == 'O3  ') then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_1 = LBOUND(o3_offline_diag, dim=1), UBOUND(o3_offline_diag, dim=1), 1
        o3_offline_diag(idx_1,nlev) = field(idx_1 + (LBOUND(field, dim=1) - LBOUND(o3_offline_diag, dim=1)))
      enddo
      !NO!$acc end kernels
    else
      if (species(1:4) == 'OH  ') then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_2 = LBOUND(oh_offline_diag, dim=1), UBOUND(oh_offline_diag, dim=1), 1
          oh_offline_diag(idx_2,nlev) = field(idx_2 + (LBOUND(field, dim=1) - LBOUND(oh_offline_diag, dim=1)))
        enddo
        !NO!$acc end kernels
      else
        if (species(1:4) == 'NO3 ') then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_3 = LBOUND(no3_offline_diag, dim=1), UBOUND(no3_offline_diag, dim=1), 1
            no3_offline_diag(idx_3,nlev) = field(idx_3 + (LBOUND(field, dim=1) - LBOUND(no3_offline_diag, dim=1)))
          enddo
          !NO!$acc end kernels
        else
          if (species(1:4) == 'HO2 ') then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_4 = LBOUND(ho2_offline_diag, dim=1), UBOUND(ho2_offline_diag, dim=1), 1
              ho2_offline_diag(idx_4,nlev) = field(idx_4 + (LBOUND(field, dim=1) - LBOUND(ho2_offline_diag, dim=1)))
            enddo
            !NO!$acc end kernels
          end if
        end if
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_set_diurnal_ox
  subroutine ukca_set_diurnal_ox_col(species, field, n_points, ix, jy)
    use ukca_config_specification_mod, only : ukca_config
    use ukca_chem_offline, only : ho2_offline_diag, no3_offline_diag, o3_offline_diag, oh_offline_diag
    use ukca_constants, only : pi, rsec_per_day
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    real, parameter :: daylen_cutoff = 3600.0
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    CHARACTER(LEN = 10), INTENT(IN) :: species
    real, dimension(n_points), intent(inout) :: field
    real :: daylen_secs
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    integer :: location
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SET_DIURNAL_OX_COL'
    integer :: widx1
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    location = ix + (jy - 1) * ukca_config%row_length
    if (species(1:4) == 'OH  ' .OR. species(1:4) == 'HO2 ') then
      if (cosza(location) > 0.0 .AND. intcosza(location) > 1.0e-1) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx = LBOUND(field, dim=1), UBOUND(field, dim=1), 1
          field(idx) = field(idx) * rsec_per_day * (cosza(location) / intcosza(location))
        enddo
        !NO!$acc end kernels
      else
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_1 = LBOUND(field, dim=1), UBOUND(field, dim=1), 1
          field(idx_1) = 0.0
        enddo
        !NO!$acc end kernels
      end if
    else
      if (species(1:4) == 'NO3 ') then
        daylen_secs = daylength1d(location) * rsec_per_day / (2.0 * pi)
        if (daylen_secs < daylen_cutoff) then
          daylen_secs = rsec_per_day
        end if
        if ((cosza(location) > 0.0 .AND. intcosza(location) > 1.0e-1) .AND. daylen_secs < rsec_per_day) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_2 = LBOUND(field, dim=1), UBOUND(field, dim=1), 1
            field(idx_2) = 0.0
          enddo
          !NO!$acc end kernels
        end if
      else
        if (species(1:4) == 'O3  ') then
        else
          cmessage(:) = 'ukca_set_diurnal_ox: Species: ' // species // ' not found'
          errcode = 1
          call ereport('ukca_set_diurnal_ox', errcode, cmessage)
        end if
      end if
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do widx1 = 1, n_points, 1
      if (field(widx1) < 0.0) then
        field(widx1) = 0.0
      end if
    enddo
    !NO!$acc end kernels
    if (species(1:4) == 'O3  ') then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_3 = LBOUND(o3_offline_diag, dim=2), UBOUND(o3_offline_diag, dim=2), 1
        o3_offline_diag(location,idx_3) = field(idx_3 + (LBOUND(field, dim=1) - LBOUND(o3_offline_diag, dim=2)))
      enddo
      !NO!$acc end kernels
    else
      if (species(1:4) == 'OH  ') then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_4 = LBOUND(oh_offline_diag, dim=2), UBOUND(oh_offline_diag, dim=2), 1
          oh_offline_diag(location,idx_4) = field(idx_4 + (LBOUND(field, dim=1) - LBOUND(oh_offline_diag, dim=2)))
        enddo
        !NO!$acc end kernels
      else
        if (species(1:4) == 'NO3 ') then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_5 = LBOUND(no3_offline_diag, dim=2), UBOUND(no3_offline_diag, dim=2), 1
            no3_offline_diag(location,idx_5) = field(idx_5 + (LBOUND(field, dim=1) - LBOUND(no3_offline_diag, dim=2)))
          enddo
          !NO!$acc end kernels
        else
          if (species(1:4) == 'HO2 ') then
            !NO!$acc kernels
            !NO!$acc loop gang vector independent
            do idx_6 = LBOUND(ho2_offline_diag, dim=2), UBOUND(ho2_offline_diag, dim=2), 1
              ho2_offline_diag(location,idx_6) = field(idx_6 + (LBOUND(field, dim=1) - LBOUND(ho2_offline_diag, dim=2)))
            enddo
            !NO!$acc end kernels
          end if
        end if
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_set_diurnal_ox_col
  subroutine ukca_int_cosz(row_length, rows, sin_latitude, cos_latitude, tan_latitude, sindec, cos_zenith_angle, int_zenith_angle)
    use ukca_constants, only : pi, rsec_per_day
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    real, dimension(row_length,rows), intent(in) :: sin_latitude
    real, dimension(row_length,rows), intent(in) :: cos_latitude
    real, dimension(row_length,rows), intent(in) :: tan_latitude
    real, intent(in) :: sindec
    real, dimension(row_length,rows), intent(in) :: cos_zenith_angle
    real, dimension(row_length,rows), intent(out) :: int_zenith_angle
    real, dimension(row_length,rows) :: daylength
    real :: decls
    real :: tan_decls
    real, dimension(row_length,rows) :: hs
    real, dimension(row_length,rows) :: cos_hs
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_INT_COSZ'
    integer :: widx2
    integer :: widx1
    integer :: widx2_1
    integer :: widx1_1
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    decls = ASIN(sindec)
    tan_decls = TAN(decls)
    if (ABS(tan_decls) < EPSILON(0.0)) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx = 1, rows, 1
        do idx_1 = 1, row_length, 1
          cos_hs(idx_1,idx) = 0.0
        enddo
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx_2 = 1, rows, 1
        do idx_3 = 1, row_length, 1
          cos_hs(idx_3,idx_2) = -tan_latitude(idx_3,idx_2) * tan_decls
        enddo
      enddo
      !NO!$acc end kernels
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(2)
    do widx2 = 1, rows, 1
      do widx1 = 1, row_length, 1
        if (ABS(cos_hs(widx1,widx2)) > 1.0 .AND. sin_latitude(widx1,widx2) * decls > 0.0) then
          cos_hs(widx1,widx2) = -0.99999999
        end if
      enddo
    enddo
    !NO!$acc end kernels
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (ABS(cos_hs) < 1.0)
      hs = ACOS(cos_hs)
      int_zenith_angle = (rsec_per_day / pi) * (hs * sin_latitude * sindec + cos_latitude * COS(decls) * SIN(hs))
      daylength = 2.0 * hs
    ELSEWHERE
      int_zenith_angle = 0.0
      daylength = 0.0
    END WHERE
    if (.NOT.ALLOCATED(intcosza)) then
      ALLOCATE(intcosza(1:row_length * rows))
    end if
    if (.NOT.ALLOCATED(cosza)) then
      ALLOCATE(cosza(1:row_length * rows))
    end if
    if (.NOT.ALLOCATED(daylength1d)) then
      ALLOCATE(daylength1d(1:row_length * rows))
    end if
    intcosza(:) = RESHAPE(int_zenith_angle, [row_length * rows])
    cosza(:) = RESHAPE(cos_zenith_angle, [row_length * rows])
    daylength1d(:) = RESHAPE(daylength, [row_length * rows])
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_int_cosz
  subroutine dealloc_diurnal_oxidant()
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'DEALLOC_DIURNAL_OXIDANT'
    real(kind=jprb) :: zhook_handle

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (ALLOCATED(daylength1d)) then
      DEALLOCATE(daylength1d)
    end if
    if (ALLOCATED(cosza)) then
      DEALLOCATE(cosza)
    end if
    if (ALLOCATED(intcosza)) then
      DEALLOCATE(intcosza)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine dealloc_diurnal_oxidant

end module ukca_diurnal_oxidant
