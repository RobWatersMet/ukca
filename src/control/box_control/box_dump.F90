module box_dump_mod

  !! Write every value passed into, or returned from, the UKCA library to a
  !! plain-text dump file, one value per line, so that two box models can be
  !! compared with `diff`.
  !!
  !! Each line is fully keyed and self-contained:
  !!
  !!     <step> <section> <name> <index> <value>
  !!
  !! e.g.
  !!
  !!     0 ukca_setup:in i_ukca_chem - 51
  !!     3 ukca_step:out tracers (1,1,1,2) 9.90333269774813254E-008
  !!
  !! so dumps can also be compared after `sort` if the two models make calls
  !! in a different order. Reals are written with 17 significant digits, which
  !! round-trips a double exactly: identical text means identical bits.
  !!
  !! This module deliberately depends only on Fortran intrinsics so the same
  !! file can be compiled into the old (UM/FCM) box model unchanged.
  !!
  !! Usage:
  !!
  !!     call dump_open("ukca_dump.txt")
  !!     call dump_section("ukca_setup:in", 0)
  !!     call dump("i_ukca_chem", i_ukca_chem)
  !!     call dump("tracers", tracer_ukca)
  !!     call dump_close()
  !!
  !! All dump calls are no-ops while no dump file is open. Instead of calling
  !! dump_open directly, a driver (or wrapper) can call dump_open_from_env,
  !! which opens the file named by the environment variable UKCA_DUMP_FILE,
  !! if set, so dumping can be switched on without changing driver code.
  !!
  !! Optional arguments that are not present, and allocatable arguments that
  !! are not allocated, are recorded with dump_marker as "<absent>" and
  !! "<unallocated>" so that differences in what is passed are also visible.

  use, intrinsic :: iso_fortran_env, only: error_unit, int64

  implicit none

  private

  public :: dump_open
  public :: dump_open_from_env
  public :: dump_marker
  public :: dump_close
  public :: dump_is_open
  public :: dump_section
  public :: dump

  !> write one variable (any supported type and rank) to the dump file
  interface dump
    module procedure dump_real_0d, dump_real_1d, dump_real_2d, &
      dump_real_3d, dump_real_4d
    module procedure dump_int_0d, dump_int_1d, dump_int_2d, dump_int_3d
    module procedure dump_logical_0d, dump_logical_1d, dump_logical_2d, &
      dump_logical_3d
    module procedure dump_char_0d, dump_char_1d
  end interface dump

  !> minimum width the variable name is padded to, for readability only
  integer, parameter :: name_width = 32
  !> format for real values: 17 significant digits round-trips a double
  character(len=*), parameter :: real_fmt = '(ES25.17E3)'

  !> environment variable naming the dump file, see dump_open_from_env
  character(len=*), parameter :: dump_env_var = 'UKCA_DUMP_FILE'

  integer :: dump_unit = -1
  logical :: env_checked = .false.
  integer :: current_step = 0
  character(len=64) :: current_section = 'unset'

contains

  ! ---------------------------------------------------------------------------
  ! File and context control
  ! ---------------------------------------------------------------------------

  subroutine dump_open(filename)
    !! Open (replacing) the dump file. Dumping is active until dump_close.
    character(len=*), intent(in) :: filename
    integer :: iostat
    character(len=256) :: iomsg

    if (dump_unit /= -1) call dump_close()

    open(newunit=dump_unit, file=filename, status='replace', &
      action='write', iostat=iostat, iomsg=iomsg)
    if (iostat /= 0) then
      write(error_unit, '(A)') 'box_dump_mod: cannot open dump file ' // &
        trim(filename) // ': ' // trim(iomsg)
      error stop 1
    end if

    write(dump_unit, '(A)') '# step section name index value'

  end subroutine dump_open

  subroutine dump_open_from_env()
    !! Open the dump file named by the environment variable UKCA_DUMP_FILE,
    !! if it is set and non-empty. Only the first call has any effect, so it
    !! is cheap to call at the start of every wrapped API routine.
    integer :: length, status
    character(len=:), allocatable :: filename

    if (env_checked) return
    env_checked = .true.

    call get_environment_variable(dump_env_var, length=length, status=status)
    if (status /= 0 .or. length == 0) return

    allocate(character(len=length) :: filename)
    call get_environment_variable(dump_env_var, value=filename)
    call dump_open(filename)

  end subroutine dump_open_from_env

  subroutine dump_close()
    !! Close the dump file, if open.
    if (dump_unit == -1) return
    close(dump_unit)
    dump_unit = -1
  end subroutine dump_close

  function dump_is_open() result(is_open)
    !! True if a dump file is open, i.e. dump calls will write output.
    logical :: is_open
    is_open = (dump_unit /= -1)
  end function dump_is_open

  subroutine dump_section(section, step)
    !! Set the section label (e.g. "ukca_step:in") and timestep used to key
    !! subsequent lines. Flushes the file so a crash leaves a usable dump.
    character(len=*), intent(in) :: section
    integer, intent(in) :: step

    current_section = section
    current_step = step
    if (dump_unit /= -1) flush(dump_unit)

  end subroutine dump_section

  subroutine dump_marker(name, marker)
    !! Write a marker in place of a value, e.g. "<absent>" for an optional
    !! argument that was not passed.
    character(len=*), intent(in) :: name, marker
    if (dump_unit == -1) return
    call write_line(name, '-', marker)
  end subroutine dump_marker

  ! ---------------------------------------------------------------------------
  ! Line writing and formatting helpers
  ! ---------------------------------------------------------------------------

  subroutine write_line(name, idx, value)
    character(len=*), intent(in) :: name, idx, value
    character(len=24) :: step_str

    write(step_str, '(I0)') current_step
    write(dump_unit, '(A)') trim(step_str) // ' ' // &
      trim(current_section) // ' ' // pad(name) // ' ' // idx // ' ' // &
      trim(adjustl(value))

  end subroutine write_line

  function pad(name) result(padded)
    character(len=*), intent(in) :: name
    character(len=max(name_width, len(name))) :: padded
    padded = name
  end function pad

  function real_str(value) result(str)
    !! Full-precision text for a real. NaNs are detected from the bit pattern
    !! (comparisons on a signalling NaN would trip -ffpe-trap=invalid).
    real, intent(in) :: value
    character(len=32) :: str
    integer(int64) :: bits
    integer(int64), parameter :: exp_mask = int(z'7FF0000000000000', int64)
    integer(int64), parameter :: man_mask = int(z'000FFFFFFFFFFFFF', int64)

    if (storage_size(value) == 64) then
      bits = transfer(value, bits)
      if (iand(bits, exp_mask) == exp_mask .and. &
        iand(bits, man_mask) /= 0_int64) then
        write(str, '(A,Z16.16)') 'NaN:', bits
        return
      end if
    end if
    write(str, real_fmt) value

  end function real_str

  function int_str(value) result(str)
    integer, intent(in) :: value
    character(len=24) :: str
    write(str, '(I0)') value
  end function int_str

  function logical_str(value) result(str)
    logical, intent(in) :: value
    character(len=1) :: str
    str = merge('T', 'F', value)
  end function logical_str

  function idx1(i) result(str)
    integer, intent(in) :: i
    character(len=64) :: str
    str = '(' // trim(int_str(i)) // ')'
  end function idx1

  function idx2(i, j) result(str)
    integer, intent(in) :: i, j
    character(len=64) :: str
    str = '(' // trim(int_str(i)) // ',' // trim(int_str(j)) // ')'
  end function idx2

  function idx3(i, j, k) result(str)
    integer, intent(in) :: i, j, k
    character(len=64) :: str
    str = '(' // trim(int_str(i)) // ',' // trim(int_str(j)) // ',' // &
      trim(int_str(k)) // ')'
  end function idx3

  function idx4(i, j, k, l) result(str)
    integer, intent(in) :: i, j, k, l
    character(len=100) :: str
    str = '(' // trim(int_str(i)) // ',' // trim(int_str(j)) // ',' // &
      trim(int_str(k)) // ',' // trim(int_str(l)) // ')'
  end function idx4

  ! ---------------------------------------------------------------------------
  ! Specific procedures for the generic dump interface. Array elements are
  ! written in Fortran (column-major) order with 1-based indices.
  ! ---------------------------------------------------------------------------

  subroutine dump_real_0d(name, value)
    character(len=*), intent(in) :: name
    real, intent(in) :: value
    if (dump_unit == -1) return
    call write_line(name, '-', real_str(value))
  end subroutine dump_real_0d

  subroutine dump_real_1d(name, value)
    character(len=*), intent(in) :: name
    real, intent(in) :: value(:)
    integer :: i
    if (dump_unit == -1) return
    do i = 1, size(value, 1)
      call write_line(name, trim(idx1(i)), real_str(value(i)))
    end do
  end subroutine dump_real_1d

  subroutine dump_real_2d(name, value)
    character(len=*), intent(in) :: name
    real, intent(in) :: value(:,:)
    integer :: i, j
    if (dump_unit == -1) return
    do j = 1, size(value, 2)
      do i = 1, size(value, 1)
        call write_line(name, trim(idx2(i, j)), real_str(value(i, j)))
      end do
    end do
  end subroutine dump_real_2d

  subroutine dump_real_3d(name, value)
    character(len=*), intent(in) :: name
    real, intent(in) :: value(:,:,:)
    integer :: i, j, k
    if (dump_unit == -1) return
    do k = 1, size(value, 3)
      do j = 1, size(value, 2)
        do i = 1, size(value, 1)
          call write_line(name, trim(idx3(i, j, k)), &
            real_str(value(i, j, k)))
        end do
      end do
    end do
  end subroutine dump_real_3d

  subroutine dump_real_4d(name, value)
    character(len=*), intent(in) :: name
    real, intent(in) :: value(:,:,:,:)
    integer :: i, j, k, l
    if (dump_unit == -1) return
    do l = 1, size(value, 4)
      do k = 1, size(value, 3)
        do j = 1, size(value, 2)
          do i = 1, size(value, 1)
            call write_line(name, trim(idx4(i, j, k, l)), &
              real_str(value(i, j, k, l)))
          end do
        end do
      end do
    end do
  end subroutine dump_real_4d

  subroutine dump_int_0d(name, value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    if (dump_unit == -1) return
    call write_line(name, '-', int_str(value))
  end subroutine dump_int_0d

  subroutine dump_int_1d(name, value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:)
    integer :: i
    if (dump_unit == -1) return
    do i = 1, size(value, 1)
      call write_line(name, trim(idx1(i)), int_str(value(i)))
    end do
  end subroutine dump_int_1d

  subroutine dump_int_2d(name, value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:,:)
    integer :: i, j
    if (dump_unit == -1) return
    do j = 1, size(value, 2)
      do i = 1, size(value, 1)
        call write_line(name, trim(idx2(i, j)), int_str(value(i, j)))
      end do
    end do
  end subroutine dump_int_2d

  subroutine dump_int_3d(name, value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:,:,:)
    integer :: i, j, k
    if (dump_unit == -1) return
    do k = 1, size(value, 3)
      do j = 1, size(value, 2)
        do i = 1, size(value, 1)
          call write_line(name, trim(idx3(i, j, k)), &
            int_str(value(i, j, k)))
        end do
      end do
    end do
  end subroutine dump_int_3d

  subroutine dump_logical_0d(name, value)
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    if (dump_unit == -1) return
    call write_line(name, '-', logical_str(value))
  end subroutine dump_logical_0d

  subroutine dump_logical_1d(name, value)
    character(len=*), intent(in) :: name
    logical, intent(in) :: value(:)
    integer :: i
    if (dump_unit == -1) return
    do i = 1, size(value, 1)
      call write_line(name, trim(idx1(i)), logical_str(value(i)))
    end do
  end subroutine dump_logical_1d

  subroutine dump_logical_2d(name, value)
    character(len=*), intent(in) :: name
    logical, intent(in) :: value(:,:)
    integer :: i, j
    if (dump_unit == -1) return
    do j = 1, size(value, 2)
      do i = 1, size(value, 1)
        call write_line(name, trim(idx2(i, j)), logical_str(value(i, j)))
      end do
    end do
  end subroutine dump_logical_2d

  subroutine dump_logical_3d(name, value)
    character(len=*), intent(in) :: name
    logical, intent(in) :: value(:,:,:)
    integer :: i, j, k
    if (dump_unit == -1) return
    do k = 1, size(value, 3)
      do j = 1, size(value, 2)
        do i = 1, size(value, 1)
          call write_line(name, trim(idx3(i, j, k)), &
            logical_str(value(i, j, k)))
        end do
      end do
    end do
  end subroutine dump_logical_3d

  subroutine dump_char_0d(name, value)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value
    if (dump_unit == -1) return
    call write_line(name, '-', '"' // trim(value) // '"')
  end subroutine dump_char_0d

  subroutine dump_char_1d(name, value)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value(:)
    integer :: i
    if (dump_unit == -1) return
    do i = 1, size(value, 1)
      call write_line(name, trim(idx1(i)), '"' // trim(value(i)) // '"')
    end do
  end subroutine dump_char_1d

end module box_dump_mod
