module ukca_chem_defs_mod
  implicit none
  integer, parameter, public :: photol_varname_len = 10
  type, public :: chch_t
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: speci
  INTEGER :: nodd
  CHARACTER(LEN = 10) :: ctype
  CHARACTER(LEN = 10) :: family
  INTEGER :: switch1
  INTEGER :: switch2
END TYPE chch_t
  type, public :: ratb_t
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  CHARACTER(LEN = 10) :: prod3
  CHARACTER(LEN = 10) :: prod4
  REAL :: k0
  REAL :: alpha
  REAL :: beta
  REAL :: pyield1
  REAL :: pyield2
  REAL :: pyield3
  REAL :: pyield4
END TYPE ratb_t
  type, public :: rath_t
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  CHARACTER(LEN = 10) :: prod3
  CHARACTER(LEN = 10) :: prod4
  REAL :: pyield1
  REAL :: pyield2
  REAL :: pyield3
  REAL :: pyield4
END TYPE rath_t
  type, public :: ratj_t
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  CHARACTER(LEN = 10) :: prod3
  CHARACTER(LEN = 10) :: prod4
  REAL :: pyield1
  REAL :: pyield2
  REAL :: pyield3
  REAL :: pyield4
  REAL :: jfacta
  CHARACTER(LEN = photol_varname_len) :: fname
END TYPE ratj_t
  type, public :: ratt_t
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  REAL :: f
  REAL :: k1
  REAL :: alpha1
  REAL :: beta1
  REAL :: k2
  REAL :: alpha2
  REAL :: beta2
  REAL :: pyield1
  REAL :: pyield2
END TYPE ratt_t
  type, public :: chch_t1
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: speci
  INTEGER :: nodd
  CHARACTER(LEN = 10) :: ctype
  CHARACTER(LEN = 10) :: family
  INTEGER :: switch1
  INTEGER :: switch2
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE chch_t1
  type, public :: ratb_t1
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  CHARACTER(LEN = 10) :: prod3
  CHARACTER(LEN = 10) :: prod4
  REAL :: k0
  REAL :: alpha
  REAL :: beta
  REAL :: pyield1
  REAL :: pyield2
  REAL :: pyield3
  REAL :: pyield4
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE ratb_t1
  type, public :: rath_t1
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  CHARACTER(LEN = 10) :: prod3
  CHARACTER(LEN = 10) :: prod4
  REAL :: pyield1
  REAL :: pyield2
  REAL :: pyield3
  REAL :: pyield4
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE rath_t1
  type, public :: ratj_t1
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  CHARACTER(LEN = 10) :: prod3
  CHARACTER(LEN = 10) :: prod4
  REAL :: pyield1
  REAL :: pyield2
  REAL :: pyield3
  REAL :: pyield4
  REAL :: jfacta
  CHARACTER(LEN = 10) :: fname
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE ratj_t1
  type, public :: ratt_t1
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: react1
  CHARACTER(LEN = 10) :: react2
  CHARACTER(LEN = 10) :: prod1
  CHARACTER(LEN = 10) :: prod2
  REAL :: f
  REAL :: k1
  REAL :: alpha1
  REAL :: beta1
  REAL :: k2
  REAL :: alpha2
  REAL :: beta2
  REAL :: pyield1
  REAL :: pyield2
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE ratt_t1
  type, public :: depvel_t
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: speci
  REAL :: velocity(30)
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE depvel_t
  type, public :: wetdep
  INTEGER :: item_no
  CHARACTER(LEN = 10) :: speci
  REAL :: params(6)
  INTEGER :: scheme
  INTEGER :: qual
  INTEGER :: disqual
  INTEGER :: version
END TYPE wetdep
  type(chch_t), allocatable, dimension(:), save, public :: chch_defs
  type(ratb_t), allocatable, dimension(:), save, public :: ratb_defs
  type(rath_t), allocatable, dimension(:), save, public :: rath_defs
  type(ratj_t), allocatable, dimension(:), save, public :: ratj_defs
  type(ratt_t), allocatable, dimension(:), save, public :: ratt_defs
  real, allocatable, dimension(:,:,:), save, public :: depvel_defs
  real, allocatable, dimension(:,:), save, public :: henry_defs
  type(chch_t1), allocatable, dimension(:), save, public :: chch_defs_new
  type(ratb_t1), allocatable, dimension(:), save, public :: ratb_defs_new
  type(rath_t1), allocatable, dimension(:), save, public :: rath_defs_new
  type(ratj_t1), allocatable, dimension(:), save, public :: ratj_defs_new
  type(ratt_t1), allocatable, dimension(:), save, public :: ratt_defs_new
  real, allocatable, dimension(:,:,:), save, public :: depvel_defs_new
  real, allocatable, dimension(:,:), save, public :: henry_defs_new
  character(len = 10), allocatable, target, save, public :: ratj_varnames(:)
  character(len = 10), allocatable, target, save, public :: ratj_data(:, :)
  public

  contains
  subroutine ukca_get_photol_reaction_data(data_ptr, varnames_ptr)
    use asad_mod, only : jpspj
    CHARACTER(LEN = 10), POINTER, INTENT(OUT) :: data_ptr(:, :)
    CHARACTER(LEN = photol_varname_len), POINTER, INTENT(OUT) :: varnames_ptr(:)
    integer :: i
    integer :: jppj

    jppj = SIZE(ratj_defs)
    if (.NOT.ALLOCATED(ratj_data)) then
      ALLOCATE(ratj_data(1:jppj,1:jpspj))
      do i = 1, jppj, 1
        ratj_data(i,1) = ratj_defs(i)%react1
        ratj_data(i,2) = ratj_defs(i)%react2
        ratj_data(i,3) = ratj_defs(i)%prod1
        ratj_data(i,4) = ratj_defs(i)%prod2
        ratj_data(i,5) = ratj_defs(i)%prod3
        ratj_data(i,6) = ratj_defs(i)%prod4
      enddo
    end if
    if (.NOT.ALLOCATED(ratj_varnames)) then
      ALLOCATE(ratj_varnames(1:jppj))
      do i = 1, jppj, 1
        ratj_varnames(i) = ratj_defs(i)%fname
      enddo
    end if
    data_ptr => ratj_data(:,:)
    varnames_ptr => ratj_varnames(:)
    return

  end subroutine ukca_get_photol_reaction_data

end module ukca_chem_defs_mod
