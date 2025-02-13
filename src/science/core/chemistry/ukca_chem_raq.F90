module ukca_chem_raq
  use ukca_chem_defs_mod, only : chch_t, ratb_t, rath_t, ratj_t, ratt_t
  implicit none
  integer, parameter, public :: ndry_raq = 16
  integer, parameter, public :: nwet_raq = 19
  type(ratb_t), dimension(113), parameter, public :: ratb_defs_raq = ratb_t('          ', '          ', '          ', '          ', '          ', '          ', 0.00E+00, 0.00, 0.00, 0.000, 0.000, 0.000, 0.000)
  type(ratt_t), dimension(12), parameter, public :: ratt_defs_raq = ratt_t('          ', '          ', '          ', '          ', 0.0, 0.00E+00, 0.00, 0.0, 0.00E+00, 0.00, 0.0, 0.000, 0.000)
  integer, parameter, public :: ndd_o3 = 1
  integer, parameter, public :: ndd_no2 = 2
  integer, parameter, public :: ndd_hono2 = 3
  integer, parameter, public :: ndd_h2o2 = 4
  integer, parameter, public :: ndd_ch4 = 5
  integer, parameter, public :: ndd_co = 6
  integer, parameter, public :: ndd_meooh = 7
  integer, parameter, public :: ndd_etooh = 8
  integer, parameter, public :: ndd_pan = 9
  integer, parameter, public :: ndd_iprooh = 10
  integer, parameter, public :: ndd_o3s = 11
  integer, parameter, public :: ndd_isooh = 12
  integer, parameter, public :: ndd_mvkooh = 13
  integer, parameter, public :: ndd_orgnit = 14
  integer, parameter, public :: ndd_h2 = 15
  integer, parameter, public :: ndd_sbuooh = 16
  integer, parameter, public :: ndw_no3 = 1
  integer, parameter, public :: ndw_n2o5 = 2
  integer, parameter, public :: ndw_ho2no2 = 3
  integer, parameter, public :: ndw_hono2 = 4
  integer, parameter, public :: ndw_h2o2 = 5
  integer, parameter, public :: ndw_hcho = 6
  integer, parameter, public :: ndw_ho2 = 7
  integer, parameter, public :: ndw_meoo = 8
  integer, parameter, public :: ndw_meooh = 9
  integer, parameter, public :: ndw_etooh = 10
  integer, parameter, public :: ndw_iprooh = 11
  integer, parameter, public :: ndw_isooh = 12
  integer, parameter, public :: ndw_ison = 13
  integer, parameter, public :: ndw_mgly = 14
  integer, parameter, public :: ndw_mvkooh = 15
  integer, parameter, public :: ndw_orgnit = 16
  integer, parameter, public :: ndw_ch3oh = 17
  integer, parameter, public :: ndw_sbuooh = 18
  integer, parameter, public :: ndw_gly = 19
  type(chch_t), dimension(58), save, public :: chch_defs_raq = [chch_t(1, 'O(3P)     ', 1, 'SS        ', '          ', 0, 0), chch_t(2, 'O(1D)     ', 1, 'SS        ', '          ', 0, 0), chch_t(3, 'OH        ', 1, 'SS        ', '          ', 0, 0), chch_t(4, 'O3        ', 1, 'TR        ', '          ', 1, 0), chch_t(5, 'NO        ', 1, 'TR        ', '          ', 0, 0), chch_t(6, 'NO3       ', 1, 'TR        ', '          ', 0, 1), chch_t(7, 'NO2       ', 1, 'TR        ', '          ', 1, 0), chch_t(8, 'N2O5      ', 2, 'TR        ', '          ', 0, 1), chch_t(9, 'HO2NO2    ', 1, 'TR        ', '          ', 0, 1), chch_t(10, 'HONO2     ', 1, 'TR        ', '          ', 1, 1), chch_t(11, 'H2O2      ', 1, 'TR        ', '          ', 1, 1), chch_t(12, 'CH4       ', 1, 'TR        ', '          ', 1, 0), chch_t(13, 'CO        ', 1, 'TR        ', '          ', 1, 0), chch_t(14, 'HCHO      ', 1, 'TR        ', '          ', 0, 1), chch_t(15, 'HO2       ', 1, 'SS        ', '          ', 0, 1), chch_t(16, 'MeOO      ', 1, 'SS        ', '          ', 0, 1), chch_t(17, 'MeOOH     ', 1, 'TR        ', '          ', 1, 1), chch_t(18, 'EtOO      ', 1, 'SS        ', '          ', 0, 0), chch_t(19, 'C2H6      ', 1, 'TR        ', '          ', 0, 0), chch_t(20, 'MeCO3     ', 1, 'SS        ', '          ', 0, 0), chch_t(21, 'EtOOH     ', 1, 'TR        ', '          ', 1, 1), chch_t(22, 'MeCHO     ', 1, 'TR        ', '          ', 0, 0), chch_t(23, 'PAN       ', 1, 'TR        ', '          ', 1, 0), chch_t(24, 's-BuOO    ', 1, 'SS        ', '          ', 0, 0), chch_t(25, 'C3H8      ', 1, 'TR        ', '          ', 0, 0), chch_t(26, 'i-PrOOH   ', 1, 'TR        ', '          ', 1, 1), chch_t(27, 'Me2CO     ', 1, 'TR        ', '          ', 0, 0), chch_t(28, 'O3S       ', 1, 'TR        ', '          ', 1, 0), chch_t(29, 'C5H8      ', 1, 'TR        ', '          ', 0, 0), chch_t(30, 'i-PrOO    ', 1, 'SS        ', '          ', 0, 0), chch_t(31, 'ISOOH     ', 1, 'TR        ', '          ', 1, 1), chch_t(32, 'ISON      ', 1, 'TR        ', '          ', 0, 1), chch_t(33, 'MGLY      ', 1, 'TR        ', '          ', 0, 1), chch_t(34, 'MVK       ', 1, 'TR        ', '          ', 0, 0), chch_t(35, 'MVKOOH    ', 1, 'TR        ', '          ', 1, 1), chch_t(36, 'MeCOCH2OO ', 1, 'SS        ', '          ', 0, 0), chch_t(37, 'MEKO2     ', 1, 'SS        ', '          ', 0, 0), chch_t(38, 'HOC2H4O2  ', 1, 'SS        ', '          ', 0, 0), chch_t(39, 'ORGNIT    ', 1, 'TR        ', '          ', 1, 1), chch_t(40, 'HOC3H6O2  ', 1, 'SS        ', '          ', 0, 0), chch_t(41, 'CH3OH     ', 1, 'TR        ', '          ', 0, 1), chch_t(42, 'OXYL1     ', 1, 'SS        ', '          ', 0, 0), chch_t(43, 'H2        ', 1, 'TR        ', '          ', 1, 0), chch_t(44, 'MEMALD1   ', 1, 'SS        ', '          ', 0, 0), chch_t(45, 'RNC2H4    ', 1, 'TR        ', '          ', 0, 0), chch_t(46, 'HOIPO2    ', 1, 'SS        ', '          ', 0, 0), chch_t(47, 'RNC3H6    ', 1, 'TR        ', '          ', 0, 0), chch_t(48, 'HOMVKO2   ', 1, 'SS        ', '          ', 0, 0), chch_t(49, 'C2H4      ', 1, 'TR        ', '          ', 0, 0), chch_t(50, 'C3H6      ', 1, 'TR        ', '          ', 0, 0), chch_t(51, 'C4H10     ', 1, 'TR        ', '          ', 0, 0), chch_t(52, 's-BuOOH   ', 1, 'TR        ', '          ', 1, 1), chch_t(53, 'MEK       ', 1, 'TR        ', '          ', 0, 0), chch_t(54, 'TOLUENE   ', 1, 'TR        ', '          ', 0, 0), chch_t(55, 'TOLP1     ', 1, 'SS        ', '          ', 0, 0), chch_t(56, 'MEMALD    ', 1, 'TR        ', '          ', 0, 0), chch_t(57, 'GLY       ', 1, 'TR        ', '          ', 0, 1), chch_t(58, 'oXYLENE   ', 1, 'TR        ', '          ', 0, 0)]
  type(rath_t), dimension(0), private :: rath_defs_raq_null
  type(rath_t), dimension(2), save, private :: rath_defs_raq_n2o5 = [rath_t('N2O5      ', '          ', 'HONO2     ', '          ', '          ', '          ', 2.000, 0.000, 0.000, 0.000), rath_t('N2O5      ', '          ', 'HONO2     ', '          ', '          ', '          ', 2.000, 0.000, 0.000, 0.000)]
  type(rath_t), allocatable, dimension(:), public :: rath_defs_raq
  type(ratj_t), dimension(23), save, public :: ratj_defs_raq = [ratj_t('O3        ', 'PHOTON    ', 'O2        ', 'O(3P)     ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jo3b      '), ratj_t('O3        ', 'PHOTON    ', 'O2        ', 'O(1D)     ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jo3a      '), ratj_t('NO2       ', 'PHOTON    ', 'NO        ', 'O(3P)     ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jno2      '), ratj_t('H2O2      ', 'PHOTON    ', 'OH        ', 'OH        ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jh2o2     '), ratj_t('HONO2     ', 'PHOTON    ', 'OH        ', 'NO2       ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jhno3     '), ratj_t('HCHO      ', 'PHOTON    ', 'HO2       ', 'HO2       ', 'CO        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jhchoa    '), ratj_t('HCHO      ', 'PHOTON    ', 'H2        ', 'CO        ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jhchob    '), ratj_t('MeCHO     ', 'PHOTON    ', 'MeOO      ', 'HO2       ', 'CO        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jaceta    '), ratj_t('MEK       ', 'PHOTON    ', 'MeCO3     ', 'EtOO      ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jaceto    '), ratj_t('Me2CO     ', 'PHOTON    ', 'MeCO3     ', 'MeOO      ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jaceto    '), ratj_t('HO2NO2    ', 'PHOTON    ', 'HO2       ', 'NO2       ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jpna      '), ratj_t('MGLY      ', 'PHOTON    ', 'MeCO3     ', 'CO        ', 'HO2       ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmkal     '), ratj_t('GLY       ', 'PHOTON    ', 'CO        ', 'CO        ', 'HO2       ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmkal     '), ratj_t('NO3       ', 'PHOTON    ', 'NO        ', 'O2        ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jno3a     '), ratj_t('NO3       ', 'PHOTON    ', 'NO2       ', 'O(3P)     ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jno3b     '), ratj_t('N2O5      ', 'PHOTON    ', 'NO3       ', 'NO2       ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jn2o5     '), ratj_t('MeOOH     ', 'PHOTON    ', 'HCHO      ', 'HO2       ', 'OH        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmhp      '), ratj_t('PAN       ', 'PHOTON    ', 'MeCO3     ', 'NO2       ', '          ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jpan      '), ratj_t('EtOOH     ', 'PHOTON    ', 'MeCHO     ', 'HO2       ', 'OH        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmhp      '), ratj_t('i-PrOOH   ', 'PHOTON    ', 'Me2CO     ', 'HO2       ', 'OH        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmhp      '), ratj_t('s-BuOOH   ', 'PHOTON    ', 'MEK       ', 'HO2       ', 'OH        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmhp      '), ratj_t('ISOOH     ', 'PHOTON    ', 'MVK       ', 'HCHO      ', 'OH        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmhp      '), ratj_t('MVKOOH    ', 'PHOTON    ', 'MGLY      ', 'HCHO      ', 'OH        ', '          ', 0.0, 0.0, 0.0, 0.0, 100.0, 'jmhp      ')]
  real, dimension(6,5,ndry_raq), save, public :: depvel_defs_raq = RESHAPE([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 1.00, 0.11, 0.56, 0.26, 0.11, 0.19, 1.00, 0.37, 0.69, 0.59, 0.46, 0.53, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.83, 0.04, 0.44, 0.06, 0.04, 0.05, 0.63, 0.06, 0.35, 0.07, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 4.00, 3.00, 3.50, 2.00, 1.00, 1.50, 2.50, 1.50, 2.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.25, 0.16, 0.71, 0.28, 0.12, 0.20, 1.25, 0.53, 0.89, 0.83, 0.78, 0.81, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.83, 0.04, 0.44, 0.06, 0.05, 0.06, 0.63, 0.06, 0.35, 0.08, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.83, 0.04, 0.44, 0.06, 0.05, 0.06, 0.63, 0.06, 0.35, 0.08, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.53, 0.04, 0.29, 0.06, 0.05, 0.06, 0.42, 0.06, 0.24, 0.07, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.00, 0.01, 0.01, 0.00, 0.00, 0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.83, 0.04, 0.44, 0.06, 0.05, 0.06, 0.63, 0.06, 0.35, 0.08, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 1.00, 0.11, 0.56, 0.26, 0.11, 0.19, 1.00, 0.37, 0.69, 0.59, 0.46, 0.53, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.83, 0.04, 0.44, 0.06, 0.05, 0.06, 0.63, 0.06, 0.35, 0.08, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.36, 0.36, 0.36, 0.19, 0.19, 0.19, 0.71, 0.04, 0.38, 0.06, 0.05, 0.06, 0.53, 0.07, 0.30, 0.08, 0.07, 0.08, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.02, 0.01, 0.02, 0.02, 0.01, 0.02, 0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.83, 0.04, 0.44, 0.06, 0.05, 0.06, 0.63, 0.06, 0.35, 0.08, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.83, 0.04, 0.44, 0.06, 0.05, 0.06, 0.63, 0.06, 0.35, 0.08, 0.06, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01], [6, 5, ndry_raq])
  real, dimension(6,nwet_raq), save, public :: henry_defs_raq = RESHAPE([0.200E+01, 0.200E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.210E+06, 0.870E+04, 0.200E+02, 0.000E+00, 0.000E+00, 0.000E+00, 0.130E+05, 0.690E+04, 0.100E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.210E+06, 0.870E+04, 0.200E+02, 0.000E+00, 0.000E+00, 0.000E+00, 0.830E+05, 0.740E+04, 0.240E-11, - 0.373E+04, 0.000E+00, 0.000E+00, 0.330E+04, 0.650E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.400E+04, 0.590E+04, 0.200E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.200E+04, 0.660E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.310E+03, 0.500E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.340E+03, 0.570E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.340E+03, 0.570E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.170E+07, 0.970E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.300E+04, 0.740E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.350E+04, 0.720E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.170E+07, 0.970E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.130E+03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.220E+03, 0.520E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.340E+03, 0.570E+04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.360E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00], [6, nwet_raq])
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEM_RAQ'
  private

  public :: ukca_init_raq

  contains
  subroutine ukca_init_raq()
    use asad_mod, only : jphk
    use ukca_config_specification_mod, only : ukca_config
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_INIT_RAQ'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    ALLOCATE(rath_defs_raq(1:jphk))
    if (ukca_config%l_ukca_classic_hetchem) then
      rath_defs_raq(:) = [rath_defs_raq_n2o5]
    else
      rath_defs_raq(:) = [rath_defs_raq_null]
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_init_raq

end module ukca_chem_raq
