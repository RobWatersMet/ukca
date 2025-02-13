module ukca_chemco_raq_init_mod
  implicit none
  type, private :: rr
    real, public :: a
    real, public :: n
    real, public :: e
    integer, public :: r
  end type rr
  type, private :: pd
    integer, public :: klo
    integer, public :: khi
    integer, public :: res
    real, public :: fc
  end type pd
  type(rr), allocatable, dimension(:), save, public :: rk
  type(pd), allocatable, dimension(:), save, public :: pdep
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMCO_RAQ_INIT_MOD'
  public

  contains
  subroutine ukca_chemco_raq_init(nr, npdep)
    use ukca_config_specification_mod, only : ukca_config
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    use umprintmgr, only : ummessage, umprint
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr
    integer, intent(in) :: npdep
    integer :: icode
    integer :: i
    integer :: r
    real :: ar
    real :: n
    real :: ea
    CHARACTER(LEN = errormessagelength) :: cmessage
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMCO_RAQ_INIT'
    TYPE(rr) :: zero_rate = rr(0.0, 0.0, 0.0, 0)
    real(kind=jprb) :: zhook_handle

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (ukca_config%l_ukca_raqaero) then
      if (nr /= 197) then
        cmessage(:) = 'Invalid number of thermal reactions for RAQ-AERO chemistry'
        icode = INT(nr)
        call ereport(routinename, icode, cmessage)
      end if
    else
      if (nr /= 192) then
        cmessage(:) = 'Invalid number of thermal reactions for RAQ chemistry'
        icode = INT(nr)
        call ereport(routinename, icode, cmessage)
      end if
    end if
    if (npdep /= 10) then
      cmessage(:) = 'Invalid number of pressure-dependent reactions'
      icode = INT(npdep)
      call ereport(routinename, icode, cmessage)
    end if
    ALLOCATE(rk(1:nr))
    rk(1) = rr(6.0E-34, - 2.6, 0.0, 2)
    rk(2) = zero_rate
    rk(3) = zero_rate
    rk(4) = rr(1.0E-31, - 1.6, 0.0, 2)
    rk(5) = rr(3.0E-11, 0.3, 0.0, 2)
    rk(6) = rr(3.2E-11, 0.0, 70.0, 3)
    rk(7) = rr(1.8E-11, 0.0, 110.0, 3)
    rk(8) = rr(2.2E-10, 0.0, 0.0, 1)
    rk(9) = zero_rate
    rk(10) = zero_rate
    rk(11) = rr(1.4E-12, 0.0, - 1310.0, 3)
    rk(12) = rr(1.4E-13, 0.0, - 2470.0, 3)
    rk(13) = rr(1.7E-12, 0.0, - 940.0, 3)
    rk(14) = rr(2.03E-16, 4.57, 693.0, 4)
    rk(15) = rr(1.8E-11, 0.0, 110.0, 3)
    rk(16) = rr(5.5E-12, 0.0, 188.0, 3)
    rk(17) = rr(3.6E-12, 0.0, 270.0, 3)
    rk(18) = zero_rate
    rk(19) = rr(4.5E-14, 0.0, - 1260.0, 3)
    rk(20) = rr(3.6E-30, - 4.1, 0.0, 2)
    rk(21) = rr(3.3E-30, - 3.0, 0.0, 2)
    rk(22) = rr(1.8E-31, - 3.2, 0.0, 2)
    rk(23) = rr(4.1E-5, 0.0, - 10650.0, 3)
    rk(24) = rr(3.2E-13, 0.0, 690.0, 3)
    rk(25) = rr(4.8E15, 0.0, - 11170.0, 3)
    rk(26) = zero_rate
    rk(27) = rr(8.5E-13, 0.0, - 2450.0, 3)
    rk(28) = rr(3.5E-12, 0.0, - 925.0, 3)
    rk(29) = rr(1.3E-03, - 3.5, - 11000.0, 4)
    rk(30) = rr(4.8E-11, 0.0, 250.0, 3)
    rk(31) = rr(2.9E-12, 0.0, - 160.0, 3)
    rk(32) = rr(4.2E-12, 0.0, 0.0, 1)
    rk(33) = rr(5.5E-12, 0.0, - 2000.0, 3)
    rk(34) = rr(3.5E-12, 0.0, 0.0, 1)
    rk(35) = rr(2.4E-14, 0.0, 460.0, 3)
    rk(36) = rr(2.2E-13, 0.0, 600.0, 3)
    rk(37) = rr(1.9E-33, 0.0, 980.0, 3)
    rk(38) = rr(1.4E-21, 0.0, 2200.0, 3)
    rk(39) = zero_rate
    rk(40) = zero_rate
    rk(41) = zero_rate
    rk(42) = zero_rate
    rk(43) = zero_rate
    rk(44) = rr(2.66E-12, 0.0, 200.0, 3)
    rk(45) = rr(1.14E-12, 0.0, 200.0, 3)
    rk(46) = rr(1.9E-12, 0.2, 0.0, 2)
    rk(47) = rr(4.1E-11, 0.0, 0.0, 1)
    rk(48) = rr(4.7E-12, 0.0, 0.0, 1)
    rk(49) = rr(9.7E+14, 0.1, - 11080.0, 4)
    rk(50) = rr(2.7E-17, 0.0, 2199.0, 3)
    rk(51) = rr(6.5E-34, 0.0, 1335.0, 3)
    rk(52) = zero_rate
    rk(53) = zero_rate
    rk(54) = zero_rate
    rk(55) = zero_rate
    rk(56) = zero_rate
    rk(57) = zero_rate
    rk(58) = zero_rate
    rk(59) = rr(1.85E-12, 0.0, - 1690.0, 3)
    rk(60) = rr(2.3E-12, 0.0, 360.0, 3)
    rk(61) = rr(7.4E-13, 0.0, - 520.0, 3)
    rk(62) = rr(1.03E-13, 0.0, 365.0, 3)
    rk(63) = rr(7.3E-12, 0.0, - 620.0, 3)
    rk(64) = zero_rate
    rk(65) = rr(4.1E-13, 0.0, 750.0, 3)
    rk(66) = rr(5.4E-12, 0.0, 135.0, 3)
    rk(67) = rr(5.8E-16, 0.0, 0.0, 1)
    rk(68) = zero_rate
    rk(69) = rr(3.54E-33, 0.0, 0.0, 1)
    rk(70) = rr(1.5E-13, 0.0, 0.0, 1)
    rk(71) = rr(6.9E-12, 0.0, - 1000.0, 3)
    rk(72) = rr(2.6E-12, 0.0, 0.0, 1)
    rk(73) = rr(2.0E-13, 0.0, 0.0, 1)
    rk(74) = rr(4.4E5, 0.0, - 3910.0, 3)
    rk(75) = rr(4.4E-12, 0.0, 365.0, 3)
    rk(76) = rr(2.7E-28, - 7.1, 0.0, 2)
    rk(77) = rr(1.2E-11, - 0.9, 0.0, 2)
    rk(78) = rr(4.9E-03, 0.0, - 12100.0, 3)
    rk(79) = rr(2.0E-11, 0.0, 0.0, 1)
    rk(80) = rr(1.1E-11, 0.0, 0.0, 1)
    rk(81) = rr(7.9E-13, 2.0, 300.0, 4)
    rk(82) = rr(5.4E+16, 0.0, - 13830.0, 3)
    rk(83) = rr(2.54E-12, 0.0, 360.0, 3)
    rk(84) = rr(2.5E-13, 0.0, 0.0, 1)
    rk(85) = zero_rate
    rk(86) = rr(1.3E-12, 0.0, - 25.0, 3)
    rk(87) = zero_rate
    rk(88) = zero_rate
    rk(89) = rr(8.8E-12, 0.0, - 1320.0, 3)
    rk(90) = rr(6.4E-14, 0.0, 0.0, 1)
    rk(91) = rr(2.9E-12, 0.0, 500.0, 3)
    rk(92) = rr(7.6E-12, 0.0, - 585.0, 3)
    rk(93) = rr(2.7E-12, 0.0, 360.0, 3)
    rk(94) = rr(1.7E-14, 0.0, 423.0, 3)
    rk(95) = rr(2.45E-12, 0.0, 360.0, 3)
    rk(96) = rr(3.8E-12, 0.0, 0.0, 1)
    rk(97) = rr(4.0E-14, 0.0, 0.0, 1)
    rk(98) = rr(9.5E-13, 0.0, - 650.0, 3)
    rk(99) = rr(3.8E-13, 0.0, 900.0, 3)
    rk(100) = rr(8.0E-12, 0.0, 0.0, 1)
    rk(101) = rr(1.51E-13, 0.0, 1300.0, 3)
    rk(102) = rr(1.66E-11, 0.0, 0.0, 1)
    rk(103) = rr(1.82E-13, 0.0, 1300.0, 3)
    rk(104) = rr(2.15E-11, 0.0, 0.0, 1)
    rk(105) = rr(2.54E-12, 0.0, 360.0, 3)
    rk(106) = rr(8.8E-13, 0.0, 0.0, 1)
    rk(107) = zero_rate
    rk(108) = rr(8.6E-29, - 3.1, 0.0, 2)
    rk(109) = rr(9.0E-12, - 0.85, 0.0, 2)
    rk(110) = rr(9.0E-12, 0.0, 0.0, 1)
    rk(111) = rr(2.0E-12, 0.0, 0.0, 1)
    rk(112) = rr(1.2E-14, 0.0, - 2630.0, 3)
    rk(113) = zero_rate
    rk(114) = zero_rate
    rk(115) = zero_rate
    rk(116) = zero_rate
    rk(117) = zero_rate
    rk(118) = zero_rate
    rk(119) = zero_rate
    rk(120) = zero_rate
    rk(121) = zero_rate
    rk(122) = zero_rate
    rk(123) = rr(2.75E-15, 0.0, - 1878.0, 3)
    rk(124) = rr(2.75E-15, 0.0, - 1878.0, 3)
    rk(125) = rr(8.0E-27, - 3.5, 0.0, 2)
    rk(126) = rr(2.54E-12, 0.0, 360.0, 3)
    rk(127) = rr(6.0E-13, 0.0, 0.0, 1)
    rk(128) = rr(7.86E-15, 0.0, - 1913.0, 3)
    rk(129) = rr(7.56E-16, 0.0, - 1521.0, 3)
    rk(130) = rr(3.0E-11, - 1.0, 0.0, 2)
    rk(131) = rr(5.2E-13, 0.0, 980.0, 3)
    rk(132) = rr(5.0E-13, 0.0, 0.0, 1)
    rk(133) = rr(2.0E-12, 0.0, 0.0, 1)
    rk(134) = rr(2.45E-13, 0.0, 1250.0, 3)
    rk(135) = rr(4.2E-11, 0.0, 0.0, 1)
    rk(136) = rr(2.23E-13, 0.0, 1250.0, 3)
    rk(137) = rr(5.77E-11, 0.0, 0.0, 1)
    rk(138) = rr(1.72E-11, 0.0, 0.0, 1)
    rk(139) = rr(1.14E-11, 0.0, 0.0, 1)
    rk(140) = zero_rate
    rk(141) = rr(2.54E-11, 0.0, 410.0, 3)
    rk(142) = rr(2.08E-12, 0.0, 180.0, 3)
    rk(143) = rr(4.13E-12, 0.0, 452.0, 3)
    rk(144) = rr(2.5E-12, 0.0, 360.0, 3)
    rk(145) = zero_rate
    rk(146) = zero_rate
    rk(147) = zero_rate
    rk(148) = zero_rate
    rk(149) = zero_rate
    rk(150) = zero_rate
    rk(151) = rr(5.7E-12, 0.0, - 4426.0, 3)
    rk(152) = rr(2.8E-12, 0.0, - 3280.0, 3)
    rk(153) = rr(3.3E-12, 2.0, - 2880.0, 4)
    rk(154) = rr(4.95E-12, 0.0, 0.0, 1)
    rk(155) = rr(4.59E-13, 0.0, - 1156.0, 3)
    rk(156) = rr(5.25E-12, 0.0, 0.0, 1)
    rk(157) = zero_rate
    rk(158) = rr(1.4E-12, 0.0, - 1860.0, 3)
    rk(159) = rr(3.03E-12, 0.0, - 446.0, 3)
    rk(160) = rr(4.16E-11, 0.0, 0.0, 1)
    rk(161) = zero_rate
    rk(162) = zero_rate
    rk(163) = zero_rate
    rk(164) = zero_rate
    rk(165) = zero_rate
    rk(166) = zero_rate
    rk(167) = zero_rate
    rk(168) = zero_rate
    rk(169) = zero_rate
    rk(170) = rr(1.36E-11, 0.0, 0.0, 1)
    rk(171) = rr(1.0E-11, 0.0, 0.0, 1)
    rk(172) = rr(5.6E-11, 0.0, 0.0, 1)
    rk(173) = rr(2.54E-12, 0.0, 360.0, 3)
    rk(174) = rr(1.18E-12, 0.0, 338.0, 3)
    rk(175) = rr(3.6E-13, 0.0, 0.0, 1)
    rk(176) = rr(1.36E-11, 0.0, 0.0, 1)
    rk(177) = rr(2.7E-12, 0.0, 0.0, 1)
    rk(178) = rr(7.0E-14, 0.0, 0.0, 1)
    rk(179) = zero_rate
    rk(180) = rr(2.5E-13, 0.0, 1300.0, 3)
    rk(181) = rr(1.0E-13, 0.0, 0.0, 1)
    rk(182) = rr(1.0E-11, 0.0, 0.0, 1)
    rk(183) = rr(1.0E-11, 0.0, 0.0, 1)
    rk(184) = zero_rate
    rk(185) = zero_rate
    rk(186) = zero_rate
    rk(187) = zero_rate
    rk(188) = zero_rate
    rk(189) = zero_rate
    rk(190) = zero_rate
    rk(191) = zero_rate
    rk(192) = zero_rate
    if (ukca_config%l_ukca_raqaero) then
      rk(193:197) = zero_rate
    end if
    do i = 1, nr, 1
      r = rk(i)%r
      ar = rk(i)%a
      n = rk(i)%n
      ea = rk(i)%e
      if (r == 1 .AND. (n /= 0.0 .OR. ea /= 0.0) .OR. r == 2 .AND. (n == 0.0 .OR. ea /= 0.0) .OR. r == 3 .AND. (n /= 0.0 .OR. ea == 0.0) .OR. r == 4 .AND. (n == 0.0 .OR. ea == 0.0)) then
        cmessage(:) = 'Faulty rate constant'
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,A,I8,A,I8,A,F16.8,A,F16.8,A,F16.8)') 'Faulty rate constant', ' i=', i, ' r=', r, ' ar=', ar, ' n=', n, ' ea=', ea
        call umprint(ummessage, src=routinename)
        icode = i
        call ereport(routinename, icode, cmessage)
      end if
    enddo
    ALLOCATE(pdep(1:npdep))
    pdep(1) = pd(4, 5, 5, 0.85)
    pdep(2) = pd(20, 46, 20, 0.35)
    pdep(3) = pd(21, 47, 21, 0.4)
    pdep(4) = pd(22, 48, 22, 0.6)
    pdep(5) = pd(23, 25, 23, 0.6)
    pdep(6) = pd(29, 49, 29, 0.35)
    pdep(7) = pd(76, 77, 77, 0.3)
    pdep(8) = pd(78, 82, 78, 0.3)
    pdep(9) = pd(108, 109, 109, 0.48)
    pdep(10) = pd(125, 130, 125, 0.5)
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_chemco_raq_init

end module ukca_chemco_raq_init_mod
