module ukca_chemco_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMCO_MOD'
  public

  contains
  subroutine ukca_chemco(nr, n_pnts, tc, m, h2o, o2, clw, fcloud, fdiss, k_dms, rc, h_plus_1d_arr)
    use ukca_constants, only : h_plus, m_air
    use ukca_config_constants_mod, only : avogadro, rho_water
    use ukca_config_specification_mod, only : ukca_config
    use asad_mod, only : jpeq, jpspec, peps
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr
    integer, intent(in) :: n_pnts
    real, dimension(n_pnts), intent(in) :: tc
    real, dimension(n_pnts), intent(in) :: m
    real, dimension(n_pnts), intent(in) :: h2o
    real, dimension(n_pnts), intent(in) :: o2
    real, dimension(n_pnts), intent(in) :: clw
    real, dimension(n_pnts), intent(in) :: fcloud
    real, dimension(n_pnts,jpspec,jpeq + 1), intent(in) :: fdiss
    real, dimension(n_pnts), intent(in) :: h_plus_1d_arr
    real, dimension(n_pnts,5), intent(out) :: k_dms
    real, dimension(n_pnts,nr), intent(out) :: rc
    integer :: i
    integer :: j
    real, dimension(n_pnts) :: vr
    real, dimension(n_pnts) :: z1
    real, dimension(n_pnts) :: z2
    real, dimension(n_pnts) :: z3
    real, dimension(n_pnts) :: z4
    real, dimension(n_pnts) :: ratiob2total
    real, dimension(n_pnts) :: ratioa2b
    real, dimension(7) :: at
    real, dimension(n_pnts) :: t300
    real, dimension(n_pnts) :: zo
    real, dimension(n_pnts) :: zi
    real, dimension(n_pnts) :: arg2
    real, dimension(n_pnts) :: zfc
    real, dimension(n_pnts) :: zr
    real, dimension(n_pnts,jpspec) :: faq
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMCO'
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
    integer :: idx_36
    integer :: idx_37
    integer :: idx_38
    integer :: idx_39
    integer :: idx_40
    integer :: idx_41
    integer :: idx_42
    integer :: idx_43
    integer :: idx_44
    integer :: idx_45
    integer :: idx_46
    integer :: idx_47
    integer :: idx_48
    integer :: idx_49
    integer :: idx_50
    integer :: idx_51
    integer :: idx_52
    integer :: idx_53
    integer :: idx_54
    integer :: idx_55
    integer :: idx_56
    integer :: idx_57
    integer :: idx_58
    integer :: idx_59
    integer :: idx_60
    integer :: idx_61
    integer :: idx_62
    integer :: idx_63
    integer :: idx_64
    integer :: idx_65
    integer :: idx_66
    integer :: idx_67
    integer :: idx_68
    integer :: idx_69
    integer :: idx_70
    integer :: idx_71
    integer :: idx_72
    integer :: idx_73
    integer :: idx_74
    integer :: idx_75
    integer :: idx_76
    integer :: idx_77
    integer :: idx_78
    integer :: idx_79
    integer :: idx_80
    integer :: idx_81
    integer :: idx_82
    integer :: idx_83
    integer :: idx_84
    integer :: idx_85
    integer :: idx_86
    integer :: idx_87
    integer :: idx_88
    integer :: idx_89
    integer :: idx_90
    integer :: idx_91
    integer :: idx_92
    integer :: idx_93
    integer :: idx_94
    integer :: idx_95
    integer :: idx_96
    integer :: idx_97
    integer :: idx_98
    integer :: idx_99
    integer :: idx_100
    integer :: idx_101
    integer :: idx_102
    integer :: idx_103
    integer :: idx_104
    integer :: idx_105
    integer :: idx_106
    integer :: idx_107
    integer :: idx_108
    integer :: idx_109
    integer :: idx_110
    integer :: idx_111
    integer :: idx_112
    integer :: idx_113
    integer :: idx_114
    integer :: idx_115
    integer :: idx_116
    integer :: idx_117
    integer :: idx_118
    integer :: idx_119
    integer :: idx_120
    integer :: idx_121
    integer :: idx_122
    integer :: idx_123
    integer :: idx_124
    integer :: idx_125
    integer :: idx_126
    integer :: idx_127
    integer :: idx_128
    integer :: idx_129
    integer :: idx_130
    integer :: idx_131
    integer :: idx_132
    integer :: idx_133
    integer :: idx_134
    integer :: idx_135
    integer :: idx_136
    integer :: idx_137
    integer :: idx_138
    integer :: idx_139
    integer :: idx_140
    integer :: idx_141
    integer :: idx_142
    integer :: idx_143
    integer :: idx_144
    integer :: idx_145
    integer :: idx_146
    integer :: idx_147
    integer :: idx_148
    integer :: widx1
    integer :: widx1_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    if (ukca_config%l_ukca_aerchem) then
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do idx = LBOUND(faq, dim=2), UBOUND(faq, dim=2), 1
        do idx_1 = LBOUND(faq, dim=1), UBOUND(faq, dim=1), 1
          faq(idx_1,idx) = fdiss(idx_1,idx,1) + fdiss(idx_1,idx,2)
        enddo
      enddo
      vr(:) = clw(:) * m(:) * m_air * (1e6 / avogadro) / rho_water
      !$acc end kernels
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do j = 1, nr, 1
      do i = 1, n_pnts, 1
        rc(i,j) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_2 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_2,1) = 3.60e-12 * EXP(270.0 / tc(idx_2))
    enddo
    !$acc loop gang vector independent
    do idx_3 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_3,2) = 4.00e-12
    enddo
    !$acc loop gang vector independent
    do idx_4 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_4,3) = 2.03e-16 * (tc(idx_4) / 300.0) ** 4.57 * EXP(693.0 / tc(idx_4))
    enddo
    !$acc loop gang vector independent
    do idx_5 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_5,4) = 2.20e-13 * EXP(600.0 / tc(idx_5))
    enddo
    !$acc loop gang vector independent
    do idx_6 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_6,4) = rc(idx_6,4) * (1.0 + 1.4e-21 * h2o(idx_6) * EXP(2200.0 / tc(idx_6)))
    enddo
    !$acc loop gang vector independent
    do idx_7 = LBOUND(ratiob2total, dim=1), UBOUND(ratiob2total, dim=1), 1
      ratiob2total(idx_7) = 1.0 / (1.0 + 498.0 * EXP(-1160.0 / tc(idx_7)))
    enddo
    !$acc loop gang vector independent
    do idx_8 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_8,5) = 3.80e-13 * EXP(780.0 / tc(idx_8))
    enddo
    !$acc loop gang vector independent
    do idx_9 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_9,5) = rc(idx_9,5) * (1.0 - ratiob2total(idx_9))
    enddo
    !$acc loop gang vector independent
    do idx_10 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_10,6) = 3.80e-13 * EXP(780.0 / tc(idx_10))
    enddo
    !$acc loop gang vector independent
    do idx_11 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_11,6) = rc(idx_11,6) * ratiob2total(idx_11)
    enddo
    !$acc loop gang vector independent
    do idx_12 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_12,7) = 3.80e-13 * EXP(900.0 / tc(idx_12))
    enddo
    !$acc loop gang vector independent
    do idx_13 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_13,8) = 2.08e-13 * EXP(980.0 / tc(idx_13))
    enddo
    !$acc loop gang vector independent
    do idx_14 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_14,9) = 1.04e-13 * EXP(980.0 / tc(idx_14))
    enddo
    !$acc loop gang vector independent
    do idx_15 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_15,10) = 2.08e-13 * EXP(980.0 / tc(idx_15))
    enddo
    !$acc loop gang vector independent
    do idx_16 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_16,11) = 1.51e-13 * EXP(1300.0 / tc(idx_16))
    enddo
    !$acc loop gang vector independent
    do idx_17 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_17,12) = 1.51e-13 * EXP(1300.0 / tc(idx_17))
    enddo
    !$acc loop gang vector independent
    do idx_18 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_18,13) = 3.05e-13 * EXP(1040.0 / tc(idx_18))
    enddo
    !$acc loop gang vector independent
    do idx_19 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_19,14) = 1.25e-13 * EXP(1040.0 / tc(idx_19))
    enddo
    !$acc loop gang vector independent
    do idx_20 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_20,15) = 1.36e-13 * EXP(1250.0 / tc(idx_20))
    enddo
    !$acc loop gang vector independent
    do idx_21 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_21,16) = 2.95e-12 * EXP(285.0 / tc(idx_21))
    enddo
    !$acc loop gang vector independent
    do idx_22 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_22,17) = 2.95e-15 * EXP(285.0 / tc(idx_22))
    enddo
    !$acc loop gang vector independent
    do idx_23 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_23,18) = 1.30e-12
    enddo
    !$acc loop gang vector independent
    do idx_24 = LBOUND(ratiob2total, dim=1), UBOUND(ratiob2total, dim=1), 1
      ratiob2total(idx_24) = 1.0 / (1.0 + EXP(1300.0 / tc(idx_24)) / 33.0)
    enddo
    !$acc loop gang vector independent
    do idx_25 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_25,19) = 1.03e-13 * EXP(365.0 / tc(idx_25))
    enddo
    !$acc loop gang vector independent
    do idx_26 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_26,19) = rc(idx_26,19) * (1.0 - ratiob2total(idx_26))
    enddo
    !$acc loop gang vector independent
    do idx_27 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_27,20) = 1.03e-13 * EXP(365.0 / tc(idx_27))
    enddo
    !$acc loop gang vector independent
    do idx_28 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_28,20) = rc(idx_28,20) * ratiob2total(idx_28)
    enddo
    !$acc loop gang vector independent
    do idx_29 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_29,21) = 1.80e-12 * EXP(500.0 / tc(idx_29))
    enddo
    !$acc loop gang vector independent
    do idx_30 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_30,22) = 2.00e-13 * EXP(500.0 / tc(idx_30))
    enddo
    !$acc loop gang vector independent
    do idx_31 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_31,23) = 2.60e-12 * EXP(380.0 / tc(idx_31))
    enddo
    !$acc loop gang vector independent
    do idx_32 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_32,24) = 2.30e-12
    enddo
    !$acc loop gang vector independent
    do idx_33 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_33,25) = 4.40e-13 * EXP(1070.0 / tc(idx_33))
    enddo
    !$acc loop gang vector independent
    do idx_34 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_34,26) = 7.50e-12 * EXP(290.0 / tc(idx_34))
    enddo
    !$acc loop gang vector independent
    do idx_35 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_35,27) = 4.00e-12
    enddo
    !$acc loop gang vector independent
    do idx_36 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_36,28) = 2.90e-12 * EXP(350.0 / tc(idx_36))
    enddo
    !$acc loop gang vector independent
    do idx_37 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_37,29) = 2.50e-12
    enddo
    !$acc loop gang vector independent
    do idx_38 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_38,30) = 2.70e-12 * EXP(360.0 / tc(idx_38))
    enddo
    !$acc loop gang vector independent
    do idx_39 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_39,31) = 2.50e-12
    enddo
    !$acc loop gang vector independent
    do idx_40 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_40,32) = 6.70e-12 * EXP(340.0 / tc(idx_40))
    enddo
    !$acc loop gang vector independent
    do idx_41 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_41,33) = 4.00e-12
    enddo
    !$acc loop gang vector independent
    do idx_42 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_42,34) = 2.80e-12 * EXP(300.0 / tc(idx_42))
    enddo
    !$acc loop gang vector independent
    do idx_43 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_43,35) = 2.50e-12
    enddo
    !$acc loop gang vector independent
    do idx_44 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_44,36) = 1.80e-11 * EXP(110.0 / tc(idx_44))
    enddo
    !$acc loop gang vector independent
    do idx_45 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_45,37) = 1.40e-12 * EXP(-1310.0 / tc(idx_45))
    enddo
    !$acc loop gang vector independent
    do idx_46 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_46,38) = 1.40e-13 * EXP(-2470.0 / tc(idx_46))
    enddo
    !$acc loop gang vector independent
    do idx_47 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_47,39) = 2.00e-12 * EXP(-2440.0 / tc(idx_47))
    enddo
    !$acc loop gang vector independent
    do idx_48 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_48,40) = 1.40e-12 * EXP(-1860.0 / tc(idx_48))
    enddo
    !$acc loop gang vector independent
    do idx_49 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_49,41) = 3.46e-12 * EXP(-1862.0 / tc(idx_49))
    enddo
    !$acc loop gang vector independent
    do idx_50 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_50,42) = 3.00e-17
    enddo
    !$acc loop gang vector independent
    do idx_51 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_51,43) = 2.50e-22
    enddo
    !$acc loop gang vector independent
    do idx_52 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_52,44) = 8.00e-12 * EXP(-2060.0 / tc(idx_52))
    enddo
    !$acc loop gang vector independent
    do idx_53 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_53,45) = 1.05e-10
    enddo
    !$acc loop gang vector independent
    do idx_54 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_54,46) = 7.50e-12
    enddo
    !$acc loop gang vector independent
    do idx_55 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_55,47) = 3.45e-11
    enddo
    !$acc loop gang vector independent
    do idx_56 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_56,48) = 2.20e-10
    enddo
    !$acc loop gang vector independent
    do idx_57 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_57,49) = 2.10e-11 * EXP(115.0 / tc(idx_57))
    enddo
    !$acc loop gang vector independent
    do idx_58 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_58,50) = 3.20e-11 * EXP(67.0 / tc(idx_58))
    enddo
    !$acc loop gang vector independent
    do idx_59 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_59,51) = 1.85e-12 * EXP(-1690.0 / tc(idx_59))
    enddo
    !$acc loop gang vector independent
    do idx_60 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_60,52) = 6.90e-12 * EXP(-1000.0 / tc(idx_60))
    enddo
    !$acc loop gang vector independent
    do idx_61 = LBOUND(ratioa2b, dim=1), UBOUND(ratioa2b, dim=1), 1
      ratioa2b(idx_61) = 226.0 * tc(idx_61) ** (-0.64) * EXP(-816.0 / tc(idx_61))
    enddo
    !$acc loop gang vector independent
    do idx_62 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_62,53) = 7.60e-12 * EXP(-585.0 / tc(idx_62))
    enddo
    !$acc loop gang vector independent
    do idx_63 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_63,53) = rc(idx_63,53) * (ratioa2b(idx_63) / (ratioa2b(idx_63) + 1.0))
    enddo
    !$acc loop gang vector independent
    do idx_64 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_64,54) = 7.60e-12 * EXP(-585.0 / tc(idx_64))
    enddo
    !$acc loop gang vector independent
    do idx_65 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_65,54) = rc(idx_65,54) / (ratioa2b(idx_65) + 1.0)
    enddo
    !$acc loop gang vector independent
    do idx_66 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_66,55) = 1.44e-13 * (1.0 + m(idx_66) / 4.2e19)
    enddo
    !$acc loop gang vector independent
    do idx_67 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_67,56) = 5.10e-12 * EXP(405.0 / tc(idx_67))
    enddo
    !$acc loop gang vector independent
    do idx_68 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_68,57) = 8.01e-12
    enddo
    !$acc loop gang vector independent
    do idx_69 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_69,58) = 1.90e-12 * EXP(190.0 / tc(idx_69))
    enddo
    !$acc loop gang vector independent
    do idx_70 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_70,59) = 7.70e-12 * EXP(-2100.0 / tc(idx_70))
    enddo
    !$acc loop gang vector independent
    do idx_71 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_71,60) = 2.90e-12 * EXP(-160.0 / tc(idx_71))
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_aerchem) then
      !$acc kernels
      !$acc loop gang vector independent
      do idx_72 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_72,60) = rc(idx_72,60) * (1.0 - faq(idx_72,13) * fcloud(idx_72))
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx_73 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_73,61) = 5.40e-12 * EXP(135.0 / tc(idx_73))
    enddo
    !$acc loop gang vector independent
    do idx_74 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_74,62) = 4.80e-11 * EXP(250.0 / tc(idx_74))
    enddo
    !$acc loop gang vector independent
    do idx_75 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_75,63) = 1.90e-12 * EXP(270.0 / tc(idx_75))
    enddo
    !$acc loop gang vector independent
    do idx_76 = LBOUND(z1, dim=1), UBOUND(z1, dim=1), 1
      z1(idx_76) = 2.4e-14 * EXP(460.0 / tc(idx_76))
    enddo
    !$acc loop gang vector independent
    do idx_77 = LBOUND(z3, dim=1), UBOUND(z3, dim=1), 1
      z3(idx_77) = 6.5e-34 * EXP(1335.0 / tc(idx_77))
    enddo
    !$acc loop gang vector independent
    do idx_78 = LBOUND(z4, dim=1), UBOUND(z4, dim=1), 1
      z4(idx_78) = 2.7e-17 * EXP(2199.0 / tc(idx_78))
    enddo
    !$acc loop gang vector independent
    do idx_79 = LBOUND(z2, dim=1), UBOUND(z2, dim=1), 1
      z2(idx_79) = z3(idx_79) * m(idx_79) / (1.0 + z3(idx_79) * m(idx_79) / z4(idx_79))
    enddo
    !$acc loop gang vector independent
    do idx_80 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_80,64) = z1(idx_80) + z2(idx_80)
    enddo
    !$acc loop gang vector independent
    do idx_81 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_81,65) = 2.50e-12 * EXP(260.0 / tc(idx_81))
    enddo
    !$acc loop gang vector independent
    do idx_82 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_82,66) = 1.02e-12 * EXP(190.0 / tc(idx_82))
    enddo
    !$acc loop gang vector independent
    do idx_83 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_83,67) = 1.89e-12 * EXP(190.0 / tc(idx_83))
    enddo
    !$acc loop gang vector independent
    do idx_84 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_84,68) = 4.00e-13 * EXP(-845.0 / tc(idx_84))
    enddo
    !$acc loop gang vector independent
    do idx_85 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_85,69) = 8.80e-12 * EXP(-1320.0 / tc(idx_85))
    enddo
    !$acc loop gang vector independent
    do idx_86 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_86,70) = 1.70e-14 * EXP(420.0 / tc(idx_86))
    enddo
    !$acc loop gang vector independent
    do idx_87 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_87,71) = 1.90e-12 * EXP(190.0 / tc(idx_87))
    enddo
    !$acc loop gang vector independent
    do idx_88 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_88,72) = 8.39e-12
    enddo
    !$acc loop gang vector independent
    do idx_89 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_89,73) = 4.40e-12 * EXP(365.0 / tc(idx_89))
    enddo
    !$acc loop gang vector independent
    do idx_90 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_90,74) = 2.00e-11
    enddo
    !$acc loop gang vector independent
    do idx_91 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_91,75) = 1.70e-12 * EXP(-940.0 / tc(idx_91))
    enddo
    !$acc loop gang vector independent
    do idx_92 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_92,76) = 6.31e-14 * (tc(idx_92) / 300.0) ** 2.6 * EXP(945.0 / tc(idx_92))
    enddo
    !$acc loop gang vector independent
    do idx_93 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_93,77) = 3.00e-14
    enddo
    !$acc loop gang vector independent
    do idx_94 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_94,78) = 1.27e-12
    enddo
    !$acc loop gang vector independent
    do idx_95 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_95,79) = 1.90e-12 * EXP(190.0 / tc(idx_95))
    enddo
    !$acc loop gang vector independent
    do idx_96 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_96,80) = 1.10e-11
    enddo
    !$acc loop gang vector independent
    do idx_97 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_97,81) = 1.90e-12 * EXP(190.0 / tc(idx_97))
    enddo
    !$acc loop gang vector independent
    do idx_98 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_98,82) = 1.66e-11
    enddo
    !$acc loop gang vector independent
    do idx_99 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_99,83) = 5.50e-12 * EXP(188.0 / tc(idx_99))
    enddo
    !$acc loop gang vector independent
    do idx_100 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_100,84) = 2.03e-16 * (tc(idx_100) / 300.0) ** 4.57 * EXP(693.0 / tc(idx_100))
    enddo
    !$acc loop gang vector independent
    do idx_101 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_101,85) = 1.70e-12 * EXP(-940.0 / tc(idx_101))
    enddo
    !$acc loop gang vector independent
    do idx_102 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_102,86) = 2.20e-10
    enddo
    !$acc loop gang vector independent
    do idx_103 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_103,87) = 2.10e-11 * EXP(115.0 / tc(idx_103))
    enddo
    !$acc loop gang vector independent
    do idx_104 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_104,88) = 3.20e-11 * EXP(67.0 / tc(idx_104))
    enddo
    at(1) = 0.00
    at(2) = 1.90e-33
    at(3) = 0.00
    at(4) = -980.0
    at(5) = 0.00e+00
    at(6) = 0.00
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_105 = LBOUND(t300, dim=1), UBOUND(t300, dim=1), 1
      t300(idx_105) = tc(idx_105) / 300.0
    enddo
    !$acc loop gang vector independent
    do idx_106 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_106) = at(2) * EXP(at(3) * LOG(t300(idx_106))) * EXP(-at(4) / tc(idx_106)) * m(idx_106)
    enddo
    !$acc loop gang vector independent
    do idx_107 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_107) = at(5) * EXP(at(6) * LOG(t300(idx_107))) * EXP(-at(7) / tc(idx_107))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,89) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,89) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,89) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    !$acc loop gang vector independent
    do idx_108 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
      rc(idx_108,89) = rc(idx_108,89) * (1.0 + 1.4e-21 * h2o(idx_108) * EXP(2200.0 / tc(idx_108)))
    enddo
    at(1) = 0.60
    at(2) = 1.80e-31
    at(3) = -3.20
    at(4) = 0.0
    at(5) = 4.70e-12
    at(6) = 0.00
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_109 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_109) = at(2) * EXP(at(3) * LOG(t300(idx_109))) * EXP(-at(4) / tc(idx_109)) * m(idx_109)
    enddo
    !$acc loop gang vector independent
    do idx_110 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_110) = at(5) * EXP(at(6) * LOG(t300(idx_110))) * EXP(-at(7) / tc(idx_110))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,90) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,90) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,90) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.60
    at(2) = 4.10e-5
    at(3) = 0.0
    at(4) = 10650.0
    at(5) = 4.80e+15
    at(6) = 0.00
    at(7) = 11170.0
    !$acc loop gang vector independent
    do idx_111 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_111) = at(2) * EXP(at(3) * LOG(t300(idx_111))) * EXP(-at(4) / tc(idx_111)) * m(idx_111)
    enddo
    !$acc loop gang vector independent
    do idx_112 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_112) = at(5) * EXP(at(6) * LOG(t300(idx_112))) * EXP(-at(7) / tc(idx_112))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,91) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,91) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,91) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.30
    at(2) = 2.70e-28
    at(3) = -7.10
    at(4) = 0.0
    at(5) = 1.21e-11
    at(6) = -0.90
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_113 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_113) = at(2) * EXP(at(3) * LOG(t300(idx_113))) * EXP(-at(4) / tc(idx_113)) * m(idx_113)
    enddo
    !$acc loop gang vector independent
    do idx_114 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_114) = at(5) * EXP(at(6) * LOG(t300(idx_114))) * EXP(-at(7) / tc(idx_114))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,92) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,92) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,92) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.30
    at(2) = 4.90e-03
    at(3) = 0.00
    at(4) = 12100.0
    at(5) = 5.40e+16
    at(6) = 0.00
    at(7) = 13830.0
    !$acc loop gang vector independent
    do idx_115 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_115) = at(2) * EXP(at(3) * LOG(t300(idx_115))) * EXP(-at(4) / tc(idx_115)) * m(idx_115)
    enddo
    !$acc loop gang vector independent
    do idx_116 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_116) = at(5) * EXP(at(6) * LOG(t300(idx_116))) * EXP(-at(7) / tc(idx_116))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,93) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,93) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,93) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.35
    at(2) = 1.30e-03
    at(3) = -3.50
    at(4) = 11000.0
    at(5) = 9.70e+14
    at(6) = 0.10
    at(7) = 11080.0
    !$acc loop gang vector independent
    do idx_117 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_117) = at(2) * EXP(at(3) * LOG(t300(idx_117))) * EXP(-at(4) / tc(idx_117)) * m(idx_117)
    enddo
    !$acc loop gang vector independent
    do idx_118 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_118) = at(5) * EXP(at(6) * LOG(t300(idx_118))) * EXP(-at(7) / tc(idx_118))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,94) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,94) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,94) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.35
    at(2) = 3.60e-30
    at(3) = -4.10
    at(4) = 0.0
    at(5) = 1.90e-12
    at(6) = 0.20
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_119 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_119) = at(2) * EXP(at(3) * LOG(t300(idx_119))) * EXP(-at(4) / tc(idx_119)) * m(idx_119)
    enddo
    !$acc loop gang vector independent
    do idx_120 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_120) = at(5) * EXP(at(6) * LOG(t300(idx_120))) * EXP(-at(7) / tc(idx_120))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,95) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,95) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,95) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.00
    at(2) = 5.70e-34
    at(3) = -2.60
    at(4) = 0.0
    at(5) = 0.00e+00
    at(6) = 0.00
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_121 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_121) = at(2) * EXP(at(3) * LOG(t300(idx_121))) * EXP(-at(4) / tc(idx_121)) * m(idx_121)
    enddo
    !$acc loop gang vector independent
    do idx_122 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_122) = at(5) * EXP(at(6) * LOG(t300(idx_122))) * EXP(-at(7) / tc(idx_122))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,96) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,96) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,96) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 1420.0
    at(2) = 7.40e-31
    at(3) = -2.40
    at(4) = 0.0
    at(5) = 3.30e-11
    at(6) = -0.30
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_123 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_123) = at(2) * EXP(at(3) * LOG(t300(idx_123))) * EXP(-at(4) / tc(idx_123)) * m(idx_123)
    enddo
    !$acc loop gang vector independent
    do idx_124 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_124) = at(5) * EXP(at(6) * LOG(t300(idx_124))) * EXP(-at(7) / tc(idx_124))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,97) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,97) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,97) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.40
    at(2) = 3.30e-30
    at(3) = -3.00
    at(4) = 0.0
    at(5) = 4.10e-11
    at(6) = 0.00
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_125 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_125) = at(2) * EXP(at(3) * LOG(t300(idx_125))) * EXP(-at(4) / tc(idx_125)) * m(idx_125)
    enddo
    !$acc loop gang vector independent
    do idx_126 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_126) = at(5) * EXP(at(6) * LOG(t300(idx_126))) * EXP(-at(7) / tc(idx_126))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,98) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,98) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,98) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.50
    at(2) = 6.90e-31
    at(3) = -0.80
    at(4) = 0.0
    at(5) = 2.60e-11
    at(6) = 0.00
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_127 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_127) = at(2) * EXP(at(3) * LOG(t300(idx_127))) * EXP(-at(4) / tc(idx_127)) * m(idx_127)
    enddo
    !$acc loop gang vector independent
    do idx_128 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_128) = at(5) * EXP(at(6) * LOG(t300(idx_128))) * EXP(-at(7) / tc(idx_128))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,99) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,99) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,99) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.30
    at(2) = 2.70e-28
    at(3) = -7.10
    at(4) = 0.0
    at(5) = 1.20e-11
    at(6) = -0.90
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_129 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_129) = at(2) * EXP(at(3) * LOG(t300(idx_129))) * EXP(-at(4) / tc(idx_129)) * m(idx_129)
    enddo
    !$acc loop gang vector independent
    do idx_130 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_130) = at(5) * EXP(at(6) * LOG(t300(idx_130))) * EXP(-at(7) / tc(idx_130))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,100) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,100) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,100) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.36
    at(2) = 1.70e-03
    at(3) = 0.00
    at(4) = 11280.0
    at(5) = 8.30e+16
    at(6) = 0.00
    at(7) = 13940.0
    !$acc loop gang vector independent
    do idx_131 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_131) = at(2) * EXP(at(3) * LOG(t300(idx_131))) * EXP(-at(4) / tc(idx_131)) * m(idx_131)
    enddo
    !$acc loop gang vector independent
    do idx_132 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_132) = at(5) * EXP(at(6) * LOG(t300(idx_132))) * EXP(-at(7) / tc(idx_132))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,101) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,101) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,101) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    at(1) = 0.00
    at(2) = 5.70e-34
    at(3) = -2.60
    at(4) = 0.0
    at(5) = 0.00e+00
    at(6) = 0.00
    at(7) = 0.0
    !$acc loop gang vector independent
    do idx_133 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
      zo(idx_133) = at(2) * EXP(at(3) * LOG(t300(idx_133))) * EXP(-at(4) / tc(idx_133)) * m(idx_133)
    enddo
    !$acc loop gang vector independent
    do idx_134 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
      zi(idx_134) = at(5) * EXP(at(6) * LOG(t300(idx_134))) * EXP(-at(7) / tc(idx_134))
    enddo
    !$acc loop gang vector independent
    do i = 1, n_pnts, 1
      if (zo(i) < peps) then
        rc(i,102) = zi(i)
      else
        if (zi(i) < peps) then
          rc(i,102) = zo(i)
        else
          if (at(1) <= 1.0) then
            zfc(i) = at(1)
          else
            zfc(i) = EXP(-tc(i) / at(1))
          end if
          zr(i) = zo(i) / zi(i)
          arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
          rc(i,102) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
        end if
      end if
    enddo
    !$acc end kernels
    if (ukca_config%l_ukca_aerchem) then
      !$acc kernels
      at(1) = 0.60
      at(2) = 3.00e-31
      at(3) = -3.30
      at(4) = 0.0
      at(5) = 1.50e-12
      at(6) = 0.00
      at(7) = 0.0
      !$acc loop gang vector independent
      do idx_135 = LBOUND(zo, dim=1), UBOUND(zo, dim=1), 1
        zo(idx_135) = at(2) * EXP(at(3) * LOG(t300(idx_135))) * EXP(-at(4) / tc(idx_135)) * m(idx_135)
      enddo
      !$acc loop gang vector independent
      do idx_136 = LBOUND(zi, dim=1), UBOUND(zi, dim=1), 1
        zi(idx_136) = at(5) * EXP(at(6) * LOG(t300(idx_136))) * EXP(-at(7) / tc(idx_136))
      enddo
      !$acc loop gang vector independent
      do i = 1, n_pnts, 1
        if (zo(i) < peps) then
          rc(i,103) = zi(i)
        else
          if (zi(i) < peps) then
            rc(i,103) = zo(i)
          else
            if (at(1) <= 1.0) then
              zfc(i) = at(1)
            else
              zfc(i) = EXP(-tc(i) / at(1))
            end if
            zr(i) = zo(i) / zi(i)
            arg2(i) = 1.0 / (1.0 + LOG10(zr(i)) ** 2)
            rc(i,103) = zo(i) / (1.0 + zr(i)) * EXP(arg2(i) * LOG(zfc(i)))
            rc(i,103) = rc(i,103) * (1.0 - faq(i,48) * fcloud(i))
          end if
        end if
      enddo
      !$acc end kernels
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (clw(:) < 1E-10)
        rc(:, 104) = 0.0
      ELSEWHERE
        rc(:, 104) = 1.10E+12 * (EXP(- 3655.3 / tc(:))) * (H_plus_1d_arr(:) / (0.1 + H_plus_1d_arr(:))) * fcloud(:) * fdiss(:, 48, 2) * fdiss(:, 13, 1) * 1000.0 / (avogadro * vr(:))
      END WHERE
      !$acc kernels
      !$acc loop gang vector independent
      do idx_137 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_137,105) = 1.12e-11 * EXP(-250.0 / tc(idx_137))
      enddo
      !$acc loop gang vector independent
      do idx_138 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_138,106) = 9.3e-39 * o2(idx_138) * EXP(5270.0 / tc(idx_138)) / (1.0 + 7.4e-29 * o2(idx_138) * EXP(5610.0 / tc(idx_138)))
      enddo
      !$acc loop gang vector independent
      do idx_139 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_139,107) = 1.9e-13 * EXP(520.0 / tc(idx_139))
      enddo
      !$acc loop gang vector independent
      do idx_140 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_140,111) = 3.5e-12 * EXP(-925.0 / tc(idx_140))
      enddo
      !$acc loop gang vector independent
      do idx_141 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_141,112) = 1.2e-11 * EXP(444.0 / tc(idx_141))
      enddo
      !$acc loop gang vector independent
      do idx_142 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_142,113) = 1.01e-15 * EXP(-732.0 / tc(idx_142))
      enddo
      !$acc loop gang vector independent
      do idx_143 = LBOUND(rc, dim=1), UBOUND(rc, dim=1), 1
        rc(idx_143,114) = 1.19e-12 * EXP(490.0 / tc(idx_143))
      enddo
      !$acc end kernels
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (clw(:) < 1E-10)
        rc(:, 115) = 0.0
      ELSEWHERE
        rc(:, 115) = (7.43E+16 * (EXP(- 5280.0 / tc(:))) * fdiss(:, 48, 3)) * fcloud(:) * fdiss(:, 3, 1) * 1000.0 / (avogadro * vr(:))
      END WHERE
      !$acc kernels
      !$acc loop gang vector independent
      do idx_144 = LBOUND(k_dms, dim=1), UBOUND(k_dms, dim=1), 1
        k_dms(idx_144,1) = 100.0
      enddo
      !$acc loop gang vector independent
      do idx_145 = LBOUND(k_dms, dim=1), UBOUND(k_dms, dim=1), 1
        k_dms(idx_145,2) = 6.3e-13
      enddo
      !$acc loop gang vector independent
      do idx_146 = LBOUND(k_dms, dim=1), UBOUND(k_dms, dim=1), 1
        k_dms(idx_146,3) = 2.2e-11
      enddo
      !$acc loop gang vector independent
      do idx_147 = LBOUND(k_dms, dim=1), UBOUND(k_dms, dim=1), 1
        k_dms(idx_147,4) = 5.0e-11
      enddo
      !$acc loop gang vector independent
      do idx_148 = LBOUND(k_dms, dim=1), UBOUND(k_dms, dim=1), 1
        k_dms(idx_148,5) = 1.2e-3
      enddo
      !$acc end kernels
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_chemco

end module ukca_chemco_mod
