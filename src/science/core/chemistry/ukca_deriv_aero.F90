module ukca_deriv_aero_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DERIV_AERO_MOD'
  public

  contains
  subroutine ukca_deriv_aero(nr_therm, n_be_calls, theta_field_size, rc, wdep, ddep, dj, k_dms, h2o, m, o2, dts, y, delta_so2_wetox_h2o2, delta_so2_wetox_o3, delta_so2_dryox_oh)
    use asad_mod, only : jppj, jpspec
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr_therm
    integer, intent(in) :: n_be_calls
    integer, intent(in) :: theta_field_size
    real, dimension(theta_field_size,nr_therm), intent(in) :: rc
    real, dimension(theta_field_size,jpspec), intent(in) :: wdep
    real, dimension(theta_field_size,jpspec), intent(in) :: ddep
    real, dimension(theta_field_size,jppj), intent(in) :: dj
    real, dimension(theta_field_size,5), intent(in) :: k_dms
    real, dimension(theta_field_size), intent(in) :: h2o
    real, dimension(theta_field_size), intent(in) :: m
    real, dimension(theta_field_size), intent(in) :: o2
    real, intent(in) :: dts
    real, dimension(theta_field_size,jpspec), intent(inout) :: y
    real, dimension(theta_field_size), intent(out) :: delta_so2_wetox_h2o2
    real, dimension(theta_field_size), intent(out) :: delta_so2_wetox_o3
    real, dimension(theta_field_size), intent(out) :: delta_so2_dryox_oh
    real, dimension(theta_field_size) :: p1
    real, dimension(theta_field_size) :: p2
    real, dimension(theta_field_size) :: p
    real, dimension(theta_field_size) :: l
    real, dimension(theta_field_size) :: l1
    real, dimension(theta_field_size) :: l2
    real, dimension(theta_field_size) :: l3
    real, dimension(theta_field_size) :: r1
    real, dimension(theta_field_size) :: r2
    real, dimension(theta_field_size) :: f_so2
    real, dimension(theta_field_size) :: f_so4
    real, dimension(theta_field_size) :: f_msa
    real, dimension(theta_field_size,jpspec) :: yp
    integer :: i
    integer :: j
    integer :: n
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DERIV_AERO'
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
    integer :: idx_149
    integer :: idx_150
    integer :: idx_151
    integer :: idx_152
    integer :: idx_153
    integer :: idx_154
    integer :: idx_155
    integer :: idx_156
    integer :: idx_157
    integer :: idx_158
    integer :: idx_159
    integer :: idx_160
    integer :: idx_161
    integer :: idx_162
    integer :: idx_163
    integer :: idx_164
    integer :: idx_165
    integer :: idx_166
    integer :: idx_167
    integer :: idx_168
    integer :: idx_169
    integer :: idx_170
    integer :: idx_171

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(delta_so2_wetox_h2o2, dim=1), UBOUND(delta_so2_wetox_h2o2, dim=1), 1
      delta_so2_wetox_h2o2(idx) = 0.0
    enddo
    !$acc loop gang vector independent
    do idx_1 = LBOUND(delta_so2_wetox_o3, dim=1), UBOUND(delta_so2_wetox_o3, dim=1), 1
      delta_so2_wetox_o3(idx_1) = 0.0
    enddo
    !$acc loop gang vector independent
    do idx_2 = LBOUND(delta_so2_dryox_oh, dim=1), UBOUND(delta_so2_dryox_oh, dim=1), 1
      delta_so2_dryox_oh(idx_2) = 0.0
    enddo
    !$acc loop seq
    do n = 1, n_be_calls, 1
      !$acc loop gang vector independent collapse(2)
      do j = 1, jpspec, 1
        do i = 1, theta_field_size, 1
          yp(i,j) = y(i,j)
        enddo
      enddo
      !$acc loop gang vector independent
      do idx_3 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_3,12) = 5.0e-7 * m(idx_3)
      enddo
      !$acc loop gang vector independent
      do idx_4 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_4,16) = 350.0e-6 * m(idx_4)
      enddo
      !$acc loop gang vector independent
      do idx_5 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_5,19) = h2o(idx_5)
      enddo
      !$acc loop gang vector independent
      do idx_6 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_6,22) = o2(idx_6)
      enddo
      !$acc loop gang vector independent
      do idx_7 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_7,23) = m(idx_7) - o2(idx_7)
      enddo
      !$acc loop gang vector independent
      do idx_8 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_8,12) = yp(idx_8,12)
      enddo
      !$acc loop gang vector independent
      do idx_9 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_9,16) = yp(idx_9,16)
      enddo
      !$acc loop gang vector independent
      do idx_10 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_10,19) = yp(idx_10,19)
      enddo
      !$acc loop gang vector independent
      do idx_11 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_11,22) = yp(idx_11,22)
      enddo
      !$acc loop gang vector independent
      do idx_12 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_12,23) = yp(idx_12,23)
      enddo
      !$acc loop seq
      do i = 1, ukca_config%nit, 1
        !$acc loop gang vector independent
        do idx_13 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_13) = +dj(idx_13,15) * y(idx_13,3) + dj(idx_13,12) * y(idx_13,5) + dj(idx_13,13) * y(idx_13,22) * 2.00 + rc(idx_13,76) * y(idx_13,10) * y(idx_13,10) + dj(idx_13,10) * y(idx_13,6) + rc(idx_13,49) * y(idx_13,2) * y(idx_13,23) + rc(idx_13,50) * y(idx_13,2) * y(idx_13,22)
        enddo
        !$acc loop gang vector independent
        do idx_14 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_14) = ddep(idx_14,1) + wdep(idx_14,1) + rc(idx_14,44) * y(idx_14,3) + rc(idx_14,83) * y(idx_14,6) + rc(idx_14,96) * y(idx_14,22)
        enddo
        !$acc loop gang vector independent
        do idx_15 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_15,1) = p(idx_15) / l(idx_15)
        enddo
        !$acc loop gang vector independent
        do idx_16 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_16) = +dj(idx_16,14) * y(idx_16,3)
        enddo
        !$acc loop gang vector independent
        do idx_17 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_17) = ddep(idx_17,2) + wdep(idx_17,2) + rc(idx_17,48) * y(idx_17,19) + rc(idx_17,49) * y(idx_17,23) + rc(idx_17,50) * y(idx_17,22) + rc(idx_17,45) * y(idx_17,14) + rc(idx_17,46) * y(idx_17,14) + rc(idx_17,47) * y(idx_17,14)
        enddo
        !$acc loop gang vector independent
        do idx_18 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_18,2) = p(idx_18) / l(idx_18)
        enddo
        !$acc loop gang vector independent
        do idx_19 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_19) = +dj(idx_19,1) * y(idx_19,34) + dj(idx_19,1) * y(idx_19,39) + dj(idx_19,17) * y(idx_19,21) + dj(idx_19,1) * y(idx_19,33) + dj(idx_19,6) * y(idx_19,9) + dj(idx_19,1) * y(idx_19,20) + dj(idx_19,1) * y(idx_19,26) + dj(idx_19,2) * y(idx_19,13) * 2.00 + rc(idx_19,82) * y(idx_19,10) * y(idx_19,34) + rc(idx_19,85) * y(idx_19,45) * y(idx_19,44) + rc(idx_19,72) * y(idx_19,10) * y(idx_19,39) + rc(idx_19,80) * y(idx_19,10) * y(idx_19,33) + rc(idx_19,57) * y(idx_19,10) * y(idx_19,26) + rc(idx_19,66) * y(idx_19,10) * y(idx_19,20) + rc(idx_19,45) * y(idx_19,2) * y(idx_19,14) + rc(idx_19,48) * y(idx_19,2) * y(idx_19,19) * 2.00 + rc(idx_19,3) * y(idx_19,11) * y(idx_19,3) + rc(idx_19,10) * y(idx_19,11) * y(idx_19,28) + rc(idx_19,1) * y(idx_19,11) * y(idx_19,4) + rc(idx_19,2) * y(idx_19,11) * y(idx_19,5)
        enddo
        !$acc loop gang vector independent
        do idx_20 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_20) = ddep(idx_20,10) + wdep(idx_20,10) + rc(idx_20,99) * y(idx_20,10) + rc(idx_20,97) * y(idx_20,4) + rc(idx_20,98) * y(idx_20,6) + rc(idx_20,99) * y(idx_20,10) + rc(idx_20,80) * y(idx_20,33) + rc(idx_20,81) * y(idx_20,34) + rc(idx_20,82) * y(idx_20,34) + rc(idx_20,77) * y(idx_20,29) + rc(idx_20,78) * y(idx_20,40) + rc(idx_20,79) * y(idx_20,33) + rc(idx_20,75) * y(idx_20,3) + rc(idx_20,76) * y(idx_20,10) + rc(idx_20,76) * y(idx_20,10) + rc(idx_20,72) * y(idx_20,39) + rc(idx_20,73) * y(idx_20,27) + rc(idx_20,74) * y(idx_20,5) + rc(idx_20,69) * y(idx_20,37) + rc(idx_20,70) * y(idx_20,37) + rc(idx_20,71) * y(idx_20,39) + rc(idx_20,66) * y(idx_20,20) + rc(idx_20,67) * y(idx_20,20) + rc(idx_20,68) * y(idx_20,41) + rc(idx_20,63) * y(idx_20,8) + rc(idx_20,64) * y(idx_20,9) + rc(idx_20,65) * y(idx_20,21) + rc(idx_20,60) * y(idx_20,13) + rc(idx_20,61) * y(idx_20,17) + rc(idx_20,62) * y(idx_20,11) + rc(idx_20,57) * y(idx_20,26) + rc(idx_20,58) * y(idx_20,26) + rc(idx_20,59) * y(idx_20,12) + rc(idx_20,54) * y(idx_20,30) + rc(idx_20,55) * y(idx_20,15) + rc(idx_20,56) * y(idx_20,35) + rc(idx_20,51) * y(idx_20,14) + rc(idx_20,52) * y(idx_20,24) + rc(idx_20,53) * y(idx_20,30) + rc(idx_20,103) * y(idx_20,48) + rc(idx_20,111) * y(idx_20,51) + rc(idx_20,112) * y(idx_20,52) + rc(idx_20,105) * y(idx_20,47) + rc(idx_20,106) * y(idx_20,47)
        enddo
        !$acc loop gang vector independent
        do idx_21 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_21,10) = p(idx_21) / l(idx_21)
        enddo
        !$acc loop gang vector independent
        do idx_22 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_22) = +rc(idx_22,96) * y(idx_22,1) * y(idx_22,22) + rc(idx_22,9) * y(idx_22,11) * y(idx_22,28) + rc(idx_22,14) * y(idx_22,11) * y(idx_22,36)
        enddo
        !$acc loop gang vector independent
        do idx_23 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_23) = ddep(idx_23,3) + wdep(idx_23,3) + dj(idx_23,15) + rc(idx_23,44) * y(idx_23,1) + rc(idx_23,75) * y(idx_23,10) + dj(idx_23,14) + rc(idx_23,3) * y(idx_23,11) + rc(idx_23,37) * y(idx_23,4) + rc(idx_23,38) * y(idx_23,6) + rc(idx_23,115) * y(idx_23,48) + rc(idx_23,113) * y(idx_23,52)
        enddo
        !$acc loop gang vector independent
        do idx_24 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_24,3) = (yp(idx_24,3) + dts * p(idx_24)) / (1.0 + dts * l(idx_24))
        enddo
        !$acc loop gang vector independent
        do idx_25 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_25) = +dj(idx_25,11) * y(idx_25,5) + dj(idx_25,17) * y(idx_25,21) + rc(idx_25,83) * y(idx_25,1) * y(idx_25,6) + dj(idx_25,10) * y(idx_25,6) - rc(idx_25,111) * y(idx_25,51) * y(idx_25,10)
        enddo
        !$acc loop gang vector independent
        do idx_26 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_26) = ddep(idx_26,4) + wdep(idx_26,4) + rc(idx_26,36) * y(idx_26,5) + rc(idx_26,37) * y(idx_26,3) + rc(idx_26,97) * y(idx_26,10) + rc(idx_26,30) * y(idx_26,32) + rc(idx_26,32) * y(idx_26,36) + rc(idx_26,34) * y(idx_26,38) + rc(idx_26,23) * y(idx_26,25) + rc(idx_26,26) * y(idx_26,28) + rc(idx_26,28) * y(idx_26,31) + rc(idx_26,1) * y(idx_26,11) + rc(idx_26,16) * y(idx_26,18) + rc(idx_26,17) * y(idx_26,18)
        enddo
        !$acc loop gang vector independent
        do idx_27 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_27,4) = (yp(idx_27,4) + dts * p(idx_27)) / (1.0 + dts * l(idx_27))
        enddo
        !$acc loop gang vector independent
        do idx_28 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_28) = +dj(idx_28,20) * y(idx_28,41) + dj(idx_28,16) * y(idx_28,29) + dj(idx_28,16) * y(idx_28,40) + dj(idx_28,9) * y(idx_28,7) + dj(idx_28,12) * y(idx_28,5) + dj(idx_28,5) * y(idx_28,8) + dj(idx_28,6) * y(idx_28,9) + rc(idx_28,94) * y(idx_28,7) + rc(idx_28,101) * y(idx_28,40) + rc(idx_28,91) * y(idx_28,8) + rc(idx_28,93) * y(idx_28,29) + rc(idx_28,77) * y(idx_28,10) * y(idx_28,29) + rc(idx_28,78) * y(idx_28,10) * y(idx_28,40) + rc(idx_28,68) * y(idx_28,10) * y(idx_28,41) + rc(idx_28,74) * y(idx_28,10) * y(idx_28,5) + rc(idx_28,63) * y(idx_28,10) * y(idx_28,8) + rc(idx_28,65) * y(idx_28,10) * y(idx_28,21) + rc(idx_28,36) * y(idx_28,4) * y(idx_28,5) * 2.00 + rc(idx_28,37) * y(idx_28,4) * y(idx_28,3) + rc(idx_28,34) * y(idx_28,38) * y(idx_28,4) + rc(idx_28,35) * y(idx_28,38) * y(idx_28,5) + rc(idx_28,32) * y(idx_28,36) * y(idx_28,4) + rc(idx_28,33) * y(idx_28,36) * y(idx_28,5) + rc(idx_28,30) * y(idx_28,32) * y(idx_28,4) + rc(idx_28,31) * y(idx_28,32) * y(idx_28,5) + rc(idx_28,28) * y(idx_28,31) * y(idx_28,4) + rc(idx_28,29) * y(idx_28,31) * y(idx_28,5) + rc(idx_28,26) * y(idx_28,28) * y(idx_28,4) + rc(idx_28,27) * y(idx_28,28) * y(idx_28,5) + rc(idx_28,23) * y(idx_28,25) * y(idx_28,4) + rc(idx_28,24) * y(idx_28,25) * y(idx_28,5) + rc(idx_28,16) * y(idx_28,18) * y(idx_28,4) + rc(idx_28,18) * y(idx_28,18) * y(idx_28,5) + rc(idx_28,1) * y(idx_28,11) * y(idx_28,4) + rc(idx_28,2) * y(idx_28,11) * y(idx_28,5)
        enddo
        !$acc loop gang vector independent
        do idx_29 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_29) = ddep(idx_29,6) + wdep(idx_29,6) + rc(idx_29,100) * y(idx_29,36) + dj(idx_29,10) + rc(idx_29,92) * y(idx_29,28) + rc(idx_29,95) * y(idx_29,5) + rc(idx_29,98) * y(idx_29,10) + rc(idx_29,38) * y(idx_29,3) + rc(idx_29,83) * y(idx_29,1) + rc(idx_29,90) * y(idx_29,11)
        enddo
        !$acc loop gang vector independent
        do idx_30 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_30,6) = (yp(idx_30,6) + dts * p(idx_30)) / (1.0 + dts * l(idx_30))
        enddo
        !$acc loop gang vector independent
        do idx_31 = LBOUND(p1, dim=1), UBOUND(p1, dim=1), 1
          p1(idx_31) = +rc(idx_31,38) * y(idx_31,6) * y(idx_31,3) + rc(idx_31,64) * y(idx_31,10) * y(idx_31,9)
        enddo
        !$acc loop gang vector independent
        do idx_32 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_32) = ddep(idx_32,5) + wdep(idx_32,5) + dj(idx_32,11) + dj(idx_32,12) + rc(idx_32,42) * y(idx_32,37) + rc(idx_32,74) * y(idx_32,10) + rc(idx_32,95) * y(idx_32,6) + rc(idx_32,39) * y(idx_32,17) + rc(idx_32,40) * y(idx_32,27) + rc(idx_32,41) * y(idx_32,35) + rc(idx_32,33) * y(idx_32,36) + rc(idx_32,35) * y(idx_32,38) + rc(idx_32,36) * y(idx_32,4) + rc(idx_32,27) * y(idx_32,28) + rc(idx_32,29) * y(idx_32,31) + rc(idx_32,31) * y(idx_32,32) + rc(idx_32,2) * y(idx_32,11) + rc(idx_32,18) * y(idx_32,18) + rc(idx_32,24) * y(idx_32,25) + rc(idx_32,107) * y(idx_32,47) + rc(idx_32,114) * y(idx_32,52)
        enddo
        !$acc loop gang vector independent
        do idx_33 = LBOUND(p2, dim=1), UBOUND(p2, dim=1), 1
          p2(idx_33) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_34 = LBOUND(r1, dim=1), UBOUND(r1, dim=1), 1
          r1(idx_34) = rc(idx_34,94) + dj(idx_34,9)
        enddo
        !$acc loop gang vector independent
        do idx_35 = LBOUND(r2, dim=1), UBOUND(r2, dim=1), 1
          r2(idx_35) = rc(idx_35,95) * y(idx_35,6)
        enddo
        !$acc loop gang vector independent
        do idx_36 = LBOUND(l1, dim=1), UBOUND(l1, dim=1), 1
          l1(idx_36) = ddep(idx_36,7) + wdep(idx_36,7) + rc(idx_36,43) * y(idx_36,19) + rc(idx_36,94) + dj(idx_36,9)
        enddo
        !$acc loop gang vector independent
        do idx_37 = LBOUND(l2, dim=1), UBOUND(l2, dim=1), 1
          l2(idx_37) = 1.0 + l(idx_37) * dts
        enddo
        !$acc loop gang vector independent
        do idx_38 = LBOUND(l3, dim=1), UBOUND(l3, dim=1), 1
          l3(idx_38) = 1.0 + l1(idx_38) * dts
        enddo
        !$acc loop gang vector independent
        do idx_39 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_39,5) = (l3(idx_39) * (yp(idx_39,5) + p1(idx_39) * dts) + r1(idx_39) * dts * yp(idx_39,7)) / (l3(idx_39) * l2(idx_39) - r1(idx_39) * r2(idx_39) * dts ** 2)
        enddo
        !$acc loop gang vector independent
        do idx_40 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_40,7) = (yp(idx_40,7) + p2(idx_40) * dts + r2(idx_40) * dts * y(idx_40,5)) / l3(idx_40)
        enddo
        !$acc loop gang vector independent
        do idx_41 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_41) = +rc(idx_41,90) * y(idx_41,11) * y(idx_41,6)
        enddo
        !$acc loop gang vector independent
        do idx_42 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_42) = ddep(idx_42,8) + wdep(idx_42,8) + rc(idx_42,63) * y(idx_42,10) + rc(idx_42,91) + dj(idx_42,5)
        enddo
        !$acc loop gang vector independent
        do idx_43 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_43,8) = (yp(idx_43,8) + dts * p(idx_43)) / (1.0 + dts * l(idx_43))
        enddo
        !$acc loop gang vector independent
        do idx_44 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_44) = +rc(idx_44,43) * y(idx_44,7) * y(idx_44,19) * 2.00 + rc(idx_44,98) * y(idx_44,10) * y(idx_44,6) + rc(idx_44,41) * y(idx_44,5) * y(idx_44,35) + rc(idx_44,42) * y(idx_44,5) * y(idx_44,37) + rc(idx_44,39) * y(idx_44,5) * y(idx_44,17) + rc(idx_44,40) * y(idx_44,5) * y(idx_44,27) + rc(idx_44,107) * y(idx_44,5) * y(idx_44,47)
        enddo
        !$acc loop gang vector independent
        do idx_45 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_45) = ddep(idx_45,9) + wdep(idx_45,9) + rc(idx_45,64) * y(idx_45,10) + dj(idx_45,6)
        enddo
        !$acc loop gang vector independent
        do idx_46 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_46,9) = (yp(idx_46,9) + dts * p(idx_46)) / (1.0 + dts * l(idx_46))
        enddo
        !$acc loop gang vector independent
        do idx_47 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_47) = +dj(idx_47,1) * y(idx_47,34) + dj(idx_47,20) * y(idx_47,41) + dj(idx_47,18) * y(idx_47,35) + dj(idx_47,1) * y(idx_47,33) + dj(idx_47,7) * y(idx_47,27) + dj(idx_47,1) * y(idx_47,20) + dj(idx_47,3) * y(idx_47,17) * 2.00 + dj(idx_47,5) * y(idx_47,8) + rc(idx_47,91) * y(idx_47,8) + dj(idx_47,1) * y(idx_47,26) + rc(idx_47,75) * y(idx_47,10) * y(idx_47,3) + rc(idx_47,84) * y(idx_47,46) * y(idx_47,44) + rc(idx_47,61) * y(idx_47,10) * y(idx_47,17) + rc(idx_47,74) * y(idx_47,10) * y(idx_47,5) + rc(idx_47,59) * y(idx_47,10) * y(idx_47,12) + rc(idx_47,60) * y(idx_47,10) * y(idx_47,13) + rc(idx_47,47) * y(idx_47,2) * y(idx_47,14) * 2.00 + rc(idx_47,55) * y(idx_47,10) * y(idx_47,15) + rc(idx_47,31) * y(idx_47,32) * y(idx_47,5) + rc(idx_47,39) * y(idx_47,5) * y(idx_47,17) + rc(idx_47,29) * y(idx_47,31) * y(idx_47,5) + rc(idx_47,30) * y(idx_47,32) * y(idx_47,4) + rc(idx_47,25) * y(idx_47,25) * y(idx_47,28) + rc(idx_47,28) * y(idx_47,31) * y(idx_47,4) + rc(idx_47,23) * y(idx_47,25) * y(idx_47,4) + rc(idx_47,24) * y(idx_47,25) * y(idx_47,5) + rc(idx_47,20) * y(idx_47,18) * y(idx_47,18) * 2.00 + rc(idx_47,21) * y(idx_47,18) * y(idx_47,28) + rc(idx_47,16) * y(idx_47,18) * y(idx_47,4) + rc(idx_47,18) * y(idx_47,18) * y(idx_47,5)
        enddo
        !$acc loop gang vector independent
        do idx_48 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_48) = ddep(idx_48,11) + wdep(idx_48,11) + rc(idx_48,89) * y(idx_48,11) + rc(idx_48,90) * y(idx_48,6) + rc(idx_48,15) * y(idx_48,38) + rc(idx_48,62) * y(idx_48,10) + rc(idx_48,89) * y(idx_48,11) + rc(idx_48,12) * y(idx_48,32) + rc(idx_48,13) * y(idx_48,36) + rc(idx_48,14) * y(idx_48,36) + rc(idx_48,9) * y(idx_48,28) + rc(idx_48,10) * y(idx_48,28) + rc(idx_48,11) * y(idx_48,31) + rc(idx_48,6) * y(idx_48,18) + rc(idx_48,7) * y(idx_48,25) + rc(idx_48,8) * y(idx_48,28) + rc(idx_48,4) * y(idx_48,11) + rc(idx_48,4) * y(idx_48,11) + rc(idx_48,5) * y(idx_48,18) + rc(idx_48,1) * y(idx_48,4) + rc(idx_48,2) * y(idx_48,5) + rc(idx_48,3) * y(idx_48,3)
        enddo
        !$acc loop gang vector independent
        do idx_49 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_49,11) = (yp(idx_49,11) + dts * p(idx_49)) / (1.0 + dts * l(idx_49))
        enddo
        !$acc loop gang vector independent
        do idx_50 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_50,12) = yp(idx_50,12)
        enddo
        !$acc loop gang vector independent
        do idx_51 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_51) = +rc(idx_51,99) * y(idx_51,10) * y(idx_51,10) + rc(idx_51,4) * y(idx_51,11) * y(idx_51,11) + rc(idx_51,89) * y(idx_51,11) * y(idx_51,11)
        enddo
        !$acc loop gang vector independent
        do idx_52 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_52) = ddep(idx_52,13) + wdep(idx_52,13) + rc(idx_52,60) * y(idx_52,10) + dj(idx_52,2) + rc(idx_52,104) * y(idx_52,48)
        enddo
        !$acc loop gang vector independent
        do idx_53 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_53,13) = (yp(idx_53,13) + dts * p(idx_53)) / (1.0 + dts * l(idx_53))
        enddo
        !$acc loop gang vector independent
        do idx_54 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_54) = +dj(idx_54,8) * y(idx_54,27)
        enddo
        !$acc loop gang vector independent
        do idx_55 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_55) = ddep(idx_55,14) + wdep(idx_55,14) + rc(idx_55,51) * y(idx_55,10) + rc(idx_55,45) * y(idx_55,2) + rc(idx_55,46) * y(idx_55,2) + rc(idx_55,47) * y(idx_55,2)
        enddo
        !$acc loop gang vector independent
        do idx_56 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_56,14) = (yp(idx_56,14) + dts * p(idx_56)) / (1.0 + dts * l(idx_56))
        enddo
        !$acc loop gang vector independent
        do idx_57 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_57) = +dj(idx_57,18) * y(idx_57,35) + dj(idx_57,7) * y(idx_57,27) + dj(idx_57,8) * y(idx_57,27) + dj(idx_57,3) * y(idx_57,17) + dj(idx_57,4) * y(idx_57,17) + rc(idx_57,39) * y(idx_57,5) * y(idx_57,17) + rc(idx_57,61) * y(idx_57,10) * y(idx_57,17)
        enddo
        !$acc loop gang vector independent
        do idx_58 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_58) = ddep(idx_58,15) + wdep(idx_58,15) + rc(idx_58,55) * y(idx_58,10)
        enddo
        !$acc loop gang vector independent
        do idx_59 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_59,15) = (yp(idx_59,15) + dts * p(idx_59)) / (1.0 + dts * l(idx_59))
        enddo
        !$acc loop gang vector independent
        do idx_60 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_60,16) = yp(idx_60,16)
        enddo
        !$acc loop gang vector independent
        do idx_61 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_61) = +dj(idx_61,20) * y(idx_61,41) + dj(idx_61,1) * y(idx_61,20) + dj(idx_61,1) * y(idx_61,39) + rc(idx_61,68) * y(idx_61,10) * y(idx_61,41) + rc(idx_61,77) * y(idx_61,10) * y(idx_61,29) + rc(idx_61,47) * y(idx_61,2) * y(idx_61,14) + rc(idx_61,66) * y(idx_61,10) * y(idx_61,20) + rc(idx_61,35) * y(idx_61,38) * y(idx_61,5) + rc(idx_61,46) * y(idx_61,2) * y(idx_61,14) + rc(idx_61,22) * y(idx_61,18) * y(idx_61,28) + rc(idx_61,34) * y(idx_61,38) * y(idx_61,4) + rc(idx_61,20) * y(idx_61,18) * y(idx_61,18) * 2.00 + rc(idx_61,21) * y(idx_61,18) * y(idx_61,28) + rc(idx_61,18) * y(idx_61,18) * y(idx_61,5) + rc(idx_61,19) * y(idx_61,18) * y(idx_61,18) + rc(idx_61,6) * y(idx_61,11) * y(idx_61,18) + rc(idx_61,16) * y(idx_61,18) * y(idx_61,4)
        enddo
        !$acc loop gang vector independent
        do idx_62 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_62) = ddep(idx_62,17) + wdep(idx_62,17) + dj(idx_62,4) + rc(idx_62,39) * y(idx_62,5) + rc(idx_62,61) * y(idx_62,10) + dj(idx_62,3)
        enddo
        !$acc loop gang vector independent
        do idx_63 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_63,17) = (yp(idx_63,17) + dts * p(idx_63)) / (1.0 + dts * l(idx_63))
        enddo
        !$acc loop gang vector independent
        do idx_64 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_64) = +dj(idx_64,7) * y(idx_64,27) + dj(idx_64,19) * y(idx_64,37) + rc(idx_64,51) * y(idx_64,10) * y(idx_64,14) + rc(idx_64,67) * y(idx_64,10) * y(idx_64,20) + rc(idx_64,27) * y(idx_64,28) * y(idx_64,5) + rc(idx_64,45) * y(idx_64,2) * y(idx_64,14) + rc(idx_64,25) * y(idx_64,25) * y(idx_64,28) + rc(idx_64,26) * y(idx_64,28) * y(idx_64,4) + rc(idx_64,10) * y(idx_64,11) * y(idx_64,28) + rc(idx_64,21) * y(idx_64,18) * y(idx_64,28)
        enddo
        !$acc loop gang vector independent
        do idx_65 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_65) = ddep(idx_65,18) + wdep(idx_65,18) + rc(idx_65,21) * y(idx_65,28) + rc(idx_65,22) * y(idx_65,28) + rc(idx_65,19) * y(idx_65,18) + rc(idx_65,20) * y(idx_65,18) + rc(idx_65,20) * y(idx_65,18) + rc(idx_65,17) * y(idx_65,4) + rc(idx_65,18) * y(idx_65,5) + rc(idx_65,19) * y(idx_65,18) + rc(idx_65,5) * y(idx_65,11) + rc(idx_65,6) * y(idx_65,11) + rc(idx_65,16) * y(idx_65,4)
        enddo
        !$acc loop gang vector independent
        do idx_66 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_66,18) = (yp(idx_66,18) + dts * p(idx_66)) / (1.0 + dts * l(idx_66))
        enddo
        !$acc loop gang vector independent
        do idx_67 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_67,19) = yp(idx_67,19)
        enddo
        !$acc loop gang vector independent
        do idx_68 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_68) = +rc(idx_68,5) * y(idx_68,11) * y(idx_68,18)
        enddo
        !$acc loop gang vector independent
        do idx_69 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_69) = ddep(idx_69,20) + wdep(idx_69,20) + rc(idx_69,66) * y(idx_69,10) + rc(idx_69,67) * y(idx_69,10) + dj(idx_69,1)
        enddo
        !$acc loop gang vector independent
        do idx_70 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_70,20) = (yp(idx_70,20) + dts * p(idx_70)) / (1.0 + dts * l(idx_70))
        enddo
        !$acc loop gang vector independent
        do idx_71 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_71) = +rc(idx_71,97) * y(idx_71,10) * y(idx_71,4)
        enddo
        !$acc loop gang vector independent
        do idx_72 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_72) = ddep(idx_72,21) + wdep(idx_72,21) + rc(idx_72,65) * y(idx_72,10) + dj(idx_72,17)
        enddo
        !$acc loop gang vector independent
        do idx_73 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_73,21) = (yp(idx_73,21) + dts * p(idx_73)) / (1.0 + dts * l(idx_73))
        enddo
        !$acc loop gang vector independent
        do idx_74 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_74,22) = yp(idx_74,22)
        enddo
        !$acc loop gang vector independent
        do idx_75 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_75,23) = yp(idx_75,23)
        enddo
        !$acc loop gang vector independent
        do idx_76 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_76) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_77 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_77) = ddep(idx_77,24) + wdep(idx_77,24) + rc(idx_77,52) * y(idx_77,10)
        enddo
        !$acc loop gang vector independent
        do idx_78 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_78,24) = (yp(idx_78,24) + dts * p(idx_78)) / (1.0 + dts * l(idx_78))
        enddo
        !$acc loop gang vector independent
        do idx_79 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_79) = +dj(idx_79,18) * y(idx_79,35) + rc(idx_79,52) * y(idx_79,10) * y(idx_79,24) + rc(idx_79,58) * y(idx_79,10) * y(idx_79,26) + rc(idx_79,32) * y(idx_79,36) * y(idx_79,4) + rc(idx_79,33) * y(idx_79,36) * y(idx_79,5)
        enddo
        !$acc loop gang vector independent
        do idx_80 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_80) = ddep(idx_80,25) + wdep(idx_80,25) + rc(idx_80,25) * y(idx_80,28) + rc(idx_80,7) * y(idx_80,11) + rc(idx_80,23) * y(idx_80,4) + rc(idx_80,24) * y(idx_80,5)
        enddo
        !$acc loop gang vector independent
        do idx_81 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_81,25) = (yp(idx_81,25) + dts * p(idx_81)) / (1.0 + dts * l(idx_81))
        enddo
        !$acc loop gang vector independent
        do idx_82 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_82) = +rc(idx_82,7) * y(idx_82,11) * y(idx_82,25)
        enddo
        !$acc loop gang vector independent
        do idx_83 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_83) = ddep(idx_83,26) + wdep(idx_83,26) + rc(idx_83,57) * y(idx_83,10) + rc(idx_83,58) * y(idx_83,10) + dj(idx_83,1)
        enddo
        !$acc loop gang vector independent
        do idx_84 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_84,26) = (yp(idx_84,26) + dts * p(idx_84)) / (1.0 + dts * l(idx_84))
        enddo
        !$acc loop gang vector independent
        do idx_85 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_85) = +rc(idx_85,78) * y(idx_85,10) * y(idx_85,40) + dj(idx_85,1) * y(idx_85,26) + rc(idx_85,25) * y(idx_85,25) * y(idx_85,28) + rc(idx_85,57) * y(idx_85,10) * y(idx_85,26) + rc(idx_85,23) * y(idx_85,25) * y(idx_85,4) + rc(idx_85,24) * y(idx_85,25) * y(idx_85,5)
        enddo
        !$acc loop gang vector independent
        do idx_86 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_86) = ddep(idx_86,27) + wdep(idx_86,27) + dj(idx_86,8) + rc(idx_86,40) * y(idx_86,5) + rc(idx_86,73) * y(idx_86,10) + dj(idx_86,7)
        enddo
        !$acc loop gang vector independent
        do idx_87 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_87,27) = (yp(idx_87,27) + dts * p(idx_87)) / (1.0 + dts * l(idx_87))
        enddo
        !$acc loop gang vector independent
        do idx_88 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_88) = +dj(idx_88,19) * y(idx_88,37) + dj(idx_88,1) * y(idx_88,39) + rc(idx_88,93) * y(idx_88,29) + dj(idx_88,16) * y(idx_88,29) + rc(idx_88,40) * y(idx_88,5) * y(idx_88,27) + rc(idx_88,73) * y(idx_88,10) * y(idx_88,27) + rc(idx_88,34) * y(idx_88,38) * y(idx_88,4) + rc(idx_88,35) * y(idx_88,38) * y(idx_88,5)
        enddo
        !$acc loop gang vector independent
        do idx_89 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_89) = ddep(idx_89,28) + wdep(idx_89,28) + rc(idx_89,26) * y(idx_89,4) + rc(idx_89,27) * y(idx_89,5) + rc(idx_89,92) * y(idx_89,6) + rc(idx_89,21) * y(idx_89,18) + rc(idx_89,22) * y(idx_89,18) + rc(idx_89,25) * y(idx_89,25) + rc(idx_89,8) * y(idx_89,11) + rc(idx_89,9) * y(idx_89,11) + rc(idx_89,10) * y(idx_89,11)
        enddo
        !$acc loop gang vector independent
        do idx_90 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_90,28) = (yp(idx_90,28) + dts * p(idx_90)) / (1.0 + dts * l(idx_90))
        enddo
        !$acc loop gang vector independent
        do idx_91 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_91) = +rc(idx_91,92) * y(idx_91,28) * y(idx_91,6)
        enddo
        !$acc loop gang vector independent
        do idx_92 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_92) = ddep(idx_92,29) + wdep(idx_92,29) + rc(idx_92,77) * y(idx_92,10) + rc(idx_92,93) + dj(idx_92,16)
        enddo
        !$acc loop gang vector independent
        do idx_93 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_93,29) = (yp(idx_93,29) + dts * p(idx_93)) / (1.0 + dts * l(idx_93))
        enddo
        !$acc loop gang vector independent
        do idx_94 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_94) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_95 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_95) = ddep(idx_95,30) + wdep(idx_95,30) + rc(idx_95,53) * y(idx_95,10) + rc(idx_95,54) * y(idx_95,10)
        enddo
        !$acc loop gang vector independent
        do idx_96 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_96,30) = (yp(idx_96,30) + dts * p(idx_96)) / (1.0 + dts * l(idx_96))
        enddo
        !$acc loop gang vector independent
        do idx_97 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_97) = +rc(idx_97,53) * y(idx_97,10) * y(idx_97,30) + rc(idx_97,79) * y(idx_97,10) * y(idx_97,33)
        enddo
        !$acc loop gang vector independent
        do idx_98 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_98) = ddep(idx_98,31) + wdep(idx_98,31) + rc(idx_98,11) * y(idx_98,11) + rc(idx_98,28) * y(idx_98,4) + rc(idx_98,29) * y(idx_98,5)
        enddo
        !$acc loop gang vector independent
        do idx_99 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_99,31) = (yp(idx_99,31) + dts * p(idx_99)) / (1.0 + dts * l(idx_99))
        enddo
        !$acc loop gang vector independent
        do idx_100 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_100) = +rc(idx_100,54) * y(idx_100,10) * y(idx_100,30) + rc(idx_100,81) * y(idx_100,10) * y(idx_100,34)
        enddo
        !$acc loop gang vector independent
        do idx_101 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_101) = ddep(idx_101,32) + wdep(idx_101,32) + rc(idx_101,12) * y(idx_101,11) + rc(idx_101,30) * y(idx_101,4) + rc(idx_101,31) * y(idx_101,5)
        enddo
        !$acc loop gang vector independent
        do idx_102 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_102,32) = (yp(idx_102,32) + dts * p(idx_102)) / (1.0 + dts * l(idx_102))
        enddo
        !$acc loop gang vector independent
        do idx_103 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_103) = +rc(idx_103,11) * y(idx_103,11) * y(idx_103,31)
        enddo
        !$acc loop gang vector independent
        do idx_104 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_104) = ddep(idx_104,33) + wdep(idx_104,33) + rc(idx_104,79) * y(idx_104,10) + rc(idx_104,80) * y(idx_104,10) + dj(idx_104,1)
        enddo
        !$acc loop gang vector independent
        do idx_105 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_105,33) = (yp(idx_105,33) + dts * p(idx_105)) / (1.0 + dts * l(idx_105))
        enddo
        !$acc loop gang vector independent
        do idx_106 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_106) = +rc(idx_106,12) * y(idx_106,11) * y(idx_106,32)
        enddo
        !$acc loop gang vector independent
        do idx_107 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_107) = ddep(idx_107,34) + wdep(idx_107,34) + rc(idx_107,81) * y(idx_107,10) + rc(idx_107,82) * y(idx_107,10) + dj(idx_107,1)
        enddo
        !$acc loop gang vector independent
        do idx_108 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_108,34) = (yp(idx_108,34) + dts * p(idx_108)) / (1.0 + dts * l(idx_108))
        enddo
        !$acc loop gang vector independent
        do idx_109 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_109) = +rc(idx_109,80) * y(idx_109,10) * y(idx_109,33) + dj(idx_109,1) * y(idx_109,33) + rc(idx_109,28) * y(idx_109,31) * y(idx_109,4) + rc(idx_109,29) * y(idx_109,31) * y(idx_109,5)
        enddo
        !$acc loop gang vector independent
        do idx_110 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_110) = ddep(idx_110,35) + wdep(idx_110,35) + rc(idx_110,41) * y(idx_110,5) + rc(idx_110,56) * y(idx_110,10) + dj(idx_110,18)
        enddo
        !$acc loop gang vector independent
        do idx_111 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_111,35) = (yp(idx_111,35) + dts * p(idx_111)) / (1.0 + dts * l(idx_111))
        enddo
        !$acc loop gang vector independent
        do idx_112 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_112) = +rc(idx_112,101) * y(idx_112,40) + dj(idx_112,16) * y(idx_112,40) + rc(idx_112,41) * y(idx_112,5) * y(idx_112,35) + rc(idx_112,56) * y(idx_112,10) * y(idx_112,35)
        enddo
        !$acc loop gang vector independent
        do idx_113 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_113) = ddep(idx_113,36) + wdep(idx_113,36) + rc(idx_113,33) * y(idx_113,5) + rc(idx_113,100) * y(idx_113,6) + rc(idx_113,13) * y(idx_113,11) + rc(idx_113,14) * y(idx_113,11) + rc(idx_113,32) * y(idx_113,4)
        enddo
        !$acc loop gang vector independent
        do idx_114 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_114,36) = (yp(idx_114,36) + dts * p(idx_114)) / (1.0 + dts * l(idx_114))
        enddo
        !$acc loop gang vector independent
        do idx_115 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_115) = +rc(idx_115,82) * y(idx_115,10) * y(idx_115,34) + dj(idx_115,1) * y(idx_115,34) + rc(idx_115,30) * y(idx_115,32) * y(idx_115,4) + rc(idx_115,31) * y(idx_115,32) * y(idx_115,5)
        enddo
        !$acc loop gang vector independent
        do idx_116 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_116) = ddep(idx_116,37) + wdep(idx_116,37) + dj(idx_116,19) + rc(idx_116,42) * y(idx_116,5) + rc(idx_116,69) * y(idx_116,10) + rc(idx_116,70) * y(idx_116,10)
        enddo
        !$acc loop gang vector independent
        do idx_117 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_117,37) = (yp(idx_117,37) + dts * p(idx_117)) / (1.0 + dts * l(idx_117))
        enddo
        !$acc loop gang vector independent
        do idx_118 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_118) = +rc(idx_118,70) * y(idx_118,10) * y(idx_118,37) + rc(idx_118,71) * y(idx_118,10) * y(idx_118,39) + rc(idx_118,42) * y(idx_118,5) * y(idx_118,37) + rc(idx_118,69) * y(idx_118,10) * y(idx_118,37)
        enddo
        !$acc loop gang vector independent
        do idx_119 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_119) = ddep(idx_119,38) + wdep(idx_119,38) + rc(idx_119,15) * y(idx_119,11) + rc(idx_119,34) * y(idx_119,4) + rc(idx_119,35) * y(idx_119,5)
        enddo
        !$acc loop gang vector independent
        do idx_120 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_120,38) = (yp(idx_120,38) + dts * p(idx_120)) / (1.0 + dts * l(idx_120))
        enddo
        !$acc loop gang vector independent
        do idx_121 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_121) = +rc(idx_121,15) * y(idx_121,11) * y(idx_121,38)
        enddo
        !$acc loop gang vector independent
        do idx_122 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_122) = ddep(idx_122,39) + wdep(idx_122,39) + rc(idx_122,71) * y(idx_122,10) + rc(idx_122,72) * y(idx_122,10) + dj(idx_122,1)
        enddo
        !$acc loop gang vector independent
        do idx_123 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_123,39) = (yp(idx_123,39) + dts * p(idx_123)) / (1.0 + dts * l(idx_123))
        enddo
        !$acc loop gang vector independent
        do idx_124 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_124) = +rc(idx_124,100) * y(idx_124,36) * y(idx_124,6)
        enddo
        !$acc loop gang vector independent
        do idx_125 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_125) = ddep(idx_125,40) + wdep(idx_125,40) + rc(idx_125,78) * y(idx_125,10) + rc(idx_125,101) + dj(idx_125,16)
        enddo
        !$acc loop gang vector independent
        do idx_126 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_126,40) = (yp(idx_126,40) + dts * p(idx_126)) / (1.0 + dts * l(idx_126))
        enddo
        !$acc loop gang vector independent
        do idx_127 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_127) = +rc(idx_127,17) * y(idx_127,18) * y(idx_127,4)
        enddo
        !$acc loop gang vector independent
        do idx_128 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_128) = ddep(idx_128,41) + wdep(idx_128,41) + rc(idx_128,68) * y(idx_128,10) + dj(idx_128,20)
        enddo
        !$acc loop gang vector independent
        do idx_129 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_129,41) = (yp(idx_129,41) + dts * p(idx_129)) / (1.0 + dts * l(idx_129))
        enddo
        !$acc loop gang vector independent
        do idx_130 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_130) = +dj(idx_130,15) * y(idx_130,44) + rc(idx_130,87) * y(idx_130,43) * y(idx_130,23) + rc(idx_130,88) * y(idx_130,43) * y(idx_130,22)
        enddo
        !$acc loop gang vector independent
        do idx_131 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_131) = ddep(idx_131,42) + wdep(idx_131,42) + rc(idx_131,102) * y(idx_131,22)
        enddo
        !$acc loop gang vector independent
        do idx_132 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_132,42) = (yp(idx_132,42) + dts * p(idx_132)) / (1.0 + dts * l(idx_132))
        enddo
        !$acc loop gang vector independent
        do idx_133 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_133) = +dj(idx_133,14) * y(idx_133,44)
        enddo
        !$acc loop gang vector independent
        do idx_134 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_134) = ddep(idx_134,43) + wdep(idx_134,43) + rc(idx_134,86) * y(idx_134,19) + rc(idx_134,87) * y(idx_134,23) + rc(idx_134,88) * y(idx_134,22)
        enddo
        !$acc loop gang vector independent
        do idx_135 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_135,43) = (yp(idx_135,43) + dts * p(idx_135)) / (1.0 + dts * l(idx_135))
        enddo
        !$acc loop gang vector independent
        do idx_136 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_136) = +rc(idx_136,102) * y(idx_136,42) * y(idx_136,22)
        enddo
        !$acc loop gang vector independent
        do idx_137 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_137) = ddep(idx_137,44) + wdep(idx_137,44) + dj(idx_137,15) + rc(idx_137,84) * y(idx_137,46) + rc(idx_137,85) * y(idx_137,45) + dj(idx_137,14)
        enddo
        !$acc loop gang vector independent
        do idx_138 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_138,44) = (yp(idx_138,44) + dts * p(idx_138)) / (1.0 + dts * l(idx_138))
        enddo
        !$acc loop gang vector independent
        do idx_139 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_139) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_140 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_140) = ddep(idx_140,45) + wdep(idx_140,45) + rc(idx_140,85) * y(idx_140,44)
        enddo
        !$acc loop gang vector independent
        do idx_141 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_141,45) = (yp(idx_141,45) + dts * p(idx_141)) / (1.0 + dts * l(idx_141))
        enddo
        !$acc loop gang vector independent
        do idx_142 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_142) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_143 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_143) = ddep(idx_143,46) + wdep(idx_143,46) + rc(idx_143,84) * y(idx_143,44)
        enddo
        !$acc loop gang vector independent
        do idx_144 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_144,46) = (yp(idx_144,46) + dts * p(idx_144)) / (1.0 + dts * l(idx_144))
        enddo
        !$acc loop gang vector independent
        do idx_145 = LBOUND(f_so2, dim=1), UBOUND(f_so2, dim=1), 1
          f_so2(idx_145) = k_dms(idx_145,1) / (k_dms(idx_145,1) + k_dms(idx_145,2) * y(idx_145,3) + k_dms(idx_145,3) * y(idx_145,6))
        enddo
        !$acc loop gang vector independent
        do idx_146 = LBOUND(f_so4, dim=1), UBOUND(f_so4, dim=1), 1
          f_so4(idx_146) = (1.0 - f_so2(idx_146)) * k_dms(idx_146,5) / (k_dms(idx_146,5) + k_dms(idx_146,4) * y(idx_146,11))
        enddo
        !$acc loop gang vector independent
        do idx_147 = LBOUND(f_msa, dim=1), UBOUND(f_msa, dim=1), 1
          f_msa(idx_147) = (1.0 - f_so2(idx_147)) * k_dms(idx_147,4) * y(idx_147,11) / (k_dms(idx_147,5) + k_dms(idx_147,4) * y(idx_147,11))
        enddo
        !$acc loop gang vector independent
        do idx_148 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_148) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_149 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_149) = ddep(idx_149,47) + wdep(idx_149,47) + rc(idx_149,105) * y(idx_149,10) + rc(idx_149,106) * y(idx_149,10) + rc(idx_149,107) * y(idx_149,5)
        enddo
        !$acc loop gang vector independent
        do idx_150 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_150,47) = (yp(idx_150,47) + dts * p(idx_150)) / (1.0 + dts * l(idx_150))
        enddo
        !$acc loop gang vector independent
        do idx_151 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_151) = 0.0 + (rc(idx_151,105) * y(idx_151,10) + rc(idx_151,106) * y(idx_151,10) + rc(idx_151,107) * y(idx_151,5)) * y(idx_151,47) * f_so2(idx_151)
        enddo
        !$acc loop gang vector independent
        do idx_152 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_152) = 0.0 + ddep(idx_152,48) + wdep(idx_152,48) + rc(idx_152,103) * y(idx_152,10) + rc(idx_152,104) * y(idx_152,13) + rc(idx_152,115) * y(idx_152,3)
        enddo
        !$acc loop gang vector independent
        do idx_153 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_153,48) = (yp(idx_153,48) + dts * p(idx_153)) / (1.0 + dts * l(idx_153))
        enddo
        !$acc loop gang vector independent
        do idx_154 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_154) = 0.0 + rc(idx_154,103) * y(idx_154,48) * y(idx_154,10) + (rc(idx_154,105) * y(idx_154,10) + rc(idx_154,106) * y(idx_154,10) + rc(idx_154,107) * y(idx_154,5)) * y(idx_154,47) * f_so4(idx_154)
        enddo
        !$acc loop gang vector independent
        do idx_155 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_155) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_156 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_156,49) = (yp(idx_156,49) + dts * p(idx_156)) / (1.0 + dts * l(idx_156))
        enddo
        !$acc loop gang vector independent
        do idx_157 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_157) = 0.0 + (rc(idx_157,105) * y(idx_157,10) + rc(idx_157,106) * y(idx_157,10) + rc(idx_157,107) * y(idx_157,5)) * y(idx_157,47) * f_msa(idx_157)
        enddo
        !$acc loop gang vector independent
        do idx_158 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_158) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_159 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_159,50) = (yp(idx_159,50) + dts * p(idx_159)) / (1.0 + dts * l(idx_159))
        enddo
        !$acc loop gang vector independent
        do idx_160 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_160) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_161 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_161) = ddep(idx_161,51) + wdep(idx_161,51) + rc(idx_161,111) * y(idx_161,10)
        enddo
        !$acc loop gang vector independent
        do idx_162 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_162,51) = (yp(idx_162,51) + dts * p(idx_162)) / (1.0 + dts * l(idx_162))
        enddo
        !$acc loop gang vector independent
        do idx_163 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_163) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_164 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_164) = ddep(idx_164,52) + rc(idx_164,112) * y(idx_164,10) + rc(idx_164,113) * y(idx_164,3) + rc(idx_164,114) * y(idx_164,5)
        enddo
        !$acc loop gang vector independent
        do idx_165 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_165,52) = (yp(idx_165,52) + dts * p(idx_165)) / (1.0 + dts * l(idx_165))
        enddo
        !$acc loop gang vector independent
        do idx_166 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_166) = 0.0 + (rc(idx_166,112) * y(idx_166,52) * y(idx_166,10) * 0.13 + rc(idx_166,113) * y(idx_166,52) * y(idx_166,3) * 0.13) + rc(idx_166,114) * y(idx_166,52) * y(idx_166,5) * 0.13
        enddo
        !$acc loop gang vector independent
        do idx_167 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_167) = ddep(idx_167,53) + wdep(idx_167,53)
        enddo
        !$acc loop gang vector independent
        do idx_168 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_168,53) = (yp(idx_168,53) + dts * p(idx_168)) / (1.0 + dts * l(idx_168))
        enddo
      enddo
      !$acc loop gang vector independent
      do idx_169 = LBOUND(delta_so2_wetox_h2o2, dim=1), UBOUND(delta_so2_wetox_h2o2, dim=1), 1
        delta_so2_wetox_h2o2(idx_169) = delta_so2_wetox_h2o2(idx_169) + rc(idx_169,104) * y(idx_169,13) * y(idx_169,48) * dts
      enddo
      !$acc loop gang vector independent
      do idx_170 = LBOUND(delta_so2_wetox_o3, dim=1), UBOUND(delta_so2_wetox_o3, dim=1), 1
        delta_so2_wetox_o3(idx_170) = delta_so2_wetox_o3(idx_170) + rc(idx_170,115) * y(idx_170,3) * y(idx_170,48) * dts
      enddo
      !$acc loop gang vector independent
      do idx_171 = LBOUND(delta_so2_dryox_oh, dim=1), UBOUND(delta_so2_dryox_oh, dim=1), 1
        delta_so2_dryox_oh(idx_171) = delta_so2_dryox_oh(idx_171) + rc(idx_171,103) * y(idx_171,10) * y(idx_171,48) * dts
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_deriv_aero

end module ukca_deriv_aero_mod
