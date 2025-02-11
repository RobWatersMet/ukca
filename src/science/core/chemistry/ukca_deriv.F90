module ukca_deriv_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DERIV_MOD'
  public

  contains
  subroutine ukca_deriv(nr, n_be_calls, n_pnts, rc, wdep, ddep, dj, h2o, m, o2, dts, y)
    use asad_mod, only : jppj, jpspec
    use ukca_config_specification_mod, only : ukca_config
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr
    integer, intent(in) :: n_be_calls
    integer, intent(in) :: n_pnts
    real, intent(in) :: dts
    real, dimension(n_pnts,nr), intent(in) :: rc
    real, dimension(n_pnts,jppj), intent(in) :: dj
    real, dimension(n_pnts,jpspec), intent(in) :: wdep
    real, dimension(n_pnts,jpspec), intent(in) :: ddep
    real, dimension(n_pnts), intent(in) :: h2o
    real, dimension(n_pnts), intent(in) :: m
    real, dimension(n_pnts), intent(in) :: o2
    real, dimension(n_pnts,jpspec), intent(inout) :: y
    integer :: i
    integer :: j
    integer :: n
    real, dimension(n_pnts) :: p
    real, dimension(n_pnts) :: p1
    real, dimension(n_pnts) :: p2
    real, dimension(n_pnts) :: r1
    real, dimension(n_pnts) :: r2
    real, dimension(n_pnts) :: l
    real, dimension(n_pnts) :: l1
    real, dimension(n_pnts) :: l2
    real, dimension(n_pnts) :: l3
    real, dimension(n_pnts,jpspec) :: yp
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DERIV'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    !$acc loop seq
    do n = 1, n_be_calls, 1
      !$acc loop gang vector independent collapse(2)
      do j = 1, jpspec, 1
        do idx = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
          yp(idx,j) = y(idx,j)
        enddo
      enddo
      !$acc loop gang vector independent
      do idx_1 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_1,12) = 5.0e-7 * m(idx_1)
      enddo
      !$acc loop gang vector independent
      do idx_2 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_2,16) = 350.0e-6 * m(idx_2)
      enddo
      !$acc loop gang vector independent
      do idx_3 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_3,19) = h2o(idx_3)
      enddo
      !$acc loop gang vector independent
      do idx_4 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_4,22) = o2(idx_4)
      enddo
      !$acc loop gang vector independent
      do idx_5 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
        yp(idx_5,23) = m(idx_5) - o2(idx_5)
      enddo
      !$acc loop gang vector independent
      do idx_6 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_6,12) = yp(idx_6,12)
      enddo
      !$acc loop gang vector independent
      do idx_7 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_7,16) = yp(idx_7,16)
      enddo
      !$acc loop gang vector independent
      do idx_8 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_8,19) = yp(idx_8,19)
      enddo
      !$acc loop gang vector independent
      do idx_9 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_9,22) = yp(idx_9,22)
      enddo
      !$acc loop gang vector independent
      do idx_10 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_10,23) = yp(idx_10,23)
      enddo
      !$acc loop seq
      do i = 1, ukca_config%nit, 1
        !$acc loop gang vector independent
        do idx_11 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_11) = 0.0 + dj(idx_11,15) * y(idx_11,3) + dj(idx_11,12) * y(idx_11,5) + dj(idx_11,13) * y(idx_11,22) * 2.00 + rc(idx_11,76) * y(idx_11,10) * y(idx_11,10) + dj(idx_11,10) * y(idx_11,6) + rc(idx_11,49) * y(idx_11,2) * y(idx_11,23) + rc(idx_11,50) * y(idx_11,2) * y(idx_11,22)
        enddo
        !$acc loop gang vector independent
        do idx_12 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_12) = ddep(idx_12,1) + wdep(idx_12,1) + rc(idx_12,44) * y(idx_12,3) + rc(idx_12,83) * y(idx_12,6) + rc(idx_12,96) * y(idx_12,22)
        enddo
        !$acc loop gang vector independent
        do idx_13 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_13,1) = p(idx_13) / l(idx_13)
        enddo
        !$acc loop gang vector independent
        do idx_14 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_14) = 0.0 + dj(idx_14,14) * y(idx_14,3)
        enddo
        !$acc loop gang vector independent
        do idx_15 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_15) = ddep(idx_15,2) + wdep(idx_15,2) + rc(idx_15,48) * y(idx_15,19) + rc(idx_15,49) * y(idx_15,23) + rc(idx_15,50) * y(idx_15,22) + rc(idx_15,45) * y(idx_15,14) + rc(idx_15,46) * y(idx_15,14) + rc(idx_15,47) * y(idx_15,14)
        enddo
        !$acc loop gang vector independent
        do idx_16 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_16,2) = p(idx_16) / l(idx_16)
        enddo
        !$acc loop gang vector independent
        do idx_17 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_17) = 0.0 + dj(idx_17,1) * y(idx_17,34) + dj(idx_17,1) * y(idx_17,39) + dj(idx_17,17) * y(idx_17,21) + dj(idx_17,1) * y(idx_17,33) + dj(idx_17,6) * y(idx_17,9) + dj(idx_17,1) * y(idx_17,20) + dj(idx_17,1) * y(idx_17,26) + dj(idx_17,2) * y(idx_17,13) * 2.00 + rc(idx_17,82) * y(idx_17,10) * y(idx_17,34) + rc(idx_17,85) * y(idx_17,45) * y(idx_17,44) + rc(idx_17,72) * y(idx_17,10) * y(idx_17,39) + rc(idx_17,80) * y(idx_17,10) * y(idx_17,33) + rc(idx_17,57) * y(idx_17,10) * y(idx_17,26) + rc(idx_17,66) * y(idx_17,10) * y(idx_17,20) + rc(idx_17,45) * y(idx_17,2) * y(idx_17,14) + rc(idx_17,48) * y(idx_17,2) * y(idx_17,19) * 2.00 + rc(idx_17,3) * y(idx_17,11) * y(idx_17,3) + rc(idx_17,10) * y(idx_17,11) * y(idx_17,28) + rc(idx_17,1) * y(idx_17,11) * y(idx_17,4) + rc(idx_17,2) * y(idx_17,11) * y(idx_17,5)
        enddo
        !$acc loop gang vector independent
        do idx_18 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_18) = ddep(idx_18,10) + wdep(idx_18,10) + rc(idx_18,99) * y(idx_18,10) + rc(idx_18,97) * y(idx_18,4) + rc(idx_18,98) * y(idx_18,6) + rc(idx_18,99) * y(idx_18,10) + rc(idx_18,80) * y(idx_18,33) + rc(idx_18,81) * y(idx_18,34) + rc(idx_18,82) * y(idx_18,34) + rc(idx_18,77) * y(idx_18,29) + rc(idx_18,78) * y(idx_18,40) + rc(idx_18,79) * y(idx_18,33) + rc(idx_18,75) * y(idx_18,3) + rc(idx_18,76) * y(idx_18,10) + rc(idx_18,76) * y(idx_18,10) + rc(idx_18,72) * y(idx_18,39) + rc(idx_18,73) * y(idx_18,27) + rc(idx_18,74) * y(idx_18,5) + rc(idx_18,69) * y(idx_18,37) + rc(idx_18,70) * y(idx_18,37) + rc(idx_18,71) * y(idx_18,39) + rc(idx_18,66) * y(idx_18,20) + rc(idx_18,67) * y(idx_18,20) + rc(idx_18,68) * y(idx_18,41) + rc(idx_18,63) * y(idx_18,8) + rc(idx_18,64) * y(idx_18,9) + rc(idx_18,65) * y(idx_18,21) + rc(idx_18,60) * y(idx_18,13) + rc(idx_18,61) * y(idx_18,17) + rc(idx_18,62) * y(idx_18,11) + rc(idx_18,57) * y(idx_18,26) + rc(idx_18,58) * y(idx_18,26) + rc(idx_18,59) * y(idx_18,12) + rc(idx_18,54) * y(idx_18,30) + rc(idx_18,55) * y(idx_18,15) + rc(idx_18,56) * y(idx_18,35) + rc(idx_18,51) * y(idx_18,14) + rc(idx_18,52) * y(idx_18,24) + rc(idx_18,53) * y(idx_18,30)
        enddo
        !$acc loop gang vector independent
        do idx_19 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_19,10) = p(idx_19) / l(idx_19)
        enddo
        !$acc loop gang vector independent
        do idx_20 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_20) = 0.0 + rc(idx_20,96) * y(idx_20,1) * y(idx_20,22) + rc(idx_20,9) * y(idx_20,11) * y(idx_20,28) + rc(idx_20,14) * y(idx_20,11) * y(idx_20,36)
        enddo
        !$acc loop gang vector independent
        do idx_21 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_21) = ddep(idx_21,3) + wdep(idx_21,3) + dj(idx_21,15) + rc(idx_21,44) * y(idx_21,1) + rc(idx_21,75) * y(idx_21,10) + dj(idx_21,14) + rc(idx_21,3) * y(idx_21,11) + rc(idx_21,37) * y(idx_21,4) + rc(idx_21,38) * y(idx_21,6)
        enddo
        !$acc loop gang vector independent
        do idx_22 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_22,3) = (yp(idx_22,3) + dts * p(idx_22)) / (1.0 + dts * l(idx_22))
        enddo
        !$acc loop gang vector independent
        do idx_23 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_23) = 0.0 + dj(idx_23,11) * y(idx_23,5) + dj(idx_23,17) * y(idx_23,21) + rc(idx_23,83) * y(idx_23,1) * y(idx_23,6) + dj(idx_23,10) * y(idx_23,6)
        enddo
        !$acc loop gang vector independent
        do idx_24 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_24) = ddep(idx_24,4) + wdep(idx_24,4) + rc(idx_24,36) * y(idx_24,5) + rc(idx_24,37) * y(idx_24,3) + rc(idx_24,97) * y(idx_24,10) + rc(idx_24,30) * y(idx_24,32) + rc(idx_24,32) * y(idx_24,36) + rc(idx_24,34) * y(idx_24,38) + rc(idx_24,23) * y(idx_24,25) + rc(idx_24,26) * y(idx_24,28) + rc(idx_24,28) * y(idx_24,31) + rc(idx_24,1) * y(idx_24,11) + rc(idx_24,16) * y(idx_24,18) + rc(idx_24,17) * y(idx_24,18)
        enddo
        !$acc loop gang vector independent
        do idx_25 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_25,4) = (yp(idx_25,4) + dts * p(idx_25)) / (1.0 + dts * l(idx_25))
        enddo
        !$acc loop gang vector independent
        do idx_26 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_26) = 0.0 + dj(idx_26,20) * y(idx_26,41) + dj(idx_26,16) * y(idx_26,29) + dj(idx_26,16) * y(idx_26,40) + dj(idx_26,9) * y(idx_26,7) + dj(idx_26,12) * y(idx_26,5) + dj(idx_26,5) * y(idx_26,8) + dj(idx_26,6) * y(idx_26,9) + rc(idx_26,94) * y(idx_26,7) + rc(idx_26,101) * y(idx_26,40) + rc(idx_26,91) * y(idx_26,8) + rc(idx_26,93) * y(idx_26,29) + rc(idx_26,77) * y(idx_26,10) * y(idx_26,29) + rc(idx_26,78) * y(idx_26,10) * y(idx_26,40) + rc(idx_26,68) * y(idx_26,10) * y(idx_26,41) + rc(idx_26,74) * y(idx_26,10) * y(idx_26,5) + rc(idx_26,63) * y(idx_26,10) * y(idx_26,8) + rc(idx_26,65) * y(idx_26,10) * y(idx_26,21) + rc(idx_26,36) * y(idx_26,4) * y(idx_26,5) * 2.00 + rc(idx_26,37) * y(idx_26,4) * y(idx_26,3) + rc(idx_26,34) * y(idx_26,38) * y(idx_26,4) + rc(idx_26,35) * y(idx_26,38) * y(idx_26,5) + rc(idx_26,32) * y(idx_26,36) * y(idx_26,4) + rc(idx_26,33) * y(idx_26,36) * y(idx_26,5) + rc(idx_26,30) * y(idx_26,32) * y(idx_26,4) + rc(idx_26,31) * y(idx_26,32) * y(idx_26,5) + rc(idx_26,28) * y(idx_26,31) * y(idx_26,4) + rc(idx_26,29) * y(idx_26,31) * y(idx_26,5) + rc(idx_26,26) * y(idx_26,28) * y(idx_26,4) + rc(idx_26,27) * y(idx_26,28) * y(idx_26,5) + rc(idx_26,23) * y(idx_26,25) * y(idx_26,4) + rc(idx_26,24) * y(idx_26,25) * y(idx_26,5) + rc(idx_26,16) * y(idx_26,18) * y(idx_26,4) + rc(idx_26,18) * y(idx_26,18) * y(idx_26,5) + rc(idx_26,1) * y(idx_26,11) * y(idx_26,4) + rc(idx_26,2) * y(idx_26,11) * y(idx_26,5)
        enddo
        !$acc loop gang vector independent
        do idx_27 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_27) = ddep(idx_27,6) + wdep(idx_27,6) + rc(idx_27,100) * y(idx_27,36) + dj(idx_27,10) + rc(idx_27,92) * y(idx_27,28) + rc(idx_27,95) * y(idx_27,5) + rc(idx_27,98) * y(idx_27,10) + rc(idx_27,38) * y(idx_27,3) + rc(idx_27,83) * y(idx_27,1) + rc(idx_27,90) * y(idx_27,11)
        enddo
        !$acc loop gang vector independent
        do idx_28 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_28,6) = (yp(idx_28,6) + dts * p(idx_28)) / (1.0 + dts * l(idx_28))
        enddo
        !$acc loop gang vector independent
        do idx_29 = LBOUND(p1, dim=1), UBOUND(p1, dim=1), 1
          p1(idx_29) = 0.0 + rc(idx_29,38) * y(idx_29,6) * y(idx_29,3) + rc(idx_29,64) * y(idx_29,10) * y(idx_29,9)
        enddo
        !$acc loop gang vector independent
        do idx_30 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_30) = ddep(idx_30,5) + wdep(idx_30,5) + dj(idx_30,11) + dj(idx_30,12) + rc(idx_30,42) * y(idx_30,37) + rc(idx_30,74) * y(idx_30,10) + rc(idx_30,95) * y(idx_30,6) + rc(idx_30,39) * y(idx_30,17) + rc(idx_30,40) * y(idx_30,27) + rc(idx_30,41) * y(idx_30,35) + rc(idx_30,33) * y(idx_30,36) + rc(idx_30,35) * y(idx_30,38) + rc(idx_30,36) * y(idx_30,4) + rc(idx_30,27) * y(idx_30,28) + rc(idx_30,29) * y(idx_30,31) + rc(idx_30,31) * y(idx_30,32) + rc(idx_30,2) * y(idx_30,11) + rc(idx_30,18) * y(idx_30,18) + rc(idx_30,24) * y(idx_30,25)
        enddo
        !$acc loop gang vector independent
        do idx_31 = LBOUND(p2, dim=1), UBOUND(p2, dim=1), 1
          p2(idx_31) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_32 = LBOUND(r1, dim=1), UBOUND(r1, dim=1), 1
          r1(idx_32) = rc(idx_32,94) + dj(idx_32,9)
        enddo
        !$acc loop gang vector independent
        do idx_33 = LBOUND(r2, dim=1), UBOUND(r2, dim=1), 1
          r2(idx_33) = rc(idx_33,95) * y(idx_33,6)
        enddo
        !$acc loop gang vector independent
        do idx_34 = LBOUND(l1, dim=1), UBOUND(l1, dim=1), 1
          l1(idx_34) = ddep(idx_34,7) + wdep(idx_34,7) + rc(idx_34,43) * y(idx_34,19) + rc(idx_34,94) + dj(idx_34,9)
        enddo
        !$acc loop gang vector independent
        do idx_35 = LBOUND(l2, dim=1), UBOUND(l2, dim=1), 1
          l2(idx_35) = 1.0 + l(idx_35) * dts
        enddo
        !$acc loop gang vector independent
        do idx_36 = LBOUND(l3, dim=1), UBOUND(l3, dim=1), 1
          l3(idx_36) = 1.0 + l1(idx_36) * dts
        enddo
        !$acc loop gang vector independent
        do idx_37 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_37,5) = (l3(idx_37) * (yp(idx_37,5) + p1(idx_37) * dts) + r1(idx_37) * dts * yp(idx_37,7)) / (l3(idx_37) * l2(idx_37) - r1(idx_37) * r2(idx_37) * dts ** 2)
        enddo
        !$acc loop gang vector independent
        do idx_38 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_38,7) = (yp(idx_38,7) + p2(idx_38) * dts + r2(idx_38) * dts * y(idx_38,5)) / l3(idx_38)
        enddo
        !$acc loop gang vector independent
        do idx_39 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_39) = 0.0 + rc(idx_39,90) * y(idx_39,11) * y(idx_39,6)
        enddo
        !$acc loop gang vector independent
        do idx_40 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_40) = ddep(idx_40,8) + wdep(idx_40,8) + rc(idx_40,63) * y(idx_40,10) + rc(idx_40,91) + dj(idx_40,5)
        enddo
        !$acc loop gang vector independent
        do idx_41 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_41,8) = (yp(idx_41,8) + dts * p(idx_41)) / (1.0 + dts * l(idx_41))
        enddo
        !$acc loop gang vector independent
        do idx_42 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_42) = 0.0 + rc(idx_42,43) * y(idx_42,7) * y(idx_42,19) * 2.00 + rc(idx_42,98) * y(idx_42,10) * y(idx_42,6) + rc(idx_42,41) * y(idx_42,5) * y(idx_42,35) + rc(idx_42,42) * y(idx_42,5) * y(idx_42,37) + rc(idx_42,39) * y(idx_42,5) * y(idx_42,17) + rc(idx_42,40) * y(idx_42,5) * y(idx_42,27)
        enddo
        !$acc loop gang vector independent
        do idx_43 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_43) = ddep(idx_43,9) + wdep(idx_43,9) + rc(idx_43,64) * y(idx_43,10) + dj(idx_43,6)
        enddo
        !$acc loop gang vector independent
        do idx_44 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_44,9) = (yp(idx_44,9) + dts * p(idx_44)) / (1.0 + dts * l(idx_44))
        enddo
        !$acc loop gang vector independent
        do idx_45 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_45) = 0.0 + dj(idx_45,1) * y(idx_45,34) + dj(idx_45,20) * y(idx_45,41) + dj(idx_45,18) * y(idx_45,35) + dj(idx_45,1) * y(idx_45,33) + dj(idx_45,7) * y(idx_45,27) + dj(idx_45,1) * y(idx_45,20) + dj(idx_45,3) * y(idx_45,17) * 2.00 + dj(idx_45,5) * y(idx_45,8) + rc(idx_45,91) * y(idx_45,8) + dj(idx_45,1) * y(idx_45,26) + rc(idx_45,75) * y(idx_45,10) * y(idx_45,3) + rc(idx_45,84) * y(idx_45,46) * y(idx_45,44) + rc(idx_45,61) * y(idx_45,10) * y(idx_45,17) + rc(idx_45,74) * y(idx_45,10) * y(idx_45,5) + rc(idx_45,59) * y(idx_45,10) * y(idx_45,12) + rc(idx_45,60) * y(idx_45,10) * y(idx_45,13) + rc(idx_45,47) * y(idx_45,2) * y(idx_45,14) * 2.00 + rc(idx_45,55) * y(idx_45,10) * y(idx_45,15) + rc(idx_45,31) * y(idx_45,32) * y(idx_45,5) + rc(idx_45,39) * y(idx_45,5) * y(idx_45,17) + rc(idx_45,29) * y(idx_45,31) * y(idx_45,5) + rc(idx_45,30) * y(idx_45,32) * y(idx_45,4) + rc(idx_45,25) * y(idx_45,25) * y(idx_45,28) + rc(idx_45,28) * y(idx_45,31) * y(idx_45,4) + rc(idx_45,23) * y(idx_45,25) * y(idx_45,4) + rc(idx_45,24) * y(idx_45,25) * y(idx_45,5) + rc(idx_45,20) * y(idx_45,18) * y(idx_45,18) * 2.00 + rc(idx_45,21) * y(idx_45,18) * y(idx_45,28) + rc(idx_45,16) * y(idx_45,18) * y(idx_45,4) + rc(idx_45,18) * y(idx_45,18) * y(idx_45,5)
        enddo
        !$acc loop gang vector independent
        do idx_46 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_46) = ddep(idx_46,11) + wdep(idx_46,11) + rc(idx_46,89) * y(idx_46,11) + rc(idx_46,90) * y(idx_46,6) + rc(idx_46,15) * y(idx_46,38) + rc(idx_46,62) * y(idx_46,10) + rc(idx_46,89) * y(idx_46,11) + rc(idx_46,12) * y(idx_46,32) + rc(idx_46,13) * y(idx_46,36) + rc(idx_46,14) * y(idx_46,36) + rc(idx_46,9) * y(idx_46,28) + rc(idx_46,10) * y(idx_46,28) + rc(idx_46,11) * y(idx_46,31) + rc(idx_46,6) * y(idx_46,18) + rc(idx_46,7) * y(idx_46,25) + rc(idx_46,8) * y(idx_46,28) + rc(idx_46,4) * y(idx_46,11) + rc(idx_46,4) * y(idx_46,11) + rc(idx_46,5) * y(idx_46,18) + rc(idx_46,1) * y(idx_46,4) + rc(idx_46,2) * y(idx_46,5) + rc(idx_46,3) * y(idx_46,3)
        enddo
        !$acc loop gang vector independent
        do idx_47 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_47,11) = (yp(idx_47,11) + dts * p(idx_47)) / (1.0 + dts * l(idx_47))
        enddo
        !$acc loop gang vector independent
        do idx_48 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_48) = 0.0 + rc(idx_48,99) * y(idx_48,10) * y(idx_48,10) + rc(idx_48,4) * y(idx_48,11) * y(idx_48,11) + rc(idx_48,89) * y(idx_48,11) * y(idx_48,11)
        enddo
        !$acc loop gang vector independent
        do idx_49 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_49) = ddep(idx_49,13) + wdep(idx_49,13) + rc(idx_49,60) * y(idx_49,10) + dj(idx_49,2)
        enddo
        !$acc loop gang vector independent
        do idx_50 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_50,13) = (yp(idx_50,13) + dts * p(idx_50)) / (1.0 + dts * l(idx_50))
        enddo
        !$acc loop gang vector independent
        do idx_51 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_51) = 0.0 + dj(idx_51,8) * y(idx_51,27)
        enddo
        !$acc loop gang vector independent
        do idx_52 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_52) = ddep(idx_52,14) + wdep(idx_52,14) + rc(idx_52,51) * y(idx_52,10) + rc(idx_52,45) * y(idx_52,2) + rc(idx_52,46) * y(idx_52,2) + rc(idx_52,47) * y(idx_52,2)
        enddo
        !$acc loop gang vector independent
        do idx_53 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_53,14) = (yp(idx_53,14) + dts * p(idx_53)) / (1.0 + dts * l(idx_53))
        enddo
        !$acc loop gang vector independent
        do idx_54 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_54) = 0.0 + dj(idx_54,18) * y(idx_54,35) + dj(idx_54,7) * y(idx_54,27) + dj(idx_54,8) * y(idx_54,27) + dj(idx_54,3) * y(idx_54,17) + dj(idx_54,4) * y(idx_54,17) + rc(idx_54,39) * y(idx_54,5) * y(idx_54,17) + rc(idx_54,61) * y(idx_54,10) * y(idx_54,17)
        enddo
        !$acc loop gang vector independent
        do idx_55 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_55) = ddep(idx_55,15) + wdep(idx_55,15) + rc(idx_55,55) * y(idx_55,10)
        enddo
        !$acc loop gang vector independent
        do idx_56 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_56,15) = (yp(idx_56,15) + dts * p(idx_56)) / (1.0 + dts * l(idx_56))
        enddo
        !$acc loop gang vector independent
        do idx_57 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_57) = 0.0 + dj(idx_57,20) * y(idx_57,41) + dj(idx_57,1) * y(idx_57,20) + dj(idx_57,1) * y(idx_57,39) + rc(idx_57,68) * y(idx_57,10) * y(idx_57,41) + rc(idx_57,77) * y(idx_57,10) * y(idx_57,29) + rc(idx_57,47) * y(idx_57,2) * y(idx_57,14) + rc(idx_57,66) * y(idx_57,10) * y(idx_57,20) + rc(idx_57,35) * y(idx_57,38) * y(idx_57,5) + rc(idx_57,46) * y(idx_57,2) * y(idx_57,14) + rc(idx_57,22) * y(idx_57,18) * y(idx_57,28) + rc(idx_57,34) * y(idx_57,38) * y(idx_57,4) + rc(idx_57,20) * y(idx_57,18) * y(idx_57,18) * 2.00 + rc(idx_57,21) * y(idx_57,18) * y(idx_57,28) + rc(idx_57,18) * y(idx_57,18) * y(idx_57,5) + rc(idx_57,19) * y(idx_57,18) * y(idx_57,18) + rc(idx_57,6) * y(idx_57,11) * y(idx_57,18) + rc(idx_57,16) * y(idx_57,18) * y(idx_57,4)
        enddo
        !$acc loop gang vector independent
        do idx_58 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_58) = ddep(idx_58,17) + wdep(idx_58,17) + dj(idx_58,4) + rc(idx_58,39) * y(idx_58,5) + rc(idx_58,61) * y(idx_58,10) + dj(idx_58,3)
        enddo
        !$acc loop gang vector independent
        do idx_59 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_59,17) = (yp(idx_59,17) + dts * p(idx_59)) / (1.0 + dts * l(idx_59))
        enddo
        !$acc loop gang vector independent
        do idx_60 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_60) = 0.0 + dj(idx_60,7) * y(idx_60,27) + dj(idx_60,19) * y(idx_60,37) + rc(idx_60,51) * y(idx_60,10) * y(idx_60,14) + rc(idx_60,67) * y(idx_60,10) * y(idx_60,20) + rc(idx_60,27) * y(idx_60,28) * y(idx_60,5) + rc(idx_60,45) * y(idx_60,2) * y(idx_60,14) + rc(idx_60,25) * y(idx_60,25) * y(idx_60,28) + rc(idx_60,26) * y(idx_60,28) * y(idx_60,4) + rc(idx_60,10) * y(idx_60,11) * y(idx_60,28) + rc(idx_60,21) * y(idx_60,18) * y(idx_60,28)
        enddo
        !$acc loop gang vector independent
        do idx_61 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_61) = ddep(idx_61,18) + wdep(idx_61,18) + rc(idx_61,21) * y(idx_61,28) + rc(idx_61,22) * y(idx_61,28) + rc(idx_61,19) * y(idx_61,18) + rc(idx_61,20) * y(idx_61,18) + rc(idx_61,20) * y(idx_61,18) + rc(idx_61,17) * y(idx_61,4) + rc(idx_61,18) * y(idx_61,5) + rc(idx_61,19) * y(idx_61,18) + rc(idx_61,5) * y(idx_61,11) + rc(idx_61,6) * y(idx_61,11) + rc(idx_61,16) * y(idx_61,4)
        enddo
        !$acc loop gang vector independent
        do idx_62 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_62,18) = (yp(idx_62,18) + dts * p(idx_62)) / (1.0 + dts * l(idx_62))
        enddo
        !$acc loop gang vector independent
        do idx_63 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_63) = 0.0 + rc(idx_63,5) * y(idx_63,11) * y(idx_63,18)
        enddo
        !$acc loop gang vector independent
        do idx_64 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_64) = ddep(idx_64,20) + wdep(idx_64,20) + rc(idx_64,66) * y(idx_64,10) + rc(idx_64,67) * y(idx_64,10) + dj(idx_64,1)
        enddo
        !$acc loop gang vector independent
        do idx_65 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_65,20) = (yp(idx_65,20) + dts * p(idx_65)) / (1.0 + dts * l(idx_65))
        enddo
        !$acc loop gang vector independent
        do idx_66 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_66) = 0.0 + rc(idx_66,97) * y(idx_66,10) * y(idx_66,4)
        enddo
        !$acc loop gang vector independent
        do idx_67 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_67) = ddep(idx_67,21) + wdep(idx_67,21) + rc(idx_67,65) * y(idx_67,10) + dj(idx_67,17)
        enddo
        !$acc loop gang vector independent
        do idx_68 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_68,21) = (yp(idx_68,21) + dts * p(idx_68)) / (1.0 + dts * l(idx_68))
        enddo
        !$acc loop gang vector independent
        do idx_69 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_69) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_70 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_70) = ddep(idx_70,24) + wdep(idx_70,24) + rc(idx_70,52) * y(idx_70,10)
        enddo
        !$acc loop gang vector independent
        do idx_71 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_71,24) = (yp(idx_71,24) + dts * p(idx_71)) / (1.0 + dts * l(idx_71))
        enddo
        !$acc loop gang vector independent
        do idx_72 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_72) = 0.0 + dj(idx_72,18) * y(idx_72,35) + rc(idx_72,52) * y(idx_72,10) * y(idx_72,24) + rc(idx_72,58) * y(idx_72,10) * y(idx_72,26) + rc(idx_72,32) * y(idx_72,36) * y(idx_72,4) + rc(idx_72,33) * y(idx_72,36) * y(idx_72,5)
        enddo
        !$acc loop gang vector independent
        do idx_73 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_73) = ddep(idx_73,25) + wdep(idx_73,25) + rc(idx_73,25) * y(idx_73,28) + rc(idx_73,7) * y(idx_73,11) + rc(idx_73,23) * y(idx_73,4) + rc(idx_73,24) * y(idx_73,5)
        enddo
        !$acc loop gang vector independent
        do idx_74 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_74,25) = (yp(idx_74,25) + dts * p(idx_74)) / (1.0 + dts * l(idx_74))
        enddo
        !$acc loop gang vector independent
        do idx_75 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_75) = 0.0 + rc(idx_75,7) * y(idx_75,11) * y(idx_75,25)
        enddo
        !$acc loop gang vector independent
        do idx_76 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_76) = ddep(idx_76,26) + wdep(idx_76,26) + rc(idx_76,57) * y(idx_76,10) + rc(idx_76,58) * y(idx_76,10) + dj(idx_76,1)
        enddo
        !$acc loop gang vector independent
        do idx_77 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_77,26) = (yp(idx_77,26) + dts * p(idx_77)) / (1.0 + dts * l(idx_77))
        enddo
        !$acc loop gang vector independent
        do idx_78 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_78) = 0.0 + rc(idx_78,78) * y(idx_78,10) * y(idx_78,40) + dj(idx_78,1) * y(idx_78,26) + rc(idx_78,25) * y(idx_78,25) * y(idx_78,28) + rc(idx_78,57) * y(idx_78,10) * y(idx_78,26) + rc(idx_78,23) * y(idx_78,25) * y(idx_78,4) + rc(idx_78,24) * y(idx_78,25) * y(idx_78,5)
        enddo
        !$acc loop gang vector independent
        do idx_79 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_79) = ddep(idx_79,27) + wdep(idx_79,27) + dj(idx_79,8) + rc(idx_79,40) * y(idx_79,5) + rc(idx_79,73) * y(idx_79,10) + dj(idx_79,7)
        enddo
        !$acc loop gang vector independent
        do idx_80 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_80,27) = (yp(idx_80,27) + dts * p(idx_80)) / (1.0 + dts * l(idx_80))
        enddo
        !$acc loop gang vector independent
        do idx_81 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_81) = 0.0 + dj(idx_81,19) * y(idx_81,37) + dj(idx_81,1) * y(idx_81,39) + rc(idx_81,93) * y(idx_81,29) + dj(idx_81,16) * y(idx_81,29) + rc(idx_81,40) * y(idx_81,5) * y(idx_81,27) + rc(idx_81,73) * y(idx_81,10) * y(idx_81,27) + rc(idx_81,34) * y(idx_81,38) * y(idx_81,4) + rc(idx_81,35) * y(idx_81,38) * y(idx_81,5)
        enddo
        !$acc loop gang vector independent
        do idx_82 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_82) = ddep(idx_82,28) + wdep(idx_82,28) + rc(idx_82,26) * y(idx_82,4) + rc(idx_82,27) * y(idx_82,5) + rc(idx_82,92) * y(idx_82,6) + rc(idx_82,21) * y(idx_82,18) + rc(idx_82,22) * y(idx_82,18) + rc(idx_82,25) * y(idx_82,25) + rc(idx_82,8) * y(idx_82,11) + rc(idx_82,9) * y(idx_82,11) + rc(idx_82,10) * y(idx_82,11)
        enddo
        !$acc loop gang vector independent
        do idx_83 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_83,28) = (yp(idx_83,28) + dts * p(idx_83)) / (1.0 + dts * l(idx_83))
        enddo
        !$acc loop gang vector independent
        do idx_84 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_84) = 0.0 + rc(idx_84,92) * y(idx_84,28) * y(idx_84,6)
        enddo
        !$acc loop gang vector independent
        do idx_85 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_85) = ddep(idx_85,29) + wdep(idx_85,29) + rc(idx_85,77) * y(idx_85,10) + rc(idx_85,93) + dj(idx_85,16)
        enddo
        !$acc loop gang vector independent
        do idx_86 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_86,29) = (yp(idx_86,29) + dts * p(idx_86)) / (1.0 + dts * l(idx_86))
        enddo
        !$acc loop gang vector independent
        do idx_87 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_87) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_88 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_88) = ddep(idx_88,30) + wdep(idx_88,30) + rc(idx_88,53) * y(idx_88,10) + rc(idx_88,54) * y(idx_88,10)
        enddo
        !$acc loop gang vector independent
        do idx_89 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_89,30) = (yp(idx_89,30) + dts * p(idx_89)) / (1.0 + dts * l(idx_89))
        enddo
        !$acc loop gang vector independent
        do idx_90 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_90) = 0.0 + rc(idx_90,53) * y(idx_90,10) * y(idx_90,30) + rc(idx_90,79) * y(idx_90,10) * y(idx_90,33)
        enddo
        !$acc loop gang vector independent
        do idx_91 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_91) = ddep(idx_91,31) + wdep(idx_91,31) + rc(idx_91,11) * y(idx_91,11) + rc(idx_91,28) * y(idx_91,4) + rc(idx_91,29) * y(idx_91,5)
        enddo
        !$acc loop gang vector independent
        do idx_92 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_92,31) = (yp(idx_92,31) + dts * p(idx_92)) / (1.0 + dts * l(idx_92))
        enddo
        !$acc loop gang vector independent
        do idx_93 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_93) = 0.0 + rc(idx_93,54) * y(idx_93,10) * y(idx_93,30) + rc(idx_93,81) * y(idx_93,10) * y(idx_93,34)
        enddo
        !$acc loop gang vector independent
        do idx_94 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_94) = ddep(idx_94,32) + wdep(idx_94,32) + rc(idx_94,12) * y(idx_94,11) + rc(idx_94,30) * y(idx_94,4) + rc(idx_94,31) * y(idx_94,5)
        enddo
        !$acc loop gang vector independent
        do idx_95 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_95,32) = (yp(idx_95,32) + dts * p(idx_95)) / (1.0 + dts * l(idx_95))
        enddo
        !$acc loop gang vector independent
        do idx_96 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_96) = 0.0 + rc(idx_96,11) * y(idx_96,11) * y(idx_96,31)
        enddo
        !$acc loop gang vector independent
        do idx_97 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_97) = ddep(idx_97,33) + wdep(idx_97,33) + rc(idx_97,79) * y(idx_97,10) + rc(idx_97,80) * y(idx_97,10) + dj(idx_97,1)
        enddo
        !$acc loop gang vector independent
        do idx_98 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_98,33) = (yp(idx_98,33) + dts * p(idx_98)) / (1.0 + dts * l(idx_98))
        enddo
        !$acc loop gang vector independent
        do idx_99 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_99) = 0.0 + rc(idx_99,12) * y(idx_99,11) * y(idx_99,32)
        enddo
        !$acc loop gang vector independent
        do idx_100 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_100) = ddep(idx_100,34) + wdep(idx_100,34) + rc(idx_100,81) * y(idx_100,10) + rc(idx_100,82) * y(idx_100,10) + dj(idx_100,1)
        enddo
        !$acc loop gang vector independent
        do idx_101 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_101,34) = (yp(idx_101,34) + dts * p(idx_101)) / (1.0 + dts * l(idx_101))
        enddo
        !$acc loop gang vector independent
        do idx_102 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_102) = 0.0 + rc(idx_102,80) * y(idx_102,10) * y(idx_102,33) + dj(idx_102,1) * y(idx_102,33) + rc(idx_102,28) * y(idx_102,31) * y(idx_102,4) + rc(idx_102,29) * y(idx_102,31) * y(idx_102,5)
        enddo
        !$acc loop gang vector independent
        do idx_103 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_103) = ddep(idx_103,35) + wdep(idx_103,35) + rc(idx_103,41) * y(idx_103,5) + rc(idx_103,56) * y(idx_103,10) + dj(idx_103,18)
        enddo
        !$acc loop gang vector independent
        do idx_104 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_104,35) = (yp(idx_104,35) + dts * p(idx_104)) / (1.0 + dts * l(idx_104))
        enddo
        !$acc loop gang vector independent
        do idx_105 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_105) = 0.0 + rc(idx_105,101) * y(idx_105,40) + dj(idx_105,16) * y(idx_105,40) + rc(idx_105,41) * y(idx_105,5) * y(idx_105,35) + rc(idx_105,56) * y(idx_105,10) * y(idx_105,35)
        enddo
        !$acc loop gang vector independent
        do idx_106 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_106) = ddep(idx_106,36) + wdep(idx_106,36) + rc(idx_106,33) * y(idx_106,5) + rc(idx_106,100) * y(idx_106,6) + rc(idx_106,13) * y(idx_106,11) + rc(idx_106,14) * y(idx_106,11) + rc(idx_106,32) * y(idx_106,4)
        enddo
        !$acc loop gang vector independent
        do idx_107 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_107,36) = (yp(idx_107,36) + dts * p(idx_107)) / (1.0 + dts * l(idx_107))
        enddo
        !$acc loop gang vector independent
        do idx_108 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_108) = 0.0 + rc(idx_108,82) * y(idx_108,10) * y(idx_108,34) + dj(idx_108,1) * y(idx_108,34) + rc(idx_108,30) * y(idx_108,32) * y(idx_108,4) + rc(idx_108,31) * y(idx_108,32) * y(idx_108,5)
        enddo
        !$acc loop gang vector independent
        do idx_109 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_109) = ddep(idx_109,37) + wdep(idx_109,37) + dj(idx_109,19) + rc(idx_109,42) * y(idx_109,5) + rc(idx_109,69) * y(idx_109,10) + rc(idx_109,70) * y(idx_109,10)
        enddo
        !$acc loop gang vector independent
        do idx_110 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_110,37) = (yp(idx_110,37) + dts * p(idx_110)) / (1.0 + dts * l(idx_110))
        enddo
        !$acc loop gang vector independent
        do idx_111 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_111) = 0.0 + rc(idx_111,70) * y(idx_111,10) * y(idx_111,37) + rc(idx_111,71) * y(idx_111,10) * y(idx_111,39) + rc(idx_111,42) * y(idx_111,5) * y(idx_111,37) + rc(idx_111,69) * y(idx_111,10) * y(idx_111,37)
        enddo
        !$acc loop gang vector independent
        do idx_112 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_112) = ddep(idx_112,38) + wdep(idx_112,38) + rc(idx_112,15) * y(idx_112,11) + rc(idx_112,34) * y(idx_112,4) + rc(idx_112,35) * y(idx_112,5)
        enddo
        !$acc loop gang vector independent
        do idx_113 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_113,38) = (yp(idx_113,38) + dts * p(idx_113)) / (1.0 + dts * l(idx_113))
        enddo
        !$acc loop gang vector independent
        do idx_114 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_114) = 0.0 + rc(idx_114,15) * y(idx_114,11) * y(idx_114,38)
        enddo
        !$acc loop gang vector independent
        do idx_115 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_115) = ddep(idx_115,39) + wdep(idx_115,39) + rc(idx_115,71) * y(idx_115,10) + rc(idx_115,72) * y(idx_115,10) + dj(idx_115,1)
        enddo
        !$acc loop gang vector independent
        do idx_116 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_116,39) = (yp(idx_116,39) + dts * p(idx_116)) / (1.0 + dts * l(idx_116))
        enddo
        !$acc loop gang vector independent
        do idx_117 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_117) = 0.0 + rc(idx_117,100) * y(idx_117,36) * y(idx_117,6)
        enddo
        !$acc loop gang vector independent
        do idx_118 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_118) = ddep(idx_118,40) + wdep(idx_118,40) + rc(idx_118,78) * y(idx_118,10) + rc(idx_118,101) + dj(idx_118,16)
        enddo
        !$acc loop gang vector independent
        do idx_119 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_119,40) = (yp(idx_119,40) + dts * p(idx_119)) / (1.0 + dts * l(idx_119))
        enddo
        !$acc loop gang vector independent
        do idx_120 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_120) = 0.0 + rc(idx_120,17) * y(idx_120,18) * y(idx_120,4)
        enddo
        !$acc loop gang vector independent
        do idx_121 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_121) = ddep(idx_121,41) + wdep(idx_121,41) + rc(idx_121,68) * y(idx_121,10) + dj(idx_121,20)
        enddo
        !$acc loop gang vector independent
        do idx_122 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_122,41) = (yp(idx_122,41) + dts * p(idx_122)) / (1.0 + dts * l(idx_122))
        enddo
        !$acc loop gang vector independent
        do idx_123 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_123) = 0.0 + dj(idx_123,15) * y(idx_123,44) + rc(idx_123,87) * y(idx_123,43) * y(idx_123,23) + rc(idx_123,88) * y(idx_123,43) * y(idx_123,22)
        enddo
        !$acc loop gang vector independent
        do idx_124 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_124) = ddep(idx_124,42) + wdep(idx_124,42) + rc(idx_124,102) * y(idx_124,22)
        enddo
        !$acc loop gang vector independent
        do idx_125 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_125,42) = (yp(idx_125,42) + dts * p(idx_125)) / (1.0 + dts * l(idx_125))
        enddo
        !$acc loop gang vector independent
        do idx_126 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_126) = 0.0 + dj(idx_126,14) * y(idx_126,44)
        enddo
        !$acc loop gang vector independent
        do idx_127 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_127) = ddep(idx_127,43) + wdep(idx_127,43) + rc(idx_127,86) * y(idx_127,19) + rc(idx_127,87) * y(idx_127,23) + rc(idx_127,88) * y(idx_127,22)
        enddo
        !$acc loop gang vector independent
        do idx_128 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_128,43) = (yp(idx_128,43) + dts * p(idx_128)) / (1.0 + dts * l(idx_128))
        enddo
        !$acc loop gang vector independent
        do idx_129 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_129) = 0.0 + rc(idx_129,102) * y(idx_129,42) * y(idx_129,22)
        enddo
        !$acc loop gang vector independent
        do idx_130 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_130) = ddep(idx_130,44) + wdep(idx_130,44) + dj(idx_130,15) + rc(idx_130,84) * y(idx_130,46) + rc(idx_130,85) * y(idx_130,45) + dj(idx_130,14)
        enddo
        !$acc loop gang vector independent
        do idx_131 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_131,44) = (yp(idx_131,44) + dts * p(idx_131)) / (1.0 + dts * l(idx_131))
        enddo
        !$acc loop gang vector independent
        do idx_132 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_132) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_133 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_133) = ddep(idx_133,45) + wdep(idx_133,45) + rc(idx_133,85) * y(idx_133,44)
        enddo
        !$acc loop gang vector independent
        do idx_134 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_134,45) = (yp(idx_134,45) + dts * p(idx_134)) / (1.0 + dts * l(idx_134))
        enddo
        !$acc loop gang vector independent
        do idx_135 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
          p(idx_135) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_136 = LBOUND(l, dim=1), UBOUND(l, dim=1), 1
          l(idx_136) = ddep(idx_136,46) + wdep(idx_136,46) + rc(idx_136,84) * y(idx_136,44)
        enddo
        !$acc loop gang vector independent
        do idx_137 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_137,46) = (yp(idx_137,46) + dts * p(idx_137)) / (1.0 + dts * l(idx_137))
        enddo
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_deriv

end module ukca_deriv_mod
