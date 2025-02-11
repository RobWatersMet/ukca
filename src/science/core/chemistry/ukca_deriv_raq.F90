module ukca_deriv_raq_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DERIV_RAQ_MOD'
  public

  contains
  subroutine ukca_deriv_raq(nr_therm, n_be_calls, n_pnts, dts, h2o, vol, rc, hrc, dd, dw, dj, ldepd, ldepw, y, ddflux, wdflux)
    use ukca_chemco_raq_mod, only : in2o5, in2o5_ss
    use asad_mod, only : jpdd, jpdw, jphk, jppj, jpspec
    use ukca_config_specification_mod, only : ukca_config
    use ukca_config_constants_mod, only : avogadro
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    use ukca_um_legacy_mod, only : oneover_v
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr_therm
    integer, intent(in) :: n_be_calls
    integer, intent(in) :: n_pnts
    real, intent(in) :: dts
    real, dimension(n_pnts), intent(in) :: h2o
    real, dimension(n_pnts), intent(in) :: vol
    real, dimension(n_pnts,nr_therm), intent(in) :: rc
    real, dimension(n_pnts,jphk), intent(in) :: hrc
    real, dimension(n_pnts,jpspec), intent(in) :: dd
    real, dimension(n_pnts,jpspec), intent(in) :: dw
    real, dimension(n_pnts,jppj), intent(in) :: dj
    logical, dimension(jpspec), intent(in) :: ldepd
    logical, dimension(jpspec), intent(in) :: ldepw
    real, dimension(n_pnts,jpspec), intent(inout) :: y
    real, dimension(n_pnts,jpdd), intent(out) :: ddflux
    real, dimension(n_pnts,jpdw), intent(out) :: wdflux
    integer :: i
    integer :: j
    integer :: n
    integer :: icnt_dd
    integer :: icnt_wd
    real :: recip_avogad
    real, dimension(n_pnts) :: p
    real, dimension(n_pnts) :: l
    real, dimension(n_pnts,jpspec) :: yp
    real, dimension(n_pnts) :: tmp1
    real, dimension(n_pnts) :: tmp2
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DERIV_RAQ'
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
    integer :: idx_172
    integer :: idx_173
    integer :: idx_174
    integer :: idx_175
    integer :: idx_176
    integer :: idx_177
    integer :: idx_178
    integer :: idx_179
    integer :: idx_180
    integer :: idx_181
    integer :: idx_182
    integer :: idx_183
    integer :: idx_184
    integer :: idx_185
    integer :: idx_186
    integer :: idx_187
    integer :: idx_188
    integer :: idx_189
    integer :: idx_190
    integer :: idx_191
    integer :: idx_192
    integer :: idx_193
    integer :: idx_194
    integer :: idx_195
    integer :: idx_196
    integer :: idx_197
    integer :: idx_198
    integer :: idx_199
    integer :: idx_200
    integer :: idx_201
    integer :: idx_202
    integer :: idx_203
    integer :: idx_204
    integer :: idx_205
    integer :: idx_206
    integer :: idx_207
    integer :: idx_208
    integer :: idx_209
    integer :: idx_210
    integer :: idx_211
    integer :: idx_212
    integer :: idx_213
    integer :: idx_214
    integer :: idx_215
    integer :: idx_216
    integer :: idx_217
    integer :: idx_218
    integer :: idx_219
    integer :: idx_220
    integer :: idx_221
    integer :: idx_222
    integer :: idx_223
    integer :: idx_224
    integer :: idx_225
    integer :: idx_226
    integer :: idx_227
    integer :: idx_228
    integer :: idx_229
    integer :: idx_230
    integer :: idx_231
    integer :: idx_232
    integer :: idx_233
    integer :: idx_234
    integer :: idx_235
    integer :: idx_236
    integer :: idx_237
    integer :: idx_238
    integer :: idx_239
    integer :: idx_240
    integer :: idx_241
    integer :: idx_242
    integer :: idx_243
    integer :: idx_244
    integer :: idx_245
    integer :: idx_246
    integer :: idx_247
    integer :: idx_248
    integer :: idx_249
    integer :: idx_250
    integer :: idx_251
    integer :: idx_252
    integer :: idx_253
    integer :: idx_254
    integer :: idx_255
    integer :: idx_256
    integer :: idx_257
    integer :: idx_258
    integer :: idx_259
    integer :: idx_260
    integer :: idx_261
    integer :: idx_262
    integer :: idx_263
    integer :: idx_264
    integer :: idx_265
    integer :: idx_266
    integer :: idx_267
    integer :: idx_268
    integer :: idx_269
    integer :: idx_270
    integer :: idx_271
    integer :: idx_272
    integer :: idx_273
    integer :: idx_274
    integer :: idx_275
    integer :: idx_276
    integer :: idx_277
    integer :: idx_278
    integer :: idx_279
    integer :: idx_280
    integer :: idx_281
    integer :: idx_282
    integer :: idx_283
    integer :: idx_284
    integer :: idx_285
    integer :: idx_286
    integer :: idx_287
    integer :: idx_288
    integer :: idx_289
    integer :: idx_290

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    recip_avogad = 1.0 / avogadro
    !$acc loop gang vector independent
    do idx = 1, n_pnts, 1
      tmp1(idx) = 0.0
    enddo
    !$acc loop gang vector independent
    do idx_1 = 1, n_pnts, 1
      tmp2(idx_1) = 0.0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_2 = LBOUND(ddflux, dim=2), UBOUND(ddflux, dim=2), 1
      do idx_3 = LBOUND(ddflux, dim=1), UBOUND(ddflux, dim=1), 1
        ddflux(idx_3,idx_2) = 0.0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_4 = LBOUND(wdflux, dim=2), UBOUND(wdflux, dim=2), 1
      do idx_5 = LBOUND(wdflux, dim=1), UBOUND(wdflux, dim=1), 1
        wdflux(idx_5,idx_4) = 0.0
      enddo
    enddo
    !$acc end kernels
    do n = 1, n_be_calls, 1
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do j = 1, jpspec, 1
        do idx_6 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
          yp(idx_6,j) = y(idx_6,j)
        enddo
      enddo
      !$acc end kernels
      do i = 1, ukca_config%nit, 1
        !$acc kernels
        !$acc loop gang vector independent
        do idx_7 = 1, n_pnts, 1
          p(idx_7) = 0.0 + dj(idx_7,3) * y(idx_7,7) + dj(idx_7,15) * y(idx_7,6) + rc(idx_7,7) * y(idx_7,2) + dj(idx_7,1) * y(idx_7,4)
        enddo
        !$acc loop gang vector independent
        do idx_8 = 1, n_pnts, 1
          l(idx_8) = 0.0 + rc(idx_8,1) + rc(idx_8,5) * y(idx_8,5) + rc(idx_8,16) * y(idx_8,7)
        enddo
        !$acc loop gang vector independent
        do idx_9 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_9) = 1.0 + dts * l(idx_9)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_10 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_10,1) = (yp(idx_10,1) + dts * p(idx_10)) * tmp1(idx_10)
        enddo
        !$acc loop gang vector independent
        do idx_11 = 1, n_pnts, 1
          p(idx_11) = 0.0 + dj(idx_11,2) * y(idx_11,4)
        enddo
        !$acc loop gang vector independent
        do idx_12 = 1, n_pnts, 1
          l(idx_12) = 0.0 + rc(idx_12,7) + rc(idx_12,8) * h2o(idx_12)
        enddo
        !$acc loop gang vector independent
        do idx_13 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_13) = 1.0 + dts * l(idx_13)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_14 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_14,2) = (yp(idx_14,2) + dts * p(idx_14)) * tmp1(idx_14)
        enddo
        !$acc loop gang vector independent
        do idx_15 = 1, n_pnts, 1
          p(idx_15) = 0.0 + dj(idx_15,23) * y(idx_15,35) + dj(idx_15,21) * y(idx_15,52) + dj(idx_15,22) * y(idx_15,31) + dj(idx_15,19) * y(idx_15,21) + dj(idx_15,20) * y(idx_15,26) + dj(idx_15,5) * y(idx_15,10)
        enddo
        !$acc loop gang vector independent
        do idx_16 = 1, n_pnts, 1
          p(idx_16) = p(idx_16) + dj(idx_16,17) * y(idx_16,17) + dj(idx_16,4) * y(idx_16,11) + dj(idx_16,4) * y(idx_16,11) + rc(idx_16,137) * y(idx_16,35) * y(idx_16,3) + rc(idx_16,182) * y(idx_16,15) * y(idx_16,55) + rc(idx_16,129) * y(idx_16,4) * y(idx_16,34) * 0.36
        enddo
        !$acc loop gang vector independent
        do idx_17 = 1, n_pnts, 1
          p(idx_17) = p(idx_17) + rc(idx_17,135) * y(idx_17,31) * y(idx_17,3) + rc(idx_17,123) * y(idx_17,4) * y(idx_17,50) * 0.28 + rc(idx_17,128) * y(idx_17,4) * y(idx_17,29) * 0.27 + rc(idx_17,102) * y(idx_17,3) * y(idx_17,26) + rc(idx_17,104) * y(idx_17,3) * y(idx_17,52) + rc(idx_17,45) * y(idx_17,3) * y(idx_17,17)
        enddo
        !$acc loop gang vector independent
        do idx_18 = 1, n_pnts, 1
          p(idx_18) = p(idx_18) + rc(idx_18,100) * y(idx_18,3) * y(idx_18,21) + rc(idx_18,17) * y(idx_18,5) * y(idx_18,15) + rc(idx_18,34) * y(idx_18,6) * y(idx_18,15) + rc(idx_18,8) * y(idx_18,2) * h2o(idx_18) * 2.00 + rc(idx_18,14) * y(idx_18,15) * y(idx_18,4)
        enddo
        !$acc loop gang vector independent
        do idx_19 = 1, n_pnts, 1
          l(idx_19) = 0.0 + rc(idx_19,177) * y(idx_19,39) + rc(idx_19,174) * y(idx_19,54) + rc(idx_19,175) * y(idx_19,54) + rc(idx_19,176) * y(idx_19,58) + rc(idx_19,160) * y(idx_19,32) + rc(idx_19,170) * y(idx_19,58)
        enddo
        !$acc loop gang vector independent
        do idx_20 = 1, n_pnts, 1
          l(idx_20) = l(idx_20) + rc(idx_20,172) * y(idx_20,56) + rc(idx_20,143) * y(idx_20,34) + rc(idx_20,154) * y(idx_20,45) + rc(idx_20,156) * y(idx_20,47) + rc(idx_20,138) * y(idx_20,33) + rc(idx_20,139) * y(idx_20,57)
        enddo
        !$acc loop gang vector independent
        do idx_21 = 1, n_pnts, 1
          l(idx_21) = l(idx_21) + rc(idx_21,141) * y(idx_21,29) + rc(idx_21,125) * y(idx_21,50) + rc(idx_21,135) * y(idx_21,31) + rc(idx_21,137) * y(idx_21,35) + rc(idx_21,102) * y(idx_21,26) + rc(idx_21,104) * y(idx_21,52)
        enddo
        !$acc loop gang vector independent
        do idx_22 = 1, n_pnts, 1
          l(idx_22) = l(idx_22) + rc(idx_22,109) * y(idx_22,49) + rc(idx_22,94) * y(idx_22,27) + rc(idx_22,98) * y(idx_22,23) + rc(idx_22,100) * y(idx_22,21) + rc(idx_22,81) * y(idx_22,51) + rc(idx_22,86) * y(idx_22,53)
        enddo
        !$acc loop gang vector independent
        do idx_23 = 1, n_pnts, 1
          l(idx_23) = l(idx_23) + rc(idx_23,92) * y(idx_23,25) + rc(idx_23,70) * y(idx_23,13) + rc(idx_23,71) * y(idx_23,19) + rc(idx_23,75) * y(idx_23,22) + rc(idx_23,59) * y(idx_23,12) + rc(idx_23,63) * y(idx_23,41)
        enddo
        !$acc loop gang vector independent
        do idx_24 = 1, n_pnts, 1
          l(idx_24) = l(idx_24) + rc(idx_24,66) * y(idx_24,14) + rc(idx_24,44) * y(idx_24,21) + rc(idx_24,44) * y(idx_24,26) + rc(idx_24,44) * y(idx_24,52) + rc(idx_24,35) * y(idx_24,10) + rc(idx_24,44) * y(idx_24,17)
        enddo
        !$acc loop gang vector independent
        do idx_25 = 1, n_pnts, 1
          l(idx_25) = l(idx_25) + rc(idx_25,45) * y(idx_25,17) + rc(idx_25,30) * y(idx_25,15) + rc(idx_25,31) * y(idx_25,11)
        enddo
        !$acc loop gang vector independent
        do idx_26 = 1, n_pnts, 1
          l(idx_26) = l(idx_26) + rc(idx_26,33) * y(idx_26,43) + rc(idx_26,13) * y(idx_26,4) + rc(idx_26,21) * y(idx_26,7) + rc(idx_26,24) * y(idx_26,9)
        enddo
        !$acc loop gang vector independent
        do idx_27 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_27) = 1.0 + dts * l(idx_27)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_28 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_28,3) = (yp(idx_28,3) + dts * p(idx_28)) * tmp1(idx_28)
        enddo
        !$acc loop gang vector independent
        do idx_29 = 1, n_pnts, 1
          p(idx_29) = 0.0 + rc(idx_29,1) * y(idx_29,1) + rc(idx_29,131) * y(idx_29,20) * y(idx_29,15) * 0.3
        enddo
        !$acc loop gang vector independent
        do idx_30 = 1, n_pnts, 1
          l(idx_30) = 0.0 + dj(idx_30,1) + dj(idx_30,2) + dd(idx_30,4) + rc(idx_30,124) * y(idx_30,50) + rc(idx_30,128) * y(idx_30,29) + rc(idx_30,129) * y(idx_30,34)
        enddo
        !$acc loop gang vector independent
        do idx_31 = 1, n_pnts, 1
          l(idx_31) = l(idx_31) + rc(idx_31,14) * y(idx_31,15) + rc(idx_31,112) * y(idx_31,49) + rc(idx_31,123) * y(idx_31,50) + rc(idx_31,11) * y(idx_31,5) + rc(idx_31,12) * y(idx_31,7) + rc(idx_31,13) * y(idx_31,3)
        enddo
        !$acc loop gang vector independent
        do idx_32 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_32) = 1.0 + dts * l(idx_32)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_33 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_33,4) = (yp(idx_33,4) + dts * p(idx_33)) * tmp1(idx_33)
        enddo
        !$acc loop gang vector independent
        do idx_34 = 1, n_pnts, 1
          p(idx_34) = 0.0 + dj(idx_34,3) * y(idx_34,7) + dj(idx_34,14) * y(idx_34,6) + rc(idx_34,16) * y(idx_34,7) * y(idx_34,1) + rc(idx_34,19) * y(idx_34,7) * y(idx_34,6)
        enddo
        !$acc loop gang vector independent
        do idx_35 = 1, n_pnts, 1
          l(idx_35) = 0.0 + rc(idx_35,173) * y(idx_35,44) + rc(idx_35,126) * y(idx_35,40) + rc(idx_35,142) * y(idx_35,46) + rc(idx_35,144) * y(idx_35,48) + rc(idx_35,95) * y(idx_35,36) + rc(idx_35,105) * y(idx_35,37)
        enddo
        !$acc loop gang vector independent
        do idx_36 = 1, n_pnts, 1
          l(idx_36) = l(idx_36) + rc(idx_36,110) * y(idx_36,38) + rc(idx_36,79) * y(idx_36,20) + rc(idx_36,83) * y(idx_36,24) + rc(idx_36,93) * y(idx_36,30) + rc(idx_36,17) * y(idx_36,15)
        enddo
        !$acc loop gang vector independent
        do idx_37 = 1, n_pnts, 1
          l(idx_37) = l(idx_37) + rc(idx_37,60) * y(idx_37,16) + rc(idx_37,72) * y(idx_37,18) + rc(idx_37,5) * y(idx_37,1) + rc(idx_37,11) * y(idx_37,4) + rc(idx_37,15) * y(idx_37,6)
        enddo
        !$acc loop gang vector independent
        do idx_38 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_38) = 1.0 + dts * l(idx_38)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_39 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_39,5) = (yp(idx_39,5) + dts * p(idx_39)) * tmp1(idx_39)
        enddo
        !$acc loop gang vector independent
        do idx_40 = 1, n_pnts, 1
          p(idx_40) = 0.0 + dj(idx_40,16) * y(idx_40,8) + rc(idx_40,35) * y(idx_40,3) * y(idx_40,10) + rc(idx_40,98) * y(idx_40,3) * y(idx_40,23) + rc(idx_40,12) * y(idx_40,7) * y(idx_40,4) + rc(idx_40,29) * y(idx_40,8)
        enddo
        !$acc loop gang vector independent
        do idx_41 = 1, n_pnts, 1
          l(idx_41) = 0.0 + dj(idx_41,14) + dj(idx_41,15) + dw(idx_41,6) + rc(idx_41,158) * y(idx_41,22) + rc(idx_41,159) * y(idx_41,29) + rc(idx_41,178) * y(idx_41,39)
        enddo
        !$acc loop gang vector independent
        do idx_42 = 1, n_pnts, 1
          l(idx_42) = l(idx_42) + rc(idx_42,152) * y(idx_42,51) + rc(idx_42,153) * y(idx_42,49) + rc(idx_42,155) * y(idx_42,50) + rc(idx_42,34) * y(idx_42,15) + rc(idx_42,67) * y(idx_42,14) + rc(idx_42,151) * y(idx_42,19)
        enddo
        !$acc loop gang vector independent
        do idx_43 = 1, n_pnts, 1
          l(idx_43) = l(idx_43) + rc(idx_43,27) * y(idx_43,6) + rc(idx_43,27) * y(idx_43,6) + rc(idx_43,32) * y(idx_43,15) + rc(idx_43,15) * y(idx_43,5) + rc(idx_43,19) * y(idx_43,7) + rc(idx_43,20) * y(idx_43,7)
        enddo
        !$acc loop gang vector independent
        do idx_44 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_44) = 1.0 + dts * l(idx_44)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_45 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_45,6) = (yp(idx_45,6) + dts * p(idx_45)) * tmp1(idx_45)
        enddo
        !$acc loop gang vector independent
        do idx_46 = 1, n_pnts, 1
          p(idx_46) = 0.0 + dj(idx_46,18) * y(idx_46,23) + dj(idx_46,15) * y(idx_46,6) + dj(idx_46,16) * y(idx_46,8) + dj(idx_46,5) * y(idx_46,10) + dj(idx_46,11) * y(idx_46,9) + rc(idx_46,177) * y(idx_46,3) * y(idx_46,39)
        enddo
        !$acc loop gang vector independent
        do idx_47 = 1, n_pnts, 1
          p(idx_47) = p(idx_47) + rc(idx_47,178) * y(idx_47,6) * y(idx_47,39) * 2.00 + rc(idx_47,160) * y(idx_47,32) * y(idx_47,3) + rc(idx_47,173) * y(idx_47,44) * y(idx_47,5) + rc(idx_47,154) * y(idx_47,45) * y(idx_47,3) + rc(idx_47,156) * y(idx_47,47) * y(idx_47,3) + rc(idx_47,142) * y(idx_47,46) * y(idx_47,5)
        enddo
        !$acc loop gang vector independent
        do idx_48 = 1, n_pnts, 1
          p(idx_48) = p(idx_48) + rc(idx_48,144) * y(idx_48,48) * y(idx_48,5) + rc(idx_48,110) * y(idx_48,38) * y(idx_48,5) + rc(idx_48,126) * y(idx_48,40) * y(idx_48,5) + rc(idx_48,95) * y(idx_48,36) * y(idx_48,5) + rc(idx_48,105) * y(idx_48,5) * y(idx_48,37) + rc(idx_48,83) * y(idx_48,24) * y(idx_48,5)
        enddo
        !$acc loop gang vector independent
        do idx_49 = 1, n_pnts, 1
          p(idx_49) = p(idx_49) + rc(idx_49,93) * y(idx_49,30) * y(idx_49,5) + rc(idx_49,78) * y(idx_49,23) + rc(idx_49,79) * y(idx_49,20) * y(idx_49,5) + rc(idx_49,60) * y(idx_49,5) * y(idx_49,16) + rc(idx_49,72) * y(idx_49,18) * y(idx_49,5) + rc(idx_49,29) * y(idx_49,8)
        enddo
        !$acc loop gang vector independent
        do idx_50 = 1, n_pnts, 1
          p(idx_50) = p(idx_50) + rc(idx_50,34) * y(idx_50,6) * y(idx_50,15) + rc(idx_50,24) * y(idx_50,9) * y(idx_50,3) + rc(idx_50,27) * y(idx_50,6) * y(idx_50,6) * 2.00 + rc(idx_50,19) * y(idx_50,7) * y(idx_50,6) + rc(idx_50,23) * y(idx_50,9)
        enddo
        !$acc loop gang vector independent
        do idx_51 = 1, n_pnts, 1
          p(idx_51) = p(idx_51) + rc(idx_51,15) * y(idx_51,5) * y(idx_51,6) * 2.00 + rc(idx_51,17) * y(idx_51,5) * y(idx_51,15) + rc(idx_51,5) * y(idx_51,1) * y(idx_51,5) + rc(idx_51,11) * y(idx_51,5) * y(idx_51,4)
        enddo
        !$acc loop gang vector independent
        do idx_52 = 1, n_pnts, 1
          l(idx_52) = 0.0 + dj(idx_52,3) + dd(idx_52,7) + rc(idx_52,77) * y(idx_52,20) + rc(idx_52,171) * y(idx_52,42) + rc(idx_52,183) * y(idx_52,55)
        enddo
        !$acc loop gang vector independent
        do idx_53 = 1, n_pnts, 1
          l(idx_53) = l(idx_53) + rc(idx_53,20) * y(idx_53,6) + rc(idx_53,21) * y(idx_53,3) + rc(idx_53,22) * y(idx_53,15) + rc(idx_53,12) * y(idx_53,4) + rc(idx_53,16) * y(idx_53,1) + rc(idx_53,19) * y(idx_53,6)
        enddo
        !$acc loop gang vector independent
        do idx_54 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_54) = 1.0 + dts * l(idx_54)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_55 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_55,7) = (yp(idx_55,7) + dts * p(idx_55)) * tmp1(idx_55)
        enddo
        !$acc loop gang vector independent
        do idx_56 = 1, n_pnts, 1
          p(idx_56) = 0.0 + rc(idx_56,20) * y(idx_56,7) * y(idx_56,6)
        enddo
        !$acc loop gang vector independent
        do idx_57 = 1, n_pnts, 1
          l(idx_57) = 0.0 + rc(idx_57,29) + dj(idx_57,16) + dw(idx_57,8)
        enddo
        !$acc end kernels
        if (in2o5 > 0) then
          !$acc kernels
          !$acc loop gang vector independent
          do idx_58 = 1, n_pnts, 1
            l(idx_58) = l(idx_58) + hrc(idx_58,in2o5)
          enddo
          !$acc end kernels
        end if
        if (in2o5_ss > 0) then
          !$acc kernels
          !$acc loop gang vector independent
          do idx_59 = 1, n_pnts, 1
            l(idx_59) = l(idx_59) + hrc(idx_59,in2o5_ss)
          enddo
          !$acc end kernels
        end if
        !$acc kernels
        !$acc loop gang vector independent
        do idx_60 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_60) = 1.0 + dts * l(idx_60)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_61 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_61,8) = (yp(idx_61,8) + dts * p(idx_61)) * tmp1(idx_61)
        enddo
        !$acc loop gang vector independent
        do idx_62 = 1, n_pnts, 1
          p(idx_62) = 0.0 + rc(idx_62,22) * y(idx_62,7) * y(idx_62,15)
        enddo
        !$acc loop gang vector independent
        do idx_63 = 1, n_pnts, 1
          l(idx_63) = 0.0 + dw(idx_63,9) + rc(idx_63,23) + rc(idx_63,24) * y(idx_63,3) + dj(idx_63,11)
        enddo
        !$acc loop gang vector independent
        do idx_64 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_64) = 1.0 + dts * l(idx_64)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_65 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_65,9) = (yp(idx_65,9) + dts * p(idx_65)) * tmp1(idx_65)
        enddo
        !$acc loop gang vector independent
        do idx_66 = 1, n_pnts, 1
          p(idx_66) = 0.0 + rc(idx_66,152) * y(idx_66,6) * y(idx_66,51) + rc(idx_66,158) * y(idx_66,6) * y(idx_66,22) + rc(idx_66,67) * y(idx_66,6) * y(idx_66,14) + rc(idx_66,151) * y(idx_66,6) * y(idx_66,19) + rc(idx_66,21) * y(idx_66,7) * y(idx_66,3) + rc(idx_66,32) * y(idx_66,6) * y(idx_66,15)
        enddo
        !$acc end kernels
        if (in2o5 > 0) then
          !$acc kernels
          !$acc loop gang vector independent
          do idx_67 = 1, n_pnts, 1
            p(idx_67) = p(idx_67) + hrc(idx_67,in2o5) * y(idx_67,8) * 2.0
          enddo
          !$acc end kernels
        end if
        !$acc kernels
        !$acc loop gang vector independent
        do idx_68 = 1, n_pnts, 1
          l(idx_68) = 0.0 + dd(idx_68,10) + rc(idx_68,35) * y(idx_68,3) + dj(idx_68,5) + dw(idx_68,10)
        enddo
        !$acc loop gang vector independent
        do idx_69 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_69) = 1.0 + dts * l(idx_69)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_70 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_70,10) = (yp(idx_70,10) + dts * p(idx_70)) * tmp1(idx_70)
        enddo
        !$acc loop gang vector independent
        do idx_71 = 1, n_pnts, 1
          p(idx_71) = 0.0 + rc(idx_71,36) * y(idx_71,15) * y(idx_71,15)
        enddo
        !$acc loop gang vector independent
        do idx_72 = 1, n_pnts, 1
          l(idx_72) = 0.0 + dd(idx_72,11) + rc(idx_72,31) * y(idx_72,3) + dj(idx_72,4) + dw(idx_72,11)
        enddo
        !$acc loop gang vector independent
        do idx_73 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_73) = 1.0 + dts * l(idx_73)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_74 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_74,11) = (yp(idx_74,11) + dts * p(idx_74)) * tmp1(idx_74)
        enddo
        !$acc loop gang vector independent
        do idx_75 = 1, n_pnts, 1
          p(idx_75) = 0.0 + rc(idx_75,123) * y(idx_75,4) * y(idx_75,50) * 0.30
        enddo
        !$acc loop gang vector independent
        do idx_76 = 1, n_pnts, 1
          l(idx_76) = 0.0 + rc(idx_76,59) * y(idx_76,3) + dd(idx_76,12)
        enddo
        !$acc loop gang vector independent
        do idx_77 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_77) = 1.0 + dts * l(idx_77)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_78 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_78,12) = (yp(idx_78,12) + dts * p(idx_78)) * tmp1(idx_78)
        enddo
        !$acc loop gang vector independent
        do idx_79 = 1, n_pnts, 1
          p(idx_79) = 0.0 + dj(idx_79,13) * y(idx_79,57) + dj(idx_79,12) * y(idx_79,33) + dj(idx_79,13) * y(idx_79,57) + dj(idx_79,7) * y(idx_79,14) + dj(idx_79,8) * y(idx_79,22) + rc(idx_79,139) * y(idx_79,57) * y(idx_79,3) * 2.0
        enddo
        !$acc loop gang vector independent
        do idx_80 = 1, n_pnts, 1
          p(idx_80) = p(idx_80) + dj(idx_80,6) * y(idx_80,14) + rc(idx_80,129) * y(idx_80,4) * y(idx_80,34) * 0.76 + rc(idx_80,138) * y(idx_80,33) * y(idx_80,3) + rc(idx_80,124) * y(idx_80,4) * y(idx_80,50) * 0.58 + rc(idx_80,128) * y(idx_80,4) * y(idx_80,29) * 0.78
        enddo
        !$acc loop gang vector independent
        do idx_81 = 1, n_pnts, 1
          p(idx_81) = p(idx_81) + rc(idx_81,112) * y(idx_81,4) * y(idx_81,49) * 0.31 + rc(idx_81,123) * y(idx_81,4) * y(idx_81,50) * 0.40 + rc(idx_81,66) * y(idx_81,3) * y(idx_81,14) + rc(idx_81,67) * y(idx_81,6) * y(idx_81,14)
        enddo
        !$acc loop gang vector independent
        do idx_82 = 1, n_pnts, 1
          l(idx_82) = 0.0 + rc(idx_82,70) * y(idx_82,3) + dd(idx_82,13)
        enddo
        !$acc loop gang vector independent
        do idx_83 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_83) = 1.0 + dts * l(idx_83)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_84 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_84,13) = (yp(idx_84,13) + dts * p(idx_84)) * tmp1(idx_84)
        enddo
        !$acc loop gang vector independent
        do idx_85 = 1, n_pnts, 1
          p(idx_85) = 0.0 + dj(idx_85,22) * y(idx_85,31) + dj(idx_85,23) * y(idx_85,35) + rc(idx_85,181) * y(idx_85,44) * y(idx_85,16) + dj(idx_85,17) * y(idx_85,17) + rc(idx_85,144) * y(idx_85,48) * y(idx_85,5) + rc(idx_85,154) * y(idx_85,45) * y(idx_85,3)
        enddo
        !$acc loop gang vector independent
        do idx_86 = 1, n_pnts, 1
          p(idx_86) = p(idx_86) + rc(idx_86,137) * y(idx_86,35) * y(idx_86,3) + rc(idx_86,142) * y(idx_86,46) * y(idx_86,5) + rc(idx_86,133) * y(idx_86,48) * y(idx_86,16) * 1.0 + rc(idx_86,135) * y(idx_86,31) * y(idx_86,3) + rc(idx_86,129) * y(idx_86,4) * y(idx_86,34) * 0.24 + rc(idx_86,132) * y(idx_86,46) * y(idx_86,16) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_87 = 1, n_pnts, 1
          p(idx_87) = p(idx_87) + rc(idx_87,127) * y(idx_87,16) * y(idx_87,40) * 2.00 + rc(idx_87,128) * y(idx_87,4) * y(idx_87,29) * 0.22 + rc(idx_87,123) * y(idx_87,4) * y(idx_87,50) + rc(idx_87,126) * y(idx_87,40) * y(idx_87,5) + rc(idx_87,112) * y(idx_87,4) * y(idx_87,49) + rc(idx_87,112) * y(idx_87,4) * y(idx_87,49) * 0.47
        enddo
        !$acc loop gang vector independent
        do idx_88 = 1, n_pnts, 1
          p(idx_88) = p(idx_88) + rc(idx_88,110) * y(idx_88,38) * y(idx_88,5) * 2.00 + rc(idx_88,111) * y(idx_88,16) * y(idx_88,38) * 3.00 + rc(idx_88,98) * y(idx_88,3) * y(idx_88,23) + rc(idx_88,106) * y(idx_88,16) * y(idx_88,37) + rc(idx_88,96) * y(idx_88,36) * y(idx_88,16) * 2.00 + rc(idx_88,97) * y(idx_88,30) * y(idx_88,16)
        enddo
        !$acc loop gang vector independent
        do idx_89 = 1, n_pnts, 1
          p(idx_89) = p(idx_89) + rc(idx_89,84) * y(idx_89,24) * y(idx_89,16) + rc(idx_89,95) * y(idx_89,36) * y(idx_89,5) + rc(idx_89,74) * y(idx_89,16) * y(idx_89,20) * 2.00 + rc(idx_89,80) * y(idx_89,16) * y(idx_89,20) + rc(idx_89,63) * y(idx_89,41) * y(idx_89,3)
        enddo
        !$acc loop gang vector independent
        do idx_90 = 1, n_pnts, 1
          p(idx_90) = p(idx_90) + rc(idx_90,73) * y(idx_90,18) * y(idx_90,16) + rc(idx_90,61) * y(idx_90,16) * y(idx_90,16) * 2.00 + rc(idx_90,62) * y(idx_90,16) * y(idx_90,16) + rc(idx_90,45) * y(idx_90,3) * y(idx_90,17) + rc(idx_90,60) * y(idx_90,5) * y(idx_90,16)
        enddo
        !$acc loop gang vector independent
        do idx_91 = 1, n_pnts, 1
          l(idx_91) = 0.0 + dj(idx_91,7) + dw(idx_91,14) + rc(idx_91,66) * y(idx_91,3) + rc(idx_91,67) * y(idx_91,6) + dj(idx_91,6)
        enddo
        !$acc loop gang vector independent
        do idx_92 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_92) = 1.0 + dts * l(idx_92)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_93 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_93,14) = (yp(idx_93,14) + dts * p(idx_93)) * tmp1(idx_93)
        enddo
        !$acc loop gang vector independent
        do idx_94 = 1, n_pnts, 1
          p(idx_94) = 0.0 + dj(idx_94,22) * y(idx_94,31) + dj(idx_94,23) * y(idx_94,35) + dj(idx_94,20) * y(idx_94,26) + dj(idx_94,21) * y(idx_94,52) + dj(idx_94,17) * y(idx_94,17) + dj(idx_94,19) * y(idx_94,21)
        enddo
        !$acc loop gang vector independent
        do idx_95 = 1, n_pnts, 1
          p(idx_95) = p(idx_95) + dj(idx_95,13) * y(idx_95,57) + dj(idx_95,13) * y(idx_95,57) + dj(idx_95,11) * y(idx_95,9) + dj(idx_95,12) * y(idx_95,33) + dj(idx_95,6) * y(idx_95,14) + dj(idx_95,8) * y(idx_95,22)
        enddo
        !$acc loop gang vector independent
        do idx_96 = 1, n_pnts, 1
          p(idx_96) = p(idx_96) + rc(idx_96,181) * y(idx_96,44) * y(idx_96,16) * 2.00 + dj(idx_96,6) * y(idx_96,14) + rc(idx_96,173) * y(idx_96,44) * y(idx_96,5) + rc(idx_96,174) * y(idx_96,3) * y(idx_96,54) + rc(idx_96,159) * y(idx_96,6) * y(idx_96,29) + rc(idx_96,170) * y(idx_96,3) * y(idx_96,58)
        enddo
        !$acc loop gang vector independent
        do idx_97 = 1, n_pnts, 1
          p(idx_97) = p(idx_97) + rc(idx_97,153) * y(idx_97,6) * y(idx_97,49) + rc(idx_97,155) * y(idx_97,6) * y(idx_97,50) + rc(idx_97,142) * y(idx_97,46) * y(idx_97,5) + rc(idx_97,144) * y(idx_97,48) * y(idx_97,5) + rc(idx_97,133) * y(idx_97,48) * y(idx_97,16) * 2.0 + rc(idx_97,139) * y(idx_97,57) * y(idx_97,3) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_98 = 1, n_pnts, 1
          p(idx_98) = p(idx_98) + rc(idx_98,129) * y(idx_98,4) * y(idx_98,34) * 0.36 + rc(idx_98,132) * y(idx_98,46) * y(idx_98,16) * 2.0 + rc(idx_98,127) * y(idx_98,16) * y(idx_98,40) * 2.00 + rc(idx_98,128) * y(idx_98,4) * y(idx_98,29) * 0.27 + rc(idx_98,124) * y(idx_98,4) * y(idx_98,50) * 0.18 + rc(idx_98,126) * y(idx_98,40) * y(idx_98,5)
        enddo
        !$acc loop gang vector independent
        do idx_99 = 1, n_pnts, 1
          p(idx_99) = p(idx_99) + rc(idx_99,112) * y(idx_99,4) * y(idx_99,49) * 0.20 + rc(idx_99,123) * y(idx_99,4) * y(idx_99,50) * 0.30 + rc(idx_99,110) * y(idx_99,38) * y(idx_99,5) + rc(idx_99,111) * y(idx_99,16) * y(idx_99,38) * 2.00 + rc(idx_99,97) * y(idx_99,30) * y(idx_99,16) * 2.00 + rc(idx_99,106) * y(idx_99,16) * y(idx_99,37)
        enddo
        !$acc loop gang vector independent
        do idx_100 = 1, n_pnts, 1
          p(idx_100) = p(idx_100) + rc(idx_100,93) * y(idx_100,30) * y(idx_100,5) + rc(idx_100,96) * y(idx_100,36) * y(idx_100,16) + rc(idx_100,84) * y(idx_100,24) * y(idx_100,16) * 2.00 + rc(idx_100,90) * y(idx_100,18) * y(idx_100,18) * 2.00 + rc(idx_100,80) * y(idx_100,16) * y(idx_100,20) + rc(idx_100,83) * y(idx_100,24) * y(idx_100,5)
        enddo
        !$acc loop gang vector independent
        do idx_101 = 1, n_pnts, 1
          p(idx_101) = p(idx_101) + rc(idx_101,72) * y(idx_101,18) * y(idx_101,5) + rc(idx_101,73) * y(idx_101,18) * y(idx_101,16) * 2.00 + rc(idx_101,67) * y(idx_101,6) * y(idx_101,14) + rc(idx_101,70) * y(idx_101,3) * y(idx_101,13) + rc(idx_101,63) * y(idx_101,41) * y(idx_101,3) + rc(idx_101,66) * y(idx_101,3) * y(idx_101,14)
        enddo
        !$acc loop gang vector independent
        do idx_102 = 1, n_pnts, 1
          p(idx_102) = p(idx_102) + rc(idx_102,60) * y(idx_102,5) * y(idx_102,16) + rc(idx_102,61) * y(idx_102,16) * y(idx_102,16) * 2.00 + rc(idx_102,31) * y(idx_102,3) * y(idx_102,11) + rc(idx_102,33) * y(idx_102,3) * y(idx_102,43) + rc(idx_102,13) * y(idx_102,3) * y(idx_102,4) + rc(idx_102,23) * y(idx_102,9)
        enddo
        !$acc loop gang vector independent
        do idx_103 = 1, n_pnts, 1
          l(idx_103) = 0.0 + rc(idx_103,180) * y(idx_103,42) + rc(idx_103,182) * y(idx_103,55) + dw(idx_103,15) + rc(idx_103,131) * y(idx_103,20) + rc(idx_103,134) * y(idx_103,46) + rc(idx_103,136) * y(idx_103,48)
        enddo
        !$acc loop gang vector independent
        do idx_104 = 1, n_pnts, 1
          l(idx_104) = l(idx_104) + rc(idx_104,99) * y(idx_104,18) + rc(idx_104,101) * y(idx_104,30) + rc(idx_104,103) * y(idx_104,24) + rc(idx_104,36) * y(idx_104,15) + rc(idx_104,36) * y(idx_104,15) + rc(idx_104,65) * y(idx_104,16)
        enddo
        !$acc loop gang vector independent
        do idx_105 = 1, n_pnts, 1
          l(idx_105) = l(idx_105) + rc(idx_105,30) * y(idx_105,3) + rc(idx_105,32) * y(idx_105,6) + rc(idx_105,34) * y(idx_105,6) + rc(idx_105,14) * y(idx_105,4) + rc(idx_105,17) * y(idx_105,5) + rc(idx_105,22) * y(idx_105,7)
        enddo
        !$acc loop gang vector independent
        do idx_106 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_106) = 1.0 + dts * l(idx_106)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_107 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_107,15) = (yp(idx_107,15) + dts * p(idx_107)) * tmp1(idx_107)
        enddo
        !$acc loop gang vector independent
        do idx_108 = 1, n_pnts, 1
          p(idx_108) = 0.0 + dj(idx_108,10) * y(idx_108,27) + rc(idx_108,131) * y(idx_108,20) * y(idx_108,15) * 0.8 + dj(idx_108,8) * y(idx_108,22) + rc(idx_108,91) * y(idx_108,20) * y(idx_108,20) * 2.00 + rc(idx_108,123) * y(idx_108,4) * y(idx_108,50) * 0.58
        enddo
        !$acc loop gang vector independent
        do idx_109 = 1, n_pnts, 1
          p(idx_109) = p(idx_109) + rc(idx_109,79) * y(idx_109,20) * y(idx_109,5) + rc(idx_109,80) * y(idx_109,16) * y(idx_109,20) + rc(idx_109,44) * y(idx_109,3) * y(idx_109,17) + rc(idx_109,59) * y(idx_109,3) * y(idx_109,12)
        enddo
        !$acc loop gang vector independent
        do idx_110 = 1, n_pnts, 1
          l(idx_110) = 0.0 + dw(idx_110,16) + rc(idx_110,132) * y(idx_110,46) + rc(idx_110,133) * y(idx_110,48) + rc(idx_110,181) * y(idx_110,44)
        enddo
        !$acc loop gang vector independent
        do idx_111 = 1, n_pnts, 1
          l(idx_111) = l(idx_111) + rc(idx_111,106) * y(idx_111,37) + rc(idx_111,111) * y(idx_111,38) + rc(idx_111,127) * y(idx_111,40) + rc(idx_111,84) * y(idx_111,24) + rc(idx_111,96) * y(idx_111,36) + rc(idx_111,97) * y(idx_111,30)
        enddo
        !$acc loop gang vector independent
        do idx_112 = 1, n_pnts, 1
          l(idx_112) = l(idx_112) + rc(idx_112,73) * y(idx_112,18) + rc(idx_112,74) * y(idx_112,20) + rc(idx_112,80) * y(idx_112,20) + rc(idx_112,62) * y(idx_112,16)
        enddo
        !$acc loop gang vector independent
        do idx_113 = 1, n_pnts, 1
          l(idx_113) = l(idx_113) + rc(idx_113,62) * y(idx_113,16) + rc(idx_113,65) * y(idx_113,15) + rc(idx_113,60) * y(idx_113,5) + rc(idx_113,61) * y(idx_113,16) + rc(idx_113,61) * y(idx_113,16)
        enddo
        !$acc loop gang vector independent
        do idx_114 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_114) = 1.0 + dts * l(idx_114)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_115 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_115,16) = (yp(idx_115,16) + dts * p(idx_115)) * tmp1(idx_115)
        enddo
        !$acc loop gang vector independent
        do idx_116 = 1, n_pnts, 1
          p(idx_116) = 0.0 + rc(idx_116,65) * y(idx_116,16) * y(idx_116,15)
        enddo
        !$acc loop gang vector independent
        do idx_117 = 1, n_pnts, 1
          l(idx_117) = 0.0 + dw(idx_117,17) + dd(idx_117,17) + rc(idx_117,44) * y(idx_117,3) + rc(idx_117,45) * y(idx_117,3) + dj(idx_117,17)
        enddo
        !$acc loop gang vector independent
        do idx_118 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_118) = 1.0 + dts * l(idx_118)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_119 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_119,17) = (yp(idx_119,17) + dts * p(idx_119)) * tmp1(idx_119)
        enddo
        !$acc loop gang vector independent
        do idx_120 = 1, n_pnts, 1
          p(idx_120) = 0.0 + rc(idx_120,151) * y(idx_120,6) * y(idx_120,19) + dj(idx_120,9) * y(idx_120,53) + rc(idx_120,44) * y(idx_120,3) * y(idx_120,21) + rc(idx_120,71) * y(idx_120,3) * y(idx_120,19)
        enddo
        !$acc loop gang vector independent
        do idx_121 = 1, n_pnts, 1
          l(idx_121) = 0.0 + rc(idx_121,90) * y(idx_121,18) + rc(idx_121,99) * y(idx_121,15) + rc(idx_121,72) * y(idx_121,5) + rc(idx_121,73) * y(idx_121,16) + rc(idx_121,90) * y(idx_121,18)
        enddo
        !$acc loop gang vector independent
        do idx_122 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_122) = 1.0 + dts * l(idx_122)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_123 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_123,18) = (yp(idx_123,18) + dts * p(idx_123)) * tmp1(idx_123)
        enddo
        !$acc loop gang vector independent
        do idx_124 = 1, n_pnts, 1
          p(idx_124) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_125 = 1, n_pnts, 1
          l(idx_125) = 0.0 + rc(idx_125,71) * y(idx_125,3) + rc(idx_125,151) * y(idx_125,6)
        enddo
        !$acc loop gang vector independent
        do idx_126 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_126) = 1.0 + dts * l(idx_126)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_127 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_127,19) = (yp(idx_127,19) + dts * p(idx_127)) * tmp1(idx_127)
        enddo
        !$acc loop gang vector independent
        do idx_128 = 1, n_pnts, 1
          p(idx_128) = 0.0 + dj(idx_128,18) * y(idx_128,23) + dj(idx_128,10) * y(idx_128,27) + dj(idx_128,12) * y(idx_128,33) + rc(idx_128,158) * y(idx_128,6) * y(idx_128,22) + dj(idx_128,9) * y(idx_128,53) + rc(idx_128,131) * y(idx_128,20) * y(idx_128,15) * 0.2
        enddo
        !$acc loop gang vector independent
        do idx_129 = 1, n_pnts, 1
          p(idx_129) = p(idx_129) + rc(idx_129,138) * y(idx_129,33) * y(idx_129,3) + rc(idx_129,105) * y(idx_129,5) * y(idx_129,37) + rc(idx_129,106) * y(idx_129,16) * y(idx_129,37)
        enddo
        !$acc loop gang vector independent
        do idx_130 = 1, n_pnts, 1
          p(idx_130) = p(idx_130) + rc(idx_130,95) * y(idx_130,36) * y(idx_130,5) + rc(idx_130,96) * y(idx_130,36) * y(idx_130,16) + rc(idx_130,75) * y(idx_130,3) * y(idx_130,22) + rc(idx_130,78) * y(idx_130,23)
        enddo
        !$acc loop gang vector independent
        do idx_131 = 1, n_pnts, 1
          l(idx_131) = 0.0 + rc(idx_131,131) * y(idx_131,15) + rc(idx_131,80) * y(idx_131,16) + rc(idx_131,91) * y(idx_131,20) + rc(idx_131,91) * y(idx_131,20) + rc(idx_131,74) * y(idx_131,16) + rc(idx_131,77) * y(idx_131,7) + rc(idx_131,79) * y(idx_131,5)
        enddo
        !$acc loop gang vector independent
        do idx_132 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_132) = 1.0 + dts * l(idx_132)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_133 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_133,20) = (yp(idx_133,20) + dts * p(idx_133)) * tmp1(idx_133)
        enddo
        !$acc loop gang vector independent
        do idx_134 = 1, n_pnts, 1
          p(idx_134) = 0.0 + rc(idx_134,99) * y(idx_134,15) * y(idx_134,18)
        enddo
        !$acc loop gang vector independent
        do idx_135 = 1, n_pnts, 1
          l(idx_135) = 0.0 + dw(idx_135,21) + dd(idx_135,21) + rc(idx_135,44) * y(idx_135,3) + rc(idx_135,100) * y(idx_135,3) + dj(idx_135,19)
        enddo
        !$acc loop gang vector independent
        do idx_136 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_136) = 1.0 + dts * l(idx_136)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_137 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_137,21) = (yp(idx_137,21) + dts * p(idx_137)) * tmp1(idx_137)
        enddo
        !$acc loop gang vector independent
        do idx_138 = 1, n_pnts, 1
          p(idx_138) = 0.0 + dj(idx_138,19) * y(idx_138,21) + rc(idx_138,127) * y(idx_138,16) * y(idx_138,40) + rc(idx_138,156) * y(idx_138,47) * y(idx_138,3) + rc(idx_138,124) * y(idx_138,4) * y(idx_138,50) + rc(idx_138,126) * y(idx_138,40) * y(idx_138,5)
        enddo
        !$acc loop gang vector independent
        do idx_139 = 1, n_pnts, 1
          p(idx_139) = p(idx_139) + rc(idx_139,105) * y(idx_139,5) * y(idx_139,37) + rc(idx_139,106) * y(idx_139,16) * y(idx_139,37) + rc(idx_139,90) * y(idx_139,18) * y(idx_139,18) * 2.00 + rc(idx_139,100) * y(idx_139,3) * y(idx_139,21) + rc(idx_139,72) * y(idx_139,18) * y(idx_139,5) + rc(idx_139,73) * y(idx_139,18) * y(idx_139,16)
        enddo
        !$acc loop gang vector independent
        do idx_140 = 1, n_pnts, 1
          l(idx_140) = 0.0 + rc(idx_140,75) * y(idx_140,3) + rc(idx_140,158) * y(idx_140,6) + dj(idx_140,8)
        enddo
        !$acc loop gang vector independent
        do idx_141 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_141) = 1.0 + dts * l(idx_141)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_142 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_142,22) = (yp(idx_142,22) + dts * p(idx_142)) * tmp1(idx_142)
        enddo
        !$acc loop gang vector independent
        do idx_143 = 1, n_pnts, 1
          p(idx_143) = 0.0 + rc(idx_143,77) * y(idx_143,20) * y(idx_143,7)
        enddo
        !$acc loop gang vector independent
        do idx_144 = 1, n_pnts, 1
          l(idx_144) = 0.0 + dd(idx_144,23) + rc(idx_144,78) + rc(idx_144,98) * y(idx_144,3) + dj(idx_144,18)
        enddo
        !$acc loop gang vector independent
        do idx_145 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_145) = 1.0 + dts * l(idx_145)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_146 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_146,23) = (yp(idx_146,23) + dts * p(idx_146)) * tmp1(idx_146)
        enddo
        !$acc loop gang vector independent
        do idx_147 = 1, n_pnts, 1
          p(idx_147) = 0.0 + rc(idx_147,152) * y(idx_147,6) * y(idx_147,51) + rc(idx_147,44) * y(idx_147,3) * y(idx_147,52) + rc(idx_147,81) * y(idx_147,3) * y(idx_147,51)
        enddo
        !$acc loop gang vector independent
        do idx_148 = 1, n_pnts, 1
          l(idx_148) = 0.0 + rc(idx_148,83) * y(idx_148,5) + rc(idx_148,84) * y(idx_148,16) + rc(idx_148,103) * y(idx_148,15)
        enddo
        !$acc loop gang vector independent
        do idx_149 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_149) = 1.0 + dts * l(idx_149)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_150 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_150,24) = (yp(idx_150,24) + dts * p(idx_150)) * tmp1(idx_150)
        enddo
        !$acc loop gang vector independent
        do idx_151 = 1, n_pnts, 1
          p(idx_151) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_152 = 1, n_pnts, 1
          l(idx_152) = 0.0 + rc(idx_152,92) * y(idx_152,3)
        enddo
        !$acc loop gang vector independent
        do idx_153 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_153) = 1.0 + dts * l(idx_153)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_154 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_154,25) = (yp(idx_154,25) + dts * p(idx_154)) * tmp1(idx_154)
        enddo
        !$acc loop gang vector independent
        do idx_155 = 1, n_pnts, 1
          p(idx_155) = 0.0 + rc(idx_155,101) * y(idx_155,15) * y(idx_155,30)
        enddo
        !$acc loop gang vector independent
        do idx_156 = 1, n_pnts, 1
          l(idx_156) = 0.0 + dw(idx_156,26) + dd(idx_156,26) + rc(idx_156,44) * y(idx_156,3) + rc(idx_156,102) * y(idx_156,3) + dj(idx_156,20)
        enddo
        !$acc loop gang vector independent
        do idx_157 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_157) = 1.0 + dts * l(idx_157)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_158 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_158,26) = (yp(idx_158,26) + dts * p(idx_158)) * tmp1(idx_158)
        enddo
        !$acc loop gang vector independent
        do idx_159 = 1, n_pnts, 1
          p(idx_159) = 0.0 + rc(idx_159,102) * y(idx_159,3) * y(idx_159,26) + dj(idx_159,20) * y(idx_159,26) + rc(idx_159,93) * y(idx_159,30) * y(idx_159,5) + rc(idx_159,97) * y(idx_159,30) * y(idx_159,16)
        enddo
        !$acc loop gang vector independent
        do idx_160 = 1, n_pnts, 1
          l(idx_160) = 0.0 + rc(idx_160,94) * y(idx_160,3) + dj(idx_160,10)
        enddo
        !$acc loop gang vector independent
        do idx_161 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_161) = 1.0 + dts * l(idx_161)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_162 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_162,27) = (yp(idx_162,27) + dts * p(idx_162)) * tmp1(idx_162)
        enddo
        !$acc loop gang vector independent
        do idx_163 = 1, n_pnts, 1
          p(idx_163) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_164 = 1, n_pnts, 1
          l(idx_164) = 0.0 + dj(idx_164,1) + dj(idx_164,2) + dd(idx_164,28) + rc(idx_164,124) * y(idx_164,50) + rc(idx_164,128) * y(idx_164,29) + rc(idx_164,129) * y(idx_164,34)
        enddo
        !$acc loop gang vector independent
        do idx_165 = 1, n_pnts, 1
          l(idx_165) = l(idx_165) + rc(idx_165,14) * y(idx_165,15) + rc(idx_165,112) * y(idx_165,49) + rc(idx_165,123) * y(idx_165,50) + rc(idx_165,11) * y(idx_165,5) + rc(idx_165,12) * y(idx_165,7) + rc(idx_165,13) * y(idx_165,3)
        enddo
        !$acc loop gang vector independent
        do idx_166 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_166) = 1.0 + dts * l(idx_166)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_167 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_167,28) = (yp(idx_167,28) + dts * p(idx_167)) * tmp1(idx_167)
        enddo
        !$acc loop gang vector independent
        do idx_168 = 1, n_pnts, 1
          p(idx_168) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_169 = 1, n_pnts, 1
          l(idx_169) = 0.0 + rc(idx_169,128) * y(idx_169,4) + rc(idx_169,141) * y(idx_169,3) + rc(idx_169,159) * y(idx_169,6)
        enddo
        !$acc loop gang vector independent
        do idx_170 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_170) = 1.0 + dts * l(idx_170)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_171 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_171,29) = (yp(idx_171,29) + dts * p(idx_171)) * tmp1(idx_171)
        enddo
        !$acc loop gang vector independent
        do idx_172 = 1, n_pnts, 1
          p(idx_172) = 0.0 + rc(idx_172,44) * y(idx_172,3) * y(idx_172,26) + rc(idx_172,92) * y(idx_172,3) * y(idx_172,25)
        enddo
        !$acc loop gang vector independent
        do idx_173 = 1, n_pnts, 1
          l(idx_173) = 0.0 + rc(idx_173,93) * y(idx_173,5) + rc(idx_173,97) * y(idx_173,16) + rc(idx_173,101) * y(idx_173,15)
        enddo
        !$acc loop gang vector independent
        do idx_174 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_174) = 1.0 + dts * l(idx_174)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_175 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_175,30) = (yp(idx_175,30) + dts * p(idx_175)) * tmp1(idx_175)
        enddo
        !$acc loop gang vector independent
        do idx_176 = 1, n_pnts, 1
          p(idx_176) = 0.0 + rc(idx_176,134) * y(idx_176,46) * y(idx_176,15)
        enddo
        !$acc loop gang vector independent
        do idx_177 = 1, n_pnts, 1
          l(idx_177) = 0.0 + dd(idx_177,31) + rc(idx_177,135) * y(idx_177,3) + dj(idx_177,22) + dw(idx_177,31)
        enddo
        !$acc loop gang vector independent
        do idx_178 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_178) = 1.0 + dts * l(idx_178)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_179 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_179,31) = (yp(idx_179,31) + dts * p(idx_179)) * tmp1(idx_179)
        enddo
        !$acc loop gang vector independent
        do idx_180 = 1, n_pnts, 1
          p(idx_180) = 0.0 + rc(idx_180,159) * y(idx_180,6) * y(idx_180,29)
        enddo
        !$acc loop gang vector independent
        do idx_181 = 1, n_pnts, 1
          l(idx_181) = 0.0 + rc(idx_181,160) * y(idx_181,3) + dw(idx_181,32)
        enddo
        !$acc loop gang vector independent
        do idx_182 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_182) = 1.0 + dts * l(idx_182)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_183 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_183,32) = (yp(idx_183,32) + dts * p(idx_183)) * tmp1(idx_183)
        enddo
        !$acc loop gang vector independent
        do idx_184 = 1, n_pnts, 1
          p(idx_184) = 0.0 + dj(idx_184,23) * y(idx_184,35) + rc(idx_184,180) * y(idx_184,42) * y(idx_184,15) + rc(idx_184,181) * y(idx_184,44) * y(idx_184,16) + rc(idx_184,170) * y(idx_184,3) * y(idx_184,58) * 0.80 + rc(idx_184,173) * y(idx_184,44) * y(idx_184,5)
        enddo
        !$acc loop gang vector independent
        do idx_185 = 1, n_pnts, 1
          p(idx_185) = p(idx_185) + rc(idx_185,137) * y(idx_185,35) * y(idx_185,3) + rc(idx_185,144) * y(idx_185,48) * y(idx_185,5) + rc(idx_185,129) * y(idx_185,4) * y(idx_185,34) + rc(idx_185,133) * y(idx_185,48) * y(idx_185,16) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_186 = 1, n_pnts, 1
          l(idx_186) = 0.0 + rc(idx_186,138) * y(idx_186,3) + dj(idx_186,12) + dw(idx_186,33)
        enddo
        !$acc loop gang vector independent
        do idx_187 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_187) = 1.0 + dts * l(idx_187)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_188 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_188,33) = (yp(idx_188,33) + dts * p(idx_188)) * tmp1(idx_188)
        enddo
        !$acc loop gang vector independent
        do idx_189 = 1, n_pnts, 1
          p(idx_189) = 0.0 + rc(idx_189,160) * y(idx_189,32) * y(idx_189,3) + dj(idx_189,22) * y(idx_189,31) + rc(idx_189,135) * y(idx_189,31) * y(idx_189,3) + rc(idx_189,142) * y(idx_189,46) * y(idx_189,5) + rc(idx_189,128) * y(idx_189,4) * y(idx_189,29) + rc(idx_189,132) * y(idx_189,46) * y(idx_189,16) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_190 = 1, n_pnts, 1
          l(idx_190) = 0.0 + rc(idx_190,129) * y(idx_190,4) + rc(idx_190,143) * y(idx_190,3)
        enddo
        !$acc loop gang vector independent
        do idx_191 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_191) = 1.0 + dts * l(idx_191)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_192 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_192,34) = (yp(idx_192,34) + dts * p(idx_192)) * tmp1(idx_192)
        enddo
        !$acc loop gang vector independent
        do idx_193 = 1, n_pnts, 1
          p(idx_193) = 0.0 + rc(idx_193,136) * y(idx_193,48) * y(idx_193,15)
        enddo
        !$acc loop gang vector independent
        do idx_194 = 1, n_pnts, 1
          l(idx_194) = 0.0 + dd(idx_194,35) + rc(idx_194,137) * y(idx_194,3) + dj(idx_194,23) + dw(idx_194,35)
        enddo
        !$acc loop gang vector independent
        do idx_195 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_195) = 1.0 + dts * l(idx_195)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_196 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_196,35) = (yp(idx_196,35) + dts * p(idx_196)) * tmp1(idx_196)
        enddo
        !$acc loop gang vector independent
        do idx_197 = 1, n_pnts, 1
          p(idx_197) = 0.0 + rc(idx_197,94) * y(idx_197,3) * y(idx_197,27)
        enddo
        !$acc loop gang vector independent
        do idx_198 = 1, n_pnts, 1
          l(idx_198) = 0.0 + rc(idx_198,95) * y(idx_198,5) + rc(idx_198,96) * y(idx_198,16)
        enddo
        !$acc loop gang vector independent
        do idx_199 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_199) = 1.0 + dts * l(idx_199)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_200 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_200,36) = (yp(idx_200,36) + dts * p(idx_200)) * tmp1(idx_200)
        enddo
        !$acc loop gang vector independent
        do idx_201 = 1, n_pnts, 1
          p(idx_201) = 0.0 + rc(idx_201,86) * y(idx_201,3) * y(idx_201,53)
        enddo
        !$acc loop gang vector independent
        do idx_202 = 1, n_pnts, 1
          l(idx_202) = 0.0 + rc(idx_202,105) * y(idx_202,5) + rc(idx_202,106) * y(idx_202,16)
        enddo
        !$acc loop gang vector independent
        do idx_203 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_203) = 1.0 + dts * l(idx_203)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_204 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_204,37) = (yp(idx_204,37) + dts * p(idx_204)) * tmp1(idx_204)
        enddo
        !$acc loop gang vector independent
        do idx_205 = 1, n_pnts, 1
          p(idx_205) = 0.0 + rc(idx_205,109) * y(idx_205,3) * y(idx_205,49)
        enddo
        !$acc loop gang vector independent
        do idx_206 = 1, n_pnts, 1
          l(idx_206) = 0.0 + rc(idx_206,110) * y(idx_206,5) + rc(idx_206,111) * y(idx_206,16)
        enddo
        !$acc loop gang vector independent
        do idx_207 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_207) = 1.0 + dts * l(idx_207)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_208 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_208,38) = (yp(idx_208,38) + dts * p(idx_208)) * tmp1(idx_208)
        enddo
        !$acc loop gang vector independent
        do idx_209 = 1, n_pnts, 1
          p(idx_209) = 0.0 + rc(idx_209,171) * y(idx_209,7) * y(idx_209,42) + rc(idx_209,183) * y(idx_209,7) * y(idx_209,55)
        enddo
        !$acc loop gang vector independent
        do idx_210 = 1, n_pnts, 1
          l(idx_210) = 0.0 + dd(idx_210,39) + rc(idx_210,177) * y(idx_210,3) + rc(idx_210,178) * y(idx_210,6) + dw(idx_210,39)
        enddo
        !$acc loop gang vector independent
        do idx_211 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_211) = 1.0 + dts * l(idx_211)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_212 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_212,39) = (yp(idx_212,39) + dts * p(idx_212)) * tmp1(idx_212)
        enddo
        !$acc loop gang vector independent
        do idx_213 = 1, n_pnts, 1
          p(idx_213) = 0.0 + rc(idx_213,125) * y(idx_213,3) * y(idx_213,50)
        enddo
        !$acc loop gang vector independent
        do idx_214 = 1, n_pnts, 1
          l(idx_214) = 0.0 + rc(idx_214,126) * y(idx_214,5) + rc(idx_214,127) * y(idx_214,16)
        enddo
        !$acc loop gang vector independent
        do idx_215 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_215) = 1.0 + dts * l(idx_215)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_216 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_216,40) = (yp(idx_216,40) + dts * p(idx_216)) * tmp1(idx_216)
        enddo
        !$acc loop gang vector independent
        do idx_217 = 1, n_pnts, 1
          p(idx_217) = 0.0 + rc(idx_217,62) * y(idx_217,16) * y(idx_217,16) + rc(idx_217,123) * y(idx_217,4) * y(idx_217,50) * 0.12
        enddo
        !$acc loop gang vector independent
        do idx_218 = 1, n_pnts, 1
          l(idx_218) = 0.0 + rc(idx_218,63) * y(idx_218,3) + dw(idx_218,41)
        enddo
        !$acc loop gang vector independent
        do idx_219 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_219) = 1.0 + dts * l(idx_219)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_220 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_220,41) = (yp(idx_220,41) + dts * p(idx_220)) * tmp1(idx_220)
        enddo
        !$acc loop gang vector independent
        do idx_221 = 1, n_pnts, 1
          p(idx_221) = 0.0 + rc(idx_221,176) * y(idx_221,3) * y(idx_221,58)
        enddo
        !$acc loop gang vector independent
        do idx_222 = 1, n_pnts, 1
          l(idx_222) = 0.0 + rc(idx_222,171) * y(idx_222,7) + rc(idx_222,180) * y(idx_222,15)
        enddo
        !$acc loop gang vector independent
        do idx_223 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_223) = 1.0 + dts * l(idx_223)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_224 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_224,42) = (yp(idx_224,42) + dts * p(idx_224)) * tmp1(idx_224)
        enddo
        !$acc loop gang vector independent
        do idx_225 = 1, n_pnts, 1
          p(idx_225) = 0.0 + dj(idx_225,7) * y(idx_225,14) + rc(idx_225,112) * y(idx_225,4) * y(idx_225,49) * 0.13 + rc(idx_225,124) * y(idx_225,4) * y(idx_225,50) * 0.24
        enddo
        !$acc loop gang vector independent
        do idx_226 = 1, n_pnts, 1
          l(idx_226) = 0.0 + rc(idx_226,33) * y(idx_226,3) + dd(idx_226,43)
        enddo
        !$acc loop gang vector independent
        do idx_227 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_227) = 1.0 + dts * l(idx_227)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_228 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_228,43) = (yp(idx_228,43) + dts * p(idx_228)) * tmp1(idx_228)
        enddo
        !$acc loop gang vector independent
        do idx_229 = 1, n_pnts, 1
          p(idx_229) = 0.0 + rc(idx_229,172) * y(idx_229,3) * y(idx_229,56)
        enddo
        !$acc loop gang vector independent
        do idx_230 = 1, n_pnts, 1
          l(idx_230) = 0.0 + rc(idx_230,173) * y(idx_230,5) + rc(idx_230,181) * y(idx_230,16)
        enddo
        !$acc loop gang vector independent
        do idx_231 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_231) = 1.0 + dts * l(idx_231)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_232 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_232,44) = (yp(idx_232,44) + dts * p(idx_232)) * tmp1(idx_232)
        enddo
        !$acc loop gang vector independent
        do idx_233 = 1, n_pnts, 1
          p(idx_233) = 0.0 + rc(idx_233,153) * y(idx_233,6) * y(idx_233,49)
        enddo
        !$acc loop gang vector independent
        do idx_234 = 1, n_pnts, 1
          l(idx_234) = 0.0 + rc(idx_234,154) * y(idx_234,3)
        enddo
        !$acc loop gang vector independent
        do idx_235 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_235) = 1.0 + dts * l(idx_235)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_236 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_236,45) = (yp(idx_236,45) + dts * p(idx_236)) * tmp1(idx_236)
        enddo
        !$acc loop gang vector independent
        do idx_237 = 1, n_pnts, 1
          p(idx_237) = 0.0 + rc(idx_237,141) * y(idx_237,3) * y(idx_237,29)
        enddo
        !$acc loop gang vector independent
        do idx_238 = 1, n_pnts, 1
          l(idx_238) = 0.0 + rc(idx_238,132) * y(idx_238,16) + rc(idx_238,134) * y(idx_238,15) + rc(idx_238,142) * y(idx_238,5)
        enddo
        !$acc loop gang vector independent
        do idx_239 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_239) = 1.0 + dts * l(idx_239)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_240 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_240,46) = (yp(idx_240,46) + dts * p(idx_240)) * tmp1(idx_240)
        enddo
        !$acc loop gang vector independent
        do idx_241 = 1, n_pnts, 1
          p(idx_241) = 0.0 + rc(idx_241,155) * y(idx_241,6) * y(idx_241,50)
        enddo
        !$acc loop gang vector independent
        do idx_242 = 1, n_pnts, 1
          l(idx_242) = 0.0 + rc(idx_242,156) * y(idx_242,3)
        enddo
        !$acc loop gang vector independent
        do idx_243 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_243) = 1.0 + dts * l(idx_243)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_244 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_244,47) = (yp(idx_244,47) + dts * p(idx_244)) * tmp1(idx_244)
        enddo
        !$acc loop gang vector independent
        do idx_245 = 1, n_pnts, 1
          p(idx_245) = 0.0 + rc(idx_245,143) * y(idx_245,3) * y(idx_245,34)
        enddo
        !$acc loop gang vector independent
        do idx_246 = 1, n_pnts, 1
          l(idx_246) = 0.0 + rc(idx_246,133) * y(idx_246,16) + rc(idx_246,136) * y(idx_246,15) + rc(idx_246,144) * y(idx_246,5)
        enddo
        !$acc loop gang vector independent
        do idx_247 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_247) = 1.0 + dts * l(idx_247)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_248 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_248,48) = (yp(idx_248,48) + dts * p(idx_248)) * tmp1(idx_248)
        enddo
        !$acc loop gang vector independent
        do idx_249 = 1, n_pnts, 1
          p(idx_249) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_250 = 1, n_pnts, 1
          l(idx_250) = 0.0 + rc(idx_250,109) * y(idx_250,3) + rc(idx_250,112) * y(idx_250,4) + rc(idx_250,153) * y(idx_250,6)
        enddo
        !$acc loop gang vector independent
        do idx_251 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_251) = 1.0 + dts * l(idx_251)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_252 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_252,49) = (yp(idx_252,49) + dts * p(idx_252)) * tmp1(idx_252)
        enddo
        !$acc loop gang vector independent
        do idx_253 = 1, n_pnts, 1
          p(idx_253) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_254 = 1, n_pnts, 1
          l(idx_254) = 0.0 + rc(idx_254,155) * y(idx_254,6) + rc(idx_254,123) * y(idx_254,4) + rc(idx_254,124) * y(idx_254,4) + rc(idx_254,125) * y(idx_254,3)
        enddo
        !$acc loop gang vector independent
        do idx_255 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_255) = 1.0 + dts * l(idx_255)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_256 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_256,50) = (yp(idx_256,50) + dts * p(idx_256)) * tmp1(idx_256)
        enddo
        !$acc loop gang vector independent
        do idx_257 = 1, n_pnts, 1
          p(idx_257) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_258 = 1, n_pnts, 1
          l(idx_258) = 0.0 + rc(idx_258,81) * y(idx_258,3) + rc(idx_258,152) * y(idx_258,6)
        enddo
        !$acc loop gang vector independent
        do idx_259 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_259) = 1.0 + dts * l(idx_259)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_260 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_260,51) = (yp(idx_260,51) + dts * p(idx_260)) * tmp1(idx_260)
        enddo
        !$acc loop gang vector independent
        do idx_261 = 1, n_pnts, 1
          p(idx_261) = 0.0 + rc(idx_261,103) * y(idx_261,15) * y(idx_261,24)
        enddo
        !$acc loop gang vector independent
        do idx_262 = 1, n_pnts, 1
          l(idx_262) = 0.0 + dw(idx_262,52) + dd(idx_262,52) + rc(idx_262,44) * y(idx_262,3) + rc(idx_262,104) * y(idx_262,3) + dj(idx_262,21)
        enddo
        !$acc loop gang vector independent
        do idx_263 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_263) = 1.0 + dts * l(idx_263)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_264 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_264,52) = (yp(idx_264,52) + dts * p(idx_264)) * tmp1(idx_264)
        enddo
        !$acc loop gang vector independent
        do idx_265 = 1, n_pnts, 1
          p(idx_265) = 0.0 + rc(idx_265,104) * y(idx_265,3) * y(idx_265,52) + dj(idx_265,21) * y(idx_265,52) + rc(idx_265,83) * y(idx_265,24) * y(idx_265,5) + rc(idx_265,84) * y(idx_265,24) * y(idx_265,16)
        enddo
        !$acc loop gang vector independent
        do idx_266 = 1, n_pnts, 1
          l(idx_266) = 0.0 + rc(idx_266,86) * y(idx_266,3) + dj(idx_266,9)
        enddo
        !$acc loop gang vector independent
        do idx_267 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_267) = 1.0 + dts * l(idx_267)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_268 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_268,53) = (yp(idx_268,53) + dts * p(idx_268)) * tmp1(idx_268)
        enddo
        !$acc loop gang vector independent
        do idx_269 = 1, n_pnts, 1
          p(idx_269) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_270 = 1, n_pnts, 1
          l(idx_270) = 0.0 + rc(idx_270,174) * y(idx_270,3) + rc(idx_270,175) * y(idx_270,3)
        enddo
        !$acc loop gang vector independent
        do idx_271 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_271) = 1.0 + dts * l(idx_271)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_272 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_272,54) = (yp(idx_272,54) + dts * p(idx_272)) * tmp1(idx_272)
        enddo
        !$acc loop gang vector independent
        do idx_273 = 1, n_pnts, 1
          p(idx_273) = 0.0 + rc(idx_273,175) * y(idx_273,3) * y(idx_273,54)
        enddo
        !$acc loop gang vector independent
        do idx_274 = 1, n_pnts, 1
          l(idx_274) = 0.0 + rc(idx_274,182) * y(idx_274,15) + rc(idx_274,183) * y(idx_274,7)
        enddo
        !$acc loop gang vector independent
        do idx_275 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_275) = 1.0 + dts * l(idx_275)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_276 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_276,55) = (yp(idx_276,55) + dts * p(idx_276)) * tmp1(idx_276)
        enddo
        !$acc loop gang vector independent
        do idx_277 = 1, n_pnts, 1
          p(idx_277) = 0.0 + rc(idx_277,180) * y(idx_277,42) * y(idx_277,15) + rc(idx_277,182) * y(idx_277,15) * y(idx_277,55) + rc(idx_277,177) * y(idx_277,3) * y(idx_277,39) + rc(idx_277,178) * y(idx_277,6) * y(idx_277,39) + rc(idx_277,170) * y(idx_277,3) * y(idx_277,58) * 0.80 + rc(idx_277,174) * y(idx_277,3) * y(idx_277,54)
        enddo
        !$acc loop gang vector independent
        do idx_278 = 1, n_pnts, 1
          l(idx_278) = 0.0 + rc(idx_278,172) * y(idx_278,3)
        enddo
        !$acc loop gang vector independent
        do idx_279 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_279) = 1.0 + dts * l(idx_279)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_280 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_280,56) = (yp(idx_280,56) + dts * p(idx_280)) * tmp1(idx_280)
        enddo
        !$acc loop gang vector independent
        do idx_281 = 1, n_pnts, 1
          p(idx_281) = 0.0 + rc(idx_281,181) * y(idx_281,44) * y(idx_281,16) + rc(idx_281,182) * y(idx_281,15) * y(idx_281,55) + rc(idx_281,177) * y(idx_281,3) * y(idx_281,39) + rc(idx_281,178) * y(idx_281,6) * y(idx_281,39) + rc(idx_281,173) * y(idx_281,44) * y(idx_281,5) + rc(idx_281,174) * y(idx_281,3) * y(idx_281,54)
        enddo
        !$acc loop gang vector independent
        do idx_282 = 1, n_pnts, 1
          l(idx_282) = 0.0 + rc(idx_282,139) * y(idx_282,3) + dj(idx_282,13) + dw(idx_282,57)
        enddo
        !$acc loop gang vector independent
        do idx_283 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_283) = 1.0 + dts * l(idx_283)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_284 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_284,57) = (yp(idx_284,57) + dts * p(idx_284)) * tmp1(idx_284)
        enddo
        !$acc loop gang vector independent
        do idx_285 = 1, n_pnts, 1
          p(idx_285) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_286 = 1, n_pnts, 1
          l(idx_286) = 0.0 + rc(idx_286,170) * y(idx_286,3) + rc(idx_286,176) * y(idx_286,3)
        enddo
        !$acc loop gang vector independent
        do idx_287 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_287) = 1.0 + dts * l(idx_287)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_288 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_288,58) = (yp(idx_288,58) + dts * p(idx_288)) * tmp1(idx_288)
        enddo
        !$acc end kernels
      enddo
      !$acc kernels
      icnt_dd = 0
      !$acc loop seq
      do j = 1, jpspec, 1
        if (ldepd(j)) then
          icnt_dd = icnt_dd + 1
          !$acc loop gang vector independent
          do i = 1, n_pnts, 1
            ddflux(i,icnt_dd) = ddflux(i,icnt_dd) + dd(i,j) * y(i,j) * vol(i) * recip_avogad
          enddo
        end if
      enddo
      icnt_wd = 0
      !$acc loop seq
      do j = 1, jpspec, 1
        if (ldepw(j)) then
          icnt_wd = icnt_wd + 1
          !$acc loop gang vector independent
          do i = 1, n_pnts, 1
            wdflux(i,icnt_wd) = wdflux(i,icnt_wd) + dw(i,j) * y(i,j) * vol(i) * recip_avogad
          enddo
        end if
      enddo
      !$acc end kernels
    enddo
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do j = 1, icnt_dd, 1
      do idx_289 = LBOUND(ddflux, dim=1), UBOUND(ddflux, dim=1), 1
        ddflux(idx_289,j) = ddflux(idx_289,j) / n_be_calls
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do j = 1, icnt_wd, 1
      do idx_290 = LBOUND(wdflux, dim=1), UBOUND(wdflux, dim=1), 1
        wdflux(idx_290,j) = wdflux(idx_290,j) / n_be_calls
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_deriv_raq

end module ukca_deriv_raq_mod
