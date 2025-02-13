module ukca_deriv_raqaero_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DERIV_RAQAERO_MOD'
  public

  contains
  subroutine ukca_deriv_raqaero(nr_therm, n_be_calls, n_pnts, dts, rc, dw, dd, dj, h2o, y, delta_so2_wetox_h2o2, delta_so2_wetox_o3, delta_so2_dryox_oh)
    use asad_mod, only : jppj, jpspec
    use ukca_config_specification_mod, only : ukca_config
    use yomhook, only : dr_hook, lhook
    use nvtx
    use parkind1, only : jpim, jprb
    use ukca_um_legacy_mod, only : oneover_v
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: nr_therm
    integer, intent(in) :: n_be_calls
    integer, intent(in) :: n_pnts
    real, intent(in) :: dts
    real, dimension(n_pnts,nr_therm), intent(in) :: rc
    real, dimension(n_pnts,jpspec), intent(in) :: dw
    real, dimension(n_pnts,jpspec), intent(in) :: dd
    real, dimension(n_pnts,jppj), intent(in) :: dj
    real, dimension(n_pnts), intent(in) :: h2o
    real, dimension(n_pnts,jpspec), intent(inout) :: y
    real, dimension(n_pnts), intent(out) :: delta_so2_wetox_h2o2
    real, dimension(n_pnts), intent(out) :: delta_so2_wetox_o3
    real, dimension(n_pnts), intent(out) :: delta_so2_dryox_oh
    integer :: i
    integer :: j
    integer :: n
    real, dimension(n_pnts) :: p
    real, dimension(n_pnts) :: l
    real, dimension(n_pnts,jpspec) :: yp
    real, dimension(n_pnts) :: tmp1
    real, dimension(n_pnts) :: tmp2
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DERIV_RAQAERO'
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
    integer :: idx_291
    integer :: idx_292
    integer :: idx_293
    integer :: idx_294
    integer :: idx_295
    integer :: idx_296
    integer :: idx_297
    integer :: idx_298
    integer :: idx_299
    integer :: idx_300
    integer :: idx_301
    integer :: idx_302
    integer :: idx_303
    integer :: idx_304
    integer :: idx_305
    integer :: idx_306
    integer :: idx_307

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
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
    !$acc loop gang vector independent
    do idx_3 = LBOUND(tmp1, dim=1), UBOUND(tmp1, dim=1), 1
      tmp1(idx_3) = 0.0
    enddo
    !$acc loop gang vector independent
    do idx_4 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
      tmp2(idx_4) = 0.0
    enddo
    !$acc end kernels
    do n = 1, n_be_calls, 1
      !$acc kernels
      !$acc loop gang vector independent collapse(2)
      do j = 1, jpspec, 1
        do idx_5 = LBOUND(yp, dim=1), UBOUND(yp, dim=1), 1
          yp(idx_5,j) = y(idx_5,j)
        enddo
      enddo
      !$acc end kernels
      do i = 1, ukca_config%nit, 1
        !$acc kernels
        !$acc loop gang vector independent
        do idx_6 = 1, n_pnts, 1
          p(idx_6) = 0.0 + dj(idx_6,3) * y(idx_6,7) + dj(idx_6,15) * y(idx_6,6) + rc(idx_6,7) * y(idx_6,2) + dj(idx_6,1) * y(idx_6,4)
        enddo
        !$acc loop gang vector independent
        do idx_7 = 1, n_pnts, 1
          l(idx_7) = 0.0 + rc(idx_7,1) + rc(idx_7,5) * y(idx_7,5) + rc(idx_7,16) * y(idx_7,7)
        enddo
        !$acc loop gang vector independent
        do idx_8 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_8) = 1.0 + dts * l(idx_8)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_9 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_9,1) = (yp(idx_9,1) + dts * p(idx_9)) * tmp1(idx_9)
        enddo
        !$acc loop gang vector independent
        do idx_10 = 1, n_pnts, 1
          p(idx_10) = 0.0 + dj(idx_10,2) * y(idx_10,4)
        enddo
        !$acc loop gang vector independent
        do idx_11 = 1, n_pnts, 1
          l(idx_11) = 0.0 + rc(idx_11,7) + rc(idx_11,8) * h2o(idx_11)
        enddo
        !$acc loop gang vector independent
        do idx_12 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_12) = 1.0 + dts * l(idx_12)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_13 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_13,2) = (yp(idx_13,2) + dts * p(idx_13)) * tmp1(idx_13)
        enddo
        !$acc loop gang vector independent
        do idx_14 = 1, n_pnts, 1
          p(idx_14) = 0.0 + dj(idx_14,23) * y(idx_14,35) + dj(idx_14,21) * y(idx_14,60) + dj(idx_14,22) * y(idx_14,31) + dj(idx_14,19) * y(idx_14,21) + dj(idx_14,20) * y(idx_14,26)
        enddo
        !$acc loop gang vector independent
        do idx_15 = 1, n_pnts, 1
          p(idx_15) = p(idx_15) + dj(idx_15,5) * y(idx_15,10) + dj(idx_15,17) * y(idx_15,17) + dj(idx_15,4) * y(idx_15,11) + dj(idx_15,4) * y(idx_15,11) + rc(idx_15,137) * y(idx_15,35) * y(idx_15,3) + rc(idx_15,182) * y(idx_15,15) * y(idx_15,63) + rc(idx_15,129) * y(idx_15,4) * y(idx_15,34) * 0.36 + rc(idx_15,135) * y(idx_15,31) * y(idx_15,3)
        enddo
        !$acc loop gang vector independent
        do idx_16 = 1, n_pnts, 1
          p(idx_16) = p(idx_16) + rc(idx_16,123) * y(idx_16,4) * y(idx_16,58) * 0.28 + rc(idx_16,128) * y(idx_16,4) * y(idx_16,29) * 0.27 + rc(idx_16,102) * y(idx_16,3) * y(idx_16,26) + rc(idx_16,104) * y(idx_16,3) * y(idx_16,60) + rc(idx_16,45) * y(idx_16,3) * y(idx_16,17) + rc(idx_16,100) * y(idx_16,3) * y(idx_16,21) + rc(idx_16,17) * y(idx_16,5) * y(idx_16,15) + rc(idx_16,34) * y(idx_16,6) * y(idx_16,15)
        enddo
        !$acc loop gang vector independent
        do idx_17 = 1, n_pnts, 1
          p(idx_17) = p(idx_17) + rc(idx_17,8) * y(idx_17,2) * h2o(idx_17) * 2.00 + rc(idx_17,14) * y(idx_17,15) * y(idx_17,4)
        enddo
        !$acc loop gang vector independent
        do idx_18 = 1, n_pnts, 1
          l(idx_18) = 0.0 + rc(idx_18,192) * y(idx_18,56) + rc(idx_18,185) * y(idx_18,44) + rc(idx_18,187) * y(idx_18,48) + rc(idx_18,188) * y(idx_18,45) + rc(idx_18,176) * y(idx_18,66) + rc(idx_18,177) * y(idx_18,39) + rc(idx_18,184) * y(idx_18,44)
        enddo
        !$acc loop gang vector independent
        do idx_19 = 1, n_pnts, 1
          l(idx_19) = l(idx_19) + rc(idx_19,172) * y(idx_19,64) + rc(idx_19,174) * y(idx_19,62) + rc(idx_19,175) * y(idx_19,62) + rc(idx_19,156) * y(idx_19,50) + rc(idx_19,160) * y(idx_19,32) + rc(idx_19,170) * y(idx_19,66) + rc(idx_19,141) * y(idx_19,29) + rc(idx_19,143) * y(idx_19,34) + rc(idx_19,154) * y(idx_19,53) + rc(idx_19,137) * y(idx_19,35) + rc(idx_19,138) * y(idx_19,33) + rc(idx_19,139) * y(idx_19,65)
        enddo
        !$acc loop gang vector independent
        do idx_20 = 1, n_pnts, 1
          l(idx_20) = l(idx_20) + rc(idx_20,109) * y(idx_20,51) + rc(idx_20,125) * y(idx_20,58) + rc(idx_20,135) * y(idx_20,31) + rc(idx_20,100) * y(idx_20,21) + rc(idx_20,102) * y(idx_20,26) + rc(idx_20,104) * y(idx_20,60) + rc(idx_20,92) * y(idx_20,25) + rc(idx_20,94) * y(idx_20,27) + rc(idx_20,98) * y(idx_20,23) + rc(idx_20,75) * y(idx_20,22) + rc(idx_20,81) * y(idx_20,59) + rc(idx_20,86) * y(idx_20,61)
        enddo
        !$acc loop gang vector independent
        do idx_21 = 1, n_pnts, 1
          l(idx_21) = l(idx_21) + rc(idx_21,66) * y(idx_21,14) + rc(idx_21,70) * y(idx_21,13) + rc(idx_21,71) * y(idx_21,19) + rc(idx_21,44) * y(idx_21,60) + rc(idx_21,59) * y(idx_21,12) + rc(idx_21,63) * y(idx_21,41) + rc(idx_21,45) * y(idx_21,17) + rc(idx_21,44) * y(idx_21,21) + rc(idx_21,44) * y(idx_21,26) + rc(idx_21,33) * y(idx_21,43) + rc(idx_21,35) * y(idx_21,10) + rc(idx_21,44) * y(idx_21,17)
        enddo
        !$acc loop gang vector independent
        do idx_22 = 1, n_pnts, 1
          l(idx_22) = l(idx_22) + rc(idx_22,28) * y(idx_22,49) + rc(idx_22,30) * y(idx_22,15) + rc(idx_22,31) * y(idx_22,11) + rc(idx_22,13) * y(idx_22,4) + rc(idx_22,21) * y(idx_22,7) + rc(idx_22,24) * y(idx_22,9)
        enddo
        !$acc loop gang vector independent
        do idx_23 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_23) = 1.0 + dts * l(idx_23)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_24 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_24,3) = (yp(idx_24,3) + dts * p(idx_24)) * tmp1(idx_24)
        enddo
        !$acc loop gang vector independent
        do idx_25 = 1, n_pnts, 1
          p(idx_25) = 0.0 + rc(idx_25,1) * y(idx_25,1) + rc(idx_25,131) * y(idx_25,20) * y(idx_25,15) * 0.3
        enddo
        !$acc loop gang vector independent
        do idx_26 = 1, n_pnts, 1
          l(idx_26) = 0.0 + dd(idx_26,4) + dj(idx_26,1) + dj(idx_26,2) + dw(idx_26,4) + rc(idx_26,190) * y(idx_26,45) + rc(idx_26,193) * y(idx_26,56) + rc(idx_26,197) * y(idx_26,45)
        enddo
        !$acc loop gang vector independent
        do idx_27 = 1, n_pnts, 1
          l(idx_27) = l(idx_27) + rc(idx_27,124) * y(idx_27,58) + rc(idx_27,128) * y(idx_27,29) + rc(idx_27,129) * y(idx_27,34) + rc(idx_27,14) * y(idx_27,15) + rc(idx_27,112) * y(idx_27,51) + rc(idx_27,123) * y(idx_27,58) + rc(idx_27,11) * y(idx_27,5) + rc(idx_27,12) * y(idx_27,7) + rc(idx_27,13) * y(idx_27,3)
        enddo
        !$acc loop gang vector independent
        do idx_28 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_28) = 1.0 + dts * l(idx_28)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_29 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_29,4) = (yp(idx_29,4) + dts * p(idx_29)) * tmp1(idx_29)
        enddo
        !$acc loop gang vector independent
        do idx_30 = 1, n_pnts, 1
          p(idx_30) = 0.0 + dj(idx_30,3) * y(idx_30,7) + dj(idx_30,14) * y(idx_30,6) + rc(idx_30,16) * y(idx_30,7) * y(idx_30,1) + rc(idx_30,19) * y(idx_30,7) * y(idx_30,6)
        enddo
        !$acc loop gang vector independent
        do idx_31 = 1, n_pnts, 1
          l(idx_31) = 0.0 + rc(idx_31,173) * y(idx_31,52) + rc(idx_31,126) * y(idx_31,40) + rc(idx_31,142) * y(idx_31,54) + rc(idx_31,144) * y(idx_31,55) + rc(idx_31,95) * y(idx_31,36) + rc(idx_31,105) * y(idx_31,37) + rc(idx_31,110) * y(idx_31,38)
        enddo
        !$acc loop gang vector independent
        do idx_32 = 1, n_pnts, 1
          l(idx_32) = l(idx_32) + rc(idx_32,79) * y(idx_32,20) + rc(idx_32,83) * y(idx_32,24) + rc(idx_32,93) * y(idx_32,30) + rc(idx_32,17) * y(idx_32,15) + rc(idx_32,60) * y(idx_32,16) + rc(idx_32,72) * y(idx_32,18) + rc(idx_32,5) * y(idx_32,1) + rc(idx_32,11) * y(idx_32,4) + rc(idx_32,15) * y(idx_32,6)
        enddo
        !$acc loop gang vector independent
        do idx_33 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_33) = 1.0 + dts * l(idx_33)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_34 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_34,5) = (yp(idx_34,5) + dts * p(idx_34)) * tmp1(idx_34)
        enddo
        !$acc loop gang vector independent
        do idx_35 = 1, n_pnts, 1
          p(idx_35) = 0.0 + dj(idx_35,16) * y(idx_35,8) + rc(idx_35,35) * y(idx_35,3) * y(idx_35,10) + rc(idx_35,98) * y(idx_35,3) * y(idx_35,23) + rc(idx_35,12) * y(idx_35,7) * y(idx_35,4) + rc(idx_35,29) * y(idx_35,8)
        enddo
        !$acc loop gang vector independent
        do idx_36 = 1, n_pnts, 1
          l(idx_36) = 0.0 + dj(idx_36,15) + dw(idx_36,6) + rc(idx_36,186) * y(idx_36,44) + rc(idx_36,194) * y(idx_36,56) + dj(idx_36,14) + rc(idx_36,158) * y(idx_36,22) + rc(idx_36,159) * y(idx_36,29) + rc(idx_36,178) * y(idx_36,39)
        enddo
        !$acc loop gang vector independent
        do idx_37 = 1, n_pnts, 1
          l(idx_37) = l(idx_37) + rc(idx_37,152) * y(idx_37,59) + rc(idx_37,153) * y(idx_37,51) + rc(idx_37,155) * y(idx_37,58) + rc(idx_37,34) * y(idx_37,15) + rc(idx_37,67) * y(idx_37,14) + rc(idx_37,151) * y(idx_37,19) + rc(idx_37,27) * y(idx_37,6) + rc(idx_37,27) * y(idx_37,6) + rc(idx_37,32) * y(idx_37,15) + rc(idx_37,15) * y(idx_37,5) + rc(idx_37,19) * y(idx_37,7) + rc(idx_37,20) * y(idx_37,7)
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
          y(idx_39,6) = (yp(idx_39,6) + dts * p(idx_39)) * tmp1(idx_39)
        enddo
        !$acc loop gang vector independent
        do idx_40 = 1, n_pnts, 1
          p(idx_40) = 0.0 + dj(idx_40,18) * y(idx_40,23) + dj(idx_40,15) * y(idx_40,6) + dj(idx_40,16) * y(idx_40,8) + dj(idx_40,5) * y(idx_40,10) + dj(idx_40,11) * y(idx_40,9)
        enddo
        !$acc loop gang vector independent
        do idx_41 = 1, n_pnts, 1
          p(idx_41) = p(idx_41) + rc(idx_41,177) * y(idx_41,3) * y(idx_41,39) + rc(idx_41,178) * y(idx_41,6) * y(idx_41,39) * 2.00 + rc(idx_41,160) * y(idx_41,32) * y(idx_41,3) + rc(idx_41,173) * y(idx_41,52) * y(idx_41,5) + rc(idx_41,154) * y(idx_41,53) * y(idx_41,3) + rc(idx_41,156) * y(idx_41,50) * y(idx_41,3) + rc(idx_41,142) * y(idx_41,54) * y(idx_41,5) + rc(idx_41,144) * y(idx_41,55) * y(idx_41,5)
        enddo
        !$acc loop gang vector independent
        do idx_42 = 1, n_pnts, 1
          p(idx_42) = p(idx_42) + rc(idx_42,110) * y(idx_42,38) * y(idx_42,5) + rc(idx_42,126) * y(idx_42,40) * y(idx_42,5) + rc(idx_42,95) * y(idx_42,36) * y(idx_42,5) + rc(idx_42,105) * y(idx_42,5) * y(idx_42,37) + rc(idx_42,83) * y(idx_42,24) * y(idx_42,5) + rc(idx_42,93) * y(idx_42,30) * y(idx_42,5) + rc(idx_42,78) * y(idx_42,23) + rc(idx_42,79) * y(idx_42,20) * y(idx_42,5)
        enddo
        !$acc loop gang vector independent
        do idx_43 = 1, n_pnts, 1
          p(idx_43) = p(idx_43) + rc(idx_43,60) * y(idx_43,5) * y(idx_43,16) + rc(idx_43,72) * y(idx_43,18) * y(idx_43,5) + rc(idx_43,29) * y(idx_43,8) + rc(idx_43,34) * y(idx_43,6) * y(idx_43,15) + rc(idx_43,24) * y(idx_43,9) * y(idx_43,3) + rc(idx_43,27) * y(idx_43,6) * y(idx_43,6) * 2.00 + rc(idx_43,19) * y(idx_43,7) * y(idx_43,6) + rc(idx_43,23) * y(idx_43,9)
        enddo
        !$acc loop gang vector independent
        do idx_44 = 1, n_pnts, 1
          p(idx_44) = p(idx_44) + rc(idx_44,15) * y(idx_44,5) * y(idx_44,6) * 2.00 + rc(idx_44,17) * y(idx_44,5) * y(idx_44,15) + rc(idx_44,5) * y(idx_44,1) * y(idx_44,5) + rc(idx_44,11) * y(idx_44,5) * y(idx_44,4)
        enddo
        !$acc loop gang vector independent
        do idx_45 = 1, n_pnts, 1
          l(idx_45) = 0.0 + dj(idx_45,3) + dd(idx_45,7) + rc(idx_45,77) * y(idx_45,20) + rc(idx_45,171) * y(idx_45,42) + rc(idx_45,183) * y(idx_45,63) + rc(idx_45,20) * y(idx_45,6) + rc(idx_45,21) * y(idx_45,3) + rc(idx_45,22) * y(idx_45,15)
        enddo
        !$acc loop gang vector independent
        do idx_46 = 1, n_pnts, 1
          l(idx_46) = l(idx_46) + rc(idx_46,12) * y(idx_46,4) + rc(idx_46,16) * y(idx_46,1) + rc(idx_46,19) * y(idx_46,6)
        enddo
        !$acc loop gang vector independent
        do idx_47 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_47) = 1.0 + dts * l(idx_47)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_48 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_48,7) = (yp(idx_48,7) + dts * p(idx_48)) * tmp1(idx_48)
        enddo
        !$acc loop gang vector independent
        do idx_49 = 1, n_pnts, 1
          p(idx_49) = 0.0 + rc(idx_49,20) * y(idx_49,7) * y(idx_49,6)
        enddo
        !$acc loop gang vector independent
        do idx_50 = 1, n_pnts, 1
          l(idx_50) = 0.0 + dw(idx_50,8) + rc(idx_50,29) + rc(idx_50,195) + dj(idx_50,16)
        enddo
        !$acc loop gang vector independent
        do idx_51 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_51) = 1.0 + dts * l(idx_51)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_52 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_52,8) = (yp(idx_52,8) + dts * p(idx_52)) * tmp1(idx_52)
        enddo
        !$acc loop gang vector independent
        do idx_53 = 1, n_pnts, 1
          p(idx_53) = 0.0 + rc(idx_53,22) * y(idx_53,7) * y(idx_53,15)
        enddo
        !$acc loop gang vector independent
        do idx_54 = 1, n_pnts, 1
          l(idx_54) = 0.0 + dw(idx_54,9) + rc(idx_54,23) + rc(idx_54,24) * y(idx_54,3) + dj(idx_54,11)
        enddo
        !$acc loop gang vector independent
        do idx_55 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_55) = 1.0 + dts * l(idx_55)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_56 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_56,9) = (yp(idx_56,9) + dts * p(idx_56)) * tmp1(idx_56)
        enddo
        !$acc loop gang vector independent
        do idx_57 = 1, n_pnts, 1
          p(idx_57) = 0.0 + rc(idx_57,186) * y(idx_57,44) * y(idx_57,6) + rc(idx_57,195) * y(idx_57,8) * 2.0 + rc(idx_57,152) * y(idx_57,6) * y(idx_57,59) + rc(idx_57,158) * y(idx_57,6) * y(idx_57,22) + rc(idx_57,67) * y(idx_57,6) * y(idx_57,14) + rc(idx_57,151) * y(idx_57,6) * y(idx_57,19)
        enddo
        !$acc loop gang vector independent
        do idx_58 = 1, n_pnts, 1
          p(idx_58) = p(idx_58) + rc(idx_58,21) * y(idx_58,7) * y(idx_58,3) + rc(idx_58,32) * y(idx_58,6) * y(idx_58,15)
        enddo
        !$acc loop gang vector independent
        do idx_59 = 1, n_pnts, 1
          l(idx_59) = 0.0 + dd(idx_59,10) + rc(idx_59,35) * y(idx_59,3) + dj(idx_59,5) + dw(idx_59,10)
        enddo
        !$acc loop gang vector independent
        do idx_60 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_60) = 1.0 + dts * l(idx_60)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_61 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_61,10) = (yp(idx_61,10) + dts * p(idx_61)) * tmp1(idx_61)
        enddo
        !$acc loop gang vector independent
        do idx_62 = 1, n_pnts, 1
          p(idx_62) = 0.0 + rc(idx_62,36) * y(idx_62,15) * y(idx_62,15) + rc(idx_62,196) * y(idx_62,15) * y(idx_62,15)
        enddo
        !$acc loop gang vector independent
        do idx_63 = 1, n_pnts, 1
          l(idx_63) = 0.0 + dw(idx_63,11) + dd(idx_63,11) + rc(idx_63,31) * y(idx_63,3) + rc(idx_63,191) * y(idx_63,45) + dj(idx_63,4)
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
          y(idx_65,11) = (yp(idx_65,11) + dts * p(idx_65)) * tmp1(idx_65)
        enddo
        !$acc loop gang vector independent
        do idx_66 = 1, n_pnts, 1
          p(idx_66) = 0.0 + rc(idx_66,123) * y(idx_66,4) * y(idx_66,58) * 0.30
        enddo
        !$acc loop gang vector independent
        do idx_67 = 1, n_pnts, 1
          l(idx_67) = 0.0 + rc(idx_67,59) * y(idx_67,3) + dd(idx_67,12)
        enddo
        !$acc loop gang vector independent
        do idx_68 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_68) = 1.0 + dts * l(idx_68)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_69 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_69,12) = (yp(idx_69,12) + dts * p(idx_69)) * tmp1(idx_69)
        enddo
        !$acc loop gang vector independent
        do idx_70 = 1, n_pnts, 1
          p(idx_70) = 0.0 + dj(idx_70,13) * y(idx_70,65) + dj(idx_70,12) * y(idx_70,33) + dj(idx_70,13) * y(idx_70,65) + dj(idx_70,7) * y(idx_70,14) + dj(idx_70,8) * y(idx_70,22)
        enddo
        !$acc loop gang vector independent
        do idx_71 = 1, n_pnts, 1
          p(idx_71) = p(idx_71) + rc(idx_71,139) * y(idx_71,65) * y(idx_71,3) * 2.0 + dj(idx_71,6) * y(idx_71,14) + rc(idx_71,129) * y(idx_71,4) * y(idx_71,34) * 0.76 + rc(idx_71,138) * y(idx_71,33) * y(idx_71,3) + rc(idx_71,124) * y(idx_71,4) * y(idx_71,58) * 0.58 + rc(idx_71,128) * y(idx_71,4) * y(idx_71,29) * 0.78 + rc(idx_71,112) * y(idx_71,4) * y(idx_71,51) * 0.31 + rc(idx_71,123) * y(idx_71,4) * y(idx_71,58) * 0.40
        enddo
        !$acc loop gang vector independent
        do idx_72 = 1, n_pnts, 1
          p(idx_72) = p(idx_72) + rc(idx_72,66) * y(idx_72,3) * y(idx_72,14) + rc(idx_72,67) * y(idx_72,6) * y(idx_72,14)
        enddo
        !$acc loop gang vector independent
        do idx_73 = 1, n_pnts, 1
          l(idx_73) = 0.0 + rc(idx_73,70) * y(idx_73,3) + dd(idx_73,13)
        enddo
        !$acc loop gang vector independent
        do idx_74 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_74) = 1.0 + dts * l(idx_74)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_75 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_75,13) = (yp(idx_75,13) + dts * p(idx_75)) * tmp1(idx_75)
        enddo
        !$acc loop gang vector independent
        do idx_76 = 1, n_pnts, 1
          p(idx_76) = 0.0 + dj(idx_76,22) * y(idx_76,31) + dj(idx_76,23) * y(idx_76,35) + rc(idx_76,186) * y(idx_76,44) * y(idx_76,6) + dj(idx_76,17) * y(idx_76,17) + rc(idx_76,181) * y(idx_76,52) * y(idx_76,16) + rc(idx_76,184) * y(idx_76,44) * y(idx_76,3)
        enddo
        !$acc loop gang vector independent
        do idx_77 = 1, n_pnts, 1
          p(idx_77) = p(idx_77) + rc(idx_77,144) * y(idx_77,55) * y(idx_77,5) + rc(idx_77,154) * y(idx_77,53) * y(idx_77,3) + rc(idx_77,137) * y(idx_77,35) * y(idx_77,3) + rc(idx_77,142) * y(idx_77,54) * y(idx_77,5) + rc(idx_77,133) * y(idx_77,55) * y(idx_77,16) * 1.0 + rc(idx_77,135) * y(idx_77,31) * y(idx_77,3) + rc(idx_77,129) * y(idx_77,4) * y(idx_77,34) * 0.24 + rc(idx_77,132) * y(idx_77,54) * y(idx_77,16) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_78 = 1, n_pnts, 1
          p(idx_78) = p(idx_78) + rc(idx_78,127) * y(idx_78,16) * y(idx_78,40) * 2.00 + rc(idx_78,128) * y(idx_78,4) * y(idx_78,29) * 0.22 + rc(idx_78,123) * y(idx_78,4) * y(idx_78,58) + rc(idx_78,126) * y(idx_78,40) * y(idx_78,5) + rc(idx_78,112) * y(idx_78,4) * y(idx_78,51) + rc(idx_78,112) * y(idx_78,4) * y(idx_78,51) * 0.47 + rc(idx_78,110) * y(idx_78,38) * y(idx_78,5) * 2.00 + rc(idx_78,111) * y(idx_78,16) * y(idx_78,38) * 3.00
        enddo
        !$acc loop gang vector independent
        do idx_79 = 1, n_pnts, 1
          p(idx_79) = p(idx_79) + rc(idx_79,98) * y(idx_79,3) * y(idx_79,23) + rc(idx_79,106) * y(idx_79,16) * y(idx_79,37) + rc(idx_79,96) * y(idx_79,36) * y(idx_79,16) * 2.00 + rc(idx_79,97) * y(idx_79,30) * y(idx_79,16) + rc(idx_79,84) * y(idx_79,24) * y(idx_79,16) + rc(idx_79,95) * y(idx_79,36) * y(idx_79,5) + rc(idx_79,74) * y(idx_79,16) * y(idx_79,20) * 2.00 + rc(idx_79,80) * y(idx_79,16) * y(idx_79,20)
        enddo
        !$acc loop gang vector independent
        do idx_80 = 1, n_pnts, 1
          p(idx_80) = p(idx_80) + rc(idx_80,63) * y(idx_80,41) * y(idx_80,3) + rc(idx_80,73) * y(idx_80,18) * y(idx_80,16) + rc(idx_80,61) * y(idx_80,16) * y(idx_80,16) * 2.00 + rc(idx_80,62) * y(idx_80,16) * y(idx_80,16) + rc(idx_80,45) * y(idx_80,3) * y(idx_80,17) + rc(idx_80,60) * y(idx_80,5) * y(idx_80,16)
        enddo
        !$acc loop gang vector independent
        do idx_81 = 1, n_pnts, 1
          l(idx_81) = 0.0 + dj(idx_81,7) + dw(idx_81,14) + rc(idx_81,66) * y(idx_81,3) + rc(idx_81,67) * y(idx_81,6) + dj(idx_81,6)
        enddo
        !$acc loop gang vector independent
        do idx_82 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_82) = 1.0 + dts * l(idx_82)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_83 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_83,14) = (yp(idx_83,14) + dts * p(idx_83)) * tmp1(idx_83)
        enddo
        !$acc loop gang vector independent
        do idx_84 = 1, n_pnts, 1
          p(idx_84) = 0.0 + dj(idx_84,23) * y(idx_84,35) + dj(idx_84,21) * y(idx_84,60) + dj(idx_84,22) * y(idx_84,31) + dj(idx_84,19) * y(idx_84,21) + dj(idx_84,20) * y(idx_84,26)
        enddo
        !$acc loop gang vector independent
        do idx_85 = 1, n_pnts, 1
          p(idx_85) = p(idx_85) + dj(idx_85,13) * y(idx_85,65) + dj(idx_85,17) * y(idx_85,17) + dj(idx_85,12) * y(idx_85,33) + dj(idx_85,13) * y(idx_85,65) + dj(idx_85,8) * y(idx_85,22) + dj(idx_85,11) * y(idx_85,9) + dj(idx_85,6) * y(idx_85,14) + dj(idx_85,6) * y(idx_85,14)
        enddo
        !$acc loop gang vector independent
        do idx_86 = 1, n_pnts, 1
          p(idx_86) = p(idx_86) + rc(idx_86,181) * y(idx_86,52) * y(idx_86,16) * 2.00 + rc(idx_86,188) * y(idx_86,45) * y(idx_86,3) + rc(idx_86,173) * y(idx_86,52) * y(idx_86,5) + rc(idx_86,174) * y(idx_86,3) * y(idx_86,62) + rc(idx_86,159) * y(idx_86,6) * y(idx_86,29) + rc(idx_86,170) * y(idx_86,3) * y(idx_86,66) + rc(idx_86,153) * y(idx_86,6) * y(idx_86,51) + rc(idx_86,155) * y(idx_86,6) * y(idx_86,58)
        enddo
        !$acc loop gang vector independent
        do idx_87 = 1, n_pnts, 1
          p(idx_87) = p(idx_87) + rc(idx_87,142) * y(idx_87,54) * y(idx_87,5) + rc(idx_87,144) * y(idx_87,55) * y(idx_87,5) + rc(idx_87,133) * y(idx_87,55) * y(idx_87,16) * 2.0 + rc(idx_87,139) * y(idx_87,65) * y(idx_87,3) * 1.0 + rc(idx_87,129) * y(idx_87,4) * y(idx_87,34) * 0.36 + rc(idx_87,132) * y(idx_87,54) * y(idx_87,16) * 2.0 + rc(idx_87,127) * y(idx_87,16) * y(idx_87,40) * 2.00 + rc(idx_87,128) * y(idx_87,4) * y(idx_87,29) * 0.27
        enddo
        !$acc loop gang vector independent
        do idx_88 = 1, n_pnts, 1
          p(idx_88) = p(idx_88) + rc(idx_88,124) * y(idx_88,4) * y(idx_88,58) * 0.18 + rc(idx_88,126) * y(idx_88,40) * y(idx_88,5) + rc(idx_88,112) * y(idx_88,4) * y(idx_88,51) * 0.20 + rc(idx_88,123) * y(idx_88,4) * y(idx_88,58) * 0.30 + rc(idx_88,110) * y(idx_88,38) * y(idx_88,5) + rc(idx_88,111) * y(idx_88,16) * y(idx_88,38) * 2.00 + rc(idx_88,97) * y(idx_88,30) * y(idx_88,16) * 2.00 + rc(idx_88,106) * y(idx_88,16) * y(idx_88,37)
        enddo
        !$acc loop gang vector independent
        do idx_89 = 1, n_pnts, 1
          p(idx_89) = p(idx_89) + rc(idx_89,93) * y(idx_89,30) * y(idx_89,5) + rc(idx_89,96) * y(idx_89,36) * y(idx_89,16) + rc(idx_89,84) * y(idx_89,24) * y(idx_89,16) * 2.00 + rc(idx_89,90) * y(idx_89,18) * y(idx_89,18) * 2.00 + rc(idx_89,80) * y(idx_89,16) * y(idx_89,20) + rc(idx_89,83) * y(idx_89,24) * y(idx_89,5) + rc(idx_89,72) * y(idx_89,18) * y(idx_89,5) + rc(idx_89,73) * y(idx_89,18) * y(idx_89,16) * 2.00
        enddo
        !$acc loop gang vector independent
        do idx_90 = 1, n_pnts, 1
          p(idx_90) = p(idx_90) + rc(idx_90,67) * y(idx_90,6) * y(idx_90,14) + rc(idx_90,70) * y(idx_90,3) * y(idx_90,13) + rc(idx_90,63) * y(idx_90,41) * y(idx_90,3) + rc(idx_90,66) * y(idx_90,3) * y(idx_90,14) + rc(idx_90,60) * y(idx_90,5) * y(idx_90,16) + rc(idx_90,61) * y(idx_90,16) * y(idx_90,16) * 2.00 + rc(idx_90,31) * y(idx_90,3) * y(idx_90,11) + rc(idx_90,33) * y(idx_90,3) * y(idx_90,43)
        enddo
        !$acc loop gang vector independent
        do idx_91 = 1, n_pnts, 1
          p(idx_91) = p(idx_91) + rc(idx_91,13) * y(idx_91,3) * y(idx_91,4) + rc(idx_91,23) * y(idx_91,9)
        enddo
        !$acc loop gang vector independent
        do idx_92 = 1, n_pnts, 1
          l(idx_92) = 0.0 + rc(idx_92,196) * y(idx_92,15) + dw(idx_92,15) + rc(idx_92,180) * y(idx_92,42) + rc(idx_92,182) * y(idx_92,63) + rc(idx_92,196) * y(idx_92,15) + rc(idx_92,131) * y(idx_92,20) + rc(idx_92,134) * y(idx_92,54) + rc(idx_92,136) * y(idx_92,55)
        enddo
        !$acc loop gang vector independent
        do idx_93 = 1, n_pnts, 1
          l(idx_93) = l(idx_93) + rc(idx_93,99) * y(idx_93,18) + rc(idx_93,101) * y(idx_93,30) + rc(idx_93,103) * y(idx_93,24) + rc(idx_93,36) * y(idx_93,15) + rc(idx_93,36) * y(idx_93,15) + rc(idx_93,65) * y(idx_93,16) + rc(idx_93,30) * y(idx_93,3) + rc(idx_93,32) * y(idx_93,6) + rc(idx_93,34) * y(idx_93,6) + rc(idx_93,14) * y(idx_93,4) + rc(idx_93,17) * y(idx_93,5) + rc(idx_93,22) * y(idx_93,7)
        enddo
        !$acc loop gang vector independent
        do idx_94 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_94) = 1.0 + dts * l(idx_94)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_95 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_95,15) = (yp(idx_95,15) + dts * p(idx_95)) * tmp1(idx_95)
        enddo
        !$acc loop gang vector independent
        do idx_96 = 1, n_pnts, 1
          p(idx_96) = 0.0 + dj(idx_96,8) * y(idx_96,22) + dj(idx_96,10) * y(idx_96,27) + rc(idx_96,185) * y(idx_96,44) * y(idx_96,3) * 1.0 + rc(idx_96,186) * y(idx_96,44) * y(idx_96,6) + rc(idx_96,131) * y(idx_96,20) * y(idx_96,15) * 0.8 + rc(idx_96,184) * y(idx_96,44) * y(idx_96,3)
        enddo
        !$acc loop gang vector independent
        do idx_97 = 1, n_pnts, 1
          p(idx_97) = p(idx_97) + rc(idx_97,91) * y(idx_97,20) * y(idx_97,20) * 2.00 + rc(idx_97,123) * y(idx_97,4) * y(idx_97,58) * 0.58 + rc(idx_97,79) * y(idx_97,20) * y(idx_97,5) + rc(idx_97,80) * y(idx_97,16) * y(idx_97,20) + rc(idx_97,44) * y(idx_97,3) * y(idx_97,17) + rc(idx_97,59) * y(idx_97,3) * y(idx_97,12)
        enddo
        !$acc loop gang vector independent
        do idx_98 = 1, n_pnts, 1
          l(idx_98) = 0.0 + dw(idx_98,16) + rc(idx_98,132) * y(idx_98,54) + rc(idx_98,133) * y(idx_98,55) + rc(idx_98,181) * y(idx_98,52) + rc(idx_98,106) * y(idx_98,37) + rc(idx_98,111) * y(idx_98,38) + rc(idx_98,127) * y(idx_98,40)
        enddo
        !$acc loop gang vector independent
        do idx_99 = 1, n_pnts, 1
          l(idx_99) = l(idx_99) + rc(idx_99,84) * y(idx_99,24) + rc(idx_99,96) * y(idx_99,36) + rc(idx_99,97) * y(idx_99,30) + rc(idx_99,73) * y(idx_99,18) + rc(idx_99,74) * y(idx_99,20) + rc(idx_99,80) * y(idx_99,20) + rc(idx_99,62) * y(idx_99,16) + rc(idx_99,62) * y(idx_99,16) + rc(idx_99,65) * y(idx_99,15) + rc(idx_99,60) * y(idx_99,5) + rc(idx_99,61) * y(idx_99,16) + rc(idx_99,61) * y(idx_99,16)
        enddo
        !$acc loop gang vector independent
        do idx_100 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_100) = 1.0 + dts * l(idx_100)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_101 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_101,16) = (yp(idx_101,16) + dts * p(idx_101)) * tmp1(idx_101)
        enddo
        !$acc loop gang vector independent
        do idx_102 = 1, n_pnts, 1
          p(idx_102) = 0.0 + rc(idx_102,65) * y(idx_102,16) * y(idx_102,15)
        enddo
        !$acc loop gang vector independent
        do idx_103 = 1, n_pnts, 1
          l(idx_103) = 0.0 + dw(idx_103,17) + dd(idx_103,17) + rc(idx_103,44) * y(idx_103,3) + rc(idx_103,45) * y(idx_103,3) + dj(idx_103,17)
        enddo
        !$acc loop gang vector independent
        do idx_104 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_104) = 1.0 + dts * l(idx_104)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_105 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_105,17) = (yp(idx_105,17) + dts * p(idx_105)) * tmp1(idx_105)
        enddo
        !$acc loop gang vector independent
        do idx_106 = 1, n_pnts, 1
          p(idx_106) = 0.0 + rc(idx_106,151) * y(idx_106,6) * y(idx_106,19) + dj(idx_106,9) * y(idx_106,61) + rc(idx_106,44) * y(idx_106,3) * y(idx_106,21) + rc(idx_106,71) * y(idx_106,3) * y(idx_106,19)
        enddo
        !$acc loop gang vector independent
        do idx_107 = 1, n_pnts, 1
          l(idx_107) = 0.0 + rc(idx_107,90) * y(idx_107,18) + rc(idx_107,99) * y(idx_107,15) + rc(idx_107,72) * y(idx_107,5) + rc(idx_107,73) * y(idx_107,16) + rc(idx_107,90) * y(idx_107,18)
        enddo
        !$acc loop gang vector independent
        do idx_108 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_108) = 1.0 + dts * l(idx_108)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_109 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_109,18) = (yp(idx_109,18) + dts * p(idx_109)) * tmp1(idx_109)
        enddo
        !$acc loop gang vector independent
        do idx_110 = 1, n_pnts, 1
          p(idx_110) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_111 = 1, n_pnts, 1
          l(idx_111) = 0.0 + rc(idx_111,71) * y(idx_111,3) + rc(idx_111,151) * y(idx_111,6)
        enddo
        !$acc loop gang vector independent
        do idx_112 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_112) = 1.0 + dts * l(idx_112)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_113 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_113,19) = (yp(idx_113,19) + dts * p(idx_113)) * tmp1(idx_113)
        enddo
        !$acc loop gang vector independent
        do idx_114 = 1, n_pnts, 1
          p(idx_114) = 0.0 + dj(idx_114,18) * y(idx_114,23) + dj(idx_114,10) * y(idx_114,27) + dj(idx_114,12) * y(idx_114,33) + rc(idx_114,158) * y(idx_114,6) * y(idx_114,22) + dj(idx_114,9) * y(idx_114,61)
        enddo
        !$acc loop gang vector independent
        do idx_115 = 1, n_pnts, 1
          p(idx_115) = p(idx_115) + rc(idx_115,131) * y(idx_115,20) * y(idx_115,15) * 0.2 + rc(idx_115,138) * y(idx_115,33) * y(idx_115,3) + rc(idx_115,105) * y(idx_115,5) * y(idx_115,37) + rc(idx_115,106) * y(idx_115,16) * y(idx_115,37) + rc(idx_115,95) * y(idx_115,36) * y(idx_115,5) + rc(idx_115,96) * y(idx_115,36) * y(idx_115,16) + rc(idx_115,75) * y(idx_115,3) * y(idx_115,22) + rc(idx_115,78) * y(idx_115,23)
        enddo
        !$acc loop gang vector independent
        do idx_116 = 1, n_pnts, 1
          l(idx_116) = 0.0 + rc(idx_116,131) * y(idx_116,15) + rc(idx_116,80) * y(idx_116,16) + rc(idx_116,91) * y(idx_116,20) + rc(idx_116,91) * y(idx_116,20) + rc(idx_116,74) * y(idx_116,16) + rc(idx_116,77) * y(idx_116,7) + rc(idx_116,79) * y(idx_116,5)
        enddo
        !$acc loop gang vector independent
        do idx_117 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_117) = 1.0 + dts * l(idx_117)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_118 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_118,20) = (yp(idx_118,20) + dts * p(idx_118)) * tmp1(idx_118)
        enddo
        !$acc loop gang vector independent
        do idx_119 = 1, n_pnts, 1
          p(idx_119) = 0.0 + rc(idx_119,99) * y(idx_119,15) * y(idx_119,18)
        enddo
        !$acc loop gang vector independent
        do idx_120 = 1, n_pnts, 1
          l(idx_120) = 0.0 + dw(idx_120,21) + dd(idx_120,21) + rc(idx_120,44) * y(idx_120,3) + rc(idx_120,100) * y(idx_120,3) + dj(idx_120,19)
        enddo
        !$acc loop gang vector independent
        do idx_121 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_121) = 1.0 + dts * l(idx_121)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_122 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_122,21) = (yp(idx_122,21) + dts * p(idx_122)) * tmp1(idx_122)
        enddo
        !$acc loop gang vector independent
        do idx_123 = 1, n_pnts, 1
          p(idx_123) = 0.0 + dj(idx_123,19) * y(idx_123,21) + rc(idx_123,127) * y(idx_123,16) * y(idx_123,40) + rc(idx_123,156) * y(idx_123,50) * y(idx_123,3) + rc(idx_123,124) * y(idx_123,4) * y(idx_123,58) + rc(idx_123,126) * y(idx_123,40) * y(idx_123,5)
        enddo
        !$acc loop gang vector independent
        do idx_124 = 1, n_pnts, 1
          p(idx_124) = p(idx_124) + rc(idx_124,105) * y(idx_124,5) * y(idx_124,37) + rc(idx_124,106) * y(idx_124,16) * y(idx_124,37) + rc(idx_124,90) * y(idx_124,18) * y(idx_124,18) * 2.00 + rc(idx_124,100) * y(idx_124,3) * y(idx_124,21) + rc(idx_124,72) * y(idx_124,18) * y(idx_124,5) + rc(idx_124,73) * y(idx_124,18) * y(idx_124,16)
        enddo
        !$acc loop gang vector independent
        do idx_125 = 1, n_pnts, 1
          l(idx_125) = 0.0 + rc(idx_125,75) * y(idx_125,3) + rc(idx_125,158) * y(idx_125,6) + dj(idx_125,8)
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
          y(idx_127,22) = (yp(idx_127,22) + dts * p(idx_127)) * tmp1(idx_127)
        enddo
        !$acc loop gang vector independent
        do idx_128 = 1, n_pnts, 1
          p(idx_128) = 0.0 + rc(idx_128,77) * y(idx_128,20) * y(idx_128,7)
        enddo
        !$acc loop gang vector independent
        do idx_129 = 1, n_pnts, 1
          l(idx_129) = 0.0 + dd(idx_129,23) + rc(idx_129,78) + rc(idx_129,98) * y(idx_129,3) + dj(idx_129,18)
        enddo
        !$acc loop gang vector independent
        do idx_130 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_130) = 1.0 + dts * l(idx_130)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_131 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_131,23) = (yp(idx_131,23) + dts * p(idx_131)) * tmp1(idx_131)
        enddo
        !$acc loop gang vector independent
        do idx_132 = 1, n_pnts, 1
          p(idx_132) = 0.0 + rc(idx_132,152) * y(idx_132,6) * y(idx_132,59) + rc(idx_132,44) * y(idx_132,3) * y(idx_132,60) + rc(idx_132,81) * y(idx_132,3) * y(idx_132,59)
        enddo
        !$acc loop gang vector independent
        do idx_133 = 1, n_pnts, 1
          l(idx_133) = 0.0 + rc(idx_133,83) * y(idx_133,5) + rc(idx_133,84) * y(idx_133,16) + rc(idx_133,103) * y(idx_133,15)
        enddo
        !$acc loop gang vector independent
        do idx_134 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_134) = 1.0 + dts * l(idx_134)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_135 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_135,24) = (yp(idx_135,24) + dts * p(idx_135)) * tmp1(idx_135)
        enddo
        !$acc loop gang vector independent
        do idx_136 = 1, n_pnts, 1
          p(idx_136) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_137 = 1, n_pnts, 1
          l(idx_137) = 0.0 + rc(idx_137,92) * y(idx_137,3)
        enddo
        !$acc loop gang vector independent
        do idx_138 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_138) = 1.0 + dts * l(idx_138)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_139 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_139,25) = (yp(idx_139,25) + dts * p(idx_139)) * tmp1(idx_139)
        enddo
        !$acc loop gang vector independent
        do idx_140 = 1, n_pnts, 1
          p(idx_140) = 0.0 + rc(idx_140,101) * y(idx_140,15) * y(idx_140,30)
        enddo
        !$acc loop gang vector independent
        do idx_141 = 1, n_pnts, 1
          l(idx_141) = 0.0 + dw(idx_141,26) + dd(idx_141,26) + rc(idx_141,44) * y(idx_141,3) + rc(idx_141,102) * y(idx_141,3) + dj(idx_141,20)
        enddo
        !$acc loop gang vector independent
        do idx_142 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_142) = 1.0 + dts * l(idx_142)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_143 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_143,26) = (yp(idx_143,26) + dts * p(idx_143)) * tmp1(idx_143)
        enddo
        !$acc loop gang vector independent
        do idx_144 = 1, n_pnts, 1
          p(idx_144) = 0.0 + rc(idx_144,102) * y(idx_144,3) * y(idx_144,26) + dj(idx_144,20) * y(idx_144,26) + rc(idx_144,93) * y(idx_144,30) * y(idx_144,5) + rc(idx_144,97) * y(idx_144,30) * y(idx_144,16)
        enddo
        !$acc loop gang vector independent
        do idx_145 = 1, n_pnts, 1
          l(idx_145) = 0.0 + rc(idx_145,94) * y(idx_145,3) + dj(idx_145,10)
        enddo
        !$acc loop gang vector independent
        do idx_146 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_146) = 1.0 + dts * l(idx_146)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_147 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_147,27) = (yp(idx_147,27) + dts * p(idx_147)) * tmp1(idx_147)
        enddo
        !$acc loop gang vector independent
        do idx_148 = 1, n_pnts, 1
          p(idx_148) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_149 = 1, n_pnts, 1
          l(idx_149) = 0.0 + dd(idx_149,28)
        enddo
        !$acc loop gang vector independent
        do idx_150 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_150) = 1.0 + dts * l(idx_150)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_151 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_151,28) = (yp(idx_151,28) + dts * p(idx_151)) * tmp1(idx_151)
        enddo
        !$acc loop gang vector independent
        do idx_152 = 1, n_pnts, 1
          p(idx_152) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_153 = 1, n_pnts, 1
          l(idx_153) = 0.0 + rc(idx_153,128) * y(idx_153,4) + rc(idx_153,141) * y(idx_153,3) + rc(idx_153,159) * y(idx_153,6)
        enddo
        !$acc loop gang vector independent
        do idx_154 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_154) = 1.0 + dts * l(idx_154)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_155 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_155,29) = (yp(idx_155,29) + dts * p(idx_155)) * tmp1(idx_155)
        enddo
        !$acc loop gang vector independent
        do idx_156 = 1, n_pnts, 1
          p(idx_156) = 0.0 + rc(idx_156,44) * y(idx_156,3) * y(idx_156,26) + rc(idx_156,92) * y(idx_156,3) * y(idx_156,25)
        enddo
        !$acc loop gang vector independent
        do idx_157 = 1, n_pnts, 1
          l(idx_157) = 0.0 + rc(idx_157,93) * y(idx_157,5) + rc(idx_157,97) * y(idx_157,16) + rc(idx_157,101) * y(idx_157,15)
        enddo
        !$acc loop gang vector independent
        do idx_158 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_158) = 1.0 + dts * l(idx_158)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_159 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_159,30) = (yp(idx_159,30) + dts * p(idx_159)) * tmp1(idx_159)
        enddo
        !$acc loop gang vector independent
        do idx_160 = 1, n_pnts, 1
          p(idx_160) = 0.0 + rc(idx_160,134) * y(idx_160,54) * y(idx_160,15)
        enddo
        !$acc loop gang vector independent
        do idx_161 = 1, n_pnts, 1
          l(idx_161) = 0.0 + dd(idx_161,31) + rc(idx_161,135) * y(idx_161,3) + dj(idx_161,22) + dw(idx_161,31)
        enddo
        !$acc loop gang vector independent
        do idx_162 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_162) = 1.0 + dts * l(idx_162)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_163 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_163,31) = (yp(idx_163,31) + dts * p(idx_163)) * tmp1(idx_163)
        enddo
        !$acc loop gang vector independent
        do idx_164 = 1, n_pnts, 1
          p(idx_164) = 0.0 + rc(idx_164,159) * y(idx_164,6) * y(idx_164,29)
        enddo
        !$acc loop gang vector independent
        do idx_165 = 1, n_pnts, 1
          l(idx_165) = 0.0 + rc(idx_165,160) * y(idx_165,3) + dw(idx_165,32)
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
          y(idx_167,32) = (yp(idx_167,32) + dts * p(idx_167)) * tmp1(idx_167)
        enddo
        !$acc loop gang vector independent
        do idx_168 = 1, n_pnts, 1
          p(idx_168) = 0.0 + dj(idx_168,23) * y(idx_168,35) + rc(idx_168,180) * y(idx_168,42) * y(idx_168,15) + rc(idx_168,181) * y(idx_168,52) * y(idx_168,16) + rc(idx_168,170) * y(idx_168,3) * y(idx_168,66) * 0.80 + rc(idx_168,173) * y(idx_168,52) * y(idx_168,5)
        enddo
        !$acc loop gang vector independent
        do idx_169 = 1, n_pnts, 1
          p(idx_169) = p(idx_169) + rc(idx_169,137) * y(idx_169,35) * y(idx_169,3) + rc(idx_169,144) * y(idx_169,55) * y(idx_169,5) + rc(idx_169,129) * y(idx_169,4) * y(idx_169,34) + rc(idx_169,133) * y(idx_169,55) * y(idx_169,16) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_170 = 1, n_pnts, 1
          l(idx_170) = 0.0 + rc(idx_170,138) * y(idx_170,3) + dj(idx_170,12) + dw(idx_170,33)
        enddo
        !$acc loop gang vector independent
        do idx_171 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_171) = 1.0 + dts * l(idx_171)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_172 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_172,33) = (yp(idx_172,33) + dts * p(idx_172)) * tmp1(idx_172)
        enddo
        !$acc loop gang vector independent
        do idx_173 = 1, n_pnts, 1
          p(idx_173) = 0.0 + rc(idx_173,160) * y(idx_173,32) * y(idx_173,3) + dj(idx_173,22) * y(idx_173,31) + rc(idx_173,135) * y(idx_173,31) * y(idx_173,3) + rc(idx_173,142) * y(idx_173,54) * y(idx_173,5) + rc(idx_173,128) * y(idx_173,4) * y(idx_173,29) + rc(idx_173,132) * y(idx_173,54) * y(idx_173,16) * 1.0
        enddo
        !$acc loop gang vector independent
        do idx_174 = 1, n_pnts, 1
          l(idx_174) = 0.0 + rc(idx_174,129) * y(idx_174,4) + rc(idx_174,143) * y(idx_174,3)
        enddo
        !$acc loop gang vector independent
        do idx_175 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_175) = 1.0 + dts * l(idx_175)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_176 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_176,34) = (yp(idx_176,34) + dts * p(idx_176)) * tmp1(idx_176)
        enddo
        !$acc loop gang vector independent
        do idx_177 = 1, n_pnts, 1
          p(idx_177) = 0.0 + rc(idx_177,136) * y(idx_177,55) * y(idx_177,15)
        enddo
        !$acc loop gang vector independent
        do idx_178 = 1, n_pnts, 1
          l(idx_178) = 0.0 + dd(idx_178,35) + rc(idx_178,137) * y(idx_178,3) + dj(idx_178,23) + dw(idx_178,35)
        enddo
        !$acc loop gang vector independent
        do idx_179 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_179) = 1.0 + dts * l(idx_179)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_180 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_180,35) = (yp(idx_180,35) + dts * p(idx_180)) * tmp1(idx_180)
        enddo
        !$acc loop gang vector independent
        do idx_181 = 1, n_pnts, 1
          p(idx_181) = 0.0 + rc(idx_181,94) * y(idx_181,3) * y(idx_181,27)
        enddo
        !$acc loop gang vector independent
        do idx_182 = 1, n_pnts, 1
          l(idx_182) = 0.0 + rc(idx_182,95) * y(idx_182,5) + rc(idx_182,96) * y(idx_182,16)
        enddo
        !$acc loop gang vector independent
        do idx_183 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_183) = 1.0 + dts * l(idx_183)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_184 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_184,36) = (yp(idx_184,36) + dts * p(idx_184)) * tmp1(idx_184)
        enddo
        !$acc loop gang vector independent
        do idx_185 = 1, n_pnts, 1
          p(idx_185) = 0.0 + rc(idx_185,86) * y(idx_185,3) * y(idx_185,61)
        enddo
        !$acc loop gang vector independent
        do idx_186 = 1, n_pnts, 1
          l(idx_186) = 0.0 + rc(idx_186,105) * y(idx_186,5) + rc(idx_186,106) * y(idx_186,16)
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
          y(idx_188,37) = (yp(idx_188,37) + dts * p(idx_188)) * tmp1(idx_188)
        enddo
        !$acc loop gang vector independent
        do idx_189 = 1, n_pnts, 1
          p(idx_189) = 0.0 + rc(idx_189,109) * y(idx_189,3) * y(idx_189,51)
        enddo
        !$acc loop gang vector independent
        do idx_190 = 1, n_pnts, 1
          l(idx_190) = 0.0 + rc(idx_190,110) * y(idx_190,5) + rc(idx_190,111) * y(idx_190,16)
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
          y(idx_192,38) = (yp(idx_192,38) + dts * p(idx_192)) * tmp1(idx_192)
        enddo
        !$acc loop gang vector independent
        do idx_193 = 1, n_pnts, 1
          p(idx_193) = 0.0 + rc(idx_193,171) * y(idx_193,7) * y(idx_193,42) + rc(idx_193,183) * y(idx_193,7) * y(idx_193,63)
        enddo
        !$acc loop gang vector independent
        do idx_194 = 1, n_pnts, 1
          l(idx_194) = 0.0 + dd(idx_194,39) + rc(idx_194,177) * y(idx_194,3) + rc(idx_194,178) * y(idx_194,6) + dw(idx_194,39)
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
          y(idx_196,39) = (yp(idx_196,39) + dts * p(idx_196)) * tmp1(idx_196)
        enddo
        !$acc loop gang vector independent
        do idx_197 = 1, n_pnts, 1
          p(idx_197) = 0.0 + rc(idx_197,125) * y(idx_197,3) * y(idx_197,58)
        enddo
        !$acc loop gang vector independent
        do idx_198 = 1, n_pnts, 1
          l(idx_198) = 0.0 + rc(idx_198,126) * y(idx_198,5) + rc(idx_198,127) * y(idx_198,16)
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
          y(idx_200,40) = (yp(idx_200,40) + dts * p(idx_200)) * tmp1(idx_200)
        enddo
        !$acc loop gang vector independent
        do idx_201 = 1, n_pnts, 1
          p(idx_201) = 0.0 + rc(idx_201,62) * y(idx_201,16) * y(idx_201,16) + rc(idx_201,123) * y(idx_201,4) * y(idx_201,58) * 0.12
        enddo
        !$acc loop gang vector independent
        do idx_202 = 1, n_pnts, 1
          l(idx_202) = 0.0 + rc(idx_202,63) * y(idx_202,3) + dw(idx_202,41)
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
          y(idx_204,41) = (yp(idx_204,41) + dts * p(idx_204)) * tmp1(idx_204)
        enddo
        !$acc loop gang vector independent
        do idx_205 = 1, n_pnts, 1
          p(idx_205) = 0.0 + rc(idx_205,176) * y(idx_205,3) * y(idx_205,66)
        enddo
        !$acc loop gang vector independent
        do idx_206 = 1, n_pnts, 1
          l(idx_206) = 0.0 + rc(idx_206,171) * y(idx_206,7) + rc(idx_206,180) * y(idx_206,15)
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
          y(idx_208,42) = (yp(idx_208,42) + dts * p(idx_208)) * tmp1(idx_208)
        enddo
        !$acc loop gang vector independent
        do idx_209 = 1, n_pnts, 1
          p(idx_209) = 0.0 + dj(idx_209,7) * y(idx_209,14) + rc(idx_209,112) * y(idx_209,4) * y(idx_209,51) * 0.13 + rc(idx_209,124) * y(idx_209,4) * y(idx_209,58) * 0.24
        enddo
        !$acc loop gang vector independent
        do idx_210 = 1, n_pnts, 1
          l(idx_210) = 0.0 + rc(idx_210,33) * y(idx_210,3) + dd(idx_210,43)
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
          y(idx_212,43) = (yp(idx_212,43) + dts * p(idx_212)) * tmp1(idx_212)
        enddo
        !$acc loop gang vector independent
        do idx_213 = 1, n_pnts, 1
          p(idx_213) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_214 = 1, n_pnts, 1
          l(idx_214) = 0.0 + rc(idx_214,184) * y(idx_214,3) + rc(idx_214,185) * y(idx_214,3) + rc(idx_214,186) * y(idx_214,6)
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
          y(idx_216,44) = (yp(idx_216,44) + dts * p(idx_216)) * tmp1(idx_216)
        enddo
        !$acc loop gang vector independent
        do idx_217 = 1, n_pnts, 1
          p(idx_217) = 0.0 + rc(idx_217,186) * y(idx_217,44) * y(idx_217,6) + rc(idx_217,187) * y(idx_217,48) * y(idx_217,3) * 0.6 + rc(idx_217,184) * y(idx_217,44) * y(idx_217,3) + rc(idx_217,185) * y(idx_217,44) * y(idx_217,3) * 0.6
        enddo
        !$acc loop gang vector independent
        do idx_218 = 1, n_pnts, 1
          l(idx_218) = 0.0 + rc(idx_218,197) * y(idx_218,4) + dw(idx_218,45) + dd(idx_218,45) + rc(idx_218,188) * y(idx_218,3) + rc(idx_218,190) * y(idx_218,4) + rc(idx_218,191) * y(idx_218,11)
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
          y(idx_220,45) = (yp(idx_220,45) + dts * p(idx_220)) * tmp1(idx_220)
        enddo
        !$acc loop gang vector independent
        do idx_221 = 1, n_pnts, 1
          p(idx_221) = 0.0 + rc(idx_221,188) * y(idx_221,45) * y(idx_221,3)
        enddo
        !$acc loop gang vector independent
        do idx_222 = 1, n_pnts, 1
          l(idx_222) = 0.0
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
          y(idx_224,46) = (yp(idx_224,46) + dts * p(idx_224)) * tmp1(idx_224)
        enddo
        !$acc loop gang vector independent
        do idx_225 = 1, n_pnts, 1
          p(idx_225) = 0.0 + rc(idx_225,187) * y(idx_225,48) * y(idx_225,3) * 0.4
        enddo
        !$acc loop gang vector independent
        do idx_226 = 1, n_pnts, 1
          l(idx_226) = 0.0
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
          y(idx_228,47) = (yp(idx_228,47) + dts * p(idx_228)) * tmp1(idx_228)
        enddo
        !$acc loop gang vector independent
        do idx_229 = 1, n_pnts, 1
          p(idx_229) = 0.0 + rc(idx_229,185) * y(idx_229,44) * y(idx_229,3) * 0.4
        enddo
        !$acc loop gang vector independent
        do idx_230 = 1, n_pnts, 1
          l(idx_230) = 0.0 + rc(idx_230,187) * y(idx_230,3) + dw(idx_230,48)
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
          y(idx_232,48) = (yp(idx_232,48) + dts * p(idx_232)) * tmp1(idx_232)
        enddo
        !$acc loop gang vector independent
        do idx_233 = 1, n_pnts, 1
          p(idx_233) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_234 = 1, n_pnts, 1
          l(idx_234) = 0.0 + rc(idx_234,28) * y(idx_234,3) + dw(idx_234,49) + dd(idx_234,49)
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
          y(idx_236,49) = (yp(idx_236,49) + dts * p(idx_236)) * tmp1(idx_236)
        enddo
        !$acc loop gang vector independent
        do idx_237 = 1, n_pnts, 1
          p(idx_237) = 0.0 + rc(idx_237,155) * y(idx_237,6) * y(idx_237,58)
        enddo
        !$acc loop gang vector independent
        do idx_238 = 1, n_pnts, 1
          l(idx_238) = 0.0 + rc(idx_238,156) * y(idx_238,3)
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
          y(idx_240,50) = (yp(idx_240,50) + dts * p(idx_240)) * tmp1(idx_240)
        enddo
        !$acc loop gang vector independent
        do idx_241 = 1, n_pnts, 1
          p(idx_241) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_242 = 1, n_pnts, 1
          l(idx_242) = 0.0 + rc(idx_242,109) * y(idx_242,3) + rc(idx_242,112) * y(idx_242,4) + rc(idx_242,153) * y(idx_242,6)
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
          y(idx_244,51) = (yp(idx_244,51) + dts * p(idx_244)) * tmp1(idx_244)
        enddo
        !$acc loop gang vector independent
        do idx_245 = 1, n_pnts, 1
          p(idx_245) = 0.0 + rc(idx_245,172) * y(idx_245,3) * y(idx_245,64)
        enddo
        !$acc loop gang vector independent
        do idx_246 = 1, n_pnts, 1
          l(idx_246) = 0.0 + rc(idx_246,173) * y(idx_246,5) + rc(idx_246,181) * y(idx_246,16)
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
          y(idx_248,52) = (yp(idx_248,52) + dts * p(idx_248)) * tmp1(idx_248)
        enddo
        !$acc loop gang vector independent
        do idx_249 = 1, n_pnts, 1
          p(idx_249) = 0.0 + rc(idx_249,153) * y(idx_249,6) * y(idx_249,51)
        enddo
        !$acc loop gang vector independent
        do idx_250 = 1, n_pnts, 1
          l(idx_250) = 0.0 + rc(idx_250,154) * y(idx_250,3)
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
          y(idx_252,53) = (yp(idx_252,53) + dts * p(idx_252)) * tmp1(idx_252)
        enddo
        !$acc loop gang vector independent
        do idx_253 = 1, n_pnts, 1
          p(idx_253) = 0.0 + rc(idx_253,141) * y(idx_253,3) * y(idx_253,29)
        enddo
        !$acc loop gang vector independent
        do idx_254 = 1, n_pnts, 1
          l(idx_254) = 0.0 + rc(idx_254,132) * y(idx_254,16) + rc(idx_254,134) * y(idx_254,15) + rc(idx_254,142) * y(idx_254,5)
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
          y(idx_256,54) = (yp(idx_256,54) + dts * p(idx_256)) * tmp1(idx_256)
        enddo
        !$acc loop gang vector independent
        do idx_257 = 1, n_pnts, 1
          p(idx_257) = 0.0 + rc(idx_257,143) * y(idx_257,3) * y(idx_257,34)
        enddo
        !$acc loop gang vector independent
        do idx_258 = 1, n_pnts, 1
          l(idx_258) = 0.0 + rc(idx_258,133) * y(idx_258,16) + rc(idx_258,136) * y(idx_258,15) + rc(idx_258,144) * y(idx_258,5)
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
          y(idx_260,55) = (yp(idx_260,55) + dts * p(idx_260)) * tmp1(idx_260)
        enddo
        !$acc loop gang vector independent
        do idx_261 = 1, n_pnts, 1
          p(idx_261) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_262 = 1, n_pnts, 1
          l(idx_262) = 0.0 + rc(idx_262,192) * y(idx_262,3) + rc(idx_262,193) * y(idx_262,4) + rc(idx_262,194) * y(idx_262,6)
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
          y(idx_264,56) = (yp(idx_264,56) + dts * p(idx_264)) * tmp1(idx_264)
        enddo
        !$acc loop gang vector independent
        do idx_265 = 1, n_pnts, 1
          p(idx_265) = 0.0 + rc(idx_265,194) * y(idx_265,56) * y(idx_265,6) * 0.14 + rc(idx_265,192) * y(idx_265,56) * y(idx_265,3) * 0.14 + rc(idx_265,193) * y(idx_265,56) * y(idx_265,4) * 0.14
        enddo
        !$acc loop gang vector independent
        do idx_266 = 1, n_pnts, 1
          l(idx_266) = 0.0 + dw(idx_266,57)
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
          y(idx_268,57) = (yp(idx_268,57) + dts * p(idx_268)) * tmp1(idx_268)
        enddo
        !$acc loop gang vector independent
        do idx_269 = 1, n_pnts, 1
          p(idx_269) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_270 = 1, n_pnts, 1
          l(idx_270) = 0.0 + rc(idx_270,155) * y(idx_270,6) + rc(idx_270,123) * y(idx_270,4) + rc(idx_270,124) * y(idx_270,4) + rc(idx_270,125) * y(idx_270,3)
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
          y(idx_272,58) = (yp(idx_272,58) + dts * p(idx_272)) * tmp1(idx_272)
        enddo
        !$acc loop gang vector independent
        do idx_273 = 1, n_pnts, 1
          p(idx_273) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_274 = 1, n_pnts, 1
          l(idx_274) = 0.0 + rc(idx_274,81) * y(idx_274,3) + rc(idx_274,152) * y(idx_274,6)
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
          y(idx_276,59) = (yp(idx_276,59) + dts * p(idx_276)) * tmp1(idx_276)
        enddo
        !$acc loop gang vector independent
        do idx_277 = 1, n_pnts, 1
          p(idx_277) = 0.0 + rc(idx_277,103) * y(idx_277,15) * y(idx_277,24)
        enddo
        !$acc loop gang vector independent
        do idx_278 = 1, n_pnts, 1
          l(idx_278) = 0.0 + dw(idx_278,60) + dd(idx_278,60) + rc(idx_278,44) * y(idx_278,3) + rc(idx_278,104) * y(idx_278,3) + dj(idx_278,21)
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
          y(idx_280,60) = (yp(idx_280,60) + dts * p(idx_280)) * tmp1(idx_280)
        enddo
        !$acc loop gang vector independent
        do idx_281 = 1, n_pnts, 1
          p(idx_281) = 0.0 + rc(idx_281,104) * y(idx_281,3) * y(idx_281,60) + dj(idx_281,21) * y(idx_281,60) + rc(idx_281,83) * y(idx_281,24) * y(idx_281,5) + rc(idx_281,84) * y(idx_281,24) * y(idx_281,16)
        enddo
        !$acc loop gang vector independent
        do idx_282 = 1, n_pnts, 1
          l(idx_282) = 0.0 + rc(idx_282,86) * y(idx_282,3) + dj(idx_282,9)
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
          y(idx_284,61) = (yp(idx_284,61) + dts * p(idx_284)) * tmp1(idx_284)
        enddo
        !$acc loop gang vector independent
        do idx_285 = 1, n_pnts, 1
          p(idx_285) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_286 = 1, n_pnts, 1
          l(idx_286) = 0.0 + rc(idx_286,174) * y(idx_286,3) + rc(idx_286,175) * y(idx_286,3)
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
          y(idx_288,62) = (yp(idx_288,62) + dts * p(idx_288)) * tmp1(idx_288)
        enddo
        !$acc loop gang vector independent
        do idx_289 = 1, n_pnts, 1
          p(idx_289) = 0.0 + rc(idx_289,175) * y(idx_289,3) * y(idx_289,62)
        enddo
        !$acc loop gang vector independent
        do idx_290 = 1, n_pnts, 1
          l(idx_290) = 0.0 + rc(idx_290,182) * y(idx_290,15) + rc(idx_290,183) * y(idx_290,7)
        enddo
        !$acc loop gang vector independent
        do idx_291 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_291) = 1.0 + dts * l(idx_291)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_292 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_292,63) = (yp(idx_292,63) + dts * p(idx_292)) * tmp1(idx_292)
        enddo
        !$acc loop gang vector independent
        do idx_293 = 1, n_pnts, 1
          p(idx_293) = 0.0 + rc(idx_293,180) * y(idx_293,42) * y(idx_293,15) + rc(idx_293,182) * y(idx_293,15) * y(idx_293,63) + rc(idx_293,177) * y(idx_293,3) * y(idx_293,39) + rc(idx_293,178) * y(idx_293,6) * y(idx_293,39) + rc(idx_293,170) * y(idx_293,3) * y(idx_293,66) * 0.80 + rc(idx_293,174) * y(idx_293,3) * y(idx_293,62)
        enddo
        !$acc loop gang vector independent
        do idx_294 = 1, n_pnts, 1
          l(idx_294) = 0.0 + rc(idx_294,172) * y(idx_294,3)
        enddo
        !$acc loop gang vector independent
        do idx_295 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_295) = 1.0 + dts * l(idx_295)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_296 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_296,64) = (yp(idx_296,64) + dts * p(idx_296)) * tmp1(idx_296)
        enddo
        !$acc loop gang vector independent
        do idx_297 = 1, n_pnts, 1
          p(idx_297) = 0.0 + rc(idx_297,181) * y(idx_297,52) * y(idx_297,16) + rc(idx_297,182) * y(idx_297,15) * y(idx_297,63) + rc(idx_297,177) * y(idx_297,3) * y(idx_297,39) + rc(idx_297,178) * y(idx_297,6) * y(idx_297,39) + rc(idx_297,173) * y(idx_297,52) * y(idx_297,5) + rc(idx_297,174) * y(idx_297,3) * y(idx_297,62)
        enddo
        !$acc loop gang vector independent
        do idx_298 = 1, n_pnts, 1
          l(idx_298) = 0.0 + rc(idx_298,139) * y(idx_298,3) + dj(idx_298,13) + dw(idx_298,65)
        enddo
        !$acc loop gang vector independent
        do idx_299 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_299) = 1.0 + dts * l(idx_299)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_300 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_300,65) = (yp(idx_300,65) + dts * p(idx_300)) * tmp1(idx_300)
        enddo
        !$acc loop gang vector independent
        do idx_301 = 1, n_pnts, 1
          p(idx_301) = 0.0
        enddo
        !$acc loop gang vector independent
        do idx_302 = 1, n_pnts, 1
          l(idx_302) = 0.0 + rc(idx_302,170) * y(idx_302,3) + rc(idx_302,176) * y(idx_302,3)
        enddo
        !$acc loop gang vector independent
        do idx_303 = LBOUND(tmp2, dim=1), UBOUND(tmp2, dim=1), 1
          tmp2(idx_303) = 1.0 + dts * l(idx_303)
        enddo
        !$acc end kernels
        call oneover_v(n_pnts, tmp2, tmp1)
        !$acc kernels
        !$acc loop gang vector independent
        do idx_304 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
          y(idx_304,66) = (yp(idx_304,66) + dts * p(idx_304)) * tmp1(idx_304)
        enddo
        !$acc end kernels
      enddo
      !$acc kernels
      !$acc loop gang vector independent
      do idx_305 = LBOUND(delta_so2_wetox_h2o2, dim=1), UBOUND(delta_so2_wetox_h2o2, dim=1), 1
        delta_so2_wetox_h2o2(idx_305) = delta_so2_wetox_h2o2(idx_305) + rc(idx_305,191) * y(idx_305,11) * y(idx_305,45) * dts
      enddo
      !$acc loop gang vector independent
      do idx_306 = LBOUND(delta_so2_wetox_o3, dim=1), UBOUND(delta_so2_wetox_o3, dim=1), 1
        delta_so2_wetox_o3(idx_306) = delta_so2_wetox_o3(idx_306) + rc(idx_306,190) * y(idx_306,4) * y(idx_306,45) * dts + rc(idx_306,197) * y(idx_306,4) * y(idx_306,45) * dts
      enddo
      !$acc loop gang vector independent
      do idx_307 = LBOUND(delta_so2_dryox_oh, dim=1), UBOUND(delta_so2_dryox_oh, dim=1), 1
        delta_so2_dryox_oh(idx_307) = delta_so2_dryox_oh(idx_307) + rc(idx_307,188) * y(idx_307,3) * y(idx_307,45) * dts
      enddo
      !$acc end kernels
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_deriv_raqaero

end module ukca_deriv_raqaero_mod
