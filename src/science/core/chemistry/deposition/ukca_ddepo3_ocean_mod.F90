module ukca_ddepo3_ocean_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_DDEPO3_OCEAN_MOD'
  private

  public :: ukca_ddepo3_ocean

  contains
  subroutine ukca_ddepo3_ocean(p0, t0, sst, usa, rc)
    use ukca_um_legacy_mod, only : rgas=>r, vkman
    use ukca_config_constants_mod, only : rhosea
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    real, parameter :: del = 3.0e-6
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    real, intent(in) :: p0
    real, intent(in) :: t0
    real, intent(in) :: sst
    real, intent(in) :: usa
    real, intent(out) :: rc
    real(kind=jprb) :: zhook_handle
    real :: rhoa
    real :: alpha
    real :: usw
    real :: m_diff
    real :: iodide
    real :: k_rate
    real :: vdw
    real :: a1
    real :: b0
    real :: delzt1
    real :: argm
    real :: psi
    real :: trm2
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_DDEPO3_OCEAN'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    rhoa = p0 / (rgas * t0)
    usw = SQRT(rhoa / rhosea) * usa
    b0 = 2.0 / vkman / usw
    alpha = 10.0 ** (-0.25 - 0.013 * (sst - 273.15))
    m_diff = 1.1e-6 * EXP(-1896.0 / sst)
    iodide = 1.46 * 1.0e6 * EXP(-9134.0 / sst)
    k_rate = EXP(-8772.2 / sst + 51.5)
    a1 = k_rate * iodide
    delzt1 = SQRT(2 * a1 * b0 * (del + m_diff * b0 / 2.0))
    argm = SQRT(a1 / m_diff) * del
    trm2 = SQRT(a1 * m_diff)
    psi = delzt1 / b0 / trm2
    vdw = SQRT(a1 * m_diff) * (psi * b_k1(delzt1) * COSH(argm) + b_k0(delzt1) * SINH(argm)) / (psi * b_k1(delzt1) * SINH(argm) + b_k0(delzt1) * COSH(argm))
    rc = 1.0 / (alpha * vdw)
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine ukca_ddepo3_ocean
  function b_k0(x)
    real, parameter :: f1 = -0.57721566
    real, parameter :: f2 = 0.42278420
    real, parameter :: f3 = 0.23069756
    real, parameter :: f4 = 0.3488590e-1
    real, parameter :: f5 = 0.262698e-2
    real, parameter :: f6 = 0.10750e-3
    real, parameter :: f7 = 0.74e-5
    real, parameter :: g1 = 1.25331414
    real, parameter :: g2 = -0.7832358e-1
    real, parameter :: g3 = 0.2189568e-1
    real, parameter :: g4 = -0.1062446e-1
    real, parameter :: g5 = 0.587872e-2
    real, parameter :: g6 = -0.251540e-2
    real, parameter :: g7 = 0.53208e-3
    real, intent(in) :: x
    real :: y
    real :: b_k0

    if (x <= 2.0) then
      y = x * x / 4.0
      b_k0 = -LOG(x / 2.0) * b_i0(x) + (f1 + y * (f2 + y * (f3 + y * (f4 + y * (f5 + y * (f6 + y * f7))))))
    else
      y = 2.0 / x
      b_k0 = (EXP(-x) / SQRT(x)) * (g1 + y * (g2 + y * (g3 + y * (g4 + y * (g5 + y * (g6 + y * g7))))))
    end if
    return

  end function b_k0
  function b_k1(x)
    real, parameter :: f1 = 1.0000
    real, parameter :: f2 = 0.15443144
    real, parameter :: f3 = -0.67278579
    real, parameter :: f4 = -0.18156897
    real, parameter :: f5 = -0.1919402e-1
    real, parameter :: f6 = -0.110404e-2
    real, parameter :: f7 = -0.4686e-4
    real, parameter :: g1 = 1.25331414
    real, parameter :: g2 = 0.23498619
    real, parameter :: g3 = -0.3655620e-1
    real, parameter :: g4 = 0.1504268e-1
    real, parameter :: g5 = -0.780353e-2
    real, parameter :: g6 = 0.325614e-2
    real, parameter :: g7 = -0.68245e-3
    real, intent(in) :: x
    real :: y
    real :: b_k1

    if (x <= 2.0) then
      y = x * x / 4.0
      b_k1 = LOG(x / 2.0) * b_i1(x) + (1.0 / x) * (f1 + y * (f2 + y * (f3 + y * (f4 + y * (f5 + y * (f6 + y * f7))))))
    else
      y = 2.0 / x
      b_k1 = (EXP(-x) / SQRT(x)) * (g1 + y * (g2 + y * (g3 + y * (g4 + y * (g5 + y * (g6 + y * g7))))))
    end if
    return

  end function b_k1
  function b_i0(x)
    real, parameter :: f1 = 1.0000
    real, parameter :: f2 = 3.5156229
    real, parameter :: f3 = 3.0899424
    real, parameter :: f4 = 1.2067492
    real, parameter :: f5 = 0.2659732
    real, parameter :: f6 = 0.360768e-1
    real, parameter :: f7 = 0.45813e-2
    real, parameter :: g1 = 0.39894228
    real, parameter :: g2 = 0.1328592e-1
    real, parameter :: g3 = 0.225319e-2
    real, parameter :: g4 = -0.157565e-2
    real, parameter :: g5 = 0.916281e-2
    real, parameter :: g6 = -0.2057706e-1
    real, parameter :: g7 = 0.2635537e-1
    real, parameter :: g8 = -0.1647633e-1
    real, parameter :: g9 = 0.392377e-2
    real, intent(in) :: x
    real :: xm
    real :: y
    real :: b_i0

    if (ABS(x) < 3.75) then
      y = (x / 3.75) ** 2
      b_i0 = f1 + y * (f2 + y * (f3 + y * (f4 + y * (f5 + y * (f6 + y * f7)))))
    else
      xm = ABS(x)
      y = 3.75 / xm
      b_i0 = (EXP(xm) / SQRT(xm)) * (g1 + y * (g2 + y * (g3 + y * (g4 + y * (g5 + y * (g6 + y * (g7 + y * (g8 + y * g9))))))))
    end if
    return

  end function b_i0
  function b_i1(x)
    real, parameter :: f1 = 0.5000
    real, parameter :: f2 = 0.87890594
    real, parameter :: f3 = 0.51498869
    real, parameter :: f4 = 0.15084934
    real, parameter :: f5 = 0.2658733e-1
    real, parameter :: f6 = 0.301532e-2
    real, parameter :: f7 = 0.32411e-3
    real, parameter :: g1 = 0.39894228
    real, parameter :: g2 = -0.3988024e-1
    real, parameter :: g3 = -0.362018e-2
    real, parameter :: g4 = 0.163801e-2
    real, parameter :: g5 = -0.1031555e-1
    real, parameter :: g6 = 0.2282967e-1
    real, parameter :: g7 = -0.2895312e-1
    real, parameter :: g8 = 0.1787654e-1
    real, parameter :: g9 = -0.420059e-2
    real, intent(in) :: x
    real :: xm
    real :: y
    real :: b_i1

    if (ABS(x) < 3.75) then
      y = (x / 3.75) ** 2
      b_i1 = x * (f1 + y * (f2 + y * (f3 + y * (f4 + y * (f5 + y * (f6 + y * f7))))))
    else
      xm = ABS(x)
      y = 3.75 / xm
      b_i1 = (EXP(xm) / SQRT(xm)) * (g1 + y * (g2 + y * (g3 + y * (g4 + y * (g5 + y * (g6 + y * (g7 + y * (g8 + y * g9))))))))
      if (x < 0.0) then
        b_i1 = -b_i1
      end if
    end if
    return

  end function b_i1

end module ukca_ddepo3_ocean_mod
