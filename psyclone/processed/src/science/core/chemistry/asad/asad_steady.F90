module asad_steady_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_STEADY_MOD'
  public

  contains
  subroutine asad_steady(kl)
    use asad_mod, only : deriv, h_in_ss, n_in_ss, nsph, nspho2, nspi, nspn, nspno, nspo1d, nspo3, nspo3p, nspoh, nss_h, nss_n, &
&nss_o1d, nss_o3p, nssi, nsspi, nsspt, nssri, nssrt, nssrx, nsst, nuni, o3p_in_ss, peps, rk, y
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    integer, parameter :: n_o3 = 1
    integer, parameter :: n_oh = 2
    integer, parameter :: n_ho2 = 3
    integer, parameter :: n_no = 4
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kl
    integer :: jl1
    integer :: jl2
    integer :: jr
    integer :: ix
    integer :: i
    integer :: j
    real, dimension(kl) :: ssnum
    real, dimension(kl) :: ssden
    real, dimension(kl,4) :: dssnum
    real, dimension(kl,4) :: dssden
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_STEADY'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    jl1 = 1
    jl2 = kl
    if (kl == 1) then
      jl1 = jl1
      jl2 = jl1
    end if
    do idx = LBOUND(deriv, dim=3), UBOUND(deriv, dim=3), 1
      do idx_1 = LBOUND(deriv, dim=2), UBOUND(deriv, dim=2), 1
        do idx_2 = LBOUND(deriv, dim=1), UBOUND(deriv, dim=1), 1
          deriv(idx_2,idx_1,idx) = 1.0
        enddo
      enddo
    enddo
    do ix = 1, nsst, 1
      do idx_3 = 1, kl, 1
        ssnum(idx_3) = 0.0
      enddo
      do idx_4 = 1, kl, 1
        ssden(idx_4) = 0.0
      enddo
      do idx_5 = 1, 4, 1
        do idx_6 = 1, kl, 1
          dssden(idx_6,idx_5) = 0.0
        enddo
      enddo
      do idx_7 = 1, 4, 1
        do idx_8 = 1, kl, 1
          dssnum(idx_8,idx_7) = 0.0
        enddo
      enddo
      do jr = 1, nsspt(ix), 1
        i = nsspi(ix,jr)
        if (i <= nuni) then
          do idx_9 = jl1, jl2, 1
            ssnum(idx_9) = ssnum(idx_9) + rk(idx_9,i) * y(idx_9,nspi(i,1))
          enddo
          if (ix < 5 .AND. nspi(i,1) == nspo3) then
            do idx_10 = jl1, jl2, 1
              dssnum(idx_10,n_o3) = dssnum(idx_10,n_o3) + rk(idx_10,i)
            enddo
          end if
          if (ix < 5 .AND. nspi(i,1) == nspno) then
            do idx_11 = jl1, jl2, 1
              dssnum(idx_11,n_no) = dssnum(idx_11,n_no) + rk(idx_11,i)
            enddo
          end if
        else
          do idx_12 = jl1, jl2, 1
            ssnum(idx_12) = ssnum(idx_12) + rk(idx_12,i) * y(idx_12,nspi(i,1)) * y(idx_12,nspi(i,2))
          enddo
          if (ix < 5) then
            if (nspi(i,1) == nspo1d) then
              do idx_13 = jl1, jl2, 1
                dssnum(idx_13,n_o3) = dssnum(idx_13,n_o3) + rk(idx_13,i) * y(idx_13,nspi(i,2)) * deriv(idx_13,nss_o1d,n_o3)
              enddo
            end if
            if (nspi(i,2) == nspo1d) then
              do idx_14 = jl1, jl2, 1
                dssnum(idx_14,n_o3) = dssnum(idx_14,n_o3) + rk(idx_14,i) * y(idx_14,nspi(i,1)) * deriv(idx_14,nss_o1d,n_o3)
              enddo
            end if
            if (nspi(i,1) == nspo3p) then
              do idx_15 = jl1, jl2, 1
                dssnum(idx_15,n_o3) = dssnum(idx_15,n_o3) + rk(idx_15,i) * y(idx_15,nspi(i,2)) * deriv(idx_15,nss_o3p,n_o3)
              enddo
            end if
            if (nspi(i,2) == nspo3p) then
              do idx_16 = jl1, jl2, 1
                dssnum(idx_16,n_o3) = dssnum(idx_16,n_o3) + rk(idx_16,i) * y(idx_16,nspi(i,1)) * deriv(idx_16,nss_o3p,n_o3)
              enddo
            end if
            if (nspi(i,1) == nspo3) then
              do idx_17 = jl1, jl2, 1
                dssnum(idx_17,n_o3) = dssnum(idx_17,n_o3) + rk(idx_17,i) * y(idx_17,nspi(i,2))
              enddo
            end if
            if (nspi(i,2) == nspo3) then
              do idx_18 = jl1, jl2, 1
                dssnum(idx_18,n_o3) = dssnum(idx_18,n_o3) + rk(idx_18,i) * y(idx_18,nspi(i,1))
              enddo
            end if
            if (nspi(i,1) == nspo3p) then
              do idx_19 = jl1, jl2, 1
                dssnum(idx_19,n_oh) = dssnum(idx_19,n_oh) + rk(idx_19,i) * y(idx_19,nspi(i,2)) * deriv(idx_19,nss_o3p,n_oh)
              enddo
            end if
            if (nspi(i,2) == nspo3p) then
              do idx_20 = jl1, jl2, 1
                dssnum(idx_20,n_oh) = dssnum(idx_20,n_oh) + rk(idx_20,i) * y(idx_20,nspi(i,1)) * deriv(idx_20,nss_o3p,n_oh)
              enddo
            end if
            if (nspi(i,1) == nspoh) then
              do idx_21 = jl1, jl2, 1
                dssnum(idx_21,n_oh) = dssnum(idx_21,n_oh) + rk(idx_21,i) * y(idx_21,nspi(i,2))
              enddo
            end if
            if (nspi(i,2) == nspoh) then
              do idx_22 = jl1, jl2, 1
                dssnum(idx_22,n_oh) = dssnum(idx_22,n_oh) + rk(idx_22,i) * y(idx_22,nspi(i,1))
              enddo
            end if
            if (nspi(i,1) == nspo3p) then
              do idx_23 = jl1, jl2, 1
                dssnum(idx_23,n_ho2) = dssnum(idx_23,n_ho2) + rk(idx_23,i) * y(idx_23,nspi(i,2)) * deriv(idx_23,nss_o3p,n_ho2)
              enddo
            end if
            if (nspi(i,2) == nspo3p) then
              do idx_24 = jl1, jl2, 1
                dssnum(idx_24,n_ho2) = dssnum(idx_24,n_ho2) + rk(idx_24,i) * y(idx_24,nspi(i,1)) * deriv(idx_24,nss_o3p,n_ho2)
              enddo
            end if
            if (nspi(i,1) == nspho2) then
              do idx_25 = jl1, jl2, 1
                dssnum(idx_25,n_ho2) = dssnum(idx_25,n_ho2) + rk(idx_25,i) * y(idx_25,nspi(i,2))
              enddo
            end if
            if (nspi(i,2) == nspho2) then
              do idx_26 = jl1, jl2, 1
                dssnum(idx_26,n_ho2) = dssnum(idx_26,n_ho2) + rk(idx_26,i) * y(idx_26,nspi(i,1))
              enddo
            end if
            if (nspi(i,1) == nspno) then
              do idx_27 = jl1, jl2, 1
                dssnum(idx_27,n_no) = dssnum(idx_27,n_no) + rk(idx_27,i) * y(idx_27,nspi(i,2))
              enddo
            end if
            if (nspi(i,2) == nspno) then
              do idx_28 = jl1, jl2, 1
                dssnum(idx_28,n_no) = dssnum(idx_28,n_no) + rk(idx_28,i) * y(idx_28,nspi(i,1))
              enddo
            end if
          end if
        end if
      enddo
      do jr = 1, nssrt(ix), 1
        i = nssri(ix,jr)
        j = nssrx(ix,jr)
        if (i <= nuni) then
          do idx_29 = jl1, jl2, 1
            ssden(idx_29) = ssden(idx_29) + rk(idx_29,i)
          enddo
        else
          do idx_30 = jl1, jl2, 1
            ssden(idx_30) = ssden(idx_30) + rk(idx_30,i) * y(idx_30,nspi(i,j))
          enddo
          if (ix < 5) then
            if (nspi(i,j) == nspo3) then
              do idx_31 = jl1, jl2, 1
                dssden(idx_31,n_o3) = dssden(idx_31,n_o3) + rk(idx_31,i)
              enddo
            end if
            if (nspi(i,j) == nspoh) then
              do idx_32 = jl1, jl2, 1
                dssden(idx_32,n_oh) = dssden(idx_32,n_oh) + rk(idx_32,i)
              enddo
            end if
            if (nspi(i,j) == nspho2) then
              do idx_33 = jl1, jl2, 1
                dssden(idx_33,n_ho2) = dssden(idx_33,n_ho2) + rk(idx_33,i)
              enddo
            end if
            if (nspi(i,j) == nspno) then
              do idx_34 = jl1, jl2, 1
                dssden(idx_34,n_no) = dssden(idx_34,n_no) + rk(idx_34,i)
              enddo
            end if
          end if
        end if
      enddo
      do idx_35 = jl1, jl2, 1
        y(idx_35,nssi(ix)) = ssnum(idx_35) / ssden(idx_35)
      enddo
      if (ix < 5) then
        do jr = 1, 4, 1
          do idx_36 = jl1, jl2, 1
            deriv(idx_36,ix,jr) = (ssden(idx_36) * dssnum(idx_36,jr) - ssnum(idx_36) * dssden(idx_36,jr)) / (ssden(idx_36) * &
&ssden(idx_36))
          enddo
        enddo
      end if
    enddo
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (y(jl1 : jl2, nspo1d) > peps)
      deriv(jl1 : jl2, nss_o1d, n_o3) = deriv(jl1 : jl2, nss_o1d, n_o3) * y(jl1 : jl2, nspo3) / y(jl1 : jl2, nspo1d)
      deriv(jl1 : jl2, nss_o1d, n_oh) = deriv(jl1 : jl2, nss_o1d, n_oh) * y(jl1 : jl2, nspoh) / y(jl1 : jl2, nspo1d)
      deriv(jl1 : jl2, nss_o1d, n_ho2) = deriv(jl1 : jl2, nss_o1d, n_ho2) * y(jl1 : jl2, nspho2) / y(jl1 : jl2, nspo1d)
      deriv(jl1 : jl2, nss_o1d, n_no) = deriv(jl1 : jl2, nss_o1d, n_no) * y(jl1 : jl2, nspno) / y(jl1 : jl2, nspo1d)
    ELSEWHERE
      deriv(jl1 : jl2, nss_o1d, n_o3) = 1.0
      deriv(jl1 : jl2, nss_o1d, n_oh) = 1.0
      deriv(jl1 : jl2, nss_o1d, n_ho2) = 1.0
      deriv(jl1 : jl2, nss_o1d, n_no) = 1.0
    END WHERE
    if (o3p_in_ss) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (y(jl1 : jl2, nspo3p) > peps)
        deriv(jl1 : jl2, nss_o3p, n_o3) = deriv(jl1 : jl2, nss_o3p, n_o3) * y(jl1 : jl2, nspo3) / y(jl1 : jl2, nspo3p)
        deriv(jl1 : jl2, nss_o3p, n_oh) = deriv(jl1 : jl2, nss_o3p, n_oh) * y(jl1 : jl2, nspoh) / y(jl1 : jl2, nspo3p)
        deriv(jl1 : jl2, nss_o3p, n_ho2) = deriv(jl1 : jl2, nss_o3p, n_ho2) * y(jl1 : jl2, nspho2) / y(jl1 : jl2, nspo3p)
        deriv(jl1 : jl2, nss_o3p, n_no) = deriv(jl1 : jl2, nss_o3p, n_no) * y(jl1 : jl2, nspno) / y(jl1 : jl2, nspo3p)
      ELSEWHERE
        deriv(jl1 : jl2, nss_o3p, n_o3) = 1.0
        deriv(jl1 : jl2, nss_o3p, n_oh) = 1.0
        deriv(jl1 : jl2, nss_o3p, n_ho2) = 1.0
        deriv(jl1 : jl2, nss_o3p, n_no) = 1.0
      END WHERE
    end if
    if (n_in_ss) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (y(jl1 : jl2, nspn) > peps)
        deriv(jl1 : jl2, nss_n, n_o3) = deriv(jl1 : jl2, nss_n, n_o3) * y(jl1 : jl2, nspo3) / y(jl1 : jl2, nspn)
        deriv(jl1 : jl2, nss_n, n_oh) = deriv(jl1 : jl2, nss_n, n_oh) * y(jl1 : jl2, nspoh) / y(jl1 : jl2, nspn)
        deriv(jl1 : jl2, nss_n, n_ho2) = deriv(jl1 : jl2, nss_n, n_ho2) * y(jl1 : jl2, nspho2) / y(jl1 : jl2, nspn)
        deriv(jl1 : jl2, nss_n, n_no) = deriv(jl1 : jl2, nss_n, n_no) * y(jl1 : jl2, nspno) / y(jl1 : jl2, nspn)
      ELSEWHERE
        deriv(jl1 : jl2, nss_n, n_o3) = 1.0
        deriv(jl1 : jl2, nss_n, n_oh) = 1.0
        deriv(jl1 : jl2, nss_n, n_ho2) = 1.0
        deriv(jl1 : jl2, nss_n, n_no) = 1.0
      END WHERE
    end if
    if (h_in_ss) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (y(jl1 : jl2, nsph) > peps)
        deriv(jl1 : jl2, nss_h, n_o3) = deriv(jl1 : jl2, nss_h, n_o3) * y(jl1 : jl2, nspo3) / y(jl1 : jl2, nsph)
        deriv(jl1 : jl2, nss_h, n_oh) = deriv(jl1 : jl2, nss_h, n_oh) * y(jl1 : jl2, nspoh) / y(jl1 : jl2, nsph)
        deriv(jl1 : jl2, nss_h, n_ho2) = deriv(jl1 : jl2, nss_h, n_ho2) * y(jl1 : jl2, nspho2) / y(jl1 : jl2, nsph)
        deriv(jl1 : jl2, nss_h, n_no) = deriv(jl1 : jl2, nss_h, n_no) * y(jl1 : jl2, nspno) / y(jl1 : jl2, nsph)
      ELSEWHERE
        deriv(jl1 : jl2, nss_h, n_o3) = 1.0
        deriv(jl1 : jl2, nss_h, n_oh) = 1.0
        deriv(jl1 : jl2, nss_h, n_ho2) = 1.0
        deriv(jl1 : jl2, nss_h, n_no) = 1.0
      END WHERE
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine asad_steady

end module asad_steady_mod
