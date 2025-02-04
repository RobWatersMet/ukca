module asad_fuljac_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_FULJAC_MOD'
  public

  contains
  subroutine asad_fuljac(n_points)
    use asad_mod, only : cdt, ctype, deriv, dpd, dpw, f, fj, frpx, jpcspf, jpfm, jpif, jpmsp, jpspec, linfam, madvtr, moffam, ndepd, ndepw, nfrpx, njcoth, nltrf, nmsjac, nmzjac, nodd, npdfr, nsjac1, nstst, ntabpd, ntrf, ntrho2, ntrno, ntro3, ntroh, nzjac1, prk, y, ztabpd
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real :: deltt
    real :: fr
    integer :: irj
    integer :: ifamd
    integer :: itrd
    integer :: i
    integer :: j
    integer :: jc
    integer :: itrcr
    integer :: j3
    integer :: jn
    integer :: js
    integer :: jl
    integer :: i1
    integer :: i2
    integer :: iss
    CHARACTER :: ityped*2
    integer, dimension(jpmsp) :: ij
    logical, dimension(jpmsp) :: lj
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_FULJAC'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    deltt = 1.0 / cdt
    fj = 0.0
    !$acc loop gang vector independent collapse(2)
    do i = 1, jpcspf, 1
      do idx = LBOUND(fj, dim=1), UBOUND(fj, dim=1), 1
        fj(idx,i,i) = -deltt * f(idx + (LBOUND(f, dim=1) - LBOUND(fj, dim=1)),i)
      enddo
    enddo
    !$acc loop seq
    do jc = 1, ntrf, 1
      itrcr = nltrf(jc)
      !$acc loop seq
      do j3 = 1, nmzjac(itrcr), 1
        irj = nzjac1(j3,itrcr)
        !$acc loop gang vector independent
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        !$acc loop seq
        do jn = 1, 2, 1
          if (lj(jn)) then
            !$acc loop gang vector independent
            do jl = 1, n_points, 1
              fj(jl,ij(jn),itrcr) = fj(jl,ij(jn),itrcr) - prk(jl,irj)
            enddo
          end if
        enddo
        if (nfrpx(irj) == 0) then
          !$acc loop seq
          do jn = 3, jpmsp, 1
            if (lj(jn)) then
              !$acc loop gang vector independent
              do jl = 1, n_points, 1
                fj(jl,ij(jn),itrcr) = fj(jl,ij(jn),itrcr) + prk(jl,irj)
              enddo
            end if
          enddo
        else
          !$acc loop seq
          do jn = 3, jpmsp, 1
            if (lj(jn)) then
              fr = frpx(nfrpx(irj) + jn - 3)
              !$acc loop gang vector independent
              do jl = 1, n_points, 1
                fj(jl,ij(jn),itrcr) = fj(jl,ij(jn),itrcr) + fr * prk(jl,irj)
              enddo
            end if
          enddo
        end if
        if (npdfr(irj,1) /= 0) then
          i1 = npdfr(irj,1)
          i2 = npdfr(irj,2)
          !$acc loop seq
          do jn = i1, i2, 1
            iss = ntabpd(jn,1)
            fr = ztabpd(jn,1)
            !$acc loop gang vector independent
            do jl = 1, n_points, 1
              fj(jl,iss,itrcr) = fj(jl,iss,itrcr) + fr * prk(jl,irj)
            enddo
          enddo
        end if
      enddo
    enddo
    !$acc loop seq
    do jc = 1, nstst, 1
      !$acc loop seq
      do j3 = 1, nmsjac(jc), 1
        irj = nsjac1(j3,jc)
        !$acc loop gang vector independent
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        !$acc loop seq
        do jn = 1, 2, 1
          if (lj(jn)) then
            !$acc loop gang vector independent
            do jl = 1, n_points, 1
              fj(jl,ij(jn),ntro3) = fj(jl,ij(jn),ntro3) - prk(jl,irj) * deriv(jl,jc,1)
              fj(jl,ij(jn),ntroh) = fj(jl,ij(jn),ntroh) - prk(jl,irj) * deriv(jl,jc,2)
              fj(jl,ij(jn),ntrho2) = fj(jl,ij(jn),ntrho2) - prk(jl,irj) * deriv(jl,jc,3)
              fj(jl,ij(jn),ntrno) = fj(jl,ij(jn),ntrno) - prk(jl,irj) * deriv(jl,jc,4)
            enddo
          end if
        enddo
        !$acc loop seq
        do jn = 3, jpmsp, 1
          if (lj(jn)) then
            !$acc loop gang vector independent
            do idx_1 = 1, n_points, 1
              fj(idx_1,ij(jn),ntro3) = fj(idx_1,ij(jn),ntro3) + prk(idx_1,irj) * deriv(idx_1,jc,1)
            enddo
            !$acc loop gang vector independent
            do idx_2 = 1, n_points, 1
              fj(idx_2,ij(jn),ntroh) = fj(idx_2,ij(jn),ntroh) + prk(idx_2,irj) * deriv(idx_2,jc,2)
            enddo
            !$acc loop gang vector independent
            do idx_3 = 1, n_points, 1
              fj(idx_3,ij(jn),ntrho2) = fj(idx_3,ij(jn),ntrho2) + prk(idx_3,irj) * deriv(idx_3,jc,3)
            enddo
            !$acc loop gang vector independent
            do idx_4 = 1, n_points, 1
              fj(idx_4,ij(jn),ntrno) = fj(idx_4,ij(jn),ntrno) + prk(idx_4,irj) * deriv(idx_4,jc,4)
            enddo
          end if
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ndepw /= 0 .OR. ndepd /= 0) then
      !$acc kernels
      !$acc loop seq
      do js = 1, jpspec, 1
        ifamd = moffam(js)
        itrd = madvtr(js)
        ityped = ctype(js)
        if (ifamd /= 0) then
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            if (ityped == jpfm .OR. ityped == jpif .AND. linfam(jl,itrd)) then
              fj(jl,ifamd,ifamd) = fj(jl,ifamd,ifamd) - nodd(js) * (dpd(jl,js) + dpw(jl,js)) * y(jl,js)
            end if
          enddo
        end if
        if (itrd /= 0) then
          !$acc loop gang vector independent
          do idx_5 = 1, n_points, 1
            fj(idx_5,itrd,itrd) = fj(idx_5,itrd,itrd) - (dpd(idx_5,js) + dpw(idx_5,js)) * y(idx_5,js)
          enddo
        end if
      enddo
      !$acc end kernels
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(3)
    do j = 1, jpcspf, 1
      do idx_6 = LBOUND(fj, dim=3), UBOUND(fj, dim=3), 1
        do idx_7 = 1, n_points, 1
          fj(idx_7,j,idx_6) = fj(idx_7,j,idx_6) / f(idx_7,idx_6 + (LBOUND(f, dim=2) - LBOUND(fj, dim=3)))
        enddo
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_fuljac

end module asad_fuljac_mod
