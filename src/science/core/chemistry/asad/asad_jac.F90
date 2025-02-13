module asad_jac_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_JAC_MOD'
  public

  contains
  subroutine asad_jac(n_points)
    use asad_mod, only : ctype, dpd, dpw, ej, f, jpcspf, jpfm, jpif, jpmsp, jpnr, jpspec, linfam, madvtr, moffam, ndepd, ndepw, njacx1, njacx2, njacx3, njcgrp, nltr3, nltrf, nmpjac, nodd, npjac1, nspi, ntr3, ntrf, peps, prk, y
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, dimension(jpmsp) :: ifam
    integer, dimension(jpmsp) :: itr
    integer :: j
    integer :: jc
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: jl
    integer :: jp
    integer :: jr
    integer :: js
    integer :: jtr
    integer :: ip1
    integer :: ip2
    integer :: ip3
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: ir1
    integer :: ir2
    integer :: itrcr
    integer :: ks
    integer :: ifamd
    integer :: itrd
    logical, dimension(jpmsp) :: gif
    logical, dimension(jpmsp) :: gtype
    logical :: gir1
    logical :: gir2
    logical :: ginclude2
    logical :: ginclude3
    logical :: ginclude4
    logical :: ginclude5
    CHARACTER :: ityped*2
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_JAC'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do jtr = 1, jpcspf, 1
      do jl = 1, n_points, 1
        ej(jl,jtr) = 0.0
      enddo
    enddo
    !$acc loop seq
    do jc = 1, ntrf, 1
      itrcr = nltrf(jc)
      ip3 = njcgrp(itrcr,3)
      ip2 = njcgrp(itrcr,2)
      ip1 = njcgrp(itrcr,1)
      !$acc loop seq
      do j3 = 1, ip3, 1
        i1 = njacx3(1,j3,itrcr)
        i2 = njacx3(2,j3,itrcr)
        i3 = njacx3(3,j3,itrcr)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          ej(jl,itrcr) = ej(jl,itrcr) - prk(jl,i1) - prk(jl,i2) - prk(jl,i3)
        enddo
      enddo
      !$acc loop seq
      do j2 = 1, ip2, 1
        i1 = njacx2(1,itrcr)
        i2 = njacx2(2,itrcr)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          ej(jl,itrcr) = ej(jl,itrcr) - prk(jl,i1) - prk(jl,i2)
        enddo
      enddo
      !$acc loop seq
      do j1 = 1, ip1, 1
        i1 = njacx1(itrcr)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          ej(jl,itrcr) = ej(jl,itrcr) - prk(jl,i1)
        enddo
      enddo
      !$acc loop seq
      do jp = 1, nmpjac(itrcr), 1
        i1 = npjac1(jp,itrcr)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          ej(jl,itrcr) = ej(jl,itrcr) + prk(jl,i1)
        enddo
      enddo
    enddo
    !$acc end kernels
    if (ntr3 /= 0) then
      do jr = 1, jpnr, 1
        !$acc kernels
        ir1 = nspi(jr,1)
        ir2 = nspi(jr,2)
        ip1 = nspi(jr,3)
        ip2 = nspi(jr,4)
        ip3 = nspi(jr,5)
        !$acc loop gang vector independent
        do j = 1, jpmsp, 1
          gif(j) = .false.
          gtype(j) = .false.
          ifam(j) = 0
          itr(j) = 0
        enddo
        !$acc end kernels
        if (ir1 /= 0) then
          ifam(1) = moffam(ir1)
          itr(1) = madvtr(ir1)
        end if
        if (ir2 /= 0) then
          ifam(2) = moffam(ir2)
          itr(2) = madvtr(ir2)
        end if
        !$acc kernels
        gir1 = .false.
        gir2 = .false.
        !$acc loop seq
        do j3 = 1, ntr3, 1
          gir1 = gir1 .OR. ifam(1) == nltr3(j3)
          gir2 = gir2 .OR. ifam(2) == nltr3(j3)
        enddo
        !$acc end kernels
        if (gir1 .OR. gir2) then
          if (ir1 /= 0) then
            gtype(1) = ctype(ir1) == jpfm(:)
            gif(1) = ctype(ir1) == jpif(:)
          end if
          if (ir2 /= 0) then
            gtype(2) = ctype(ir2) == jpfm(:)
            gif(2) = ctype(ir2) == jpif(:)
          end if
          if (ip1 /= 0) then
            ifam(3) = moffam(ip1)
            gtype(3) = ctype(ip1) == jpfm(:)
            gif(3) = ctype(ip1) == jpif(:)
            itr(3) = madvtr(ip1)
          end if
          if (ip2 /= 0) then
            ifam(4) = moffam(ip2)
            gtype(4) = ctype(ip2) == jpfm(:)
            gif(4) = ctype(ip2) == jpif(:)
            itr(4) = madvtr(ip2)
          end if
          if (ip3 /= 0) then
            ifam(5) = moffam(ip3)
            gtype(5) = ctype(ip3) == jpfm(:)
            gif(5) = ctype(ip3) == jpif(:)
            itr(5) = madvtr(ip3)
          end if
        end if
        if (gir1) then
          !$acc kernels
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            ginclude2 = gif(2) .AND. linfam(jl,itr(2)) .OR. gtype(2)
            ginclude3 = gif(3) .AND. linfam(jl,itr(3)) .OR. gtype(3)
            ginclude4 = gif(4) .AND. linfam(jl,itr(4)) .OR. gtype(4)
            ginclude5 = gif(5) .AND. linfam(jl,itr(5)) .OR. gtype(5)
            if (gif(1) .AND. (.NOT.linfam(jl,itr(1)))) then
              ks = 0
            else
              ks = -nodd(ir1)
              if (ifam(2) == ifam(1) .AND. ginclude2) then
                ks = ks - nodd(ir2)
              end if
              if (ifam(3) == ifam(1) .AND. ginclude3) then
                ks = ks + nodd(ip1)
              end if
              if (ifam(4) == ifam(1) .AND. ginclude4) then
                ks = ks + nodd(ip2)
              end if
              if (ifam(5) == ifam(1) .AND. ginclude5) then
                ks = ks + nodd(ip3)
              end if
            end if
            if (ks /= 0) then
              ej(jl,ifam(1)) = ej(jl,ifam(1)) + prk(jl,jr) * ks
            end if
          enddo
          !$acc end kernels
        end if
        if (gir2) then
          !$acc kernels
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            ginclude2 = gif(2) .AND. linfam(jl,itr(2)) .OR. gtype(2)
            ginclude3 = gif(3) .AND. linfam(jl,itr(3)) .OR. gtype(3)
            ginclude4 = gif(4) .AND. linfam(jl,itr(4)) .OR. gtype(4)
            ginclude5 = gif(5) .AND. linfam(jl,itr(5)) .OR. gtype(5)
            if (gif(2) .AND. (.NOT.linfam(jl,itr(2)))) then
              ks = 0
            else
              ks = -nodd(ir2)
              if (ifam(1) == ifam(2) .AND. ginclude2) then
                ks = ks - nodd(ir1)
              end if
              if (ifam(3) == ifam(2) .AND. ginclude3) then
                ks = ks + nodd(ip1)
              end if
              if (ifam(4) == ifam(2) .AND. ginclude4) then
                ks = ks + nodd(ip2)
              end if
              if (ifam(5) == ifam(2) .AND. ginclude5) then
                ks = ks + nodd(ip3)
              end if
            end if
            if (ks /= 0) then
              ej(jl,ifam(2)) = ej(jl,ifam(2)) + prk(jl,jr) * ks
            end if
          enddo
          !$acc end kernels
        end if
      enddo
    end if
    if (ndepw /= 0 .OR. ndepd /= 0) then
      do js = 1, jpspec, 1
        ifamd = moffam(js)
        itrd = madvtr(js)
        ityped = ctype(js)
        if (ifamd /= 0) then
          !$acc kernels
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            if (ityped == jpfm(:) .OR. ityped == jpif(:) .AND. linfam(jl,itrd)) then
              ej(jl,ifamd) = ej(jl,ifamd) - nodd(js) * (dpd(jl,js) + dpw(jl,js)) * y(jl,js)
            end if
          enddo
          !$acc end kernels
        end if
        if (itrd /= 0) then
          !$acc kernels
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            ej(jl,itrd) = ej(jl,itrd) - (dpd(jl,js) + dpw(jl,js)) * y(jl,js)
          enddo
          !$acc end kernels
        end if
      enddo
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do jtr = 1, jpcspf, 1
      do jl = 1, n_points, 1
        if (f(jl,jtr) > peps) then
          ej(jl,jtr) = ej(jl,jtr) / f(jl,jtr)
        else
          ej(jl,jtr) = 0.0
        end if
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_jac

end module asad_jac_mod
