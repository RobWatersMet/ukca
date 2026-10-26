module asad_prls_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_PRLS_MOD'
  public

  contains
  subroutine asad_prls(kl, knspec, kspec, ldepem)
    use asad_mod, only : dpd, dpw, frpx, jpnr, jpspec, ngrp, nldepx, nnfrp, nprdx1, nprdx2, nprdx3, nspi, ntabfp, nuni, pd, prk, &
&prod, rk, slos, y
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kl
    integer, intent(in) :: knspec
    integer, dimension(knspec), intent(in) :: kspec
    logical, intent(in) :: ldepem
    integer :: j
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: jl
    integer :: jr
    integer :: js
    integer :: jskip
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: ir1
    integer :: ir2
    integer :: ip1
    integer :: ip2
    integer :: ip3
    integer :: iunip1
    integer :: iss
    integer :: ir
    integer :: ix
    integer :: istart
    integer :: iend
    real :: fr
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_PRLS'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    do js = 1, jpspec, 1
      do jl = 1, kl, 1
        prod(jl,js) = 0.0
        slos(jl,js) = 0.0
      enddo
    enddo
    do jr = 1, nuni, 1
      ir1 = nspi(jr,1)
      do jl = 1, kl, 1
        prk(jl,jr) = rk(jl,jr) * y(jl,ir1)
      enddo
    enddo
    iunip1 = nuni + 1
    do jr = iunip1, jpnr, 1
      ir1 = nspi(jr,1)
      ir2 = nspi(jr,2)
      do jl = 1, kl, 1
        prk(jl,jr) = rk(jl,jr) * y(jl,ir1) * y(jl,ir2)
      enddo
    enddo
    do jskip = 0, jpspec, jpspec
      do j = 1, knspec, 1
        js = kspec(j) + jskip
        ip3 = ngrp(js,3)
        ip2 = ngrp(js,2)
        ip1 = ngrp(js,1)
        do j3 = 1, ip3, 1
          i1 = nprdx3(1,j3,js)
          i2 = nprdx3(2,j3,js)
          i3 = nprdx3(3,j3,js)
          do jl = 1, kl, 1
            pd(jl,js) = pd(jl,js) + prk(jl,i1) + prk(jl,i2) + prk(jl,i3)
          enddo
        enddo
        do j2 = 1, ip2, 1
          i1 = nprdx2(1,js)
          i2 = nprdx2(2,js)
          do jl = 1, kl, 1
            pd(jl,js) = pd(jl,js) + prk(jl,i1) + prk(jl,i2)
          enddo
        enddo
        do j1 = 1, ip1, 1
          i1 = nprdx1(js)
          do jl = 1, kl, 1
            pd(jl,js) = pd(jl,js) + prk(jl,i1)
          enddo
        enddo
      enddo
    enddo
    do jr = 1, nnfrp, 1
      iss = ntabfp(jr,1)
      ir = ntabfp(jr,2)
      ix = ntabfp(jr,3)
      fr = frpx(ix)
      do jl = 1, kl, 1
        pd(jl,iss) = pd(jl,iss) + fr * prk(jl,ir)
      enddo
    enddo
    if (ldepem) then
      istart = nldepx(1)
      iend = nldepx(2)
      do j = istart, iend, 1
        js = nldepx(j)
        do jl = 1, kl, 1
          slos(jl,js) = slos(jl,js) + (dpd(jl,js) + dpw(jl,js)) * y(jl,js)
        enddo
      enddo
      istart = nldepx(3)
      iend = nldepx(4)
      do j = istart, iend, 1
        js = nldepx(j)
        do jl = 1, kl, 1
          slos(jl,js) = slos(jl,js) + dpd(jl,js) * y(jl,js)
        enddo
      enddo
      istart = nldepx(5)
      iend = nldepx(6)
      do j = istart, iend, 1
        js = nldepx(j)
        do jl = 1, kl, 1
          slos(jl,js) = slos(jl,js) + dpw(jl,js) * y(jl,js)
        enddo
      enddo
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine asad_prls

end module asad_prls_mod
