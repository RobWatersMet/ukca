module asad_diffun_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_DIFFUN_MOD'
  public

  contains
  subroutine asad_diffun(kl)
    use asad_mod, only : fdot, jpcspf, linfam, madvtr, moffam, nf, nlf, nodd, prod, slos, ydot
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use asad_prls_mod, only : asad_prls
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kl
    integer :: ifam
    integer :: itr
    integer :: j
    integer :: jl
    integer :: jtr
    integer :: js
    logical :: gfam
    logical :: gtr
    logical, save :: gdepem = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_DIFFUN'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(2)
    do jtr = 1, jpcspf, 1
      do jl = 1, kl, 1
        fdot(jl,jtr) = 0.0
      enddo
    enddo
    !NO!$acc end kernels
    call asad_prls(kl, nf, nlf, gdepem)
    !NO!$acc kernels
    !NO!$acc loop seq
    do j = 1, nf, 1
      js = nlf(j)
      ifam = moffam(js)
      itr = madvtr(js)
      gfam = ifam /= 0
      gtr = itr /= 0
      !NO!$acc loop gang vector independent
      do jl = 1, kl, 1
        ydot(jl,js) = prod(jl,js) - slos(jl,js)
        if (gfam .AND. (.NOT.gtr .OR. gtr .AND. linfam(jl,itr))) then
          fdot(jl,ifam) = fdot(jl,ifam) + nodd(js) * ydot(jl,js)
        end if
        if (gtr) then
          fdot(jl,itr) = ydot(jl,js)
        end if
      enddo
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_diffun

end module asad_diffun_mod
