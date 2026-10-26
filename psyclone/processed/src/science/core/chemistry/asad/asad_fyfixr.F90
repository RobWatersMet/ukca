module asad_fyfixr_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_FYFIXR_MOD'
  public

  contains
  subroutine asad_fyfixr(n_points)
    use asad_mod, only : ctype, f, jpif, linfam, madvtr, majors, moffam, nlmajmin, ratio, y
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer :: istart
    integer :: iend
    integer :: j
    integer :: jl
    integer :: js
    integer :: ifam
    integer :: imaj
    integer :: itr
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_FYFIXR'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    istart = nlmajmin(1)
    iend = nlmajmin(2)
    do j = istart, iend, 1
      js = nlmajmin(j)
      ifam = moffam(js)
      do jl = 1, n_points, 1
        y(jl,js) = f(jl,ifam) * ratio(jl,js)
      enddo
    enddo
    istart = nlmajmin(3)
    iend = nlmajmin(4)
    do j = istart, iend, 1
      js = nlmajmin(j)
      ifam = moffam(js)
      imaj = majors(ifam)
      if (ctype(js) /= jpif) then
        do jl = 1, n_points, 1
          y(jl,js) = y(jl,imaj) * ratio(jl,js)
        enddo
      else
        itr = madvtr(js)
        do jl = 1, n_points, 1
          if (linfam(jl,itr)) then
            y(jl,js) = y(jl,imaj) * ratio(jl,js)
          end if
        enddo
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine asad_fyfixr

end module asad_fyfixr_mod
