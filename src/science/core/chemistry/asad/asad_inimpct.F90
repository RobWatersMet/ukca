module asad_inimpct_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_INIMPCT_MOD'
  public

  contains
  subroutine asad_inimpct()
    use asad_mod, only : jpcspf, jppj, jpspj, madvtr, moffam, nlpdv, nltrim, nprkx, nspi
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: j
    integer :: ji
    integer :: jii
    integer :: jm
    integer :: jp
    integer :: jt
    integer :: ipr
    integer :: ip
    integer :: ientry
    integer :: ntr
    integer :: ientries
    integer :: indx
    integer :: inf
    integer :: ipos
    integer :: icount
    integer :: ir
    integer :: ir2
    integer :: irk
    integer :: irk2
    integer :: iss
    integer :: is2
    integer :: irtr
    integer :: irtr2
    integer :: inr
    integer :: iprod
    integer, dimension(jpspj - 2) :: itr
    integer, dimension(jpspj - 2) :: ifam
    integer, dimension(0:(jpspj - 2) * jppj,2) :: iwork
    logical, dimension(jpcspf) :: gtr
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INIMPCT'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !$acc kernels
    !$acc loop gang vector independent
    do jt = 1, jpcspf, 1
      gtr(jt) = .false.
    enddo
    !$acc loop gang vector independent
    do jt = 0, jpcspf, 1
      nltrim(jt,1) = 0
      nltrim(jt,2) = 0
      nltrim(jt,3) = 0
    enddo
    !$acc loop gang vector independent
    do jt = 1, (jpspj - 2) * jppj, 1
      iwork(jt,1) = 0
      iwork(jt,2) = 0
      nlpdv(jt,1) = 0
      nlpdv(jt,2) = 0
    enddo
    iwork(0,1) = 0
    iwork(0,2) = 0
    !$acc end kernels
    do jp = 1, jppj, 1
      ipr = nprkx(jp)
      ir = nspi(ipr,1)
      if (madvtr(ir) == 0 .AND. moffam(ir) == 0) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Cycle_Stmt
        CYCLE
      end if
      !$acc kernels
      !$acc loop gang vector independent
      do jm = 3, jpspj, 1
        ip = nspi(ipr,jm)
        if (ip /= 0) then
          ifam(jm - 2) = moffam(ip)
          itr(jm - 2) = madvtr(ip)
        else
          ifam(jm - 2) = 0
          itr(jm - 2) = 0
        end if
      enddo
      !$acc loop seq
      do j = 1, jpspj - 2, 1
        if (itr(j) /= 0) then
          gtr(itr(j)) = .true.
          iwork(0,1) = iwork(0,1) + 1
          ientry = iwork(0,1)
          iwork(ientry,1) = jp
          iwork(ientry,2) = itr(j)
        else
          if (ifam(j) /= 0 .AND. ifam(j) /= moffam(ir)) then
            gtr(ifam(j)) = .true.
            iwork(0,1) = iwork(0,1) + 1
            ientry = iwork(0,1)
            iwork(ientry,1) = jp
            iwork(ientry,2) = ifam(j)
          end if
        end if
      enddo
      !$acc end kernels
    enddo
    !$acc kernels
    do j = 1, jpcspf, 1
      if (gtr(j)) then
        nltrim(0,1) = nltrim(0,1) + 1
        nltrim(nltrim(0,1),1) = j
      end if
    enddo
    ntr = nltrim(0,1)
    ientries = iwork(0,1)
    indx = 1
    !$acc loop seq
    do j = 1, ntr, 1
      jt = nltrim(j,1)
      nltrim(j,3) = indx
      !$acc loop seq
      do ji = 1, ientries, 1
        iprod = iwork(ji,2)
        if (jt == iprod) then
          nlpdv(indx,1) = iwork(ji,1)
          iwork(ji,1) = 0
          iwork(ji,2) = 0
          indx = indx + 1
        end if
      enddo
      nltrim(j,2) = indx - nltrim(j,3)
    enddo
    !$acc end kernels
    do j = 1, ntr, 1
      !$acc kernels
      jt = nltrim(j,1)
      inf = nltrim(j,2)
      ipos = nltrim(j,3)
      icount = 0
      !$acc loop seq
      do ji = 1, inf, 1
        ir = nlpdv(ipos + ji - 1,1)
        if (ir /= 0) then
          irk = nprkx(ir)
          iss = nspi(irk,1)
          irtr = madvtr(iss)
          if (irtr == 0) then
            irtr = moffam(iss)
          end if
          icount = icount + 1
          iwork(icount,1) = nlpdv(ipos + ji - 1,1)
          nlpdv(ipos + ji - 1,1) = 0
          !$acc loop seq
          do jii = ji + 1, inf, 1
            ir2 = nlpdv(ipos + jii - 1,1)
            if (ir2 /= 0) then
              irk2 = nprkx(ir2)
              is2 = nspi(irk2,1)
              irtr2 = madvtr(is2)
              if (irtr2 == 0) then
                irtr2 = moffam(is2)
              end if
              if (irtr2 == irtr) then
                icount = icount + 1
                iwork(icount,1) = nlpdv(ipos + jii - 1,1)
                nlpdv(ipos + jii - 1,1) = 0
              end if
            end if
          enddo
        end if
      enddo
      !$acc end kernels
      do ji = 1, inf, 1
        indx = ipos + ji - 1
        nlpdv(indx,1) = iwork(ji,1)
        nlpdv(indx,2) = 1
        if (ji /= inf) then
          do inr = 1, inf - ji, 1
            ir = iwork(ji + inr,1)
            ir2 = iwork(ji,1)
            irk = nprkx(ir)
            irk2 = nprkx(ir2)
            iss = nspi(irk,1)
            is2 = nspi(irk2,1)
            irtr = madvtr(iss)
            irtr2 = madvtr(is2)
            if (irtr == 0) then
              irtr = moffam(iss)
            end if
            if (irtr2 == 0) then
              irtr2 = moffam(is2)
            end if
            if (irtr /= irtr2) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Exit_Stmt
              EXIT
            else
              nlpdv(indx,2) = nlpdv(indx,2) + 1
            end if
          enddo
        end if
      enddo
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_inimpct

end module asad_inimpct_mod
